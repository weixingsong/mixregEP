# ============================================================
# EM Algorithm 2 simulation for
# "Adaptive Robust Mixture Learning via Exponential Power Distribution"
#
# Same simulation models as the profile-likelihood study:
# EP errors: (alpha1,alpha2) =
#   (0.5,1.0), (0.5,1.5), (0.5,2.0),
#   (1.0,1.5), (1.0,2.0), (1.5,2.0)
# t errors: (nu1,nu2) = (3,3), (3,10), (3,20)
#
# Algorithm 2 directly updates component-specific alpha_k by
# numerically maximizing Q_k(alpha_k | S, theta^(t)).
#
# Required packages:
# install.packages(c("mvtnorm", "stabledist", "openxlsx"))
# ============================================================

library(mvtnorm)
if (!requireNamespace("stabledist", quietly=TRUE))
  stop("Install package 'stabledist'.")
if (!requireNamespace("openxlsx", quietly=TRUE))
  stop("Install package 'openxlsx'.")

# ---------- EP density and random generation ----------
logdenEp <- function(y, mu, sigma, alpha) {
  log(alpha)-log(2)-log(sigma)-lgamma(1/alpha)-
    (abs(y-mu)/sigma)^alpha
}

rpe <- function(n, mu=0, sigma=1, alpha=2) {
  G <- rgamma(n, shape=1/alpha)
  S <- sample(c(-1,1), n, replace=TRUE)
  mu + sigma*S*G^(1/alpha)
}

ep_scale_unit_variance <- function(alpha)
  sqrt(gamma(1/alpha)/gamma(3/alpha))

gen_ep_error <- function(n, alpha)
  rpe(n, 0, ep_scale_unit_variance(alpha), alpha)

gen_t_error_unit_variance <- function(n, nu) {
  if (nu <= 2) stop("nu must exceed 2.")
  rt(n, df=nu)*sqrt((nu-2)/nu)
}

# ---------- numerical helpers ----------
row_log_sum_exp <- function(A) {
  m <- apply(A,1,max)
  m + log(rowSums(exp(A-m)))
}

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x-m)))
}

# Positive stable h_gamma with Laplace transform exp(-t^gamma).
# In Nolan S1 form: beta=1, scale=cos(pi*gamma/2)^(1/gamma).
log_h_positive_stable <- function(v, gamma) {
  if (gamma <= 0 || gamma >= 1) return(rep(-Inf,length(v)))
  v <- pmax(v, .Machine$double.xmin)
  sc <- cos(pi*gamma/2)^(1/gamma)
  stabledist::dstable(v, alpha=gamma, beta=1,
                      gamma=sc, delta=0, pm=1, log=TRUE)
}

# Gauss-Laguerre quadrature for int_0^inf exp(-z) f(z) dz.
gauss_laguerre <- function(nq=20) {
  i <- seq_len(nq)
  J <- diag(2*i-1, nq, nq)
  if (nq > 1) {
    off <- seq_len(nq-1)
    J[cbind(1:(nq-1),2:nq)] <- off
    J[cbind(2:nq,1:(nq-1))] <- off
  }
  ee <- eigen(J, symmetric=TRUE)
  ord <- order(ee$values)
  list(nodes=ee$values[ord], weights=ee$vectors[1,ord]^2)
}

# ---------- Q(alpha_k | S, theta^(t)) ----------
# For candidate alpha, terms involving (3/2)log V are constant
# in alpha and are omitted from the maximization.
conditional_V_quadrature <- function(e, alpha_old, lag, e_floor=1e-4) {
  ae <- max(abs(e), e_floor)
  e2 <- ae^2
  v <- lag$nodes/e2
  lw <- log(lag$weights) + log_h_positive_stable(v,alpha_old/2) - log(e2)
  lw <- lw - log_sum_exp(lw)
  list(v=v, w=exp(lw))
}

prepare_component_quadrature <- function(e,tau,alpha_old,lag,e_floor=1e-4) {
  n <- length(e); nq <- length(lag$nodes)
  V <- W <- matrix(NA_real_,n,nq)
  for (i in seq_len(n)) {
    q <- conditional_V_quadrature(e[i],alpha_old,lag,e_floor)
    V[i,] <- q$v; W[i,] <- q$w
  }
  list(V=V,W=W,tau=tau)
}

Q_alpha_component <- function(alpha_candidate, prep) {
  if (!is.finite(alpha_candidate) || alpha_candidate <= 0 || alpha_candidate >= 2)
    return(-Inf)
  nk <- sum(prep$tau)
  const <- nk*(log(alpha_candidate)+0.5*log(pi)-lgamma(1/alpha_candidate))
  lh <- log_h_positive_stable(as.vector(prep$V),alpha_candidate/2)
  if (any(!is.finite(lh))) return(-Inf)
  LH <- matrix(lh,nrow=nrow(prep$V),ncol=ncol(prep$V))
  Ei <- rowSums(prep$W*LH)
  const + sum(prep$tau*Ei)
}

update_alpha_Q <- function(e,tau,alpha_old,lag,
                           alpha_bounds=c(0.25,1.98),
                           e_floor=1e-4,optimize_tol=1e-3) {
  prep <- prepare_component_quadrature(e,tau,alpha_old,lag,e_floor)
  f <- function(a) -Q_alpha_component(a,prep)
  oo <- optimize(f, interval=alpha_bounds, tol=optimize_tol)
  list(alpha=oo$minimum,Q=-oo$objective)
}

# ---------- label switching ----------
relabel_EP <- function(beta,prob,sigma,alpha) {
  b0 <- rbind(c(0,1,1),c(0,-1,-1))
  d1 <- sum((beta[1,]-b0[1,])^2)+sum((beta[2,]-b0[2,])^2)
  d2 <- sum((beta[2,]-b0[1,])^2)+sum((beta[1,]-b0[2,])^2)
  if (d2 < d1) {
    o <- c(2,1); beta <- beta[o,,drop=FALSE]
    prob <- prob[o]; sigma <- sigma[o]; alpha <- alpha[o]
  }
  list(beta=beta,prob=prob,sigma=sigma,alpha=alpha)
}

# ============================================================
# EM Algorithm 2
# ============================================================
EM_EP_algorithm2 <- function(y,X,alpha_init=c(1,1),max_iter=200,tol=1e-5,
                             alpha_bounds=c(0.25,1.98),nq=20,
                             resid_floor=1e-6,e_floor=1e-4,
                             optimize_tol=1e-3,sigma_floor=1e-6,
                             verbose=FALSE) {
  n <- length(y); K <- 2
  if (alpha_bounds[1] <= 0 || alpha_bounds[2] >= 2)
    stop("alpha_bounds must lie strictly inside (0,2).")

  prob <- c(.5,.5)
  beta <- rbind(c(.1,.9,.9),c(-.1,-.9,-.9))
  km <- kmeans(y,2)
  sigma <- sapply(1:2,function(k) {
    z <- y[km$cluster==k]; s <- if(length(z)>1) sd(z) else sd(y)
    if(!is.finite(s)||s<=0) s <- sd(y)
    max(s,sigma_floor)
  })
  alpha <- pmin(pmax(alpha_init,alpha_bounds[1]),alpha_bounds[2])
  lag <- gauss_laguerre(nq)
  oldll <- -Inf; converged <- FALSE; monotone <- TRUE
  llpath <- numeric(0); apath <- matrix(NA,max_iter+1,2); apath[1,] <- alpha

  for (iter in seq_len(max_iter)) {
    # E-step under theta^(t)
    logjoint <- delta <- eold <- matrix(NA_real_,n,K)
    for (k in 1:K) {
      mu <- as.vector(X%*%beta[k,]); r <- y-mu
      ar <- pmax(abs(r),resid_floor)
      logjoint[,k] <- log(prob[k])+logdenEp(y,mu,sigma[k],alpha[k])
      delta[,k] <- alpha[k]*sigma[k]^(2-alpha[k])*ar^(alpha[k]-2)
      eold[,k] <- r/sigma[k]
    }
    ld <- row_log_sum_exp(logjoint)
    tau <- exp(logjoint-ld)

    # M-step for pi, beta, sigma
    pnew <- colMeans(tau); pnew <- pmax(pnew,1e-12); pnew <- pnew/sum(pnew)
    bnew <- beta; snew <- sigma
    for (k in 1:K) {
      w <- tau[,k]*delta[,k]
      WX <- sqrt(w)*X; Wy <- sqrt(w)*y
      bk <- tryCatch(qr.solve(WX,Wy),error=function(e) NULL)
      if (is.null(bk)||any(!is.finite(bk))) {
        A <- crossprod(WX); b <- crossprod(WX,Wy)
        bk <- tryCatch(solve(A+1e-6*diag(ncol(X)),b),error=function(e) NULL)
      }
      if (is.null(bk)||any(!is.finite(bk))) bk <- beta[k,]
      bnew[k,] <- as.vector(bk)
      rnew <- y-as.vector(X%*%bnew[k,])
      s2 <- sum(tau[,k]*delta[,k]*rnew^2)/sum(tau[,k])
      snew[k] <- sqrt(max(s2,sigma_floor^2))
    }

    # M-step for alpha_k: E-step quantities remain frozen at theta^(t)
    anew <- alpha
    for (k in 1:K) {
      au <- tryCatch(update_alpha_Q(eold[,k],tau[,k],alpha[k],lag,
                                    alpha_bounds,e_floor,optimize_tol),
                     error=function(e) NULL)
      if (!is.null(au) && is.finite(au$alpha)) anew[k] <- au$alpha
    }

    prob <- pnew; beta <- bnew; sigma <- snew; alpha <- anew
    apath[iter+1,] <- alpha

    # observed likelihood
    LJ <- matrix(NA_real_,n,K)
    for (k in 1:K) {
      mu <- as.vector(X%*%beta[k,])
      LJ[,k] <- log(prob[k])+logdenEp(y,mu,sigma[k],alpha[k])
    }
    latest <- sum(row_log_sum_exp(LJ)); llpath <- c(llpath,latest)
    if (verbose)
      cat("iter",iter,"loglik",round(latest,5),"alpha",round(alpha,4),"\n")
    if (!is.finite(latest)) break
    if (is.finite(oldll) && latest < oldll-1e-6) monotone <- FALSE
    if (is.finite(oldll) && abs(latest-oldll) <= tol*(1+abs(oldll))) {
      converged <- TRUE; break
    }
    oldll <- latest
  }

  rr <- relabel_EP(beta,prob,sigma,alpha)
  list(beta=rr$beta,prob=rr$prob,sigma=rr$sigma,alpha=rr$alpha,
       loglik=latest,iter=iter,converged=converged,monotone=monotone,
       loglik_path=llpath,alpha_path=apath[seq_len(iter+1),,drop=FALSE])
}

# ---------- one replication ----------
one_rep_algorithm2 <- function(n,error_family=c("EP","t"),true_param,
                               alpha_init=c(1,1),max_iter=200,tol=1e-5,
                               alpha_bounds=c(0.25,1.98),nq=20,
                               optimize_tol=1e-3) {
  error_family <- match.arg(error_family)
  x <- rmvnorm(n,c(0,0),diag(2)); X <- cbind(1,x)
  G1 <- runif(n)<=.25
  if (error_family=="EP") {
    e1 <- gen_ep_error(n,true_param[1]); e2 <- gen_ep_error(n,true_param[2])
  } else {
    e1 <- gen_t_error_unit_variance(n,true_param[1])
    e2 <- gen_t_error_unit_variance(n,true_param[2])
  }
  y <- ifelse(G1,x[,1]+x[,2]+e1,-x[,1]-x[,2]+e2)
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(EM_EP_algorithm2(y,X,alpha_init,max_iter,tol,
                                   alpha_bounds,nq,optimize_tol=optimize_tol),
                  error=function(e) NULL)
  runtime <- proc.time()[["elapsed"]]-t0
  if (is.null(fit)||!is.finite(fit$loglik))
    return(list(success=FALSE,runtime=runtime))
  c(fit,list(success=TRUE,runtime=runtime))
}

# ============================================================
# Main simulation
# ============================================================
run_algorithm2_simulation <- function(
    n=200,total=200,
    ep_setups=rbind(c(.5,1),c(.5,1.5),c(.5,2),c(1,1.5),c(1,2),c(1.5,2)),
    t_setups=rbind(c(3,3),c(3,10),c(3,20)),
    alpha_init=c(1,1),alpha_bounds=c(.25,1.98),
    max_iter=200,tol=1e-5,nq=20,optimize_tol=1e-3,
    seed=56789,save_prefix=NULL) {

  set.seed(seed)
  setups <- rbind(data.frame(family="EP",true1=ep_setups[,1],true2=ep_setups[,2]),
                  data.frame(family="t", true1=t_setups[,1], true2=t_setups[,2]))
  btrue <- rbind(c(0,1,1),c(0,-1,-1)); ptrue <- .25
  reps <- runt <- sums <- vector("list",nrow(setups))

  for (s in seq_len(nrow(setups))) {
    fam <- setups$family[s]; tp <- c(setups$true1[s],setups$true2[s])
    cat("\n=== Setup",s,"/",nrow(setups),fam,"(",tp[1],",",tp[2],") n=",n,"===\n")
    est <- matrix(NA,total,11)
    colnames(est) <- c("beta10","beta11","beta12","beta20","beta21","beta22",
                       "pi1","sigma1","sigma2","alpha1","alpha2")
    rt <- ll <- it <- rep(NA_real_,total); conv <- success <- rep(FALSE,total)
    mono <- rep(NA,total); st <- proc.time()[["elapsed"]]

    for (b in seq_len(total)) {
      a <- one_rep_algorithm2(n,fam,tp,alpha_init,max_iter,tol,
                              alpha_bounds,nq,optimize_tol)
      rt[b] <- a$runtime; success[b] <- isTRUE(a$success)
      if (success[b]) {
        est[b,] <- c(a$beta[1,],a$beta[2,],a$prob[1],a$sigma,a$alpha)
        ll[b] <- a$loglik; it[b] <- a$iter; conv[b] <- a$converged; mono[b] <- a$monotone
      }
      if (b%%10==0 || b==total)
        cat("replication",b,"/",total," runtime",round(rt[b],2),"sec\n")
    }

    setup_runtime <- proc.time()[["elapsed"]]-st
    ok <- success & complete.cases(est)
    truth <- c(btrue[1,],btrue[2,],ptrue)
    err <- sweep(est[ok,1:7,drop=FALSE],2,truth,"-")
    sm <- data.frame(parameter=colnames(est)[1:7],truth=truth,
                     bias=colMeans(err),mse=colMeans(err^2),
                     mcse_bias=apply(err,2,sd)/sqrt(sum(ok)),
                     mcse_mse=apply(err^2,2,sd)/sqrt(sum(ok)))
    fep <- data.frame(parameter=c("sigma1","sigma2","alpha1","alpha2"),
                      mean=colMeans(est[ok,8:11,drop=FALSE]),
                      sd=apply(est[ok,8:11,drop=FALSE],2,sd))
    aa <- NULL
    if (fam=="EP") {
      ae <- sweep(est[ok,10:11,drop=FALSE],2,tp,"-")
      aa <- data.frame(parameter=c("alpha1","alpha2"),truth=tp,
                       bias=colMeans(ae),mse=colMeans(ae^2),
                       mcse_bias=apply(ae,2,sd)/sqrt(sum(ok)),
                       mcse_mse=apply(ae^2,2,sd)/sqrt(sum(ok)))
    }
    rs <- data.frame(setup=s,family=fam,true1=tp[1],true2=tp[2],n=n,
                     successful_replications=sum(ok),total_replications=total,
                     convergence_rate=mean(conv[ok]),monotone_rate=mean(mono[ok],na.rm=TRUE),
                     total_runtime_sec=setup_runtime,mean_rep_runtime_sec=mean(rt,na.rm=TRUE),
                     median_rep_runtime_sec=median(rt,na.rm=TRUE),sd_rep_runtime_sec=sd(rt,na.rm=TRUE),
                     min_rep_runtime_sec=min(rt,na.rm=TRUE),max_rep_runtime_sec=max(rt,na.rm=TRUE))

    cat("\nBias/MSE/MCSE:\n"); print(sm,row.names=FALSE,digits=5)
    cat("\nFitted EP parameters:\n"); print(fep,row.names=FALSE,digits=5)
    if (!is.null(aa)) { cat("\nAlpha accuracy:\n"); print(aa,row.names=FALSE,digits=5) }
    cat("\nRuntime/convergence:\n"); print(rs,row.names=FALSE,digits=5)

    reps[[s]] <- data.frame(setup=s,family=fam,true1=tp[1],true2=tp[2],replication=1:total,
                            n=n,success=success,converged=conv,monotone=mono,iterations=it,
                            loglik=ll,runtime_sec=rt,est)
    runt[[s]] <- rs; sums[[s]] <- list(estimation=sm,fitted_ep=fep,alpha_accuracy=aa)
  }

  out <- list(n=n,total=total,alpha_bounds=alpha_bounds,alpha_init=alpha_init,
              setup_table=setups,setup_summaries=sums,
              replication_results=do.call(rbind,reps),runtime_summary=do.call(rbind,runt))

  if (!is.null(save_prefix)) {
    saveRDS(out,paste0(save_prefix,"_full_results.rds"))
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb,"Replication_Results")
    openxlsx::writeData(wb,"Replication_Results",out$replication_results,withFilter=TRUE)
    openxlsx::addWorksheet(wb,"Runtime_Summary")
    openxlsx::writeData(wb,"Runtime_Summary",out$runtime_summary,withFilter=TRUE)
    for (s in seq_along(sums)) {
      sh <- sprintf("S%02d_Est",s); openxlsx::addWorksheet(wb,sh); openxlsx::writeData(wb,sh,sums[[s]]$estimation)
      sh <- sprintf("S%02d_EP",s); openxlsx::addWorksheet(wb,sh); openxlsx::writeData(wb,sh,sums[[s]]$fitted_ep)
      if (!is.null(sums[[s]]$alpha_accuracy)) {
        sh <- sprintf("S%02d_Alpha",s); openxlsx::addWorksheet(wb,sh); openxlsx::writeData(wb,sh,sums[[s]]$alpha_accuracy)
      }
    }
    openxlsx::saveWorkbook(wb,paste0(save_prefix,"_Algorithm2_results.xlsx"),overwrite=TRUE)
  }
  out
}

# ============================================================
# Boxplots: manuscript-style layout
# ============================================================
plot_algorithm2_alpha <- function(sim_result,family=c("EP","t"),file=NULL,width=9,height=5.5) {
  family <- match.arg(family); d <- sim_result$replication_results
  d <- d[d$family==family & d$success,,drop=FALSE]
  ids <- sort(unique(d$setup)); info <- sim_result$setup_table[ids,,drop=FALSE]
  labs <- if (family=="EP") paste0("EP(",info$true1,",",info$true2,")") else paste0("t(",info$true1,",",info$true2,")")
  if (!is.null(file)) { pdf(file,width=width,height=height); on.exit(dev.off(),add=TRUE) }
  op <- par(no.readonly=TRUE); on.exit(par(op),add=TRUE); par(mfrow=c(1,2),mar=c(5,4,2,1))
  boxplot(lapply(ids,function(s)d$alpha1[d$setup==s]),names=labs,las=2,
          ylab=expression(hat(alpha)[1]),main=expression("Algorithm 2: fitted "*alpha[1]),ylim=sim_result$alpha_bounds)
  if (family=="EP") for(j in seq_along(ids)) if(info$true1[j] <= sim_result$alpha_bounds[2])
    segments(j-.35,info$true1[j],j+.35,info$true1[j],lty=3)
  boxplot(lapply(ids,function(s)d$alpha2[d$setup==s]),names=labs,las=2,
          ylab=expression(hat(alpha)[2]),main=expression("Algorithm 2: fitted "*alpha[2]),ylim=sim_result$alpha_bounds)
  if (family=="EP") for(j in seq_along(ids)) if(info$true2[j] <= sim_result$alpha_bounds[2])
    segments(j-.35,info$true2[j],j+.35,info$true2[j],lty=3)
}

plot_algorithm2_runtime <- function(sim_result,file=NULL,width=9,height=5.5) {
  d <- sim_result$replication_results; info <- sim_result$setup_table
  labs <- ifelse(info$family=="EP",paste0("EP(",info$true1,",",info$true2,")"),paste0("t(",info$true1,",",info$true2,")"))
  if (!is.null(file)) { pdf(file,width=width,height=height); on.exit(dev.off(),add=TRUE) }
  boxplot(lapply(seq_len(nrow(info)),function(s)d$runtime_sec[d$setup==s & d$success]),names=labs,las=2,
          ylab="Running time (seconds)",main="EM Algorithm 2 running time")
}

plot_all_algorithm2 <- function(sim_result,prefix=NULL) {
  if (is.null(prefix)) {
    plot_algorithm2_alpha(sim_result,"EP"); plot_algorithm2_alpha(sim_result,"t"); plot_algorithm2_runtime(sim_result)
  } else {
    plot_algorithm2_alpha(sim_result,"EP",paste0(prefix,"_alpha_EP_errors.pdf"))
    plot_algorithm2_alpha(sim_result,"t", paste0(prefix,"_alpha_t_errors.pdf"))
    plot_algorithm2_runtime(sim_result,paste0(prefix,"_runtime.pdf"))
  }
}

# ============================================================
# Suggested quick test
# ============================================================
# quickA2 <- run_algorithm2_simulation(
#   n=100,total=2,
#   ep_setups=rbind(c(.5,1)),
#   t_setups=rbind(c(3,3)),
#   max_iter=30,nq=12,
#   save_prefix="Algorithm2_quick")
# plot_all_algorithm2(quickA2)

# Full runs:
 simA2_200 <- run_algorithm2_simulation(n=200,total=200,save_prefix="Algorithm2_n200")
# simA2_400 <- run_algorithm2_simulation(n=400,total=200,save_prefix="Algorithm2_n400")
 plot_all_algorithm2(simA2_200,prefix="Algorithm2_n200")
# plot_all_algorithm2(simA2_400,prefix="Algorithm2_n400")
