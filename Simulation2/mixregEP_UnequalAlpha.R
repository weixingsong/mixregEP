# ============================================================
# Simulation study: EP mixture regression with unequal
# component-specific shape parameters alpha_1 != alpha_2
#
# Based on the framework of mixregEP_Simulation1.R.
#
# Main features:
#   1. Two-component EP mixture regression.
#   2. Component-specific shape parameters alpha = (alpha1, alpha2).
#   3. Two-dimensional profile-likelihood grid search.
#   4. Component-specific scale parameters sigma_k.
#   5. Candidate alpha values are restricted to 0 < alpha <= 2,
#      where the normal scale-mixture EM representation is valid.
#   6. Runtime is recorded for every Monte Carlo replication and
#      summarized for every simulation setup.
#   7. Results can be saved automatically as RDS/CSV files.
# ============================================================

library(mvtnorm)

# EP density / log-density
denEp <- function(y, mu, sigma, alpha) {
  cst <- alpha / (2 * sigma * gamma(1 / alpha))
  cst * exp(-(abs(y - mu) / sigma)^alpha)
}

logdenEp <- function(y, mu, sigma, alpha) {
  log(alpha) - log(2) - log(sigma) - lgamma(1 / alpha) -
    (abs(y - mu) / sigma)^alpha
}

# EP random generator
rpe <- function(n, mu = 0, sigma = 1, alpha = 2) {
  G <- rgamma(n, shape = 1 / alpha)
  S <- sample(c(-1, 1), n, replace = TRUE)
  mu + sigma * S * G^(1 / alpha)
}

# EP variance when sigma = 1
ep_var_unit_scale <- function(alpha) {
  gamma(3 / alpha) / gamma(1 / alpha)
}

# EP scale giving variance 1
ep_scale_unit_variance <- function(alpha) {
  sqrt(gamma(1 / alpha) / gamma(3 / alpha))
}

# Generate EP errors with variance 1
gen_ep_error <- function(n, alpha) {
  sigma <- ep_scale_unit_variance(alpha)
  rpe(n, mu = 0, sigma = sigma, alpha = alpha)
}

# Student-t errors standardized to variance one
gen_t_error_unit_variance <- function(n, nu) {
  if (nu <= 2) stop("nu must exceed 2.")
  rt(n, df = nu) * sqrt((nu - 2) / nu)
}

# Stable row-wise log-sum-exp
row_log_sum_exp <- function(A) {
  m <- apply(A, 1, max)
  m + log(rowSums(exp(A - m)))
}

# Label switching: align component 1 with beta=(0,1,1)
# and component 2 with beta=(0,-1,-1). Reorder all parameters.
relabel_EP <- function(beta, prob, sigma, alpha) {
  beta0 <- rbind(c(0, 1, 1), c(0, -1, -1))

  d_identity <- sum((beta[1, ] - beta0[1, ])^2) +
                sum((beta[2, ] - beta0[2, ])^2)
  d_swap <- sum((beta[2, ] - beta0[1, ])^2) +
            sum((beta[1, ] - beta0[2, ])^2)

  if (d_swap < d_identity) {
    ord <- c(2, 1)
    beta <- beta[ord, , drop = FALSE]
    prob <- prob[ord]
    sigma <- sigma[ord]
    alpha <- alpha[ord]
  }

  list(beta = beta, prob = prob, sigma = sigma, alpha = alpha)
}

# ============================================================
# EM algorithm for FIXED unequal alpha=(alpha1,alpha2)
# ============================================================
EM_EP_unequal_alpha <- function(
    y, X, alpha,
    max_iter = 500,
    tol = 1e-6,
    resid_eps = 1e-8,
    sigma_floor = 1e-6) {

  n <- length(y)
  group <- 2

  if (length(alpha) != group) stop("alpha must be a vector of length 2.")
  if (any(alpha <= 0 | alpha > 2)) {
    stop("For the scale-mixture EM algorithm, alpha_k must satisfy 0 < alpha_k <= 2.")
  }

  # initialization
  prob <- c(0.5, 0.5)
  beta <- rbind(c(0.1, 0.9, 0.9),
                c(-0.1, -0.9, -0.9))

  clust <- kmeans(y, centers = 2)
  sigma <- sapply(1:2, function(k) {
    z <- y[clust$cluster == k]
    s <- if (length(z) > 1) sd(z) else sd(y)
    if (!is.finite(s) || s <= 0) sd(y) else s
  })
  sigma <- pmax(sigma, sigma_floor)

  old_loglik <- -Inf
  converged <- FALSE

  for (iter in 1:max_iter) {

    # E-step
    log_joint <- matrix(NA_real_, n, group)
    dik <- matrix(NA_real_, n, group)

    for (k in 1:group) {
      mu_k <- as.vector(X %*% beta[k, ])
      resid_abs <- pmax(abs(y - mu_k), resid_eps)

      log_joint[, k] <- log(prob[k]) +
        logdenEp(y, mu_k, sigma[k], alpha[k])

      dik[, k] <- alpha[k] * sigma[k]^(2 - alpha[k]) *
        resid_abs^(alpha[k] - 2)
    }

    log_denom <- row_log_sum_exp(log_joint)
    tik <- exp(log_joint - log_denom)
    wik <- tik * dik

    # M-step: mixing proportions
    prob_new <- colMeans(tik)
    prob_new <- pmax(prob_new, 1e-12)
    prob_new <- prob_new / sum(prob_new)

    beta_new <- beta
    sigma_new <- sigma

    # M-step: beta_k and sigma_k
    for (k in 1:group) {
      W <- wik[, k]
      WX <- sqrt(W) * X
      Wy <- sqrt(W) * y

      bk <- tryCatch(qr.solve(WX, Wy), error = function(e) NULL)

      if (is.null(bk) || any(!is.finite(bk))) {
        XtWX <- crossprod(WX)
        XtWy <- crossprod(WX, Wy)
        lambda <- 1e-6
        bk <- tryCatch(
          solve(XtWX + lambda * diag(ncol(X)), XtWy),
          error = function(e) NULL
        )
      }

      if (is.null(bk) || any(!is.finite(bk))) bk <- beta[k, ]

      beta_new[k, ] <- as.vector(bk)

      resid_new <- y - as.vector(X %*% beta_new[k, ])
      denom <- sum(tik[, k])
      sigma2_k <- sum(tik[, k] * dik[, k] * resid_new^2) / denom
      sigma_new[k] <- sqrt(max(sigma2_k, sigma_floor^2))
    }

    prob <- prob_new
    beta <- beta_new
    sigma <- sigma_new

    # observed-data log-likelihood
    log_joint_new <- matrix(NA_real_, n, group)
    for (k in 1:group) {
      mu_k <- as.vector(X %*% beta[k, ])
      log_joint_new[, k] <- log(prob[k]) +
        logdenEp(y, mu_k, sigma[k], alpha[k])
    }
    latest <- sum(row_log_sum_exp(log_joint_new))

    if (!is.finite(latest)) break
    if (is.finite(old_loglik) && abs(latest - old_loglik) <= tol) {
      converged <- TRUE
      break
    }
    old_loglik <- latest
  }

  tmp <- relabel_EP(beta, prob, sigma, alpha)

  list(
    beta = tmp$beta,
    prob = tmp$prob,
    sigma = tmp$sigma,
    alpha = tmp$alpha,
    loglik = latest,
    iter = iter,
    converged = converged
  )
}

# ============================================================
# Two-dimensional profile likelihood for (alpha1, alpha2)
# ============================================================
profile_EP_unequal_alpha <- function(
    y, X, alpha_grid,
    max_iter = 500,
    tol = 1e-6,
    resid_eps = 1e-8) {

  if (any(alpha_grid <= 0 | alpha_grid > 2)) {
    stop("alpha_grid must be contained in (0,2] for this EM implementation.")
  }

  grid <- expand.grid(
    alpha1 = alpha_grid,
    alpha2 = alpha_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  loglik <- rep(-Inf, nrow(grid))
  fits <- vector("list", nrow(grid))

  t0 <- proc.time()[["elapsed"]]

  for (g in seq_len(nrow(grid))) {
    a <- c(grid$alpha1[g], grid$alpha2[g])

    fit <- tryCatch(
      EM_EP_unequal_alpha(
        y = y, X = X, alpha = a,
        max_iter = max_iter, tol = tol,
        resid_eps = resid_eps
      ),
      error = function(e) NULL
    )

    if (!is.null(fit) && is.finite(fit$loglik)) {
      loglik[g] <- fit$loglik
      fits[[g]] <- fit
    }
  }

  runtime <- proc.time()[["elapsed"]] - t0

  if (all(!is.finite(loglik))) {
    return(list(
      best = NULL,
      grid = transform(grid, loglik = loglik),
      runtime = runtime
    ))
  }

  idx <- which.max(loglik)

  list(
    best = fits[[idx]],
    best_grid_row = grid[idx, , drop = FALSE],
    grid = transform(grid, loglik = loglik),
    runtime = runtime
  )
}

# ============================================================
# One Monte Carlo replication
# Data may be generated from EP errors or Student-t errors.
# The fitted model is always mixregEP with component-specific
# (alpha1, alpha2) selected by 2-D profile likelihood.
# ============================================================
one_replication <- function(
    n, error_family=c("EP","t"), true_param, alpha_grid,
    max_iter=500, tol=1e-6, resid_eps=1e-8) {

  error_family <- match.arg(error_family)
  x <- rmvnorm(n, c(0,0), diag(2))
  X <- cbind(1,x)
  G1 <- runif(n) <= 0.25

  if (error_family=="EP") {
    e1 <- gen_ep_error(n,true_param[1])
    e2 <- gen_ep_error(n,true_param[2])
  } else {
    e1 <- gen_t_error_unit_variance(n,true_param[1])
    e2 <- gen_t_error_unit_variance(n,true_param[2])
  }

  y1 <- x[,1]+x[,2]+e1
  y2 <- -x[,1]-x[,2]+e2
  y <- ifelse(G1,y1,y2)

  prof <- profile_EP_unequal_alpha(
    y,X,alpha_grid,max_iter=max_iter,tol=tol,resid_eps=resid_eps)

  if (is.null(prof$best))
    return(list(success=FALSE,runtime=prof$runtime))

  fit <- prof$best
  list(success=TRUE,beta=fit$beta,prob=fit$prob,sigma=fit$sigma,
       alpha=fit$alpha,loglik=fit$loglik,iterations=fit$iter,
       converged=fit$converged,runtime=prof$runtime)
}

# ============================================================
# Main simulation: EP-generated and t-generated errors
# ============================================================
run_EP_profile_simulation <- function(
    n=200, total=200,
    alpha_grid=seq(0.2,2,by=0.1),
    ep_setups=rbind(c(.5,1),c(.5,1.5),c(.5,2),
                    c(1,1.5),c(1,2),c(1.5,2)),
    t_setups=rbind(c(3,3),c(3,10),c(3,20)),
    seed=56789,max_iter=500,tol=1e-6,resid_eps=1e-8,
    save_prefix=NULL) {

  set.seed(seed)
  setup_table <- rbind(
    data.frame(family="EP",true1=ep_setups[,1],true2=ep_setups[,2]),
    data.frame(family="t", true1=t_setups[,1], true2=t_setups[,2])
  )
  beta_true <- rbind(c(0,1,1),c(0,-1,-1))
  pi_true <- .25
  allrep <- list(); runtab <- list(); sums <- list()

  for (s in seq_len(nrow(setup_table))) {
    fam <- setup_table$family[s]
    tp <- c(setup_table$true1[s],setup_table$true2[s])
    est <- matrix(NA,total,11)
    colnames(est)<-c("beta10","beta11","beta12","beta20","beta21","beta22",
                    "pi1","sigma1","sigma2","alpha1","alpha2")
    runtime <- rep(NA,total); success <- rep(FALSE,total)
    ll<-rep(NA,total); it<-rep(NA,total); conv<-rep(NA,total)
    t0<-proc.time()[["elapsed"]]

    for(j in seq_len(total)) {
      a<-one_replication(n,fam,tp,alpha_grid,max_iter,tol,resid_eps)
      runtime[j]<-a$runtime; success[j]<-isTRUE(a$success)
      if(success[j]) {
        est[j,]<-c(a$beta[1,],a$beta[2,],a$prob[1],a$sigma,a$alpha)
        ll[j]<-a$loglik; it[j]<-a$iterations; conv[j]<-a$converged
      }
      if(j%%10==0 || j==total)
        cat("setup",s,"replication",j,"/",total,
            "runtime",round(runtime[j],2),"sec 
")
    }

    setup_runtime<-proc.time()[["elapsed"]]-t0
    ok<-success & complete.cases(est)
    truth<-c(beta_true[1,],beta_true[2,],pi_true)
    err<-sweep(est[ok,1:7,drop=FALSE],2,truth,"-")
    summary_table<-data.frame(
      parameter=colnames(est)[1:7],truth=truth,
      bias=colMeans(err),mse=colMeans(err^2),
      mcse_bias=apply(err,2,sd)/sqrt(sum(ok)),
      mcse_mse=apply(err^2,2,sd)/sqrt(sum(ok))
    )
    fitted_alpha<-data.frame(
      parameter=c("sigma1","sigma2","alpha1","alpha2"),
      mean=colMeans(est[ok,8:11,drop=FALSE]),
      sd=apply(est[ok,8:11,drop=FALSE],2,sd)
    )
    alpha_accuracy<-NULL
    if(fam=="EP") {
      ae<-sweep(est[ok,10:11,drop=FALSE],2,tp,"-")
      alpha_accuracy<-data.frame(
        parameter=c("alpha1","alpha2"),truth=tp,
        bias=colMeans(ae),mse=colMeans(ae^2),
        mcse_bias=apply(ae,2,sd)/sqrt(sum(ok)),
        mcse_mse=apply(ae^2,2,sd)/sqrt(sum(ok)))
    }
    rt<-data.frame(setup=s,family=fam,true1=tp[1],true2=tp[2],n=n,
      successful_replications=sum(ok),total_replications=total,
      total_runtime_sec=setup_runtime,
      mean_rep_runtime_sec=mean(runtime,na.rm=TRUE),
      median_rep_runtime_sec=median(runtime,na.rm=TRUE),
      sd_rep_runtime_sec=sd(runtime,na.rm=TRUE),
      min_rep_runtime_sec=min(runtime,na.rm=TRUE),
      max_rep_runtime_sec=max(runtime,na.rm=TRUE))

    cat("
=== Setup",s,":",fam,"(",tp[1],",",tp[2],") ===
")
    print(summary_table,row.names=FALSE,digits=5)
    cat("
Fitted EP parameters:
"); print(fitted_alpha,row.names=FALSE,digits=5)
    if(!is.null(alpha_accuracy)) {
      cat("
Alpha accuracy:
"); print(alpha_accuracy,row.names=FALSE,digits=5)
    }
    cat("
Runtime:
"); print(rt,row.names=FALSE,digits=5); cat("
")

    allrep[[s]]<-data.frame(setup=s,family=fam,true1=tp[1],true2=tp[2],
      replication=1:total,n=n,success=success,converged=conv,loglik=ll,
      iterations=it,runtime_sec=runtime,est)
    runtab[[s]]<-rt
    sums[[s]]<-list(estimation=summary_table,fitted_EP=fitted_alpha,
                    alpha_accuracy=alpha_accuracy)
  }

  out<-list(n=n,total=total,alpha_grid=alpha_grid,setup_table=setup_table,
            setup_summaries=sums,
            replication_results=do.call(rbind,allrep),
            runtime_summary=do.call(rbind,runtab))

  if(!is.null(save_prefix)) {
    saveRDS(out,paste0(save_prefix,"_full_results.rds"))
    write.csv(out$replication_results,paste0(save_prefix,"_replication_results.csv"),
              row.names=FALSE)
    write.csv(out$runtime_summary,paste0(save_prefix,"_runtime_summary.csv"),
              row.names=FALSE)
  }
  out
}

# ============================================================
# Plotting functions
# ============================================================
plot_EP_alpha_boxplots <- function(sim_result,family=c("EP","t"),
                                   file=NULL,width=9,height=5.5) {
  family<-match.arg(family)
  d<-sim_result$replication_results
  d<-d[d$family==family & d$success,,drop=FALSE]
  ids<-sort(unique(d$setup))
  info<-sim_result$setup_table[ids,,drop=FALSE]
  labs<-if(family=="EP") paste0("EP(",info$true1,",",info$true2,")")
        else paste0("t(",info$true1,",",info$true2,")")
  if(!is.null(file)){pdf(file,width=width,height=height);on.exit(dev.off(),add=TRUE)}
  op<-par(no.readonly=TRUE);on.exit(par(op),add=TRUE)
  par(mfrow=c(2,1),mar=c(3,5.5,2,1))
  boxplot(lapply(ids,function(s)d$alpha1[d$setup==s]),names=labs,
          ylab=expression(hat(alpha)[1]),#xlab="Error distribution",
          #main=expression("Selected "*alpha[1]),
          ylim= c(min(sim_result$alpha_grid) - 0.1,
                  max(sim_result$alpha_grid) + 0.1),las=1)
  if(family=="EP") for(j in seq_along(ids))
    segments(j-.35,info$true1[j],j+.35,info$true1[j],lty=3)
  abline(h=c(0.5,1,1.5),lty=2)
  boxplot(lapply(ids,function(s)d$alpha2[d$setup==s]),names=labs,
          ylab=expression(hat(alpha)[2]),#xlab="Error distribution",
          #main=expression("Selected "*alpha[2]),
          ylim= c(min(sim_result$alpha_grid) - 0.1,
                  max(sim_result$alpha_grid) + 0.1),las=1)
  if(family=="EP") for(j in seq_along(ids))
    segments(j-.35,info$true2[j],j+.35,info$true2[j],lty=3)
  abline(h=c(1,1.5,2),lty=2)
}

plot_EP_runtime_boxplot <- function(sim_result,file=NULL,width=9,height=5.5) {
  d<-sim_result$replication_results;i<-sim_result$setup_table
  labs<-ifelse(i$family=="EP",paste0("EP(",i$true1,",",i$true2,")"),
               paste0("t(",i$true1,",",i$true2,")"))
  if(!is.null(file)){pdf(file,width=width,height=height);on.exit(dev.off(),add=TRUE)}
  op<-par(no.readonly=TRUE);on.exit(par(op),add=TRUE)
  par(mar=c(5.5,5.5,2,1))
  boxplot(lapply(seq_len(nrow(i)),function(s)d$runtime_sec[d$setup==s & d$success]),
          names=labs,xlab="Error distribution",ylab="Running time (seconds)",
          main="EP profile-likelihood running time",las=1)
}

plot_all_EP_figures <- function(sim_result,prefix=NULL) {
  if(is.null(prefix)) {
    plot_EP_alpha_boxplots(sim_result,"EP")
    #plot_EP_alpha_boxplots(sim_result,"t")
    #plot_EP_runtime_boxplot(sim_result)
  } else {
    plot_EP_alpha_boxplots(sim_result,"EP",paste0(prefix,"_alpha_EP_errors.pdf"))
    #plot_EP_alpha_boxplots(sim_result,"t", paste0(prefix,"_alpha_t_errors.pdf"))
    #plot_EP_runtime_boxplot(sim_result,paste0(prefix,"_runtime.pdf"))
  }
}

# ============================================================
# Recommended runs
# ============================================================
 simEP200 <- run_EP_profile_simulation(
   n=200,total=200,alpha_grid=seq(0.3,2,by=0.1),
   save_prefix="mixregEP_profile_n200")
 plot_all_EP_figures(simEP200)
 plot_all_EP_figures(simEP200,prefix="mixregEP_n200")

 