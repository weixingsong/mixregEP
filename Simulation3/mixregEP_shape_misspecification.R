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
    alpha_grid=seq(0.5,2,by=0.1),
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
  par(mfrow=c(1,2),mar=c(5,4,2,1))
  boxplot(lapply(ids,function(s)d$alpha1[d$setup==s]),names=labs,
          ylab=expression(hat(alpha)[1]),xlab="Error distribution",
          main=expression("Selected "*alpha[1]),ylim=range(sim_result$alpha_grid),las=2)
  if(family=="EP") for(j in seq_along(ids))
    segments(j-.35,info$true1[j],j+.35,info$true1[j],lty=3)
  boxplot(lapply(ids,function(s)d$alpha2[d$setup==s]),names=labs,
          ylab=expression(hat(alpha)[2]),xlab="Error distribution",
          main=expression("Selected "*alpha[2]),ylim=range(sim_result$alpha_grid),las=2)
  if(family=="EP") for(j in seq_along(ids))
    segments(j-.35,info$true2[j],j+.35,info$true2[j],lty=3)
}

plot_EP_runtime_boxplot <- function(sim_result,file=NULL,width=9,height=5.5) {
  d<-sim_result$replication_results;i<-sim_result$setup_table
  labs<-ifelse(i$family=="EP",paste0("EP(",i$true1,",",i$true2,")"),
               paste0("t(",i$true1,",",i$true2,")"))
  if(!is.null(file)){pdf(file,width=width,height=height);on.exit(dev.off(),add=TRUE)}
  boxplot(lapply(seq_len(nrow(i)),function(s)d$runtime_sec[d$setup==s & d$success]),
          names=labs,xlab="Error distribution",ylab="Running time (seconds)",
          main="EP profile-likelihood running time",las=2)
}

plot_all_EP_figures <- function(sim_result,prefix=NULL) {
  if(is.null(prefix)) {
    plot_EP_alpha_boxplots(sim_result,"EP")
    plot_EP_alpha_boxplots(sim_result,"t")
    plot_EP_runtime_boxplot(sim_result)
  } else {
    plot_EP_alpha_boxplots(sim_result,"EP",paste0(prefix,"_alpha_EP_errors.pdf"))
    plot_EP_alpha_boxplots(sim_result,"t", paste0(prefix,"_alpha_t_errors.pdf"))
    plot_EP_runtime_boxplot(sim_result,paste0(prefix,"_runtime.pdf"))
  }
}

# ============================================================
# Recommended runs
# ============================================================
# simEP200 <- run_EP_profile_simulation(
#   n=200,total=200,alpha_grid=seq(0.5,2,by=0.1),
#   save_prefix="mixregEP_profile_n200")
# plot_all_EP_figures(simEP200)
# plot_all_EP_figures(simEP200,prefix="mixregEP_n200")


# ============================================================
# 10. Comparison of component-specific shape versus common shape
# ============================================================
#
# Goal:
#   (A) Correct/flexible model:
#         alpha_1 and alpha_2 are selected separately by a
#         two-dimensional profile likelihood.
#
#   (B) Restricted/common-shape model:
#         alpha_1 = alpha_2 = alpha,
#         selected by a one-dimensional profile likelihood.
#
# IMPORTANT:
# For the EP-generated unequal-shape setups, model (A) is the
# correctly specified shape model and model (B) is misspecified.
# For t-generated data, neither EP shape model is literally the
# true error distribution; the comparison there is a comparison
# between a more flexible and a restricted EP working model.
#
# Because both models are fitted to EXACTLY THE SAME Monte Carlo
# sample in each replication, loss comparisons are paired.
# This substantially improves the precision of the comparison.
# ============================================================


# ------------------------------------------------------------
# Common-shape profile likelihood
#
# The existing EM_EP_unequal_alpha() can be reused by setting
# alpha = c(a, a).  Component-specific sigma_k are still allowed.
# ------------------------------------------------------------
profile_EP_common_alpha <- function(
    y, X, alpha_grid,
    max_iter = 500,
    tol = 1e-6,
    resid_eps = 1e-8) {

  if (any(alpha_grid <= 0 | alpha_grid > 2)) {
    stop("alpha_grid must be contained in (0,2].")
  }

  loglik <- rep(-Inf, length(alpha_grid))
  fits <- vector("list", length(alpha_grid))

  t0 <- proc.time()[["elapsed"]]

  for (g in seq_along(alpha_grid)) {
    a <- alpha_grid[g]

    fit <- tryCatch(
      EM_EP_unequal_alpha(
        y = y,
        X = X,
        alpha = c(a, a),
        max_iter = max_iter,
        tol = tol,
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
      alpha = NA_real_,
      loglik = -Inf,
      runtime = runtime
    ))
  }

  idx <- which.max(loglik)

  list(
    best = fits[[idx]],
    alpha = alpha_grid[idx],
    loglik = loglik[idx],
    profile = data.frame(
      alpha = alpha_grid,
      loglik = loglik
    ),
    runtime = runtime
  )
}


# ------------------------------------------------------------
# Generate ONE dataset, then fit BOTH models to exactly the
# same observations.
# ------------------------------------------------------------
one_rep_compare_shape_models <- function(
    n,
    error_family = c("EP", "t"),
    true_param,
    alpha_grid,
    max_iter = 500,
    tol = 1e-6,
    resid_eps = 1e-8) {

  error_family <- match.arg(error_family)

  # ---------- Generate one common dataset ----------
  x <- rmvnorm(n, mean = c(0, 0), sigma = diag(2))
  X <- cbind(1, x)

  G1 <- runif(n) <= 0.25

  if (error_family == "EP") {
    e1 <- gen_ep_error(n, true_param[1])
    e2 <- gen_ep_error(n, true_param[2])
  } else {
    e1 <- gen_t_error_unit_variance(n, true_param[1])
    e2 <- gen_t_error_unit_variance(n, true_param[2])
  }

  y1 <- x[, 1] + x[, 2] + e1
  y2 <- -x[, 1] - x[, 2] + e2
  y <- ifelse(G1, y1, y2)

  # ---------- Flexible unequal-shape fit ----------
  fit_u <- profile_EP_unequal_alpha(
    y = y,
    X = X,
    alpha_grid = alpha_grid,
    max_iter = max_iter,
    tol = tol,
    resid_eps = resid_eps
  )

  # ---------- Restricted common-shape fit ----------
  fit_c <- profile_EP_common_alpha(
    y = y,
    X = X,
    alpha_grid = alpha_grid,
    max_iter = max_iter,
    tol = tol,
    resid_eps = resid_eps
  )

  if (is.null(fit_u$best) || is.null(fit_c$best)) {
    return(list(
      success = FALSE,
      runtime_unequal = fit_u$runtime,
      runtime_common = fit_c$runtime
    ))
  }

  list(
    success = TRUE,
    unequal = fit_u$best,
    common = fit_c$best,
    alpha_common = fit_c$alpha,
    runtime_unequal = fit_u$runtime,
    runtime_common = fit_c$runtime
  )
}


# ------------------------------------------------------------
# Main paired comparison
#
# Estimation loss is quantified in several complementary ways:
#
# 1. Parameter-specific MSE:
#       E[(hat(theta_j)-theta_j)^2]
#
# 2. Paired excess squared error for each parameter:
#       D_j^(b) =
#       (hat(theta_j,common)-theta_j)^2
#       -
#       (hat(theta_j,unequal)-theta_j)^2.
#
#    Mean(D_j) > 0 means the common-shape restriction loses
#    estimation accuracy.
#
# 3. MSE inflation:
#       100 * (MSE_common/MSE_unequal - 1) percent.
#
# 4. Aggregate beta loss:
#       L_beta = (1/6) sum_{j=1}^6
#                (hat(beta_j)-beta_j)^2.
#
# 5. Aggregate structural loss:
#       L_struct = (1/7) [
#          sum_{j=1}^6 (hat(beta_j)-beta_j)^2
#          + (hat(pi_1)-pi_1)^2
#       ].
#
# 6. For EP-generated cases only, shape loss:
#       unequal:
#         L_alpha,U =
#           { (hat(alpha1)-alpha1)^2
#             +(hat(alpha2)-alpha2)^2 }/2
#
#       common:
#         L_alpha,C =
#           { (hat(alpha)-alpha1)^2
#             +(hat(alpha)-alpha2)^2 }/2.
#
# All loss differences use paired Monte Carlo replications, and
# their MCSE is
#       sd(L_common - L_unequal)/sqrt(B).
# ------------------------------------------------------------
run_shape_misspecification_study <- function(
    n = 200,
    total = 200,
    alpha_grid = seq(0.5, 2.0, by = 0.1),
    ep_setups = rbind(
      c(0.5, 1.0),
      c(0.5, 1.5),
      c(0.5, 2.0),
      c(1.0, 1.5),
      c(1.0, 2.0),
      c(1.5, 2.0)
    ),
    t_setups = rbind(
      c(3, 3),
      c(3, 10),
      c(3, 20)
    ),
    seed = 56789,
    max_iter = 500,
    tol = 1e-6,
    resid_eps = 1e-8,
    save_prefix = NULL) {

  set.seed(seed)

  ep_setups <- as.matrix(ep_setups)
  t_setups <- as.matrix(t_setups)

  setup_table <- rbind(
    data.frame(
      family = "EP",
      true1 = ep_setups[, 1],
      true2 = ep_setups[, 2]
    ),
    data.frame(
      family = "t",
      true1 = t_setups[, 1],
      true2 = t_setups[, 2]
    )
  )

  beta_true <- c(0, 1, 1, 0, -1, -1)
  pi_true <- 0.25
  theta_true <- c(beta_true, pi_true)

  parameter_names <- c(
    "beta10", "beta11", "beta12",
    "beta20", "beta21", "beta22",
    "pi1"
  )

  rep_all <- list()
  param_summary_all <- list()
  loss_summary_all <- list()
  shape_summary_all <- list()
  runtime_all <- list()

  for (s in seq_len(nrow(setup_table))) {

    family <- setup_table$family[s]
    tp <- c(setup_table$true1[s], setup_table$true2[s])

    cat(
      "\n============================================================\n",
      "Shape-model comparison: setup ", s, "/",
      nrow(setup_table), "\n",
      "Error family: ", family,
      "; true parameters = (", tp[1], ", ", tp[2], ")\n",
      "n = ", n, "; replications = ", total, "\n",
      "============================================================\n",
      sep = ""
    )

    # Structural estimates:
    # beta10,beta11,beta12,beta20,beta21,beta22,pi1
    est_u <- matrix(
      NA_real_, nrow = total, ncol = 7,
      dimnames = list(NULL, parameter_names)
    )
    est_c <- est_u

    alpha1_u <- alpha2_u <- alpha_c <- rep(NA_real_, total)

    loglik_u <- loglik_c <- rep(NA_real_, total)
    runtime_u <- runtime_c <- rep(NA_real_, total)
    success <- rep(FALSE, total)

    for (b in seq_len(total)) {

      ans <- one_rep_compare_shape_models(
        n = n,
        error_family = family,
        true_param = tp,
        alpha_grid = alpha_grid,
        max_iter = max_iter,
        tol = tol,
        resid_eps = resid_eps
      )

      runtime_u[b] <- ans$runtime_unequal
      runtime_c[b] <- ans$runtime_common
      success[b] <- isTRUE(ans$success)

      if (success[b]) {

        fu <- ans$unequal
        fc <- ans$common

        est_u[b, ] <- c(
          fu$beta[1, ],
          fu$beta[2, ],
          fu$prob[1]
        )

        est_c[b, ] <- c(
          fc$beta[1, ],
          fc$beta[2, ],
          fc$prob[1]
        )

        alpha1_u[b] <- fu$alpha[1]
        alpha2_u[b] <- fu$alpha[2]
        alpha_c[b] <- ans$alpha_common

        loglik_u[b] <- fu$loglik
        loglik_c[b] <- fc$loglik
      }

      if (b %% 10 == 0 || b == total) {
        cat(
          "replication ", b, "/", total,
          "; unequal runtime = ", round(runtime_u[b], 2),
          " sec; common runtime = ", round(runtime_c[b], 2),
          " sec\n",
          sep = ""
        )
      }
    }

    ok <- success &
      complete.cases(est_u) &
      complete.cases(est_c)

    B <- sum(ok)

    if (B == 0) {
      warning("No successful paired replications for setup ", s)
      next
    }

    # ========================================================
    # Parameter-specific estimation comparison
    # ========================================================
    err_u <- sweep(est_u[ok, , drop = FALSE], 2, theta_true, "-")
    err_c <- sweep(est_c[ok, , drop = FALSE], 2, theta_true, "-")

    sq_u <- err_u^2
    sq_c <- err_c^2
    paired_sq_diff <- sq_c - sq_u

    mse_u <- colMeans(sq_u)
    mse_c <- colMeans(sq_c)

    parameter_comparison <- data.frame(
      setup = s,
      family = family,
      true1 = tp[1],
      true2 = tp[2],
      parameter = parameter_names,

      bias_unequal = colMeans(err_u),
      bias_common = colMeans(err_c),

      mse_unequal = mse_u,
      mse_common = mse_c,

      excess_mse_common_minus_unequal = mse_c - mse_u,

      mse_ratio_common_over_unequal =
        ifelse(mse_u > 0, mse_c / mse_u, NA_real_),

      percent_mse_increase =
        ifelse(
          mse_u > 0,
          100 * (mse_c / mse_u - 1),
          NA_real_
        ),

      mcse_paired_excess_mse =
        apply(paired_sq_diff, 2, sd) / sqrt(B),

      row.names = NULL
    )

    # ========================================================
    # Aggregate estimation losses
    # ========================================================
    beta_loss_u <- rowMeans(sq_u[, 1:6, drop = FALSE])
    beta_loss_c <- rowMeans(sq_c[, 1:6, drop = FALSE])

    pi_loss_u <- sq_u[, 7]
    pi_loss_c <- sq_c[, 7]

    structural_loss_u <- rowMeans(sq_u)
    structural_loss_c <- rowMeans(sq_c)

    summarize_paired_loss <- function(name, loss_u, loss_c) {

      dd <- loss_c - loss_u

      data.frame(
        setup = s,
        family = family,
        true1 = tp[1],
        true2 = tp[2],
        loss = name,
        mean_loss_unequal = mean(loss_u),
        mean_loss_common = mean(loss_c),
        excess_loss_common_minus_unequal = mean(dd),
        loss_ratio_common_over_unequal =
          ifelse(mean(loss_u) > 0,
                 mean(loss_c) / mean(loss_u),
                 NA_real_),
        percent_loss_increase =
          ifelse(mean(loss_u) > 0,
                 100 * (mean(loss_c) / mean(loss_u) - 1),
                 NA_real_),
        mcse_paired_excess_loss =
          sd(dd) / sqrt(length(dd)),
        row.names = NULL
      )
    }

    aggregate_loss <- rbind(
      summarize_paired_loss(
        "beta_MSE_average",
        beta_loss_u,
        beta_loss_c
      ),
      summarize_paired_loss(
        "pi1_squared_error",
        pi_loss_u,
        pi_loss_c
      ),
      summarize_paired_loss(
        "structural_MSE_average",
        structural_loss_u,
        structural_loss_c
      )
    )

    # ========================================================
    # Shape-estimation loss for true EP cases
    # ========================================================
    shape_loss_summary <- NULL

    if (family == "EP") {

      shape_loss_u <- (
        (alpha1_u[ok] - tp[1])^2 +
          (alpha2_u[ok] - tp[2])^2
      ) / 2

      shape_loss_c <- (
        (alpha_c[ok] - tp[1])^2 +
          (alpha_c[ok] - tp[2])^2
      ) / 2

      shape_diff <- shape_loss_c - shape_loss_u

      shape_loss_summary <- data.frame(
        setup = s,
        family = family,
        alpha1_true = tp[1],
        alpha2_true = tp[2],
        mean_shape_loss_unequal = mean(shape_loss_u),
        mean_shape_loss_common = mean(shape_loss_c),
        excess_shape_loss_common_minus_unequal =
          mean(shape_diff),
        shape_loss_ratio_common_over_unequal =
          ifelse(
            mean(shape_loss_u) > 0,
            mean(shape_loss_c) / mean(shape_loss_u),
            NA_real_
          ),
        percent_shape_loss_increase =
          ifelse(
            mean(shape_loss_u) > 0,
            100 * (
              mean(shape_loss_c) /
                mean(shape_loss_u) - 1
            ),
            NA_real_
          ),
        mcse_paired_excess_shape_loss =
          sd(shape_diff) / sqrt(B),
        mean_alpha1_unequal = mean(alpha1_u[ok]),
        mean_alpha2_unequal = mean(alpha2_u[ok]),
        mean_alpha_common = mean(alpha_c[ok]),
        row.names = NULL
      )
    }

    # ========================================================
    # Likelihood and runtime diagnostics
    # ========================================================
    ll_gain <- loglik_u[ok] - loglik_c[ok]

    runtime_summary <- data.frame(
      setup = s,
      family = family,
      true1 = tp[1],
      true2 = tp[2],
      paired_successes = B,
      total_replications = total,

      mean_runtime_unequal_sec =
        mean(runtime_u[ok]),

      mean_runtime_common_sec =
        mean(runtime_c[ok]),

      runtime_ratio_unequal_over_common =
        mean(runtime_u[ok]) /
        mean(runtime_c[ok]),

      mean_loglik_gain_unequal_minus_common =
        mean(ll_gain),

      mcse_loglik_gain =
        sd(ll_gain) / sqrt(B),

      row.names = NULL
    )

    # ========================================================
    # Replication-level paired loss dataset
    # ========================================================
    rep_df <- data.frame(
      setup = s,
      family = family,
      true1 = tp[1],
      true2 = tp[2],
      replication = which(ok),

      beta_loss_unequal = beta_loss_u,
      beta_loss_common = beta_loss_c,
      beta_excess_loss =
        beta_loss_c - beta_loss_u,

      pi_loss_unequal = pi_loss_u,
      pi_loss_common = pi_loss_c,
      pi_excess_loss =
        pi_loss_c - pi_loss_u,

      structural_loss_unequal =
        structural_loss_u,
      structural_loss_common =
        structural_loss_c,
      structural_excess_loss =
        structural_loss_c -
        structural_loss_u,

      alpha1_unequal = alpha1_u[ok],
      alpha2_unequal = alpha2_u[ok],
      alpha_common = alpha_c[ok],

      loglik_unequal = loglik_u[ok],
      loglik_common = loglik_c[ok],
      loglik_gain =
        loglik_u[ok] - loglik_c[ok],

      runtime_unequal_sec = runtime_u[ok],
      runtime_common_sec = runtime_c[ok],

      row.names = NULL
    )

    # Add parameter estimates to paired replication file
    for (j in seq_along(parameter_names)) {
      rep_df[[paste0(parameter_names[j], "_unequal")]] <-
        est_u[ok, j]
      rep_df[[paste0(parameter_names[j], "_common")]] <-
        est_c[ok, j]
    }

    # ========================================================
    # Console output
    # ========================================================
    cat("\nParameter-specific comparison:\n")
    print(
      parameter_comparison,
      row.names = FALSE,
      digits = 5
    )

    cat("\nAggregate estimation-loss comparison:\n")
    print(
      aggregate_loss,
      row.names = FALSE,
      digits = 5
    )

    if (!is.null(shape_loss_summary)) {
      cat("\nShape-loss comparison:\n")
      print(
        shape_loss_summary,
        row.names = FALSE,
        digits = 5
      )
    }

    cat("\nLikelihood/runtime diagnostics:\n")
    print(
      runtime_summary,
      row.names = FALSE,
      digits = 5
    )

    # Store
    rep_all[[s]] <- rep_df
    param_summary_all[[s]] <- parameter_comparison
    loss_summary_all[[s]] <- aggregate_loss
    shape_summary_all[[s]] <- shape_loss_summary
    runtime_all[[s]] <- runtime_summary
  }

  replication_results <- do.call(rbind, rep_all)
  parameter_comparison <- do.call(rbind, param_summary_all)
  aggregate_loss_comparison <- do.call(rbind, loss_summary_all)

  nonnull_shape <- shape_summary_all[
    !vapply(shape_summary_all, is.null, logical(1))
  ]

  shape_loss_comparison <-
    if (length(nonnull_shape) > 0)
      do.call(rbind, nonnull_shape)
    else
      NULL

  runtime_likelihood_comparison <-
    do.call(rbind, runtime_all)

  out <- list(
    n = n,
    total = total,
    alpha_grid = alpha_grid,
    setup_table = setup_table,
    replication_results = replication_results,
    parameter_comparison = parameter_comparison,
    aggregate_loss_comparison =
      aggregate_loss_comparison,
    shape_loss_comparison =
      shape_loss_comparison,
    runtime_likelihood_comparison =
      runtime_likelihood_comparison
  )

  # ==========================================================
  # Save all results
  # ==========================================================
  if (!is.null(save_prefix)) {

    saveRDS(
      out,
      paste0(
        save_prefix,
        "_shape_model_comparison.rds"
      )
    )

    write.csv(
      replication_results,
      paste0(
        save_prefix,
        "_paired_replication_results.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      parameter_comparison,
      paste0(
        save_prefix,
        "_parameter_comparison.csv"
      ),
      row.names = FALSE
    )

    write.csv(
      aggregate_loss_comparison,
      paste0(
        save_prefix,
        "_aggregate_loss_comparison.csv"
      ),
      row.names = FALSE
    )

    if (!is.null(shape_loss_comparison)) {
      write.csv(
        shape_loss_comparison,
        paste0(
          save_prefix,
          "_shape_loss_comparison.csv"
        ),
        row.names = FALSE
      )
    }

    write.csv(
      runtime_likelihood_comparison,
      paste0(
        save_prefix,
        "_runtime_likelihood_comparison.csv"
      ),
      row.names = FALSE
    )

    # Excel workbook, if openxlsx is available.
    if (requireNamespace("openxlsx", quietly = TRUE)) {

      wb <- openxlsx::createWorkbook()

      openxlsx::addWorksheet(
        wb, "Parameter_Comparison"
      )
      openxlsx::writeData(
        wb, "Parameter_Comparison",
        parameter_comparison
      )

      openxlsx::addWorksheet(
        wb, "Aggregate_Loss"
      )
      openxlsx::writeData(
        wb, "Aggregate_Loss",
        aggregate_loss_comparison
      )

      if (!is.null(shape_loss_comparison)) {
        openxlsx::addWorksheet(
          wb, "Shape_Loss"
        )
        openxlsx::writeData(
          wb, "Shape_Loss",
          shape_loss_comparison
        )
      }

      openxlsx::addWorksheet(
        wb, "Runtime_Likelihood"
      )
      openxlsx::writeData(
        wb, "Runtime_Likelihood",
        runtime_likelihood_comparison
      )

      openxlsx::addWorksheet(
        wb, "Replication_Results"
      )
      openxlsx::writeData(
        wb, "Replication_Results",
        replication_results
      )

      openxlsx::saveWorkbook(
        wb,
        paste0(
          save_prefix,
          "_shape_model_comparison.xlsx"
        ),
        overwrite = TRUE
      )
    } else {
      message(
        "Package 'openxlsx' is not installed; ",
        "CSV and RDS files were saved, but no Excel file was created."
      )
    }
  }

  out
}


# ============================================================
# 11. Comparison plots
# ============================================================

# ------------------------------------------------------------
# Boxplots of paired structural estimation loss:
# unequal-shape versus common-shape model.
# ------------------------------------------------------------
plot_shape_model_loss <- function(
    sim_result,
    loss = c("beta", "structural"),
    file = NULL,
    width = 10,
    height = 6) {

  loss <- match.arg(loss)
  d <- sim_result$replication_results
  info <- sim_result$setup_table

  if (loss == "beta") {
    uvar <- "beta_loss_unequal"
    cvar <- "beta_loss_common"
    ylab <- "Average squared error of regression coefficients"
    ttl <- "Regression estimation loss"
  } else {
    uvar <- "structural_loss_unequal"
    cvar <- "structural_loss_common"
    ylab <- "Average squared structural-parameter error"
    ttl <- "Overall estimation loss"
  }

  labels <- ifelse(
    info$family == "EP",
    paste0(
      "EP(",
      info$true1, ",",
      info$true2, ")"
    ),
    paste0(
      "t(",
      info$true1, ",",
      info$true2, ")"
    )
  )

  if (!is.null(file)) {
    pdf(file, width = width, height = height)
    on.exit(dev.off(), add = TRUE)
  }

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  par(
    mar = c(8, 4.5, 2.5, 1),
    las = 2
  )

  vals <- list()
  nms <- character(0)

  for (s in seq_len(nrow(info))) {

    ds <- d[d$setup == s, , drop = FALSE]

    vals[[length(vals) + 1]] <- ds[[uvar]]
    vals[[length(vals) + 1]] <- ds[[cvar]]

    nms <- c(
      nms,
      paste0(labels[s], "\nUnequal"),
      paste0(labels[s], "\nCommon")
    )
  }

  boxplot(
    vals,
    names = nms,
    ylab = ylab,
    main = ttl,
    las = 2
  )
}


# ------------------------------------------------------------
# Boxplot of the paired excess loss
#   L_common - L_unequal.
#
# Zero is the reference line.
# Positive values favor the component-specific-shape model.
# ------------------------------------------------------------
plot_excess_loss <- function(
    sim_result,
    loss = c("beta", "structural"),
    file = NULL,
    width = 9,
    height = 5.5) {

  loss <- match.arg(loss)
  d <- sim_result$replication_results
  info <- sim_result$setup_table

  var <- if (loss == "beta") {
    "beta_excess_loss"
  } else {
    "structural_excess_loss"
  }

  labels <- ifelse(
    info$family == "EP",
    paste0(
      "EP(",
      info$true1, ",",
      info$true2, ")"
    ),
    paste0(
      "t(",
      info$true1, ",",
      info$true2, ")"
    )
  )

  vals <- lapply(
    seq_len(nrow(info)),
    function(s)
      d[[var]][d$setup == s]
  )

  if (!is.null(file)) {
    pdf(file, width = width, height = height)
    on.exit(dev.off(), add = TRUE)
  }

  boxplot(
    vals,
    names = labels,
    xlab = "Error distribution",
    ylab = expression(
      L[common] - L[unequal]
    ),
    main = paste(
      "Excess",
      ifelse(loss == "beta",
             "regression", "overall"),
      "estimation loss from common shape"
    ),
    las = 2
  )

  abline(h = 0, lty = 3)
}


# ------------------------------------------------------------
# Shape-loss boxplots for EP-generated cases only.
# ------------------------------------------------------------
plot_shape_loss_paired <- function(
    sim_result,
    file = NULL,
    width = 9,
    height = 5.5) {

  d <- sim_result$replication_results
  info <- sim_result$setup_table

  ep_ids <- which(info$family == "EP")

  labels <- paste0(
    "EP(",
    info$true1[ep_ids], ",",
    info$true2[ep_ids], ")"
  )

  vals <- list()
  nms <- character(0)

  for (j in seq_along(ep_ids)) {

    s <- ep_ids[j]
    ds <- d[d$setup == s, , drop = FALSE]

    a1 <- info$true1[s]
    a2 <- info$true2[s]

    loss_u <- (
      (ds$alpha1_unequal - a1)^2 +
        (ds$alpha2_unequal - a2)^2
    ) / 2

    loss_c <- (
      (ds$alpha_common - a1)^2 +
        (ds$alpha_common - a2)^2
    ) / 2

    vals[[length(vals) + 1]] <- loss_u
    vals[[length(vals) + 1]] <- loss_c

    nms <- c(
      nms,
      paste0(labels[j], "\nUnequal"),
      paste0(labels[j], "\nCommon")
    )
  }

  if (!is.null(file)) {
    pdf(file, width = width, height = height)
    on.exit(dev.off(), add = TRUE)
  }

  boxplot(
    vals,
    names = nms,
    ylab = "Average squared shape-parameter error",
    #main = "Loss from imposing a common EP shape",
    las = 2
  )
}


plot_all_shape_comparison_figures <- function(
    sim_result,
    prefix = NULL) {

  if (is.null(prefix)) {

    plot_shape_model_loss(
      sim_result,
      loss = "beta"
    )

    plot_excess_loss(
      sim_result,
      loss = "beta"
    )

    plot_shape_model_loss(
      sim_result,
      loss = "structural"
    )

    plot_shape_loss_paired(
      sim_result
    )

  } else {

    plot_shape_model_loss(
      sim_result,
      loss = "beta",
      file = paste0(
        prefix,
        "_beta_loss.pdf"
      )
    )

    plot_excess_loss(
      sim_result,
      loss = "beta",
      file = paste0(
        prefix,
        "_beta_excess_loss.pdf"
      )
    )

    plot_shape_model_loss(
      sim_result,
      loss = "structural",
      file = paste0(
        prefix,
        "_structural_loss.pdf"
      )
    )

    plot_shape_loss_paired(
      sim_result,
      file = paste0(
        prefix,
        "_shape_loss.pdf"
      )
    )
  }
}


# ============================================================
# 12. Recommended run
# ============================================================
#
# This replaces the previous one-model simulation when the goal
# is to study common-shape misspecification.
#
 simMiss200 <- run_shape_misspecification_study(
   n = 200,
   total = 200,
   alpha_grid = seq(0.5, 2.0, by = 0.1),
   save_prefix = "EP_shape_misspec_n200"
 )
#
 simMiss400 <- run_shape_misspecification_study(
   n = 400,
   total = 200,
   alpha_grid = seq(0.5, 2.0, by = 0.1),
   save_prefix = "EP_shape_misspec_n400"
 )

 print(simMiss200$parameter_comparison)
 print(simMiss200$aggregate_loss_comparison)
 print(simMiss200$shape_loss_comparison)

 plot_all_shape_comparison_figures(
   simMiss200
 )

 plot_all_shape_comparison_figures(
   simMiss400,
   prefix = "EP_shape_misspec_n200"
 )
#
# ============================================================
