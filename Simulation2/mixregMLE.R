# ============================================================
# Simulation study: EP mixture regression with unequal
# component-specific shape parameters alpha_1 != alpha_2
#
# Based on the framework of mixregEP_Simulation1.R.
#
# Main features:
#   1. Two-component EP mixture regression.
#   2. Component-specific shape parameters alpha = (alpha1, alpha2).
#   3. Joint MLE updates from mle.R (no profile grid search).
#   4. Component-specific scale parameters sigma_k.
#   5. Both alpha parameters are estimated continuously by Newton updates.
#      alpha_grid is retained as compatibility metadata, not an MLE constraint.
#   6. Runtime is recorded for every Monte Carlo replication and
#      summarized for every simulation setup.
#   7. Results can be saved automatically as RDS/CSV files.
#   8. Printed tables and CSV results use exactly three decimal places.
#      Numeric return values and RDS files retain full calculation precision.
#   9. Fitted sigma/alpha tables report bias and MSE relative to EP truth.
#      Student-t scenarios have no true EP sigma/alpha; their accuracy is NA.
# ============================================================

library(mvtnorm)

# Format a reporting copy only; preserve numeric data used in calculations.
# Identifiers, sample sizes, counts, and iteration numbers remain integers.
format_results_3dp <- function(tab) {
  integer_columns <- c("setup", "replication", "n", "iterations",
                       "successful_replications", "total_replications")
  for (name in names(tab)) {
    if (is.numeric(tab[[name]]) && !(name %in% integer_columns)) {
      x <- tab[[name]]
      # Avoid displaying small negative values as -0.000.
      x[is.finite(x) & round(x, 3) == 0] <- 0
      tab[[name]] <- trimws(formatC(x, format = "f", digits = 3,
                                    decimal.mark = "."))
      tab[[name]][is.na(x) & !is.nan(x)] <- NA_character_
    }
  }
  tab
}

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
# MLE estimator transplanted from mle.R
# ============================================================
check_loglik <- function(value, stage) {
  if (length(value) != 1L || !is.finite(value)) {
    stop(paste("Non-finite log-likelihood", stage), call. = FALSE)
  }
  invisible(value)
}

# Label only after fitting; move prob, beta, sigma and alpha together.
# Relabeling only beta/prob inside EM would mismatch the component parameters.
label_fit <- function(prob, beta, sigma, alpha) {
  truth <- c(0.25, 0, 1, 1, 0.75, 0, -1, -1)
  theta <- c(prob[1], beta[1, ], prob[2], beta[2, ])
  swapped <- theta[c(5:8, 1:4)]
  order <- if (sum((swapped - truth)^2) < sum((theta - truth)^2)) {
    c(2L, 1L)
  } else c(1L, 2L)
  list(prob = prob[order], beta = beta[order, , drop = FALSE],
       sigma = sigma[order], alpha = alpha[order])
}

# The Newton update used in the MLE script, with the component's own old alpha.
update_alpha <- function(tau, residual, sigma, alpha) {
  z <- residual / sigma
  # For alpha > 0, z^alpha*log(z) and z^alpha*log(z)^2 both tend to 0
  # as z tends to 0. Set these terms to their limits at exact zero residuals.
  term1 <- term2 <- numeric(length(z))
  nz <- z > 0
  term1[nz] <- z[nz]^alpha * log(z[nz])
  term2[nz] <- z[nz]^alpha * log(z[nz])^2
  gafun <- 1 / alpha^2 + 2 * digamma(1 / alpha) / alpha^3 +
    trigamma(1 / alpha) / alpha^4
  pd1 <- sum(tau * (1 / alpha + digamma(1 / alpha) / alpha^2)) -
    sum(tau * term1)
  pd2 <- -sum(gafun * tau) - sum(tau * term2)
  # In particular, component 2 uses alpha2 here, correcting the original
  # alp10 in its Newton denominator to alp20.
  0.001 + (alpha - 0.001) * exp(-pd1 / (pd2 * (alpha - 0.001) + pd1))
}

fit_mle_EP <- function(y, X, max_iter = 101L, tol = 1e-4) {
  y <- as.numeric(y)
  prob <- c(0.5, 0.5)
  beta <- rbind(c(0.01, 0.9, 0.9), c(0.01, -0.9, -0.9))
  sigma <- c(1, 1)
  alpha <- c(1, 1)
  mu <- X %*% t(beta)
  mix_density <- prob[1] * denEp(y, mu[, 1], sigma[1], alpha[1]) +
    prob[2] * denEp(y, mu[, 2], sigma[2], alpha[2])
  lik0 <- sum(log(mix_density))
  check_loglik(lik0, "at initialization")
  converged <- FALSE
  
  for (iteration in seq_len(max_iter)) {
    f1 <- denEp(y, mu[, 1], sigma[1], alpha[1])
    f2 <- denEp(y, mu[, 2], sigma[2], alpha[2])
    fsum <- prob[1] * f1 + prob[2] * f2
    if (any(!is.finite(fsum)) || any(fsum <= 0)) {
      stop("Invalid mixture density in E-step", call. = FALSE)
    }
    tau <- cbind(prob[1] * f1 / fsum, prob[2] * f2 / fsum)
    prob_new <- c(mean(tau[, 1]), 1 - mean(tau[, 1]))
    beta_new <- beta
    sigma_new <- alpha_new <- numeric(2)
    
    for (component in 1:2) {
      residual <- abs(y - mu[, component])
      mass <- sum(tau[, component])
      if (!is.finite(mass) || mass <= 0) {
        stop("Empty mixture component", call. = FALSE)
      }
      weights <- sqrt(tau[, component] *
                        (residual + 1e-6)^(alpha[component] - 2))
      if (any(!is.finite(weights))) {
        stop("Non-finite regression weights", call. = FALSE)
      }
      # Row scaling is equivalent to the previous diag(weights) multiplication.
      Xs <- weights * X
      ys <- weights * y
      beta_new[component, ] <- as.numeric(solve(crossprod(Xs), crossprod(Xs, ys)))
      sigma_new[component] <- sqrt(
        (alpha[component] * sum(tau[, component] * residual^alpha[component]) /
           mass)^(2 / alpha[component]) + 0.0001)
      alpha_new[component] <- update_alpha(tau[, component], residual,
                                           sigma_new[component], alpha[component])
    }
    if (any(!is.finite(c(prob_new, beta_new, sigma_new, alpha_new))) ||
        any(prob_new <= 0) || any(sigma_new <= 0) || any(alpha_new <= 0)) {
      stop("Non-finite or invalid updated parameters", call. = FALSE)
    }
    mu_new <- X %*% t(beta_new)
    lik1 <- sum(log(
      prob_new[1] * denEp(y, mu_new[, 1], sigma_new[1], alpha_new[1]) +
        prob_new[2] * denEp(y, mu_new[, 2], sigma_new[2], alpha_new[2])))
    check_loglik(lik1, paste("at iteration", iteration))
    converged <- abs(lik1 - lik0) < tol
    
    # Store the latest parameters BEFORE testing the stopping condition.
    # The previous script otherwise returned parameters from the preceding step.
    prob <- prob_new
    beta <- beta_new
    sigma <- sigma_new
    alpha <- alpha_new
    mu <- mu_new
    lik0 <- lik1
    if (converged) break
  }
  fit <- label_fit(prob, beta, sigma, alpha)
  c(fit, list(loglik = lik1, iterations = iteration, converged = converged))
}

# ============================================================
# Compatibility adapter for the original profile-fitting interface.
# There is one MLE fit per data set; no alpha grid is searched.
# alpha_grid and resid_eps are accepted for existing callers but are unused.
# The transplanted fitter uses residual + 1e-6 and scale regularizer 0.0001,
# exactly as in mle.R. max_iter and tol retain this script's original defaults.
# Finite iteration-limit fits are retained (converged = FALSE); failed fits
# remain failed replications, with no replacement data sets or extra RNG draws.
# grid/best_grid_row contain the continuous fitted alphas for compatibility;
# they are no longer profile-likelihood evaluations on alpha_grid.
# ============================================================
profile_EP_unequal_alpha <- function(
    y, X, alpha_grid,
    max_iter = 500,
    tol = 1e-6,
    resid_eps = 1e-8) {
  
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    fit_mle_EP(y, X, max_iter = max_iter, tol = tol),
    error = function(e) NULL
  )
  runtime <- proc.time()[["elapsed"]] - t0
  
  if (is.null(fit)) {
    return(list(
      best = NULL,
      grid = data.frame(alpha1 = numeric(), alpha2 = numeric(), loglik = numeric()),
      runtime = runtime
    ))
  }
  
  # Preserve the original fitter's iteration-field name for existing callers.
  fit$iter <- fit$iterations
  fit$iterations <- NULL
  fitted_alpha <- data.frame(alpha1 = fit$alpha[1], alpha2 = fit$alpha[2])
  list(
    best = fit,
    best_grid_row = fitted_alpha,
    grid = transform(fitted_alpha, loglik = fit$loglik),
    runtime = runtime
  )
}

# ============================================================
# One Monte Carlo replication
# Data may be generated from EP errors or Student-t errors.
# The fitted model is always mixregEP with component-specific
# (alpha1, alpha2) estimated jointly using the method in mle.R.
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
# Original function name, defaults, setup order, outputs, and exports retained.
# alpha_grid/resid_eps remain accepted for compatibility; the MLE ignores them.
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
            "runtime",formatC(runtime[j],format="f",digits=3,decimal.mark="."),"sec 
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
    # The EP scale is not the error standard deviation: the generators use
    # sigma_k = sqrt(gamma(1/alpha_k) / gamma(3/alpha_k)) for unit variance.
    # Student-t degrees of freedom do not supply true EP shape/scale values.
    fitted_truth<-if(fam=="EP") c(ep_scale_unit_variance(tp),tp) else rep(NA_real_,4)
    fitted_alpha<-data.frame(
      parameter=c("sigma1","sigma2","alpha1","alpha2"),
      truth=fitted_truth,bias=NA_real_,mse=NA_real_
    )
    if(fam=="EP") {
      fitted_error<-sweep(est[ok,8:11,drop=FALSE],2,fitted_truth,"-")
      fitted_alpha$bias<-colMeans(fitted_error)
      fitted_alpha$mse<-colMeans(fitted_error^2)
    }
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
=== Setup",s,":",fam,"(",formatC(tp[1],format="f",digits=3,decimal.mark="."),
        ",",formatC(tp[2],format="f",digits=3,decimal.mark="."),") ===
")
    print(format_results_3dp(summary_table),row.names=FALSE)
    cat("
Fitted EP parameter accuracy (bias and MSE):
"); print(format_results_3dp(fitted_alpha),row.names=FALSE)
    if(fam=="t")
      cat("EP sigma/alpha truth is undefined for Student-t errors; bias and MSE are NA.\n")
    if(!is.null(alpha_accuracy)) {
      cat("
Alpha accuracy:
"); print(format_results_3dp(alpha_accuracy),row.names=FALSE)
    }
    cat("
Runtime:
"); print(format_results_3dp(rt),row.names=FALSE); cat("
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
    write.csv(format_results_3dp(out$replication_results),paste0(save_prefix,"_replication_results.csv"),
              row.names=FALSE)
    write.csv(format_results_3dp(out$runtime_summary),paste0(save_prefix,"_runtime_summary.csv"),
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
  # Keep the original grid-based plotting range as a minimum, expanding it
  # to include continuous MLEs (which can exceed the former upper bound of 2).
  alpha_ylim <- range(c(sim_result$alpha_grid, d$alpha1, d$alpha2), finite=TRUE) +
    c(-0.1, 0.1)
  if(!is.null(file)){pdf(file,width=width,height=height);on.exit(dev.off(),add=TRUE)}
  op<-par(no.readonly=TRUE);on.exit(par(op),add=TRUE)
  par(mfrow=c(2,1),mar=c(3,5.5,2,1))
  boxplot(lapply(ids,function(s)d$alpha1[d$setup==s]),names=labs,
          ylab=expression(hat(alpha)[1]),#xlab="Error distribution",
          #main=expression("Selected "*alpha[1]),
          #ylim=alpha_ylim,las=1)
          ylim=c(0,2),las=1)
  if(family=="EP") for(j in seq_along(ids))
    segments(j-.35,info$true1[j],j+.35,info$true1[j],lty=3)
  abline(h=c(0.5,1,1.5),lty=2)
  boxplot(lapply(ids,function(s)d$alpha2[d$setup==s]),names=labs,
          ylab=expression(hat(alpha)[2]),#xlab="Error distribution",
          #main=expression("Selected "*alpha[2]),
          #ylim=alpha_ylim,las=1)
          ylim=c(0,2),las=1)
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
          main="EP MLE running time",las=1)
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


