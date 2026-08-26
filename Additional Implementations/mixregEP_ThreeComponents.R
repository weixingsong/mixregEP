# ============================================================
# Three-component EP mixture regression with component-specific
# shape parameters selected by three-dimensional profile likelihood
#
# Adapted from mixregEP_UnequalAlpha_ProfileSimulation.R
# ============================================================


# One computational point is important: if alpha_grid = seq(0.5, 2, by = 0.1) 
# then there are 16^3=4096 candidate \alpha_1,\alpha_2,\alpha_3) combinations 
# for each Monte Carlo replication, so the full simulation can be very expensive. 
# I included a coarse-grid quick-test example in the file.


library(mvtnorm)

denEp <- function(y, mu, sigma, alpha) {
  cst <- alpha / (2 * sigma * gamma(1 / alpha))
  cst * exp(-(abs(y - mu) / sigma)^alpha)
}

logdenEp <- function(y, mu, sigma, alpha) {
  log(alpha) - log(2) - log(sigma) - lgamma(1 / alpha) -
    (abs(y - mu) / sigma)^alpha
}

rpe <- function(n, mu = 0, sigma = 1, alpha = 2) {
  G <- rgamma(n, shape = 1 / alpha)
  S <- sample(c(-1, 1), n, replace = TRUE)
  mu + sigma * S * G^(1 / alpha)
}

ep_scale_unit_variance <- function(alpha) {
  sqrt(gamma(1 / alpha) / gamma(3 / alpha))
}

gen_ep_error <- function(n, alpha) {
  rpe(n, 0, ep_scale_unit_variance(alpha), alpha)
}

gen_t_error_unit_variance <- function(n, nu) {
  if (nu <= 2) stop("nu must exceed 2.")
  rt(n, df = nu) * sqrt((nu - 2) / nu)
}

row_log_sum_exp <- function(A) {
  m <- apply(A, 1, max)
  m + log(rowSums(exp(A - m)))
}

relabel_EP_3 <- function(beta, prob, sigma, alpha) {
  beta_true <- rbind(
    c(0,  1,  1),
    c(0, -1, -1),
    c(0,  1, -1)
  )
  
  perms <- rbind(
    c(1,2,3), c(1,3,2), c(2,1,3),
    c(2,3,1), c(3,1,2), c(3,2,1)
  )
  
  loss <- apply(perms, 1, function(ord) {
    sum((beta[ord, , drop = FALSE] - beta_true)^2)
  })
  
  ord <- perms[which.min(loss), ]
  
  list(
    beta  = beta[ord, , drop = FALSE],
    prob  = prob[ord],
    sigma = sigma[ord],
    alpha = alpha[ord]
  )
}

EM_EP_three_alpha <- function(
    y, X, alpha,
    max_iter = 500,
    tol = 1e-6,
    resid_eps = 1e-8,
    sigma_floor = 1e-6) {
  
  n <- length(y)
  K <- 3
  
  if (length(alpha) != K) stop("alpha must have length 3.")
  if (any(alpha <= 0 | alpha > 2)) {
    stop("alpha_k must lie in (0,2].")
  }
  
  prob <- c(0.30, 0.30, 0.40)
  
  beta <- rbind(
    c( 0.10,  0.90,  0.90),
    c(-0.10, -0.90, -0.90),
    c( 0.10,  0.90, -0.90)
  )
  
  clust <- kmeans(y, centers = K)
  sigma <- sapply(seq_len(K), function(k) {
    z <- y[clust$cluster == k]
    s <- if (length(z) > 1) sd(z) else sd(y)
    if (!is.finite(s) || s <= 0) s <- sd(y)
    max(s, sigma_floor)
  })
  
  old_loglik <- -Inf
  converged <- FALSE
  latest <- -Inf
  
  for (iter in seq_len(max_iter)) {
    
    log_joint <- matrix(NA_real_, n, K)
    dik <- matrix(NA_real_, n, K)
    
    for (k in seq_len(K)) {
      mu_k <- as.vector(X %*% beta[k, ])
      resid_abs <- pmax(abs(y - mu_k), resid_eps)
      
      log_joint[, k] <- log(prob[k]) +
        logdenEp(y, mu_k, sigma[k], alpha[k])
      
      dik[, k] <- alpha[k] * sigma[k]^(2 - alpha[k]) *
        resid_abs^(alpha[k] - 2)
    }
    
    log_denom <- row_log_sum_exp(log_joint)
    tau <- exp(log_joint - log_denom)
    wik <- tau * dik
    
    prob_new <- colMeans(tau)
    prob_new <- pmax(prob_new, 1e-12)
    prob_new <- prob_new / sum(prob_new)
    
    beta_new <- beta
    sigma_new <- sigma
    
    for (k in seq_len(K)) {
      W <- pmax(wik[, k], 1e-12)
      WX <- sqrt(W) * X
      Wy <- sqrt(W) * y
      
      bk <- tryCatch(qr.solve(WX, Wy), error = function(e) NULL)
      
      if (is.null(bk) || any(!is.finite(bk))) {
        XtWX <- crossprod(WX)
        XtWy <- crossprod(WX, Wy)
        bk <- tryCatch(
          solve(XtWX + 1e-6 * diag(ncol(X)), XtWy),
          error = function(e) NULL
        )
      }
      
      if (is.null(bk) || any(!is.finite(bk))) bk <- beta[k, ]
      
      beta_new[k, ] <- as.vector(bk)
      
      resid_new <- y - as.vector(X %*% beta_new[k, ])
      denom <- sum(tau[, k])
      
      if (is.finite(denom) && denom > 0) {
        sigma2_k <- sum(tau[, k] * dik[, k] * resid_new^2) / denom
        sigma_new[k] <- sqrt(max(sigma2_k, sigma_floor^2))
      }
    }
    
    prob <- prob_new
    beta <- beta_new
    sigma <- sigma_new
    
    log_joint_new <- matrix(NA_real_, n, K)
    for (k in seq_len(K)) {
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
  
  tmp <- relabel_EP_3(beta, prob, sigma, alpha)
  
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

profile_EP_three_alpha <- function(
    y, X, alpha_grid,
    max_iter = 500,
    tol = 1e-6,
    resid_eps = 1e-8) {
  
  if (any(alpha_grid <= 0 | alpha_grid > 2)) {
    stop("alpha_grid must lie in (0,2].")
  }
  
  grid <- expand.grid(
    alpha1 = alpha_grid,
    alpha2 = alpha_grid,
    alpha3 = alpha_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  loglik <- rep(-Inf, nrow(grid))
  fits <- vector("list", nrow(grid))
  
  t0 <- proc.time()[["elapsed"]]
  
  for (g in seq_len(nrow(grid))) {
    a <- c(grid$alpha1[g], grid$alpha2[g], grid$alpha3[g])
    
    fit <- tryCatch(
      EM_EP_three_alpha(
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

one_replication_3comp <- function(
    n,
    error_family = c("EP","t"),
    true_param,
    alpha_grid,
    max_iter = 500,
    tol = 1e-6,
    resid_eps = 1e-8) {
  
  error_family <- match.arg(error_family)
  
  x <- rmvnorm(n, c(0,0), diag(2))
  X <- cbind(1, x)
  
  g <- sample(1:3, n, replace = TRUE, prob = c(0.30,0.30,0.40))
  
  if (error_family == "EP") {
    e1 <- gen_ep_error(n, true_param[1])
    e2 <- gen_ep_error(n, true_param[2])
    e3 <- gen_ep_error(n, true_param[3])
  } else {
    e1 <- gen_t_error_unit_variance(n, true_param[1])
    e2 <- gen_t_error_unit_variance(n, true_param[2])
    e3 <- gen_t_error_unit_variance(n, true_param[3])
  }
  
  y1 <-  x[,1] + x[,2] + e1
  y2 <- -x[,1] - x[,2] + e2
  y3 <-  x[,1] - x[,2] + e3
  
  y <- ifelse(g == 1, y1, ifelse(g == 2, y2, y3))
  
  prof <- profile_EP_three_alpha(
    y, X, alpha_grid,
    max_iter = max_iter,
    tol = tol,
    resid_eps = resid_eps
  )
  
  if (is.null(prof$best)) {
    return(list(success = FALSE, runtime = prof$runtime))
  }
  
  fit <- prof$best
  
  list(
    success = TRUE,
    beta = fit$beta,
    prob = fit$prob,
    sigma = fit$sigma,
    alpha = fit$alpha,
    loglik = fit$loglik,
    iterations = fit$iter,
    converged = fit$converged,
    runtime = prof$runtime
  )
}

run_three_component_simulation <- function(
    n = 200,
    total = 200,
    alpha_grid = seq(0.5,2.0,by=0.1),
    ep_setups = rbind(
      c(0.5,0.5,0.5),
      c(1.0,1.0,1.0),
      c(2.0,2.0,2.0)
    ),
    t_setups = rbind(
      c(3,3,3),
      c(6,6,6),
      c(20,20,20)
    ),
    seed = 56789,
    max_iter = 500,
    tol = 1e-6,
    resid_eps = 1e-8,
    save_prefix = NULL) {
  
  set.seed(seed)
  
  setup_table <- rbind(
    data.frame(family="EP",true1=ep_setups[,1],true2=ep_setups[,2],true3=ep_setups[,3]),
    data.frame(family="t", true1=t_setups[,1], true2=t_setups[,2], true3=t_setups[,3])
  )
  
  beta_true <- rbind(
    c(0, 1, 1),
    c(0,-1,-1),
    c(0, 1,-1)
  )
  
  pi_true <- c(0.30,0.30,0.40)
  
  rep_results <- vector("list", nrow(setup_table))
  setup_results <- vector("list", nrow(setup_table))
  runtime_results <- vector("list", nrow(setup_table))
  fitted_EP_results <- vector("list", nrow(setup_table))
  alpha_accuracy_results <- replicate(nrow(setup_table), data.frame(), simplify = FALSE)
  
  for (s in seq_len(nrow(setup_table))) {
    
    fam <- setup_table$family[s]
    tp <- c(setup_table$true1[s], setup_table$true2[s], setup_table$true3[s])
    
    cat("\n=== Setup",s,":",fam,"(",tp[1],",",tp[2],",",tp[3],") ===\n")
    
    est <- matrix(NA_real_, total, 17)
    colnames(est) <- c(
      "beta10","beta11","beta12",
      "beta20","beta21","beta22",
      "beta30","beta31","beta32",
      "pi1","pi2",
      "sigma1","sigma2","sigma3",
      "alpha1","alpha2","alpha3"
    )
    
    rt <- ll <- rep(NA_real_, total)
    it <- rep(NA_integer_, total)
    conv <- rep(NA, total)
    success <- rep(FALSE, total)
    
    setup_t0 <- proc.time()[["elapsed"]]
    
    for (b in seq_len(total)) {
      ans <- one_replication_3comp(
        n, fam, tp, alpha_grid,
        max_iter, tol, resid_eps
      )
      
      rt[b] <- ans$runtime
      success[b] <- isTRUE(ans$success)
      
      if (success[b]) {
        est[b,] <- c(
          ans$beta[1,], ans$beta[2,], ans$beta[3,],
          ans$prob[1], ans$prob[2],
          ans$sigma, ans$alpha
        )
        ll[b] <- ans$loglik
        it[b] <- ans$iterations
        conv[b] <- ans$converged
      }
      
      if (b %% 10 == 0 || b == total) {
        cat("replication",b,"/",total,
            "runtime",round(rt[b],2),"sec\n")
      }
    }
    
    setup_runtime <- proc.time()[["elapsed"]] - setup_t0
    ok <- success & complete.cases(est)
    B <- sum(ok)
    
    structural_truth <- c(
      beta_true[1,], beta_true[2,], beta_true[3,],
      pi_true[1], pi_true[2]
    )
    
    if (B > 0) {
      err <- sweep(est[ok,1:11,drop=FALSE], 2, structural_truth, "-")
      
      sm <- data.frame(
        parameter = colnames(est)[1:11],
        truth = structural_truth,
        bias = colMeans(err),
        mcse_bias = apply(err,2,sd)/sqrt(B),
        mse = colMeans(err^2),
        mcse_mse = apply(err^2,2,sd)/sqrt(B),
        row.names = NULL
      )
      
      fep <- data.frame(
        parameter = c("sigma1","sigma2","sigma3","alpha1","alpha2","alpha3"),
        mean = colMeans(est[ok,12:17,drop=FALSE]),
        sd = apply(est[ok,12:17,drop=FALSE],2,sd),
        row.names = NULL
      )
    } else {
      sm <- data.frame(
        parameter = colnames(est)[1:11],
        truth = structural_truth,
        bias = NA_real_, mcse_bias = NA_real_,
        mse = NA_real_, mcse_mse = NA_real_
      )
      
      fep <- data.frame(
        parameter = c("sigma1","sigma2","sigma3","alpha1","alpha2","alpha3"),
        mean = NA_real_, sd = NA_real_
      )
    }
    
    # Keep a real data-frame object for every setup.  For t-generated
    # data there is no true EP alpha vector, so alpha accuracy is not
    # defined; use a zero-row data frame instead of NULL.  This avoids
    # R's special behavior whereby assigning NULL with [[ ]] deletes a
    # list element and can later produce "subscript out of bounds".
    aa <- data.frame(
      parameter = character(0),
      truth = numeric(0),
      bias = numeric(0),
      mcse_bias = numeric(0),
      mse = numeric(0),
      mcse_mse = numeric(0),
      stringsAsFactors = FALSE
    )
    
    if (fam == "EP" && B > 0) {
      ae <- sweep(est[ok,15:17,drop=FALSE], 2, tp, "-")
      aa <- data.frame(
        parameter = c("alpha1","alpha2","alpha3"),
        truth = tp,
        bias = colMeans(ae),
        mcse_bias = apply(ae,2,sd)/sqrt(B),
        mse = colMeans(ae^2),
        mcse_mse = apply(ae^2,2,sd)/sqrt(B),
        row.names = NULL
      )
    }
    
    rs <- data.frame(
      setup=s, family=fam,
      true1=tp[1], true2=tp[2], true3=tp[3],
      n=n,
      successful_replications=B,
      total_replications=total,
      total_runtime_sec=setup_runtime,
      mean_rep_runtime_sec=mean(rt,na.rm=TRUE),
      median_rep_runtime_sec=median(rt,na.rm=TRUE),
      sd_rep_runtime_sec=sd(rt,na.rm=TRUE),
      min_rep_runtime_sec=min(rt,na.rm=TRUE),
      max_rep_runtime_sec=max(rt,na.rm=TRUE)
    )
    
    cat("\nBias/MSE/MCSE:\n")
    print(sm,row.names=FALSE,digits=5)
    
    cat("\nFitted EP parameters:\n")
    print(fep,row.names=FALSE,digits=5)
    
    if (nrow(aa) > 0) {
      cat("\nAlpha accuracy:\n")
      print(aa,row.names=FALSE,digits=5)
    }
    
    cat("\nRuntime:\n")
    print(rs,row.names=FALSE,digits=5)
    
    rep_results[[s]] <- data.frame(
      setup=s, family=fam,
      true1=tp[1], true2=tp[2], true3=tp[3],
      replication=1:total, n=n,
      success=success, converged=conv,
      loglik=ll, iterations=it,
      runtime_sec=rt, est
    )
    
    setup_results[[s]] <- sm
    runtime_results[[s]] <- rs
    fitted_EP_results[[s]] <- fep
    # aa is always a data frame (possibly with zero rows), so ordinary
    # double-bracket list assignment is now safe.
    alpha_accuracy_results[[s]] <- aa
  }
  
  out <- list(
    n=n,
    total=total,
    alpha_grid=alpha_grid,
    setup_table=setup_table,
    setup_summaries=setup_results,
    fitted_EP_summaries=fitted_EP_results,
    alpha_accuracy_summaries=alpha_accuracy_results,
    replication_results=do.call(rbind,rep_results),
    runtime_summary=do.call(rbind,runtime_results)
  )
  
  if (!is.null(save_prefix)) {
    saveRDS(out,paste0(save_prefix,"_full_results.rds"))
    write.csv(out$replication_results,
              paste0(save_prefix,"_replication_results.csv"),
              row.names=FALSE)
    write.csv(out$runtime_summary,
              paste0(save_prefix,"_runtime_summary.csv"),
              row.names=FALSE)
    
    for (s in seq_along(setup_results)) {
      write.csv(setup_results[[s]],
                paste0(save_prefix,"_setup",s,"_structural_summary.csv"),
                row.names=FALSE)
      write.csv(fitted_EP_results[[s]],
                paste0(save_prefix,"_setup",s,"_fitted_EP_summary.csv"),
                row.names=FALSE)
      if (length(alpha_accuracy_results) >= s &&
          is.data.frame(alpha_accuracy_results[[s]]) &&
          nrow(alpha_accuracy_results[[s]]) > 0) {
        write.csv(alpha_accuracy_results[[s]],
                  paste0(save_prefix,"_setup",s,"_alpha_accuracy.csv"),
                  row.names=FALSE)
      }
    }
  }
  
  out
}

plot_three_component_alpha <- function(
    sim_result, family=c("EP","t"),
    file=NULL, width=10, height=7) {
  
  family <- match.arg(family)
  d <- sim_result$replication_results
  d <- d[d$family==family & d$success,,drop=FALSE]
  
  ids <- sort(unique(d$setup))
  info <- sim_result$setup_table[ids,,drop=FALSE]
  
  labs <- if (family=="EP") {
    paste0("EP(",info$true1,",",info$true2,",",info$true3,")")
  } else {
    paste0("t(",info$true1,",",info$true2,",",info$true3,")")
  }
  
  if (!is.null(file)) {
    pdf(file,width=width,height=height)
    on.exit(dev.off(),add=TRUE)
  }
  
  op <- par(no.readonly=TRUE)
  on.exit(par(op),add=TRUE)
  
  par(mfrow=c(3,1),mar=c(5.5,5.5,2,1))
  
  for (k in 1:3) {
    varname <- paste0("alpha",k)
    vals <- lapply(ids,function(s)d[[varname]][d$setup==s])
    
    boxplot(
      vals,
      names=labs,
      ylab=bquote(hat(alpha)[.(k)]),
      main=bquote("Selected "*alpha[.(k)]),
      ylim=range(sim_result$alpha_grid),
      las=1
    )
    
    if (family=="EP") {
      tr <- info[[paste0("true",k)]]
      for (j in seq_along(ids)) {
        segments(j-.35,tr[j],j+.35,tr[j],lty=3)
      }
    }
  }
}

plot_three_component_runtime <- function(
    sim_result,file=NULL,width=10,height=5.5) {
  
  d <- sim_result$replication_results
  info <- sim_result$setup_table
  
  labs <- ifelse(
    info$family=="EP",
    paste0("EP(",info$true1,",",info$true2,",",info$true3,")"),
    paste0("t(",info$true1,",",info$true2,",",info$true3,")")
  )
  
  if (!is.null(file)) {
    pdf(file,width=width,height=height)
    on.exit(dev.off(),add=TRUE)
  }
  
  boxplot(
    lapply(seq_len(nrow(info)),
           function(s)d$runtime_sec[d$setup==s & d$success]),
    names=labs,
    xlab="Error distribution",
    ylab="Running time (seconds)",
    main="Three-component EP profile-likelihood running time",
    las=1
  )
}

plot_all_three_component_figures <- function(sim_result,prefix=NULL) {
  if (is.null(prefix)) {
    plot_three_component_alpha(sim_result,"EP")
    plot_three_component_alpha(sim_result,"t")
    plot_three_component_runtime(sim_result)
  } else {
    plot_three_component_alpha(
      sim_result,"EP",
      paste0(prefix,"_alpha_EP_errors.pdf")
    )
    plot_three_component_alpha(
      sim_result,"t",
      paste0(prefix,"_alpha_t_errors.pdf")
    )
    plot_three_component_runtime(
      sim_result,
      paste0(prefix,"_runtime.pdf")
    )
  }
}

# ============================================================
# Suggested runs
# ============================================================

# Quick test:
quick3 <- run_three_component_simulation(
  n=100,
  total=2,
  alpha_grid=seq(0.5,2,by=0.5),
  save_prefix="three_component_quick"
)
plot_all_three_component_figures(quick3)

# Full run:
# CAUTION: seq(0.5,2,by=0.1) gives 16^3=4096 EM fits
# for EACH Monte Carlo replication.
#
# sim3_200 <- run_three_component_simulation(
#   n=200,
#   total=200,
#   alpha_grid=seq(0.5,2,by=0.1),
#   save_prefix="three_component_n200"
# )
#
# sim3_400 <- run_three_component_simulation(
#   n=400,
#   total=200,
#   alpha_grid=seq(0.5,2,by=0.1),
#   save_prefix="three_component_n400"
# )
