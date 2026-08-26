# Simulation Study for
#   Adaptive Robust Mixture Learning via Exponential Power Errors
# Figure 2,4 and Table 1, 2 in Simulation 1,  mixregEP 

library(mvtnorm)

# ===============================
# EP density
# ===============================
denEp <- function(y, mu, sigma, alpha) {
  cst <- alpha / (2 * sigma * gamma(1/alpha))
  cst * exp(-(abs(y - mu)/sigma)^alpha)
}

# ===============================
# EP random generator
# ===============================
rpe <- function(n, mu = 0, sigma = 1, alpha = 2) {
  G <- rgamma(n, shape = 1/alpha)
  S <- sample(c(-1, 1), n, replace = TRUE)
  mu + sigma * S * G^(1/alpha)
}

# ===============================
# Error generator
# ===============================
gen_error <- function(n, dist) {
  switch(as.character(dist),
         "1" = rt(n, 3) / sqrt(3),
         "2" = rt(n, 10) / sqrt(5/4),
         "3" = rt(n, 20) / sqrt(10/9),
         "4" = rpe(n, 0, 1, 0.5) / sqrt(120),
         "5" = rpe(n, 0, 1, 1.0) / sqrt(2),
         "6" = rpe(n, 0, 1, 1.5) / sqrt(0.7386),
         "7" = rpe(n, 0, 1, 1.8) / sqrt(0.5640),
         "8" = rpe(n, 0, 1, 2) / sqrt(0.5)
  )
}

# ===============================
# EM algorithm (FINAL STABLE VERSION)
# ===============================

EM_EP <- function(y, X, alpha, max_iter = 500, tol = 1e-6) {
  n <- length(y)
  group <- 2
  
  # initialization
  prob <- c(0.5, 0.5)
  beta <- rbind(c(0.1, 0.9, 0.9), c(-0.1, -0.9, -0.9))
  
  clust <- kmeans(y, 2)
  sig01 <- sd(y[clust$cluster == 1])
  sig02 <- sd(y[clust$cluster == 2])
  
  sig0 <- ifelse(sig01==0 | sig02==0,
                 max(sig01,sig02),
                 (sig01+sig02)/2)
  
  prest <- -10^10
  
  for (run in 1:max_iter) {
    
    tik <- matrix(0, n, group)
    dik <- matrix(0, n, group)
    
    # =========================
    # E-step
    # =========================
    for (k in 1:group) {
      mu_k <- X %*% beta[k, ]
      resid <- abs(y - mu_k)
      
      tik[, k] <- prob[k] * denEp(y, mu_k, sig0, alpha)
      dik[, k] <- alpha * sig0^(2 - alpha) * resid^(alpha - 2)
    }
    
    row_sum <- rowSums(tik)
    row_sum[row_sum == 0] <- 1e-12
    tik <- tik / row_sum
    
    dik[!is.finite(dik)] <- 1e6
    wik <- tik * dik
    
    # stabilize mixing weights
    prob <- colMeans(tik)
    prob <- pmax(prob, 1e-6)
    prob <- prob / sum(prob)
    
    #sse <- numeric(group)
    
    # =========================
    # M-step (robust)
    # =========================
    sse0=0;
    for (k in 1:group) {
      W <- wik[, k]
      W[W < 1e-10] <- 1e-10
      
      beta_old <- beta[k, ]
      
      WX <- sqrt(W) * X
      Wy <- sqrt(W) * y
      
      # ---- Step 1: QR solve ----
      beta_new <- tryCatch({
        qr.solve(WX, Wy)
      }, error = function(e) NULL)
      
      # ---- Step 2: Ridge fallback ----
      if (is.null(beta_new) || any(!is.finite(beta_new))) {
        XtWX <- t(WX) %*% WX
        XtWy <- t(WX) %*% Wy
        
        lambda <- 1e-6
        beta_new <- tryCatch({
          solve(XtWX + lambda * diag(ncol(X)), XtWy)
        }, error = function(e) NULL)
      }
      
      # ---- Step 3: fallback to previous ----
      if (is.null(beta_new) || any(!is.finite(beta_new))) {
        beta_new <- beta_old
      }
      
      beta[k, ] <- beta_new
      
      resid <- y - X %*% beta[k, ]
      sse0 <- sse0+sum(W * resid^2) 
    }
    
    # stabilize sigma
    sig0 <- sqrt(pmax(sse0/n, 1e-6))
    
    # =========================
    # log-likelihood
    # =========================
    gsumk <- matrix(0, n, group)
    for (k in 1:group) {
      gsumk[, k] <- prob[k] * denEp(y, X %*% beta[k, ], sig0, alpha)
    }
    
    lik_vec <- rowSums(gsumk)
    lik_vec[lik_vec <= 0] <- 1e-12
    
    latest <- sum(log(lik_vec))
    
    if (!is.finite(latest) || !is.finite(prest)) break
    if ((latest - prest) <= tol) break
    
    prest <- latest
  }
  
  # Label Swiching
  
  btemp=c(prob[1],beta[1,],prob[2],beta[2,]);
  btemp=slab(btemp);
  prob=c(btemp[1],btemp[5]);
  beta=rbind(btemp[2:4],btemp[6:8]);

  list(beta = beta, prob = prob, loglik = latest)
}

############################################
# Main simulation (as a function)
############################################
run_main_simulation <- function(
    n, total, alphaseq, dist_set)
{
  set.seed(56789)

  BV <- matrix(NA_real_, length(dist_set)*7, 2,
               dimnames=list(NULL,c("MSE","Bias")))
  AV <- matrix(NA_real_, length(dist_set), 3,
               dimnames=list(NULL,c("Mean","SD","MCSE_Mean")))
  ALP <- matrix(NA_real_, nrow=total, ncol=length(dist_set))

  parnames <- c("beta10","beta11","beta12",
                "beta20","beta21","beta22","pi1")
  truth <- c(0,1,1,0,-1,-1,0.25)
  MC_summary <- vector("list", length(dist_set))
  runtime_summary <- vector("list", length(dist_set))
  replication_runtime <- matrix(NA_real_, total, length(dist_set))

  for (ii in seq_along(dist_set)) {
    dist <- dist_set[ii]
    est <- matrix(NA_real_, total, 7, dimnames=list(NULL,parnames))
    alp <- rep(NA_real_, total)
    setup_t0 <- proc.time()[["elapsed"]]

    for (j in 1:total) {
      rep_t0 <- proc.time()[["elapsed"]]

      x <- rmvnorm(n, c(0, 0), diag(2))
      x1 <- x[,1]; x2 <- x[,2]
      X <- cbind(1, x)
      p1 <- (runif(n) <= 0.25); p2 <- 1-p1
      e1 <- gen_error(n, dist); e2 <- gen_error(n, dist)
      y <- (x1+x2+e1)*p1 + (-x1-x2+e2)*p2

      likelihood <- rep(-Inf, length(alphaseq))
      results_list <- vector("list", length(alphaseq))
      for (k in seq_along(alphaseq)) {
        fit <- tryCatch(EM_EP(y, X, alphaseq[k]), error=function(e) NULL)
        if (!is.null(fit) && is.finite(fit$loglik)) {
          likelihood[k] <- fit$loglik
          results_list[[k]] <- fit
        }
      }

      if (any(is.finite(likelihood))) {
        idx <- which.max(likelihood)
        best <- results_list[[idx]]
        alp[j] <- alphaseq[idx]
        est[j,] <- c(best$beta[1,], best$beta[2,], best$prob[1])
      }
      replication_runtime[j,ii] <- proc.time()[["elapsed"]] - rep_t0
    }

    ok <- complete.cases(est) & is.finite(alp)
    B <- sum(ok)
    if (B > 0) {
      E <- est[ok,,drop=FALSE]
      err <- sweep(E,2,truth,"-")
      meanest <- colMeans(E)
      esd <- apply(E,2,sd)
      bias <- colMeans(err)
      mse <- colMeans(err^2)
      mcsum <- data.frame(
        parameter=parnames, truth=truth,
        mean_estimate=meanest,
        empirical_sd=esd,
        mcse_mean_estimate=esd/sqrt(B),
        bias=bias,
        mcse_bias=apply(err,2,sd)/sqrt(B),
        mse=mse,
        mcse_mse=apply(err^2,2,sd)/sqrt(B),
        row.names=NULL)
      BV[(ii-1)*7+seq_len(7),] <- cbind(mse,bias)
      aok <- alp[ok]
      AV[ii,] <- c(mean(aok),sd(aok),sd(aok)/sqrt(B))
    } else {
      mcsum <- data.frame(parameter=parnames,truth=truth,
        mean_estimate=NA,empirical_sd=NA,mcse_mean_estimate=NA,
        bias=NA,mcse_bias=NA,mse=NA,mcse_mse=NA)
    }

    MC_summary[[ii]] <- mcsum
    ALP[,ii] <- alp
    rt <- replication_runtime[,ii]
    runtime_summary[[ii]] <- data.frame(
      distribution=dist,
      successful_replications=B,
      total_replications=total,
      total_runtime_sec=proc.time()[["elapsed"]]-setup_t0,
      mean_rep_runtime_sec=mean(rt,na.rm=TRUE),
      median_rep_runtime_sec=median(rt,na.rm=TRUE),
      sd_rep_runtime_sec=sd(rt,na.rm=TRUE),
      min_rep_runtime_sec=min(rt,na.rm=TRUE),
      max_rep_runtime_sec=max(rt,na.rm=TRUE))

    cat("\nDistribution index =",dist,"\n")
    cat("Parameter estimates, Bias, MSE, and Monte Carlo SEs:\n")
    print(mcsum,row.names=FALSE,digits=5)
    cat("Selected alpha: mean, SD, MCSE(mean):\n")
    print(AV[ii,,drop=FALSE],digits=5)
    cat("Running time (seconds):\n")
    print(runtime_summary[[ii]],row.names=FALSE,digits=5)
  }

  names(MC_summary) <- as.character(dist_set)
  list(BV=round(BV,3), AV=round(AV,3), ALP=ALP,
       MC_summary=MC_summary,
       runtime_summary=do.call(rbind,runtime_summary),
       replication_runtime_sec=replication_runtime)
}


distname=c("t(3)","t(10)","t(20)","EP(0.5)","EP(1)","EP(1.5)","EP(1.8)","EP(2)")
dist=seq(8)
n=200
total=200
alphaseq=seq(0.2,2,by=0.1)

myreg=run_main_simulation(n, total, alphaseq, dist)
round(myreg$BV, 3)
round(myreg$AV, 3) 
myreg$MC_summary
myreg$runtime_summary
boxplot.matrix(myreg$ALP,names=distname,ylim=c(0.2,2.1),ylab = TeX("$\\hat{\\alpha}$"))
abline(h=c(0.5,1,1.5,1.8,2),lty=3)
