# Simulation Study for
#   Adaptive Robust Mixture Learning via Exponential Power Errors
# Figure 3,5 and Table 1, 2 in Simulation 1, mixregT 


############################################
# Robust Mixture Linear Regression (t-dist)
# Optimized Version
############################################

dent_vec <- function(y, mu, sig, v) {
  c0 <- gamma((v + 1) / 2) / (sqrt(pi * v) * gamma(v / 2))
  z <- (y - mu) / sig
  c0 / sig * (1 + z^2 / v)^(-(v + 1) / 2)
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
         "8" = rpe(n, 0, 1, 2.0) / sqrt(0.5)
  )
}

############################################
# EM algorithm for fixed df (v)
############################################
EM_t_mix <- function(y, X, beta, sig, pr, v,
                     max_iter = 500, tol = 1e-6) {
  
  n <- length(y)
  m <- length(pr)
  p <- ncol(X)
  
  r  <- matrix(0, n, m)
  pk <- matrix(0, n, m)
  u  <- matrix(0, n, m)
  
  loglik_old <- -Inf
  
  for (iter in 1:max_iter) {
    
    # ---------- E-step ----------
    for (j in 1:m) {
      mu <- X %*% beta[j, ]
      r[, j] <- (y - mu) / sig
      
      dens <- dent_vec(y, mu, sig, v)
      pk[, j] <- pr[j] * dens
      u[, j]  <- (v + 1) / (v + r[, j]^2)
    }
    
    row_sum <- rowSums(pk)
    row_sum[row_sum < 1e-12] <- 1e-12
    pk <- pk / row_sum
    
    loglik_new <- sum(log(row_sum))
    
    if (abs(loglik_new - loglik_old) < tol) break
    loglik_old <- loglik_new
    
    # ---------- M-step ----------
    pr <- colMeans(pk)
    
    for (j in 1:m) {
      w <- pk[, j] * u[, j]
      
      XtWX <- crossprod(X, X * w)
      XtWy <- crossprod(X, w * y)
      
      beta[j, ] <- solve(XtWX, XtWy)
      
      r[, j] <- y - X %*% beta[j, ]
    }
    
    sig <- sqrt(sum(pk * (r^2) * u) / sum(pk))
  }
  
  # Label Swiching
  
  btemp=c(pr[1],beta[1,],pr[2],beta[2,]);
  btemp=slab(btemp);
  pr=c(btemp[1],btemp[5]);
  beta=rbind(btemp[2:4],btemp[6:8]);
  
  
  list(beta = beta, sig = sig, pr = pr, loglik = loglik_new)
}

############################################
# Main Simulation
############################################


run_simulation <- function(n, total, dist, maxv) {
  set.seed(56789)

  parnames <- c("beta10","beta11","beta12",
                "beta20","beta21","beta22","pi1")
  truth <- c(0,1,1,0,-1,-1,0.25)
  est <- matrix(NA_real_, total, 7, dimnames=list(NULL,parnames))
  dfk <- rep(NA_real_, total)
  rep_runtime <- rep(NA_real_, total)
  sim_t0 <- proc.time()[["elapsed"]]

  for (k in 1:total) {
    rep_t0 <- proc.time()[["elapsed"]]
    repeat {
      p1 <- (runif(n) <= 0.25); p2 <- 1-p1
      x1 <- rnorm(n); x2 <- rnorm(n)
      e1 <- gen_error(n, dist); e2 <- gen_error(n, dist)
      y <- (x1+x2+e1)*p1 + (-x1-x2+e2)*p2
      X <- cbind(1,x1,x2)

      km <- kmeans(y,2)
      sig01 <- sd(y[km$cluster==1]); sig02 <- sd(y[km$cluster==2])
      sig0 <- ifelse(sig01==0 | sig02==0,max(sig01,sig02),(sig01+sig02)/2)
      beta <- rbind(c(-0.1,-0.9,-0.9),c(0.1,-0.9,0.9))
      pr <- c(0.5,0.5)
      lik <- rep(-Inf,maxv); results <- vector("list",maxv)
      for (v in 1:maxv) {
        fit <- try(EM_t_mix(y,X,beta,sig0,pr,v),silent=TRUE)
        if (!inherits(fit,"try-error") && is.finite(fit$loglik)) {
          lik[v] <- fit$loglik; results[[v]] <- fit
        }
      }
      if (any(is.finite(lik))) break
    }

    fv <- which.max(lik); best <- results[[fv]]
    est[k,] <- c(best$beta[1,],best$beta[2,],best$pr[1])
    dfk[k] <- fv
    rep_runtime[k] <- proc.time()[["elapsed"]] - rep_t0
  }

  ok <- complete.cases(est) & is.finite(dfk); B <- sum(ok)
  E <- est[ok,,drop=FALSE]
  err <- sweep(E,2,truth,"-")
  meanest <- colMeans(E); esd <- apply(E,2,sd)
  bias <- colMeans(err); mse <- colMeans(err^2)

  MC_summary <- data.frame(
    parameter=parnames, truth=truth,
    mean_estimate=meanest,
    empirical_sd=esd,
    mcse_mean_estimate=esd/sqrt(B),
    bias=bias,
    mcse_bias=apply(err,2,sd)/sqrt(B),
    mse=mse,
    mcse_mse=apply(err^2,2,sd)/sqrt(B),
    row.names=NULL)

  BV <- cbind(MSE=mse,Bias=bias)
  dfo <- dfk[ok]
  Df <- c(mean=mean(dfo),sd=sd(dfo),mcse_mean=sd(dfo)/sqrt(B))
  runtime_summary <- data.frame(
    distribution=dist,successful_replications=B,total_replications=total,
    total_runtime_sec=proc.time()[["elapsed"]]-sim_t0,
    mean_rep_runtime_sec=mean(rep_runtime,na.rm=TRUE),
    median_rep_runtime_sec=median(rep_runtime,na.rm=TRUE),
    sd_rep_runtime_sec=sd(rep_runtime,na.rm=TRUE),
    min_rep_runtime_sec=min(rep_runtime,na.rm=TRUE),
    max_rep_runtime_sec=max(rep_runtime,na.rm=TRUE))

  cat("\nParameter estimates, Bias, MSE, and Monte Carlo SEs:\n")
  print(MC_summary,row.names=FALSE,digits=5)
  cat("Selected df: mean, SD, MCSE(mean):\n"); print(Df,digits=5)
  cat("Running time (seconds):\n"); print(runtime_summary,row.names=FALSE,digits=5)

  list(BV=round(BV,3),Df=round(Df,3),MC_summary=MC_summary,
       runtime_summary=runtime_summary,
       replication_runtime_sec=rep_runtime,
       raw=list(beta=est[,1:6,drop=FALSE],pi=est[,7],df=dfk))
}


distname=c("t(3)","t(10)","t(20)","EP(0.5)","EP(1)","EP(1.5)","EP(1.8)","EP(2)")


n <- 200
total <- 200
distseq <- seq(length(distname))

DFK <- matrix(0, nrow=total, ncol=length(distname))
runtime_all <- vector("list",length(distname))
MC_all <- vector("list",length(distname))

for(dist in distseq)
{
  mysim=run_simulation(n, total, dist, maxv = 25)
  cat("distribution=",distname[dist],"\n")
  print(mysim$BV)
  print(mysim$Df)
  print(mysim$MC_summary)
  print(mysim$runtime_summary)
  DFK[,dist]=mysim$raw[[3]]
  runtime_all[[dist]] <- transform(mysim$runtime_summary,distribution_name=distname[dist])
  MC_all[[dist]] <- transform(mysim$MC_summary,distribution_name=distname[dist])
}

runtime_summary_all <- do.call(rbind,runtime_all)
MC_summary_all <- do.call(rbind,MC_all)
cat("\nOverall runtime summary:\n")
print(runtime_summary_all,row.names=FALSE,digits=5)

boxplot.matrix(DFK, names = distname,ylab="Fitted df")
abline(h=c(3,10,20),lty=3)
