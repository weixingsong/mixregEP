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
# Error generator
# ===============================

# EP random generator
# ===============================
rpe <- function(n, mu = 0, sigma = 1, alpha = 2) {
  G <- rgamma(n, shape = 1/alpha)
  S <- sample(c(-1, 1), n, replace = TRUE)
  mu + sigma * S * G^(1/alpha)
}

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
  
  # storage
  b107=b117=b127=b207=b217=b227=pi17=rep(0,total)
  dfk = rep(0, total)
  
  for (k in 1:total) {
    
    repeat {
      
      # -------- Data generation --------
      u  <- runif(n)
      p1 <- (u <= 0.25)
      p2 <- 1 - p1
      
      x1 <- rnorm(n)
      x2 <- rnorm(n)
      
      e1 <- gen_error(n, dist)
      e2 <- gen_error(n, dist)
      
      y1 <- x1 + x2 + e1
      y2 <- -x1 - x2 + e2
      
      y <- y1*p1 + y2*p2
      X <- cbind(1, x1, x2)
      
      # -------- Initialization --------
      km <- kmeans(y, 2)
      
      sig01 <- sd(y[km$cluster==1])
      sig02 <- sd(y[km$cluster==2])
      
      sig0 <- ifelse(sig01==0 | sig02==0,
                     max(sig01,sig02),
                     (sig01+sig02)/2)
      
      beta <- rbind(
        c(-0.1,-0.9,-0.9),
        c(0.1,-0.9,0.9)
      )
      
      pr <- c(0.5,0.5)
      
      lik <- rep(-Inf, maxv)
      results <- vector("list", maxv)
      
      # -------- Loop over df --------
      for (v in 1:maxv) {
        
        fit <- try(
          EM_t_mix(y, X, beta, sig0, pr, v),
          silent = TRUE
        )
        
        if (!inherits(fit, "try-error")) {
          lik[v] <- fit$loglik
          results[[v]] <- fit
        }
      }
      
      if (all(is.finite(lik))) break
    }
    
    # -------- Select best df --------
    fv <- which.max(lik)
    best <- results[[fv]]
    
    beta <- best$beta
    pr   <- best$pr
    
    # -------- Store --------
    b107[k]=beta[1,1]; b117[k]=beta[1,2]; b127[k]=beta[1,3]
    b207[k]=beta[2,1]; b217[k]=beta[2,2]; b227[k]=beta[2,3]
    pi17[k]=pr[1]
    dfk[k]=fv
  }
  
  ############################################
  # TRUE PARAMETERS
  ############################################
  tb10=0; tb11=1; tb12=1
  tb20=0; tb21=-1; tb22=-1
  tpi1=0.25
  
  ############################################
  # BV Calculation
  ############################################
  RL <- cbind(b107 - tb10, b117 - tb11, b127 - tb12,
              b207 - tb20, b217 - tb21, b227 - tb22,
              pi17 - tpi1)
  
  BV <- as.matrix(cbind(
    apply(RL, 2, function(x) mean(x^2)),   # MSE
    apply(RL, 2, mean)                     # Bias
  ))
  
  ############################################
  # DF summary
  ############################################
  Df <- c(mean(dfk),sd(dfk))
  
  list(
    BV = round(BV,3),
    Df = round(Df,3),
    raw = list(beta = cbind(b107,b117,b127,b207,b217,b227),
               pi = pi17,
               df = dfk)
  )
}


distname=c("t(3)","t(10)","t(20)","EP(0.5)","EP(1)","EP(1.5)","EP(1.8)","EP(2)")


n <- 200
total <- 200
distseq <- seq(length(distname))

DFK <- matrix(0, nrow=total, ncol=length(distname))

for(dist in distseq)
{  
  mysim=run_simulation(n, total, dist, maxv = 25)
  cat("distribution=",distname[dist],"\n")
  print(mysim$BV)
  print(mysim$Df)
  DFK[,dist]=mysim$raw[[3]]
}

boxplot.matrix(DFK, names = distname)
abline(h=c(3,10,20),lty=3)
