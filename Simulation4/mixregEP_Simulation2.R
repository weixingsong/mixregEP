# Simulation Study for
#   Adaptive Robust Mixture Learning via Exponential Power Errors
# Table 3,4,5,6 in Simulation 2,  mixregEP 

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
    n, total, alphaseq, dist_set, mcd)
{
  set.seed(56789)
  
  BV <- matrix(0, length(dist_set)*7, 2)
  AV <- matrix(0, length(dist_set), 2)
  ALP <- matrix(0, nrow=total,ncol=length(dist_set))
  
  for (dist in dist_set) {
    
    b107 = b117 = b127 = b207 = b217 = b227 = pi17 = numeric(total)
    alp <- numeric(total)
    
    for (j in 1:total) {
      
      x <- rmvnorm(n, c(0, 0), diag(2))
      x1 <- x[,1]; x2 <- x[,2]

      u <- runif(n)
      p1 <- (u <= 0.25)
      p2 <- 1 - p1
      
      if (dist == 1) {
       e1 <- (rchisq(n,3) - 3)/sqrt(6)
       e2 <- (rchisq(n,3) - 3)/sqrt(6)
      } else if (dist == 2){
       e1 <- (3 - rchisq(n,3))/sqrt(6)
       e2 <- (3 - rchisq(n,3))/sqrt(6)
      } else if (dist == 3){
       u=runif(n,0,1)  
       e1 = (u <= 0.95) * (rnorm(n,0,1)) + (u > 0.95) * (rnorm(n,0,5))
       u=runif(n,0,1)
       e2 = (u <= 0.95) * (rnorm(n,0,1)) + (u > 0.95) * (rnorm(n,0,5))
      } else
      {
        e1 = rnorm(n,0,1)
        e2 = rnorm(n,0,1)
      }  
      y1 <- x1 + x2 + e1
      y2 <- -x1 - x2 + e2
      y <- y1 * p1 + y2 * p2
      
      if(dist ==4)
      {
       outn=round(n*0.05)
       x1[1:outn]=x2[1:outn]=20
       #y[1:outn]=100
      }
      
      # MCD
      if (mcd=="T")
      {
        X1 = cbind(x1,x2)
        mv=lev.p(X1,method="mcd");
        TX=cbind(X1[,1]-mv$mean[1],X1[,2]-mv$mean[2]);
        ind=which(diag(TX%*%solve(mv$covar)%*%t(TX))>=qchisq(0.975,2));
        if(length(ind)>0)
        {
          X1=X1[-ind,];
          y=y[-ind];
        }
        X <- cbind(1, X1)
      } else
      {
        X = cbind(1,x1,x2)
      }  
      
      likelihood <- numeric(length(alphaseq))
      results_list <- vector("list", length(alphaseq))
      
      for (k in seq_along(alphaseq)) {
        alpha <- alphaseq[k]
        fit <- EM_EP(y, X, alpha)
        
        if (is.finite(fit$loglik)) {
          likelihood[k] <- fit$loglik
          results_list[[k]] <- fit
        } else {
          likelihood[k] <- -Inf
        }
      }
      
      if (all(!is.finite(likelihood))) {
        next   # skip this j
      }
      
      idx <- which.max(likelihood)
      
      alp[j] <- alphaseq[idx]
      
      best <- results_list[[idx]] 
      
      b107[j] <- best$beta[1,1]
      b117[j] <- best$beta[1,2]
      b127[j] <- best$beta[1,3]

      b207[j] <- best$beta[2,1]
      b217[j] <- best$beta[2,2]
      b227[j] <- best$beta[2,3]

      pi17[j] <- best$prob[1]
    }
    
    # truth
    tb10=0; tb11=1; tb12=1
    tb20=0; tb21=-1; tb22=-1
    tpi1=0.25; 
    
    RL <- cbind(b107 - tb10, b117 - tb11, b127 - tb12, 
                b207 - tb20, b217 - tb21, b227 - tb22, 
                pi17 - tpi1)
    
    BV[seq((dist-1)*7+1, dist*7), ] <- cbind(apply(RL, 2, function(x) mean(x^2)),
            apply(RL, 2, mean))
    ALP[,dist]=alp
    AV[dist, ] <- c(mean(alp), var(alp))
  }
  
  list( BV = round(BV, 3),  AV = round(AV, 3), ALP=ALP)
}

distname=c("chi(3)","-chi(3)","y-Outlier","x-y-Outlier")


dist=seq(length(distname))
n=200
total=200
alphaseq=seq(0.2,2.8,by=0.1)

myreg=run_main_simulation(n, total, alphaseq, dist, mcd="T")
round(myreg$BV, 3)
round(myreg$AV, 3)
boxplot.matrix(myreg$ALP,names=distname,ylim=c(0.4,3.0))
abline(h=c(0.5,1,1.5,2,2.5),lty=3)
