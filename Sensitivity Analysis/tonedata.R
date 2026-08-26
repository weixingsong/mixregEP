# Sensitivity Study for
#   Adaptive Robust Mixture Learning via Exponential Power Errors
# Figure 6 and 7

library(mvtnorm)

# ===============================
# EP density
# ===============================
denEp <- function(y, mu, sigma, alpha) {
  cst <- alpha / (2 * sigma * gamma(1/alpha))
  cst * exp(-(abs(y - mu)/sigma)^alpha)
}

# ===============================
# EM algorithm (FINAL STABLE VERSION)
# ===============================

EM_EP <- function(y, X, alpha, max_iter = 500, tol = 1e-6) {
  n <- length(y)
  group <- 2
  
  # initialization
  prob <- c(0.5, 0.5)
  #beta <- rbind(c(-0.5,2), c(0.89,0.01))
  beta <- rbind(c(0.1,0.9), c(1.85,0.13))
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
      
      #beta_new = ginv(t(WX)%*%WX)%*%t(WX)%*%Wy
      
      #---- Step 1: QR solve ----
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

  btemp=c(prob[1],beta[1,],prob[2],beta[2,]);
  prob=c(btemp[1],btemp[4]);
  beta=rbind(btemp[2:3],btemp[5:6]);

  list(beta = beta, prob = prob, loglik = latest)
}

      tonedata=read.table("tonedata.txt",header=TRUE)
      xx=tonedata$stretchratio
      x=c(xx,rep(1.5,10)) #10 pairs (0,4)
      yy=tonedata$tuned
      y=c(yy,rep(5,10)) #10 pairs (0,4)
      X=cbind(1,x)
      alphaseq = seq(0.2, 1.3, by=0.1)
      
      
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

      alp <- alphaseq[idx]
      
      best <- results_list[[idx]] 
      beta1 = best$beta[1,]
      beta2 = best$beta[2,]
      round(best$beta,3)
      alp
      
      plot(x[1:150], y[1:150], xlab="strechratio", ylab="tuned",ylim=c(min(y),max(y)))
      points(x[151:160],y[151:160],pch="*",cex=1.4)
      
      abline(c(best$beta[1,1], best$beta[1,2]), col=2, lwd=2)
      abline(c(best$beta[2,1], best$beta[2,2]), col=4, lwd=2)
      
      # format numbers to fixed width
      b1 <- format(round(best$beta[1,1],3), nsmall=3, width=6)
      b2 <- format(round(best$beta[1,2],3), nsmall=3, width=6)
      b3 <- format(round(best$beta[2,1],3), nsmall=3, width=6)
      b4 <- format(round(best$beta[2,2],3), nsmall=3, width=6)
      a0 <- format(round(alp,3), nsmall=3)
      
      legend("topright",
             legend = c(
               bquote(beta[10] == .(b1) ~~ beta[11] == .(b2)),
               bquote(beta[20] == .(b3) ~~ beta[21] == .(b4)),
               bquote(alpha == .(a0))
             ),
             col = c(2, 4, NA),
             lwd = c(2, 2, NA),
             lty = c(1, 1, NA),
             bty = "n",
             cex = 0.8,        # ↓ smaller font
             y.intersp = 1.4   # ↑ more vertical spacing
      )
      
      