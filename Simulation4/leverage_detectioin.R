lev.p = function(X, method="classic")
 {
  if(method == "sd")
   {
     # Stahel-Donoho Estimator
     # Package needed: rrcov
     Temp=CovSde(X)
     mx = Temp@center
     cx = Temp@cov
   }
  else if(method=="mcd")
   {
    # Fast MCD Estimate
    # Package Needed: robustbase
    Temp=covMcd(X)
    mx=Temp$center
    cx=Temp$cov
   }
  else
   {
    # Sample mean and variance
    mx=apply(X,2,mean)
    cx=var(X)
   }
  list(mean=mx,covar=cx)
 }