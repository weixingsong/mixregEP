library(mclust)

# Use both columns
data("banknote", package = "mclust")

Xdata <- as.matrix(banknote[, c(4,7)])  # Six measurements

# Fit two Gaussian components
fit <- Mclust(Xdata, G = 2)

# Add cluster assignments
clusterid <- fit$classification

# Display results
summary(fit)
table(clusterid)

# Plot using the original two columns
plot(Xdata[,1], Xdata[,2],
     col = clusterid,
     pch = 19,
     xlab = "Right",
     ylab = "Diagonal")

# Posterior membership probabilities
head(fit$z)

x1=Xdata[,1][clusterid=="1"]
y1=Xdata[,2][clusterid=="1"]
lm(y1~x1)
var(y1)

#Coefficients:
#  (Intercept)           x1  
#      73.2705       0.5085  
# var(y1)= 0.3403236

x2=Xdata[,1][clusterid=="2"]
y2=Xdata[,2][clusterid=="2"]
lm(y2~x2)
var(y2)

#Coefficients:
#  (Intercept)           x2  
#    143.69412     -0.01653  
# var(y2) = 0.1569008