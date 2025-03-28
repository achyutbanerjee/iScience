Temporal dynamics of alien flora accumulation
================
Banerjee AK, Bhowmick AR
July 27, 2023

- [Preparation](#preparation)
- [Use](#use)
- [Naturalized range](#naturalized-range)

**Objective: To identify the naturalized range and categories of uses
that can influence minimum residence time (MRT)**

# Preparation

``` r
#Load the required libraries
library(glmnet)
library(selectiveInference)

#Set the working directory
setwd("~/regression_analysis")
```

# Use

``` r
data = read.csv("use_data.csv")
head(data)
names(data)
#For all alien species irrespective of their invasion categories
y = data$MRT
x = subset(data[,-c(1,2,3)])
x = as.matrix(x)
Ind = apply(x, 2, FUN = var) != 0 # drop predictors with no variation
x = x[,Ind]
colnames(x)
x = scale(x, center = TRUE, scale = TRUE)
# choice of the best tuning parameter lambda for L1 regularization
cv.out = cv.glmnet(x, y, alpha = 1)
plot(cv.out)
bestlam = cv.out$lambda.min
bestlam
out = glmnet(x, y, family = "gaussian", standardize = FALSE)
plot(out)
beta = coef(out, x = x, y = y, s = bestlam/nrow(x), exact = TRUE)[-1]
out_infer = fixedLassoInf(x, y, beta = beta, lambda = bestlam, alpha = 0.1)
print(out_infer)
cat("The final set of variables are\n")
colnames(x)[which(out_infer$pv < 0.1)]

#For three invasion categories separately
data = read.csv("use_data.csv")
status = "Invasive" #change 'Invasive', 'Naturalized', 'Alien'
y = data$MRT[data$Invasion.status == status]
x = subset(data[,-c(1,3)], Invasion.status == status)[,-c(1)]
x = as.matrix(x)
Ind = apply(x, 2, FUN = var) != 0 # drop predictors with no variation
x = x[,Ind]
colnames(x)
x = scale(x, center = TRUE, scale = TRUE)
# choice of the best tuning parameter lambda for L1 regularization
cv.out = cv.glmnet(x, y, alpha = 1)
plot(cv.out)
bestlam = cv.out$lambda.min
bestlam
out = glmnet(x, y, family = "gaussian", standardize = FALSE)
plot(out)
beta = coef(out, x = x, y = y, s = bestlam/nrow(x), exact = TRUE)[-1]
out_infer = fixedLassoInf(x, y, beta = beta, lambda = bestlam, alpha = 0.1)
print(out_infer)
cat("The final set of variables are\n")
colnames(x)[which(out_infer$pv < 0.1)]
```

# Naturalized range

``` r
data = read.csv("naturalized_data.csv")
head(data)
names(data)
#For all alien species irrespective of their invasion categories
y = data$MRT
x = subset(data[,-c(1,2,3)])
x = as.matrix(x)
Ind = apply(x, 2, FUN = var) != 0 # drop predictors with no variation
x = x[,Ind]
colnames(x)
x = scale(x, center = TRUE, scale = TRUE)
# choice of the best tuning parameter lambda for L1 regularization
cv.out = cv.glmnet(x, y, alpha = 1)
plot(cv.out)
bestlam = cv.out$lambda.min
bestlam
out = glmnet(x, y, family = "gaussian", standardize = FALSE)
plot(out)
beta = coef(out, x = x, y = y, s = bestlam/nrow(x), exact = TRUE)[-1]
out_infer = fixedLassoInf(x, y, beta = beta, lambda = bestlam, alpha = 0.1)
print(out_infer)
cat("The final set of variables are\n")
colnames(x)[which(out_infer$pv < 0.1)]

#For three invasion categories separately
data = read.csv("naturalized_data.csv")
status = "Invasive" #change 'Invasive', 'Naturalized', 'Alien'
y = data$MRT[data$Invasion.status == status]
x = subset(data[,-c(1,3)], Invasion.status == status)[,-c(1)]
x = as.matrix(x)
Ind = apply(x, 2, FUN = var) != 0 # drop predictors with no variation
x = x[,Ind]
colnames(x)
x = scale(x, center = TRUE, scale = TRUE)
# choice of the best tuning parameter lambda for L1 regularization
cv.out = cv.glmnet(x, y, alpha = 1)
plot(cv.out)
bestlam = cv.out$lambda.min
bestlam
out = glmnet(x, y, family = "gaussian", standardize = FALSE)
plot(out)
beta = coef(out, x = x, y = y, s = bestlam/nrow(x), exact = TRUE)[-1]
out_infer = fixedLassoInf(x, y, beta = beta, lambda = bestlam, alpha = 0.1)
print(out_infer)
cat("The final set of variables are\n")
colnames(x)[which(out_infer$pv < 0.1)]
```

**END**
