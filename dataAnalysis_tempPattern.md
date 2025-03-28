Temporal dynamics of alien flora accumulation
================
Banerjee AK, Bhowmick AR
July 27, 2023

- [Preparation](#preparation)
- [Nonlinear functions](#nonlinear-functions)
  - [Fitting](#fitting)
  - [Plotting](#plotting)
- [Breakpoint identification](#breakpoint-identification)
  - [Plotting](#plotting-1)
- [Linear regression](#linear-regression)
  - [Plotting](#plotting-2)
- [Bootstrap analysis](#bootstrap-analysis)
  - [Plotting](#plotting-3)
  - [Bootstrap estimates](#bootstrap-estimates)
- [Invasive aliens](#invasive-aliens)
- [Naturalized aliens](#naturalized-aliens)
- [Introduced aliens](#introduced-aliens)

**Objectives: 1) To analyse the temporal trend of the first record rate
of the alien flora in China; 2) To identify the break-point where the
linear relationship between the first record rate and year changes; 3)
To identify the increasing or decreasing trend of the first record rate
after the break-point; 4) To identify the temporal trend of the first
record rate in future**

**Section 1: For three invasion categories together**

# Preparation

``` r
#Load the required package
library(segmented)

#Set the working directory
setwd("~/time_series_analysis")
```

# Nonlinear functions

## Fitting

``` r
########################################################
##          Fitting nonlinear functions               ##
########################################################
#The following models have been fitted with the observed data:
#(A) Linear model: y = a + bx
#(B) Exponential model: y = a*exp(b*x)
#(C) Saturating model: y = a*(1-exp(-b*x))
#(D) Hyperbolic model: y = a*x/(x + c)
#(E) Sigmoidal model: y = a*x^2/(x^2 + c^2)

D = read.csv(file = "data.csv", header = TRUE)
#D<-subset(D,select = c('Year','Count_introduced'))#un-comment and change into 3 categories
#D = D[-c(14, 17),] #to remove outliers, make this un-comment
d = D                    # just for execution
head(d)
d$Year = 1:nrow(d)
head(d)

#linear fitting
fit_linear = lm(Count ~ Year, data = d)
fit_linear = lm(Count_introduced ~ Year, data = d)#change into 3 categories
summary(fit_linear)
AIC(fit_linear)

#exponential fitting
fit_exponential = nls(Count ~ a*exp(b*Year), data = d,
                      start = list(a = 0.5, b = 0))
fit_exponential = nls(Count_introduced ~ a*exp(b*Year), data = d,
                      start = list(a = 0.5, b = 0))#change into 3 categories
summary(fit_exponential)
AIC(fit_exponential)

#saturating fitting
fit_saturating = nls(Count ~ a*(1-exp(-b*Year)), data = d, 
                     start = list(a = 5, b = 1))
fit_saturating = nls(Count_introduced ~ a*(1-exp(-b*Year)), data = d, 
                     start = list(a = 5, b = 1))#change into 3 categories
summary(fit_saturating)
AIC(fit_saturating)

#hyperbolic fitting
fit_hyperbolic = nls(Count ~ a*Year/(Year + b),
                      data = d, start = list(a =4.5, b = 6.75))
fit_hyperbolic = nls(Count_introduced ~ a*Year/(Year + b),
                      data = d, start = list(a =4.5, b = 6.75))#change into 3 categories
summary(fit_hyperbolic)
AIC(fit_hyperbolic)

#sigmoidal fitting
fit_sigmoidal = nls(Count ~ a*Year^2/(Year^2 + b^2),
                    data = d, start = list(a =4.5, b = 6.75))
fit_sigmoidal = nls(Count_introduced ~ a*Year^2/(Year^2 + b^2),
                    data = d, start = list(a =4.5, b = 6.75))#change into 3 categories
summary(fit_sigmoidal)
AIC(fit_sigmoidal)

#Compare AIC values
AIC_vals = data.frame(AIC(fit_linear), AIC(fit_exponential), AIC(fit_saturating),
             AIC(fit_hyperbolic), AIC(fit_sigmoidal))
names(AIC_vals) = c("Linear", "Exponential", "Saturating", "Hyperbolic", "Sigmoidal")
AIC_vals

#Confidence interval for the saturation point a
confint(fit_saturating)
confint(fit_hyperbolic)
confint(fit_sigmoidal)
```

## Plotting

``` r
la_hat = coefficients(fit_linear)[1] #linear
lb_hat = coefficients(fit_linear)[2] 
ea_hat = coefficients(fit_exponential)[1] #exponential
eb_hat = coefficients(fit_exponential)[2] 
sa_hat = coefficients(fit_saturating)[1] #saturating
sb_hat = coefficients(fit_saturating)[2] 
ha_hat = coefficients(fit_hyperbolic)[1] #hyperbolic
hb_hat = coefficients(fit_hyperbolic)[2] 
sia_hat = coefficients(fit_sigmoidal)[1] #sigmoidal
sib_hat = coefficients(fit_sigmoidal)[2]

{plot(d$Year, d$Count, type = "p", pch = 19, col = "darkgrey",
     xlab = "Year", ylab = "First record rate (per year)",
     ylim = c(-0.5, 33))
  curve(la_hat + lb_hat*x, col = 2, lwd = 2, lty = 1, add = TRUE)
  curve(ea_hat*exp(eb_hat*x), col = 3, lwd = 2, lty = 1,
      add = TRUE)
  curve(sa_hat*(1-exp(-sb_hat*x)), add = TRUE, col = 4,
      lwd = 2, lty = 1)
  curve(ha_hat*x/(x + hb_hat), add = TRUE, col = 5, lty = 1,
      lwd = 2)
  curve(sia_hat*x^2/(x^2 + sib_hat^2), add = TRUE, col = 6, lty = 1,
      lwd = 2)
}
```

# Breakpoint identification

``` r
d<-read.csv("data.csv")
x = d$Year; y = d$Count
lin.mod = lm(y ~ x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=1900)
summary(segmented.mod)

#For three invasion categories
d<-read.csv("data.csv")
d<-subset(d,select = c('Year','Count_introduced'))#change categories
x = d$Year; y = d$Count_introduced #change categories
lin.mod = lm(y ~ x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=1900)
summary(segmented.mod)
```

## Plotting

``` r
x1 = segmented.mod$psi[2] - 2*segmented.mod$psi[3]
x2 = segmented.mod$psi[2] + 2*segmented.mod$psi[3]
x3 = segmented.mod$psi[2] - segmented.mod$psi[3]
x4 = segmented.mod$psi[2] + segmented.mod$psi[3]
{plot(d$Year, d$Count, pch = 19, type = "p", col = "magenta",
     xlab = "Year", ylab = "First record rate (per year)")
  polygon(c(x1,x2, x2, x1),c(32,32, -0.8, -0.8), col = "lightgrey",
        border = NA, lwd = 2, lty = 2)
  polygon(c(x3,x4, x4, x3),c(32,32, -0.8, -0.8), col = "grey",
        lwd = 2, lty = 2, border = NA)
 plot(segmented.mod, add = TRUE, col = "blue")
 points(x,y, pch = 19, type = "p", col = "magenta")
}
```

# Linear regression

``` r
#Linear regression from 1918(+/- SE) to 2020
d<-read.csv("data.csv")
d1<-subset(d,Year >= 1918)
lm_fit_1918 = lm( d1$Count ~ d1$Year)
summary(lm_fit_1918)
d2<-subset(d,Year >= 1926)
lm_fit_1926 = lm( d2$Count ~ d2$Year)
summary(lm_fit_1926)
d3<-subset(d,Year >= 1934)
lm_fit_1934 = lm( d3$Count ~ d3$Year)
summary(lm_fit_1934)
```

## Plotting

``` r
{plot(d1$Year, d1$Count,xlab = "Year",
     ylab = "First record rate (per year)", pch = 19, col = "darkgrey")
  curve(coef(lm_fit_1921)[1]+coef(lm_fit_1921)[2]*x,
      add = TRUE, lwd = 3, col = "blue")
  abline(v = 1921, lwd = 3, col = "blue", lty = 2)
  curve(coef(lm_fit_1930)[1]+coef(lm_fit_1930)[2]*x,
      add = TRUE, lwd = 3, col = "magenta", 1930, 2020)
  abline(v = 1930, lwd = 3, col = "magenta", lty = 2)
  curve(coef(lm_fit_1939)[1]+coef(lm_fit_1939)[2]*x,
      add = TRUE, lwd = 3, col = "red", 1939, 2020)
  abline(v = 1939, lwd = 3, col = "red", lty = 2)
  legend("topright", legend = c("Year 1921 (Estimated)", "Year 1930 (1 S.E.)", "Year 1939 (2 S.E.)"),
       col = c("blue", "magenta", "red"), lwd = c(2,2,2), lty = c(2,2,2),
       bty = "n")
}
```

# Bootstrap analysis

``` r
B = 1000                                 # number of bootstrap samples
new_x_vals = as.data.frame(1939:2040)    # years for prediction
names(new_x_vals) = "Year"
prediction_intervals = matrix(data = NA, # store prediction intervals
                              nrow = 2, ncol = nrow(new_x_vals))
predict_new_x = matrix(data = NA, nrow = B, 
                       ncol = nrow(new_x_vals))
for (i in 1:B) {
  b_ind = sample(1:nrow(d3), replace = TRUE)    # bootstrap indices
  b_data = d3[b_ind, ]                          # bootstrap data set
  names(b_data) = names(d3)                     # converting to data frame
  b_fit = lm(Count ~ Year, data = b_data)       # fitting on bootstrap data
  b_predict_new_x = predict(b_fit,              # predicted values
                            newdata = new_x_vals)
  predict_new_x[i,] = b_predict_new_x
}
for (j in 1:nrow(new_x_vals)) {                 # extract prediction interval
  prediction_intervals[,j] = quantile(predict_new_x[,j],c(1, 99)/100)
}
```

## Plotting

``` r
#Plot for the bootstrap confidence interval
{plot(d3$Year, d3$Count, pch = 21, bg="azure2",cex=1, type = "p",col="black",
     xlab = "Year", ylab = "First record rate (per year)",
     xlim = c(1939, 2040))
  lines(1939:2040, prediction_intervals[1,], col = "firebrick1", lwd= 2)
  lines(1939:2040, prediction_intervals[2,], col = "firebrick1", lwd = 2)
  legend("topleft", "99% prediction interval", lwd = 2,
       col = "red", bty = "n")
}
#Plot for bootstrap distribution of the predicted first record rate in 2030
{hist(predict_new_x[,93], probability = TRUE, main = "", 
     xlab = "Predicted First Record Rate (2030)",
     ylab = "Density",
     breaks = 30, xlim = c(0,6),col = "seashell1")
  lines(density(predict_new_x[,93]), lwd = 3, col = "magenta2", lty = 2)
  curve(dnorm(x, mean = mean(predict_new_x[,93]), sd = sd(predict_new_x[,93])),
      col = "blue", lwd = 3, add = TRUE)
  abline(v = quantile(predict_new_x[,93],c(2.5, 97.5)/100), lwd = 3,
       col = "firebrick2", lty = 3)
  points(mean(predict_new_x[,93]), 0, col = "firebrick2", pch = 19, cex = 1.5)
  legend("topright", legend = c("Kernel density", "Normal density", "95% CI"), lwd = c(2,2, 3),
       col = c("red", "blue", "magenta2"), lty = c(2,1, 3), bty = "n")
}
#Plot for bootstrap distribution of the predicted first record rate in 2040
{hist(predict_new_x[,102], probability = TRUE, main = "", 
     xlab = "Predicted First Record Rate (2040)",
     ylab = "Density",
     breaks = 30, xlim = c(0,6),col = "seashell1")
  lines(density(predict_new_x[,102]), lwd = 2, col = "red", lty = 2)
  curve(dnorm(x, mean = mean(predict_new_x[,102]), sd = sd(predict_new_x[,102])),
      col = "blue", lwd = 2, add = TRUE)
  abline(v = quantile(predict_new_x[,102],c(2.5,97.5)/100), lwd = 3,
       col = "magenta2", lty = 3)
  points(mean(predict_new_x[,102]), 0, col = "magenta2", pch = 19, cex = 1.5)
  legend("topright", legend = c("kernel", "normal", "95% CI"), lwd = c(2,2, 3),
       col = c("red", "blue", "magenta2"), lty = c(2,1, 3), bty = "n")
}
```

## Bootstrap estimates

``` r
mean(predict_new_x[,93]) #for 2030
quantile(predict_new_x[,93],c(2.5, 97.5)/100)
mean(predict_new_x[,102]) #for 2040
quantile(predict_new_x[,102],c(2.5, 97.5)/100) 
```

**Section 2 - For three invasion categories separately**

# Invasive aliens

``` r
###############################
##Fitting nonlinear functions##
###############################
setwd("~/time_series_analysis")
D = read.csv(file = "data.csv", header = TRUE)
D<-subset(D,select = c('Year','Count_invasive'))
d = D
head(d)
d$Year = 1:nrow(d)
head(d)
#linear fitting
fit_linear = lm(Count_invasive ~ Year, data = d)
summary(fit_linear)
AIC(fit_linear)
#exponential fitting
fit_exponential = nls(Count_invasive ~ a*exp(b*Year), data = d,
                      start = list(a = 0.5, b = 0))
summary(fit_exponential)
AIC(fit_exponential)
#saturating fitting
fit_saturating = nls(Count_invasive  ~ a*(1-exp(-b*Year)), data = d, 
                     start = list(a = 3, b = .05),
                     control = nls.control(maxiter = 100))
AIC(fit_saturating)
summary(fit_saturating)
#hyperbolic fitting
fit_hyperbolic = nls(Count_invasive ~ a*Year/(Year + b),
                      data = d, start = list(a =4.5, b = 6.75))
summary(fit_hyperbolic)
AIC(fit_hyperbolic)
#sigmoidal fitting
fit_sigmoidal = nls(Count_invasive ~ a*Year^2/(Year^2 + b^2),
                    data = d, start = list(a =4.5, b = 6.75))
summary(fit_sigmoidal)
AIC(fit_sigmoidal)
#Compare AIC values
AIC_vals = data.frame(AIC(fit_linear), AIC(fit_exponential), AIC(fit_saturating),
             AIC(fit_hyperbolic), AIC(fit_sigmoidal))
names(AIC_vals) = c("Linear", "Exponential", "Saturating", "Hyperbolic", "Sigmoidal")
AIC_vals
#Confidence interval for the saturation point a
confint(fit_saturating)
confint(fit_hyperbolic)
confint(fit_sigmoidal)

###############################
##  Break-point analysis     ##
###############################
setwd("~/time_series_analysis")
d = read.csv(file = "data.csv", header = TRUE)
d<-subset(d,select = c('Year','Count_invasive'))
x = d$Year; y = d$Count_invasive
lin.mod = lm(y ~ x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=1900)
summary(segmented.mod)

###############################
##    Linear regression      ##
###############################
#Linear regression from 1921 (+/- SE) to 2020
d1<-subset(d,Year >= 1921)
lm_fit_1921 = lm( d1$Count_invasive ~ d1$Year)
summary(lm_fit_1921)
d2<-subset(d,Year >= 1930)
lm_fit_1930 = lm( d2$Count_invasive ~ d2$Year)
summary(lm_fit_1930)
d3<-subset(d,Year >= 1939)
lm_fit_1939 = lm( d3$Count_invasive ~ d3$Year)
summary(lm_fit_1939)
#Plot function remains same as section 1

###############################
##    Bootstrap analysis     ##
###############################
D = read.csv(file = "data.csv", header = TRUE)
d = subset(D,Year >= 1932)
new_x_vals = data.frame(Year = 1932:2040)
B = 1000
names(new_x_vals) = "Year"
prediction_intervals = matrix(data = NA,
                              nrow = 2, ncol = nrow(new_x_vals))
predict_new_x = matrix(data = NA, nrow = B, 
                       ncol = nrow(new_x_vals))
for (i in 1:B) {
  b_ind = sample(1:nrow(d), replace = TRUE)          
  b_data = d[b_ind, ]                                
  names(b_data) = names(d)                           
  b_fit = lm(Count_invasive ~ Year, data = b_data)   
  b_predict_new_x = predict(b_fit,                   
                            newdata = new_x_vals)
  predict_new_x[i,] = b_predict_new_x
}
for (j in 1:nrow(new_x_vals)) {                 # extract prediction interval
  prediction_intervals[,j] = quantile(predict_new_x[,j],c(2.5, 97.5)/100)
}
#Plot for the bootstrap confidence interval
{plot(d$Year, d$Count_invasive, pch = 21, bg="azure2",cex=1, type = "p",col="black",
     xlab = "Year", ylab = "First record rate (per year)",
     xlim = c(1932, 2040))
  lines(1932:2040, prediction_intervals[1,], col = "firebrick1", lwd= 2)
  lines(1932:2040, prediction_intervals[2,], col = "firebrick1", lwd = 2)
  legend("topright", "99% confidence interval", lwd = 2,
       col = "red", bty = "n")
}
#Plot for bootstrap distribution of the predicted first record rate in 2030
{ind = which(new_x_vals == 2030)
  hist(predict_new_x[,ind], probability = TRUE, main = "", 
     xlab = "Predicted First Record Rate (2030)",
     breaks = 30, xlim = c(-0.25, 2.25),col="seashell1")
  lines(density(predict_new_x[,ind]), lwd = 2, col = "magenta2", lty = 2)
  curve(dnorm(x, mean = mean(predict_new_x[,ind]), 
            sd = sd(predict_new_x[,ind])),
      col = "blue", lwd = 2, add = TRUE)
  abline(v = quantile(predict_new_x[,ind],c(2.5, 97.5)/100), lwd = 3,
       col = "firebrick2", lty = 3)
  points(mean(predict_new_x[,ind]), 0, col = "firebrick2", 
       pch = 19, cex = 2)
  legend("topright", legend = c("kernel", "normal", "95% CI"), 
       lwd = c(2,2, 3), col = c("red", "blue", "magenta"), 
       lty = c(2,1, 3), bty = "n")
}
#Plot for bootstrap distribution of the predicted first record rate in 2040
{ind = which(new_x_vals == 2040)
  hist(predict_new_x[,ind], probability = TRUE, main = "", 
     xlab = "Predicted First Record Rate (2040)",
     breaks = 30, ylim = c(0,1.2), xlim = c(-0.5, 2.5),col = "seashell1")
  lines(density(predict_new_x[,ind]), lwd = 2, col = "magenta2", lty = 2)
  curve(dnorm(x, mean = mean(predict_new_x[,ind]), 
            sd = sd(predict_new_x[,ind])),
      col = "blue", lwd = 2, add = TRUE)
  abline(v = quantile(predict_new_x[,ind],c(2.5, 97.5)/100), 
       lwd = 3, col = "firebrick2", lty = 3)
  points(mean(predict_new_x[,ind]), 0, col = "firebrick2", 
       pch = 19, cex = 2)
  legend("topright", legend = c("kernel", "normal", "95% CI"), 
       lwd = c(2,2, 3), col = c("red", "blue", "magenta"), 
       lty = c(2,1, 3), bty = "n")
}

#Bootstrap estimates
mean(predict_new_x[,ind]) #for 2030
quantile(predict_new_x[,ind],c(2.5, 97.5)/100)
mean(predict_new_x[,ind]) #for 2040
quantile(predict_new_x[,ind],c(2.5, 97.5)/100)
```

# Naturalized aliens

``` r
###############################
##Fitting nonlinear functions##
###############################
D = read.csv(file = "data.csv", header = TRUE)
D<-subset(D,select = c('Year','Count_naturalized'))
#D = D[-c(14, 17),] #to remove outliers, make this un-comment
d = D                    # just for execution
head(d)
d$Year = 1:nrow(d)
head(d)
#linear fitting
fit_linear = lm(Count_naturalized ~ Year, data = d)
summary(fit_linear)
AIC(fit_linear)
#exponential fitting
fit_exponential = nls(Count_naturalized ~ a*exp(b*Year), data = d,
                      start = list(a = 0.5, b = 0))
summary(fit_exponential)
AIC(fit_exponential)
#saturating fitting
fit_saturating = nls(Count_naturalized ~ a*(1-exp(-b*Year)), data = d, 
                     start = list(a = 3, b = .05),
                     control = nls.control(maxiter = 100))
summary(fit_saturating)
AIC(fit_saturating)
#hyperbolic fitting
fit_hyperbolic = nls(Count_naturalized ~ a*Year/(Year + b),
                     data = d, start = list(a =4.5, b = 6.75))
summary(fit_hyperbolic)
AIC(fit_hyperbolic)
#sigmoidal fitting
fit_sigmoidal = nls(Count_naturalized ~ a*Year^2/(Year^2 + b^2),
                      data = d, start = list(a =4.5, b = 6.75))
summary(fit_sigmoidal)
AIC(fit_sigmoidal)
#Compare AIC values
AIC_vals = data.frame(AIC(fit_linear), AIC(fit_exponential), AIC(fit_saturating),
                      AIC(fit_hyperbolic), AIC(fit_sigmoidal))
names(AIC_vals) = c("Linear", "Exponential", "Saturating", "Hyperbolic", "Sigmoidal")
AIC_vals
#Confidence interval for the saturation point a
confint(fit_saturating)
confint(fit_hyperbolic)
confint(fit_sigmoidal)

###############################
##  Break-point analysis     ##
###############################
setwd("~/time_series_analysis")
d = read.csv(file = "data.csv", header = TRUE)
d<-subset(d,select = c('Year','Count_naturalized'))
x = d$Year; y = d$Count_naturalized
lin.mod = lm(y ~ x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=1900)
summary(segmented.mod)

###############################
##    Linear regression      ##
###############################
#Linear regression from 1918 (+/- SE) to 2020
d1<-subset(d,Year >= 1918)
lm_fit_1918 = lm( d1$Count_naturalized ~ d1$Year)
summary(lm_fit_1918)
d2<-subset(d,Year >= 1928)
lm_fit_1928 = lm( d2$Count_naturalized ~ d2$Year)
summary(lm_fit_1928)
d3<-subset(d,Year >= 1938)
lm_fit_1938 = lm( d3$Count_naturalized ~ d3$Year)
summary(lm_fit_1938)
#Plot function remains same as section 1

###############################
##    Bootstrap analysis     ##
###############################
D = read.csv(file = "data.csv", header = TRUE)
d = subset(D,Year >= 1942)
new_x_vals = data.frame(Year = 1942:2040)
B = 1000                          
names(new_x_vals) = "Year"
prediction_intervals = matrix(data = NA,
                              nrow = 2, ncol = nrow(new_x_vals))
predict_new_x = matrix(data = NA, nrow = B, 
                       ncol = nrow(new_x_vals))
for (i in 1:B) {
  b_ind = sample(1:nrow(d), replace = TRUE)
  b_data = d[b_ind, ]
  names(b_data) = names(d)
  b_fit = lm(Count_naturalized ~ Year, data = b_data)
  b_predict_new_x = predict(b_fit,newdata = new_x_vals)
  predict_new_x[i,] = b_predict_new_x
}
for (j in 1:nrow(new_x_vals)) {
  prediction_intervals[,j] = quantile(predict_new_x[,j],c(2.5, 97.5)/100)
}
#Plot for the bootstrap confidence interval
{plot(d$Year, d$Count_naturalized, pch = 21, bg="azure2",cex=1, type = "p",col="black", 
     xlab = "Year", 
     ylab = "First record rate (per year)",
     xlim = c(1942, 2040))
  lines(1942:2040, prediction_intervals[1,], col = "firebrick1", lwd= 2)
  lines(1942:2040, prediction_intervals[2,], col = "firebrick1", lwd = 2)
  legend("topright", "95% confidence interval", lwd = 2,
       col = "red", bty = "n")
}
#Plot for bootstrap distribution of the predicted first record rate in 2030
{ind = which(new_x_vals == 2030)
  hist(predict_new_x[,ind], probability = TRUE, main = "", 
     xlab = "Predicted First Record Rate (2030)",
     breaks = 30, xlim = c(-0.25,1.3),col = "seashell1")
  lines(density(predict_new_x[,ind]), lwd = 2, 
      col = "magenta2", lty = 2)
  curve(dnorm(x, mean = mean(predict_new_x[,ind]), 
            sd = sd(predict_new_x[,ind])),
      col = "blue", lwd = 2, add = TRUE)
  abline(v = quantile(predict_new_x[,ind],c(2.5, 97.5)/100), 
       lwd = 3, col = "firebrick2", lty = 3)
  points(mean(predict_new_x[,ind]), 0, col = "firebrick2", 
       pch = 19, cex = 2)
  legend("topright", legend = c("kernel", "normal", "95% CI"), 
       lwd = c(2,2, 3), col = c("red", "blue", "magenta"), 
       lty = c(2,1, 3), bty = "n")
}
#Plot for bootstrap distribution of the predicted first record rate in 2040
{ind = which(new_x_vals == 2040)
  hist(predict_new_x[,ind], probability = TRUE, main = "", 
     xlab = "Predicted First Record Rate (2040)",
     breaks = 30,col = "seashell1")
  lines(density(predict_new_x[,ind]), lwd = 2, col = "magenta2", lty = 2)
  curve(dnorm(x, mean = mean(predict_new_x[,ind]), 
            sd = sd(predict_new_x[,ind])),
      col = "blue", lwd = 2, add = TRUE)
  abline(v = quantile(predict_new_x[,ind],c(2.5, 97.5)/100), lwd = 3,
       col = "firebrick2", lty = 3)
  points(mean(predict_new_x[,ind]), 0, col = "firebrick2", 
       pch = 19, cex = 2)
  legend("topright", legend = c("kernel", "normal", "95% CI"), 
       lwd = c(2,2, 3), col = c("red", "blue", "magenta"), 
       lty = c(2,1, 3), bty = "n")
  points(mean(predict_new_x[,ind]), 0, col = "magenta", 
       pch = 19, cex = 1.5)
}

#Bootstrap estimates
mean(predict_new_x[,ind]) #for 2030
quantile(predict_new_x[,ind],c(2.5, 97.5)/100)
mean(predict_new_x[,ind]) #for 2040
quantile(predict_new_x[,ind],c(2.5, 97.5)/100)
```

# Introduced aliens

``` r
###############################
##Fitting nonlinear functions##
###############################
D = read.csv(file = "data.csv", header = TRUE)
D<-subset(D,select = c('Year','Count_introduced'))
d = D
head(d)
d$Year = 1:nrow(d)
head(d)
#linear fitting
fit_linear = lm(Count_introduced ~ Year, data = d)
summary(fit_linear)
AIC(fit_linear)
#exponential fitting
fit_exponential = nls(Count_introduced ~ a*exp(b*Year), data = d,
                      start = list(a = 0.5, b = 0))
summary(fit_exponential)
AIC(fit_exponential)
#saturating fitting
fit_saturating = nls(Count_introduced ~ a*(1-exp(-b*Year)), data = d, 
                     start = list(a = 3, b = .05),
                     control = nls.control(maxiter = 100))
summary(fit_saturating)
AIC(fit_saturating)
#hyperbolic fitting
fit_hyperbolic = nls(Count_introduced ~ a*Year/(Year + b),
                     data = d, start = list(a =4.5, b = 6.75))
summary(fit_hyperbolic)
AIC(fit_hyperbolic)
#sigmoidal fitting
fit_sigmoidal = nls(Count_introduced ~ a*Year^2/(Year^2 + b^2),
                      data = d, start = list(a =4.5, b = 6.75))
summary(fit_sigmoidal)
AIC(fit_sigmoidal)
#Compare AIC values
AIC_vals = data.frame(AIC(fit_linear), AIC(fit_exponential), AIC(fit_saturating),
                      AIC(fit_hyperbolic), AIC(fit_sigmoidal))
names(AIC_vals) = c("Linear", "Exponential", "Saturating", "Hyperbolic", "Sigmoidal")
AIC_vals
#Confidence interval for the saturation point a
confint(fit_saturating)
confint(fit_hyperbolic)
confint(fit_sigmoidal)

###############################
##  Break-point analysis     ##
###############################
setwd("~/time_series_analysis")
d = read.csv(file = "data.csv", header = TRUE)
d<-subset(d,select = c('Year','Count_introduced'))
x = d$Year; y = d$Count_introduced
lin.mod = lm(y ~ x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=1900)
summary(segmented.mod)

###############################
##    Linear regression      ##
###############################
#Linear regression from 1910 (+/- SE) to 2020
d1<-subset(d,Year >= 1910)
lm_fit_1910 = lm( d1$Count_introduced ~ d1$Year)
summary(lm_fit_1910)
d2<-subset(d,Year >= 1922)
lm_fit_1922 = lm( d2$Count_introduced ~ d2$Year)
summary(lm_fit_1922)
d3<-subset(d,Year >= 1934)
lm_fit_1934 = lm( d3$Count_introduced ~ d3$Year)
summary(lm_fit_1934)
#Plot function remains same as section 1

###############################
##    Bootstrap analysis     ##
###############################
D = read.csv(file = "data.csv", header = TRUE)
d = subset(D,Year >= 1925)
new_x_vals = data.frame(Year = 1925:2040)
B = 1000    
names(new_x_vals) = "Year"
prediction_intervals = matrix(data = NA,
                              nrow = 2, ncol = nrow(new_x_vals))
predict_new_x = matrix(data = NA, nrow = B, 
                       ncol = nrow(new_x_vals))
for (i in 1:B) {
  b_ind = sample(1:nrow(d), replace = TRUE)
  b_data = d[b_ind, ]
  names(b_data) = names(d)
  b_fit = lm(Count_introduced ~ Year, data = b_data)
  b_predict_new_x = predict(b_fit,newdata = new_x_vals)
  predict_new_x[i,] = b_predict_new_x
}
for (j in 1:nrow(new_x_vals)) {
  prediction_intervals[,j] = quantile(predict_new_x[,j],c(2.5, 97.5)/100)
}
#Plot for the bootstrap confidence interval
{plot(d$Year, d$Count_introduced, pch = 21, bg="azure2",cex=1, type = "p",col="black", 
     xlab = "Year", 
     ylab = "First record rate (per year)",
     xlim = c(1925, 2040))
  lines(1925:2040, prediction_intervals[1,], col = "firebrick1", lwd= 2)
  lines(1925:2040, prediction_intervals[2,], col = "firebrick1", lwd = 2)
  legend("topright", "95% confidence interval", lwd = 2,
       col = "red", bty = "n")
}
#Plot for bootstrap distribution of the predicted first record rate in 2030
{ind = which(new_x_vals == 2030)
  hist(predict_new_x[,ind], probability = TRUE, main = "", 
     xlab = "Predicted First Record Rate (2030)",
     breaks = 30,col = "seashell1")
  lines(density(predict_new_x[,ind]), lwd = 2, 
      col = "magenta2", lty = 2)
  curve(dnorm(x, mean = mean(predict_new_x[,ind]), 
            sd = sd(predict_new_x[,ind])),
      col = "blue", lwd = 2, add = TRUE)
  abline(v = quantile(predict_new_x[,ind],c(2.5, 97.5)/100), 
       lwd = 3, col = "firebrick2", lty = 3)
  points(mean(predict_new_x[,ind]), 0, col = "firebrick2", 
       pch = 19, cex = 2)
  legend("topright", legend = c("kernel", "normal", "95% CI"), 
       lwd = c(2,2, 3), col = c("red", "blue", "magenta"), 
       lty = c(2,1, 3), bty = "n")
  points(mean(predict_new_x[,ind]), 0, col = "magenta", 
       pch = 19, cex = 1.5)
}
#Plot for bootstrap distribution of the predicted first record rate in 2040
{ind = which(new_x_vals == 2040)
  hist(predict_new_x[,ind], probability = TRUE, main = "", 
     xlab = "Predicted First Record Rate (2040)",
     breaks = 30,col = "seashell1")
  lines(density(predict_new_x[,ind]), lwd = 2, col = "magenta2", lty = 2)
  curve(dnorm(x, mean = mean(predict_new_x[,ind]), 
            sd = sd(predict_new_x[,ind])),
      col = "blue", lwd = 2, add = TRUE)
  abline(v = quantile(predict_new_x[,ind],c(2.5, 97.5)/100), lwd = 3,
       col = "firebrick2", lty = 3)
  points(mean(predict_new_x[,ind]), 0, col = "firebrick2", 
       pch = 19, cex = 2)
  legend("topright", legend = c("kernel", "normal", "95% CI"), 
       lwd = c(2,2, 3), col = c("red", "blue", "magenta"), 
       lty = c(2,1, 3), bty = "n")
  points(mean(predict_new_x[,ind]), 0, col = "magenta", 
       pch = 19, cex = 1.5)
}
#Bootstrap estimates
mean(predict_new_x[,ind]) #for 2030
quantile(predict_new_x[,ind],c(2.5, 97.5)/100)
mean(predict_new_x[,ind]) #for 2040
quantile(predict_new_x[,ind],c(2.5, 97.5)/100)
```

**END**
