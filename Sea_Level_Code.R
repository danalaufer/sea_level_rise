# Set working directory
setwd('~/S21/174')

# Load packages
library(qpcR)
library(MASS)
library(forecast)
library(astsa)
library(GeneCycle)
source("plot.roots.R")

# Load in data
sea_level <- read.csv('gulf_sea_level.csv', header = TRUE)
sea_level$sea_level <- apply(sea_level[-1], 1, mean, na.rm = TRUE)
sea_level <- sea_level[,c('year', 'sea_level')]
range(sea_level$year)

# Make a dataframe for evenly spaced data
even_data <- data.frame(year = seq(1993, 2020.75, by = 1/12), sea_level = NA)

# Create the evenly spaced time series data by interpolating data points
for (i in 1:length(even_data$year)){
  year <- even_data[i, 'year']
  diff <- year - sea_level$year
  closest_year <- sea_level[which.min(abs(diff)),'year']
  
  if (closest_year == year){
    even_data[i, 'sea_level'] <- sea_level[which.min(abs(diff)), 'sea_level']
    next
  }
  
  if (closest_year > year) {
    even_data[i, 'sea_level'] <- approx(c(sea_level[which.min(abs(diff))-1, 'year'], 
                                          sea_level[which.min(abs(diff)), 'year']),
                                        c(sea_level[which.min(abs(diff))-1, 'sea_level'], 
                                          sea_level[which.min(abs(diff)), 'sea_level']),
                                     xout = even_data[i, 'year'])$y
    next
  }
  
  if (closest_year < year){
    even_data[i, 'sea_level'] <- approx(c(sea_level[which.min(abs(diff)), 'year'], 
                                          sea_level[which.min(abs(diff))+1, 'year']),
                                        c(sea_level[which.min(abs(diff)), 'sea_level'], 
                                          sea_level[which.min(abs(diff))+1, 'sea_level']),
                                        xout = even_data[i, 'year'])$y
    next
  }
}

even_data_total <- ts(even_data[,'sea_level'], start = c(1993, 1), frequency = 12)

# Compare the even and uneven datasets - should look very similar
par(mfrow=c(1,2))

plot(sea_level[,2] ~ sea_level[,1], type = 'l', main = 'Raw Time Series',
     xlab = 'Year', ylab = 'Gulf Sea Level')
abline(lm(sea_level[,'sea_level'] ~ sea_level[,1]), xlim = c(1993, 2020), col = 'purple')

plot.ts(even_data_total, main = 'Even Time Series', xlab = 'Year', 
        ylab = '', ylim = c(min(sea_level[,2]), max(sea_level[,2])))
abline(lm(even_data_total ~ as.numeric(0:(length(even_data_total)-1)/12+1993)), col = 'darkgreen')

hist(sea_level[,'sea_level'],breaks = 15, 
     main = 'Histogram of Raw Time Series', xlab = 'Gulf Sea Level', col = 'violet')
hist(even_data_total, breaks = 15, 
     main = 'Histogram of Even Time Series', xlab = 'Gulf Sea Level', col = 'lightgreen')

# Split into training and testing sets
even_data_ts <- ts(even_data[1:283,'sea_level'], start = c(1993, 1), frequency = 12)
test_data <- even_data[284:334, 'sea_level']

par(mfrow = c(1,1))
plot.ts(even_data_ts)
abline(lm(even_data_ts ~ as.numeric(0:(length(even_data_ts)-1)/12+1993)), col = 'red')

# Check Boxcox transformation
# First make data positive
even_data_pos <- even_data_ts - min(even_data_ts) + 1

# Now perform boxcox
bc <- boxcox(even_data_pos ~ as.numeric(1:length(even_data_pos)))
(lambda <- bc$x[which(bc$y == max(bc$y))])
data_bc <- (1/lambda)*(even_data_pos^lambda - 1)

# Compare the original and transformed datasets 
par(mfrow=c(1,2))
plot.ts(even_data_ts, main = 'Original Time Series', xlab = 'Year', ylab = 'Gulf Sea Level')
abline(lm(even_data_ts ~ as.numeric(0:(length(even_data_ts)-1)/12+1993)), col = 'darkgreen')

plot.ts(data_bc, main = 'Boxcox Transformed Time Series', xlab = 'Year', ylab = 'Gulf Sea Level')
abline(lm(data_bc ~ as.numeric(0:(length(data_bc)-1)/12+1993)), col = 'darkorange')

hist(even_data_ts, main = 'Histogram of Original Time Series', xlab = 'Gulf Sea Level',
     col = 'lightgreen')
hist(data_bc, main = 'Histogram of Boxcox Transformed Time Series', xlab = 'Gulf Sea Level',
     col = 'orange')

# The histogram of original data almost looks better than transformed, 
# so let's just use original data

# Find variance of data as a baseline 
var(even_data_ts)                       # 3000

# Difference data once to remove linear trend
y <- diff(even_data_ts,1)
par(mfrow=c(1,1))
plot.ts(y, main = 'Even Data Once Differenced', ylab = 'Gulf Sea Level')
abline(lm(y ~ as.numeric(0:(length(y)-1)/12+1993)), col = 'red')
abline(h = mean(y), col = 'darkgreen', lty = 'dashed')
var(y)                                  # 1182

par(mfrow=c(1,2))
acf(y, 75, main = '')
pacf(y, 75, main = '')
title(main = 'ACF and PACF for Once Differenced Data', outer = TRUE,
      line = -2)

# Looks like data needs to be differenced at lag 12
par(mfrow = c(1,1))
y12 <- diff(y, 12)   
plot.ts(y12, main = 'Even Data Differenced at Lags 1 and 12', ylab = 'Gulf Sea Level')
abline(lm(y12 ~ as.numeric(0:(length(y12)-1)/12+1993)), col = 'blue')
abline(h = mean(y), col = 'red', lty = 'dashed')
var(y12)                                # 927

var(diff(y12, 1))                       # 2598
var(diff(y12, 12))                      # 2728 

# Plot ACF and PACF
par(mfrow=c(1,2))
acf(y12, 75, main = '')
pacf(y12, 75, main = '')
title(main = 'ACF and PACF for Stationary Data', outer = TRUE,
      line = -2)
# Looks like possibly p = 1 2 or 3, q = 1, P = 5, Q = 1

# Test SARIMA models
for (p in 0:3){
  for (q in 0:1){ 
    for (P in 0:2){
      for (Q in 0:1){
        print(paste('SARIMA(', p, ',1,', q, ')(', P, ',1,', Q, ')[12]:',
                    AICc(arima(data_bc, order = c(p,1,q), seasonal = list(order = c(P, 1, Q), period = 12), method = "ML"))))
      }
    }
  }
}

# The following models have the lowest AICc

model1 <- arima(even_data_ts, order=c(0,1,1), 
                seasonal = list(order = c(0,1,1), period = 12), method="ML")

model2 <- arima(even_data_ts, order=c(1,1,1), 
                seasonal = list(order = c(0,1,1), period = 12), method="ML") 

model3 <- arima(even_data_ts, order=c(0,1,1), 
                seasonal = list(order = c(1,1,1), period = 12), method="ML") 

model4 <- arima(even_data_ts, order=c(3,1,1), 
                seasonal = list(order = c(0,1,1), period = 12), method="ML") 

# Examine model 1 - both coefficients are significant 
#       Since it is pure MA model, it is stationary
par(mfrow=c(1,1))
model1
# X_t = Z_t ( 1 - 0.6135 B )( 1 - 0.8032 B^(12) )
# Plot the roots
plot.roots(NULL, polyroot(c(1, -0.6135)), main="Nonseasonal MA Roots of Model 1")

# Although the SMA1 coefficient is close to -1, the confidence intervals for it 
# do not cover -1, so we can assume that this doesn't violate invertibility
# All roots are outside unit circle, so this model is stationary and invertible

# Examine model 2 - all coefficients are significant
model2
# X_t (1 - 0.2579 B) = Z_t ( 1 - 0.7994 B ) ( 1 - 0.8149 B^(12) )

#   Check that it is invertible and stationary by plotting roots
plot.roots(NULL, polyroot(c(1, -0.2579)), main="Nonseasonal AR Roots of Model 2")
plot.roots(NULL, polyroot(c(1, -0.7994)), main="Nonseasonal MA Roots Part of Model 2")

# Although the SMA1 coefficient is close to -1, the confidence intervals for it 
# do not cover -1, so we can assume that this doesn't violate invertibility
# All roots are outside unit circle, so this model is stationary and invertible

# Examine model 3 - SAR1 coefficient is not significant
#   Setting that to zero, we obtain model 1, which has lower AICc
model3

# Examine model 4 - 3rd AR component is not statistically significant
#   Setting that coefficient to zero, we obtain a model with a lower AICc
model4
AICc(arima(even_data_ts, order=c(2,1,1), 
           seasonal = list(order = c(0,1,1), period = 12), method="ML"))

new_model4 <- arima(even_data_ts, order=c(2,1,1), 
      seasonal = list(order = c(0,1,1), period = 12), method="ML") 
new_model4

# We can see that in the new model 4, the CI for MA1 coefficeint includes -1
# So, let's not use this model - it seems like it might not be invertible

# Analysis of Residuals
res1 <- model1$residuals
res2 <- model2$residuals

# Plot residuals for Model1
par(mfrow = c(1,1))
hist(res1,breaks=20, col="pink", xlab="", prob = TRUE,
     main = 'Histogram of Model 1 Residuals')
m <- mean(res1)
std <- sqrt(var(res1))
curve( dnorm(x,m,std), add=TRUE)
plot.ts(res1, main = 'Plot of Model 1 Residuals', ylab = 'Residuals')
abline(lm(res1 ~ as.numeric(0:(length(res1)-1)/12+1993)), col = 'red')
abline(h=mean(res1), col="blue", lty = 'dashed')

# Plot residuals for Model2
hist(res2, breaks=20, col="turquoise", xlab="", prob = TRUE,
     main = 'Histogram of Model 2 Residuals')
m <- mean(res2)
std <- sqrt(var(res2))
curve( dnorm(x,m,std), add=TRUE)
plot.ts(res2, main = 'Plot of Model 2 Residuals', ylab = 'Residuals')
abline(lm(res2 ~ as.numeric(0:(length(res2)-1)/12+1993)), col = 'red')
abline(h=mean(res2), col="blue", lty = 'dashed')


# Check ACFs and PACFs - Both look good
par(mfrow=c(2,2), oma = c(0, 0, 1, 0))
acf(res1, main = '')
pacf(res1, main = '')
title(main = 'ACF and PACF for Model 1 Residuals', outer = TRUE,
      line = -2)

acf(res2, main = '')
pacf(res2, main = '')
title(main = 'ACF and PACF for Model 2 Residuals', outer = TRUE,
      line = -23)

# Use ar() to find order - Both are order 0
ar(res1)
ar(res2)

# qqnorm plot looks okay
par(mfrow = c(1,2))
qqnorm(res1, main = 'Model 1 Normal Q-Q Plot'); qqline(res1)
qqnorm(res2, main = 'Model 2 Normal Q-Q Plot'); qqline(res2)

# Shapiro test - Both pass
shapiro.test(res1)
shapiro.test(res2)

# Check Kolmogorov-Smirnov Test - Both pass, model2 maybe looks best
cpgram(res1,main="Model 1 Residuals")
cpgram(res2,main="Model 2 Residuals")
title(main = 'Kolmogorov-Smirnov Test', outer = TRUE,
      line = -3)

# Check Fisher's test - Both pass
fisher.g.test(res1)
fisher.g.test(res2)


# Box Tests - Only model2 passes the Mcleod-Li Test
Box.test(res1, lag = 17, type = c('Box-Pierce'), fitdf = 2)
Box.test(res1, lag = 17, type = c('Ljung-Box'), fitdf = 2)
Box.test(res1^2, lag = 17, type = c('Ljung-Box'), fitdf = 0)

Box.test(res2, lag = 17, type = c('Box-Pierce'), fitdf = 3)
Box.test(res2, lag = 17, type = c('Ljung-Box'), fitdf = 3)
Box.test(res2^2, lag = 17, type = c('Ljung-Box'), fitdf = 0)


# Try predicting the test values to check accuracy with model3
test_preds <- sarima.for(xdata=even_data_ts, n.ahead=50, p=1, d=1, q=1, 
                    P=0, D=1, Q=1, S=12, plot = FALSE)
test_U.tr <- test_preds$pred + 1.96*test_preds$se
test_L.tr <- test_preds$pred - 1.96*test_preds$se

# Plot predicted values on even data
par(mfrow = c(1,1))
ts.plot(even_data_total, xlim = c(2015, 2021), ylim = c(-100, 300), 
        main='Predicted Values On Testing Set', xlab = 'Year',
        ylab = 'Gulf of Mexico Sea Level')
lines(test_U.tr, col="blue", lty="dashed")
lines(test_L.tr, col="blue", lty="dashed")
points(test_preds$pred, col = 'red', cex = 0.75)
lines(test_preds$pred, col = 'red', lty = 3)

# Plot predicted values on raw data
plot(sea_level[,2] ~ sea_level[,1], type = 'l', xlim = c(2015, 2021), ylim = c(-100, 300),
     main = 'Predicted Values on Testing Set, on Raw Data', xlab = 'Year', 
     ylab = 'Gulf of Mexico Sea Level')
lines(test_U.tr, col="blue", lty="dashed")
lines(test_L.tr, col="blue", lty="dashed")
points(test_preds$pred, col = 'red', cex = 0.75)
lines(test_preds$pred, col = 'red', lty = 3)

 # Looks pretty good, now predict future values 10 years ahead
preds <- sarima.for(xdata=even_data_total, n.ahead=120, p=1, d=1, q=1, 
                    P=0, D=1, Q=1, S=12, plot = FALSE)
U.tr <- preds$pred + 1.96*preds$se
L.tr <- preds$pred - 1.96*preds$se

# Plot predicted values on even data
ts.plot(even_data_total, xlim = c(2015, 2031), ylim = c(-150, 400),
        main = 'Predicted Values, Ten Years Ahead', ylab = 'Gulf of Mexico Sea Level',
        xlab = 'Year')
lines(U.tr, col="blue", lty="dashed")
lines(L.tr, col="blue", lty="dashed")
points(preds$pred, col = 'red', cex = 0.75)
lines(preds$pred, col = 'red', lty = 3)

# Plot predicted values on raw data
plot(sea_level[,2] ~ sea_level[,1], type = 'l', xlim = c(2015, 2031), ylim = c(-150, 400),
     main = 'Predicted Values on Testing Set, on Raw Data', xlab = 'Year', 
     ylab = 'Gulf of Mexico Sea Level')
lines(U.tr, col="blue", lty="dashed")
lines(L.tr, col="blue", lty="dashed")
points(preds$pred, col = 'red', cex = 0.75)
lines(preds$pred, col = 'red', lty = 3)
 
# Check last observed value and last predicted value
tail(sea_level$sea_level, 1)
tail(preds$pred, 1)

# 204 increase in 10 years
preds$pred[length(preds$pred)] - sea_level[nrow(sea_level), 2]

# How long until it rises 1220 mm?
preds <- sarima.for(xdata=even_data_total, n.ahead=2999, p=1, d=1, q=1, 
                    P=0, D=1, Q=1, S=12, plot = FALSE)
tail(preds$pred, 1)
