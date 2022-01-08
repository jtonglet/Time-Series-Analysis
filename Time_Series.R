#############################################################
### R code for Advanced Time Series  Analysis Homework   ###
###                    Jonathan Tonglet                  ###
###            Professor : Christophe Croux              ###
############################################################


#Load  libraries
library(ggplot2)
library(forecast)
library(CADFtest)
library(readr)

#Load Dataset
app <- read_delim("app.csv", 
                  ";",
                  escape_double = FALSE, 
                  trim_ws = TRUE)
View(app)
attach(app)
dim(app)



#####################################
### Part 1 : Univariate analysis ####
#####################################


#Create a time series (TS) object representing the number of new users  per hour
#First observation : 22th of December at 9 A.M.
new_users_ts <- ts(newusers, 
                   frequency = 24,
                   start = c(22,9))

plot.ts(new_users_ts)

#Plot the correlogram 
ggAcf(new_users_ts) #The TS seems not stationary. It shows a clear seasonal pattern.

#Perform the unit root test for stationarity
max.lag <- round(sqrt(length(new_users_ts))) 
CADFtest(new_users_ts,
         type = "drift", 
         criterion = "BIC", 
         max.lag.y = max.lag) #The TS is not stationary.


#Perform the Ljung-Box test for white noise
Box.test(new_users_ts,  
         lag = 15, 
         type = "Ljung-Box")  #Strong evidence that the TS is not white noise


#Plot the month and seasonal plots
ggmonthplot(x= new_users_ts) +
  ggtitle("Monthplot") +xlab("Hour")

ggseasonplot(x = new_users_ts,  
             year.labels = 
             TRUE, continuous = TRUE) + 
  ggtitle("Seasonal plot") + xlab("Hour")


#Apply the seasonal differences operator to make the TS stationary
#The seasonal pattern repeats itself every day (= every 24 hours)
snew_users_ts <- diff(new_users_ts, 
                      lag = 24) 
plot.ts(snew_users_ts)


#Perform the unit root test for stationarity
max.lag <- round(sqrt(length(snew_users_ts))) 

CADFtest(snew_users_ts, 
         type = "drift", 
         criterion = "BIC", 
         max.lag.y = max.lag) #We reject that the TS in seasonal differences has a unit root.


#Perform the Ljung-Box test for white noise
Box.test(snew_users_ts, 
         lag = 15,
         type = "Ljung-Box") #We do not reject that the TS is white noise.




#ARIMA models

#Model specification 
ggAcf(snew_users_ts)  #AR(1) with one seasonal repetition at lag 24
ggPacf(snew_users_ts) #MA(1) with one seasonal repetition at lag 24


#We estimate three different models based on those specification

#SARIMA(1,0,1)(1,0,1)
model1 = arima(snew_users_ts,
               order = c(1,0,1), 
               seasonal = c(1,0,1))
model1
# sar1 coefficeient is not significant

#SARIMA(1,0,1)(0,0,1)
model2 = arima(snew_users_ts,
               order = c(1,0,1), 
               seasonal = c(0,0,1))
model2
#All coefficients are significant
 
#SARIMA(1,0,1)(1,0,0)
model3 = arima(snew_users_ts, 
               order = c(1,0,1),
               seasonal = c(1,0,0))
model3
#All coefficients are significant

#Model 2 and 3 are kept for validation


#Model validation
res2<- model2$residuals #residuals

ggAcf(res2) + ggtitle("SARIMA model 2 : residuals correlogram")
plot.ts(res2)
Box.test(res2, 
         lag = 15, 
         type = "Ljung-Box") #Residuals are white noise
#Model2 is valid

res3<- model3$residuals #residuals

ggAcf(res3) + ggtitle("SARIMA model 3 : residuals correlogram")
plot.ts(res3)
Box.test(res3, 
         lag = 15, 
         type = "Ljung-Box") #Residuals are white noise
#Model3 is valid


#Models comparison (BIC)
AIC(model2, k = log(145)) #806.4
AIC(model3, k = log(145)) #819.6
#The model 2, SARIMA(1,0,1)(0,0,1) is preferred





#####################################
###    Part 2 : Forecasting      ####
#####################################


#We make forecasts of the number of new users for the next 24 hours 

forecast <- function(ts, 
                     model, 
                     n.ahead = 24) {
  
  forecast <- predict(model, n.ahead = n.ahead) 
  expected <- forecast$pred 
  
  #Forecasted values are computed in differences. 
  #We sum each forecasted value with the value on the same hour the day before
  for(i in 1:length(expected)){
    expected[i] <- expected[i] + ts[i+169-n.ahead]
  }

  #95% prediction intervals
  lower <- expected - qnorm(0.975)*forecast$se
  upper <- expected + qnorm(0.975)*forecast$se

  plot.ts(ts, 
          xlim= c(22,32),
          ylim = c(0,25))
  lines(expected_, col ="red")
  lines(upper, col = "blue")
  lines(lower, col = "blue")   
  
}


#Forecasts with model 2
forecast(new_users_ts, model2, 24)

#Forecasts with model 2
forecast(new_users_ts, model3, 24)




#Computing forecast error with expanding window approach

expanding_window <- function(ts, arima_order = c(1,0,1), arima_seasonal = c(0,0,1)) {
  y <- ts
  S = round(0.75*length(y)) 
  h = 1  #One-step ahead forecast
  
  
  error.h <- c()
  
  for(i in S:(length(y)-h))
  {
    mymodel.sub <- arima(y[1:i], 
                         order = arima_order, 
                         seasonal = arima_seasonal)
    
    predict.h <- predict(mymodel.sub, 
                         n.ahead = h)$pred[h]
    
    error.h <- c(error.h, y[i+1], predict.h)
  } 
  return(error.h)
  
  }


error2.h <- expanding_window(snew_users_ts, c(1,0,1), c(0,0,1))
error3.h <- expanding_window(snew_users_ts, c(1,0,1), c(1,0,0))


#Mean absolute error
MAE2 <- mean(abs(error2.h)) #2.15
MAE3 <- mean(abs(error3.h)) #2.16
#Results are really close to each other

#Diebold-Mariano test
dm.test(error2.h,error3.h, h = 1, power = 1)
# Model 2 and model 3 performances are not significantly different






########################################
###  Part 3 : Multivariate analysis  ###
########################################

#Create a time series for the number of active users per hour 
active_users_ts <- ts(users, 
                      frequency = 24, 
                      start = c(22,9))
plot.ts(active_users_ts)

#Plot a multiple graph of new users and active users per hour
ts.plot(new_users_ts,
        active_users_ts,
        col = c("blue","orange"))


#Plot the correlogram of the TS
ggAcf(active_users_ts) #The TS shows a seasonal pattern

#Perform the unit root test for stationarity
max.lag <- round(sqrt(length(new_users_ts))) #13
CADFtest(active_users_ts, 
         type = "drift",
         criterion = "BIC",  
         max.lag.y = max.lag)   #Strong evidence that the TS is not stationary


#Perform the Ljung-Box test for white noise
Box.test(active_users_ts,  
         lag = 15, 
         type = "Ljung-Box") #Strong evidence that the TS is not white noise




#To make the TS stationary, we go into difference and seasonal differences
#The seasonal pattern repeats itself every day (= every 24 hours)
#The TS  has a trend
dsactive_users_ts <- diff(diff(active_users_ts, lag = 24)) 
plot.ts(dsactive_users_ts)

#Perform the unit root test for stationarity
max.lag <- round(sqrt(length(dsactive_users_ts)))
CADFtest(dsactive_users_ts,
         type = "drift",
         criterion = "BIC",
         max.lag.y = max.lag)


#Perform the Ljung-Box test for white noise
Box.test(dsactive_users_ts, 
         lag = 15, 
         type = "Ljung-Box")  #The TS is not white noise




#Multivariate model : Autoregressive DLM of order 2
lag <- 2
n <- length(snew_users_ts)
sY.0 <- snew_users_ts[(lag+1):n]
dsX.0 <- dsactive_users_ts[(lag+1):n]
sY.1 <- snew_users_ts[lag:(n-1)]
dsX.1 <- dsactive_users_ts[lag:(n-1)]
sY.2 <- snew_users_ts[(lag-1):(n-2)]
dsX.2 <-  dsactive_users_ts[(lag-1):(n-2)]


fit_adlm <- lm(sY.0 ~  dsX.1 + dsX.2 + sY.1 + sY.2)
summary(fit_adlm)
#X.1 and Y.1 are significant, X.2, Y.2 are not



#Model validation
res_adlm <- fit_adlm$residuals

ggAcf(res_adlm)
plot.ts(res_adlm)
Box.test(res_adlm, 
         lag = 15,
         type = "Ljung-Box") #Model is validated


#test for Granger causality
fit_adlm_small <- lm(sY.0 ~sY.1+sY.2)
anova(fit_adlm, fit_adlm_small)   #We strongly reject that there is no Granger causality


