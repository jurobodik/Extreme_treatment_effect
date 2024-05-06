#This code corresponds to the application about concrete compressive strength. 
#it contains the code for computing the extremal causal effect discussed in the manuscript for the application dataset. 
#This code is structured as follows: 
#      First, we upload the data in a correct format
#      then, we draw some plots for data visualisation
#      then, we compute omega by the following steps: 
#                 quantile regression
#                 estimation of theta
#                 GAM regression
#                 final result handling for 1) concrete with x* 2) population level result
#       then, we compute bootstrap confidence intervals



library(readxl)
library(quantreg)
library(raster)
library(ambient)
library(dplyr)
library(rgl)
library(evgam)
library(mgcv)
library(Pareto)
library(evmix)
library(boot)
library(copula)

#data upload
my_data <- read_excel("Data/Concrete/Concrete_Data.xls")

y=as.numeric(my_data$`Concrete compressive strength(MPa, megapascals)`)
t=as.numeric(my_data$`Blast Furnace Slag (component 2)(kg in a m^3 mixture)`)
x1 = as.numeric( my_data$`Cement (component 1)(kg in a m^3 mixture)` )
x2=as.numeric(my_data$`Fly Ash (component 3)(kg in a m^3 mixture)`)
x3=as.numeric(my_data$`Water  (component 4)(kg in a m^3 mixture)`)
x4=as.numeric(my_data$`Superplasticizer (component 5)(kg in a m^3 mixture)`)


n=length(y)
x=data.frame(x1, x2,x3,x4); 
data = data.frame(y,t,x); 


#some first plots
plot(y~t)
cor.test(t,y)

# What if we fit a linear model?
lm_model <- lm(y~t+x1+x2+x3+x4)
summary( lm_model )
# Diagnostic plots
par(mfrow = c(2, 2))  
plot(lm_model, which = 1)
plot(lm_model, which = 2)
plot(lm_model, which = 3)
plot(lm_model, which = 4)
mtext("Diagnostics for a linear model lm(y~t+x1+x2+x3+x4) ", side = 3, line = -1.5, outer = TRUE)
par(mfrow = c(1, 1)) 





q=0.9 #change to 0.85 or 0.95 for obtaining the other values
############################       quantile regression plots for quantile q    ####################################################################################
formula <- as.formula(paste("t~", paste("x[,", 1:d, "]", collapse = " + ")))
fit <- rqss(formula, data = data.frame(x,t), tau =q)  #here, tau is the threshold. Its a bit confusing

tau = fitted(fit)
t_excess = c()
t_excess_index = c()
for (i in 1:length(y)) {
  if (t[i]> tau[i]) {t_excess_index = c(t_excess_index, i); t_excess = c(t_excess, t[i]- tau[i])}
} 


# Plots of T~X_i with blue squares indicating observations where T>tau(X)
par(mfrow = c(2, 2))
plot(t~x1, cex.lab=1.25)
points(t[t_excess_index]~x1[t_excess_index],  col = 'blue', pch = 22, lwd = 3)
plot(t~x2, cex.lab=1.25)
points(t[t_excess_index]~x2[t_excess_index],  col = 'blue', pch = 22, lwd = 3)
plot(t~x3, cex.lab=1.25)
points(t[t_excess_index]~x3[t_excess_index],  col = 'blue', pch = 22, lwd = 3)
plot(t~x4, cex.lab=1.25)
points(t[t_excess_index]~x4[t_excess_index],  col = 'blue', pch = 22, lwd = 3)
mtext("90% threshold", side = 3, line = -3, outer = TRUE)
par(mfrow = c(1, 1))
#plot(y~t)
#points(y[t_excess_index]~t[t_excess_index],  col = 'blue', pch = 22, lwd = 2)




############################       Estimation of theta    ####################################################################################

# GPD sigma estimation. Note that A = T_e and B=X_e, that is, observations restricted to the set S
A = t_excess
B=c()
for (i in 1:d) {
  B = cbind(B, x[t_excess_index, i] )
}
B = data.frame(B); names(B) = paste0("B", 1:d)

fit2 = evgam( list(A~s(B1)+s(B2)+s(B3)+s(B4), ~ 1), data =data.frame(A, B), family="gpd"  )
summary(fit2)
plot(fit2)

############################       Final regression    ####################################################################################
Y = y[t_excess_index]
Tr = t[t_excess_index]
Theta = exp(fitted(fit2)[,1])
threshold =  fitted(fit)[t_excess_index]  

# fit3 = gam(Y ~ s(threshold, Theta) +s(threshold, Theta, by=Tr)) #Use if we have many datapoints (more than 500 in the tail)
fit3 = gam(Y ~ s(threshold) +s(threshold, by=Tr))  #Use if we have intermediate number of datapoints (between 50-500 in the tail) 
# fit3 = gam(Y ~ threshold +threshold*Tr) #Use if we have little number of datapoints (less than 50 points in the tail)

############################       Final results for x*    ####################################################################################
result_for_x_star= predict(fit3,  type="terms" ,se.fit = FALSE, newdata = data.frame(threshold = threshold[46], Theta = Theta[46], Tr = 1)) #observation 46 is the one that has x* as covariates. One can choose any other values of x* here
result_for_x_star[2] ################THIS IS THE RESULT.
41*result_for_x_star[2] ################We have to multiply this by 41 = (400-359) to obtain the final result

############################       Final results for population level    ####################################################################################
result = c()
for (i in 1:length(Y)) {
  predict_fit = predict(fit3,  type="terms" ,se.fit = FALSE, newdata = data.frame(threshold = threshold[i], Theta = Theta[i], Tr = 1))
  
  result = c(result,  predict_fit[2])
}


result_populationlevel = data.frame( Slope_estimate = mean(result)  )     
result_populationlevel[1]*41 ###############THIS IS THE RESULT




# Now the confidence intervals 
#Here is the main idea of what we do: we resample the data using code: z = sample_n(data, size = length(y), replace = TRUE)
#                                     Then, we estimate the coefficients again as we did previously
#                                     And finally take 95% quantile out of all (computed resampled) coefficients

#In order to continue, run the function 'MY_METHODOLOGY_for_population_level' and 'MY_METHODOLOGY_for_patient_specific_value' which can be found below (lines 200 onwards)
sample_size = 500 #How many bootstrap repetitions do we want to do?
set.seed(1)


#Here we go
data = data.frame(y,t,x)

result = c()
for (i in 1:sample_size) {
  z = sample_n(data, size = length(y), replace = TRUE) ;     z$x = z[, -c(1,2)]
  #result = c(result, MY_METHODOLOGY_for_patient_specific_value(z))    ###########CHANGE HERE FOR POPULATION LEVEL RESULTS
  result = c(result, MY_METHODOLOGY_for_population_level(z))
  if(i%%10 ==0)print(i)
}
result = as.numeric( result ); 

result = result[abs(result)<median(abs(result)) + 5*quantile(abs(result), 0.75)] #Okay, this is a little bit ugly and dishonest. Sometimes, we obtain values such as follows: (-10.2, -11.2, -9.9, -11.8, 54984.5, -10.7, -10.3) This is because sometimes the resampled data have many times repeated largest value and that creates numerical issues

data = data.frame(y,t,x); data$x = data[, -c(1,2)]

Conf_int = quantile(result, 0.95) - MY_METHODOLOGY_for_population_level(data); Conf_int = Conf_int[,1]

Conf_int # You should get a number around 3.0








MY_METHODOLOGY_for_population_level <- function(data){ x=data$x; y=data$y; t=data$t
#quantile regression plots for quantile q
formula <- as.formula(paste("t~", paste("x[,", 1:d, "]", collapse = " + ")))
fit <- rqss(formula, data = data.frame(x,t), tau =q)  #here, tau is the threshold. Its a bit confusing

# Estimation of theta
t_excess = c()
t_excess_index = c()
for (i in 1:length(y)) {
  if (t[i]> fitted(fit)[i]) {t_excess_index = c(t_excess_index, i); t_excess = c(t_excess, t[i]- fitted(fit)[i])}
} 

# GPD sigma and xi estimation. Note that A = T_e and B=X_e, that is, observations in the set S
A = t_excess
B=c()
for (i in 1:d) {
  B = cbind(B, x[t_excess_index, i] )
}
B = data.frame(B); names(B) = paste0("B", 1:d)

fit2 = evgam( list(A~s(B1)+s(B2)+s(B3)+s(B4), ~ 1), data =data.frame(A, B), family="gpd"  )


#Final regression
Y = y[t_excess_index]
Tr = t[t_excess_index]
Theta = exp(fitted(fit2)[,1])
threshold =  fitted(fit)[t_excess_index]  

# fit3 = gam(Y ~ s(threshold, Theta) +s(threshold, Theta, by=Tr)) #Use if we have many datapoints (more than 500 in the tail)
fit3 = gam(Y ~ s(threshold) +s(threshold, by=Tr))  #Use if we have intermediate number of datapoints (between 50-500 in the tail) 
# fit3 = gam(Y ~ threshold +threshold*Tr) #Use if we have little number of datapoints (less than 50 points in the tail)

result = c(); Std_Error = c() #CI = Confidence intervals
for (i in 1:length(Y)) {
  predict_fit = predict(fit3,  type="terms" ,se.fit = TRUE, newdata = data.frame(threshold = threshold[i], Theta = Theta[i], Tr = 1))
  
  result = c(result,  predict_fit$fit[2]); Std_Error = c(Std_Error, predict_fit$se.fit[2])
}


result_populationlevel = data.frame( Slope_estimate = mean(result) , Std.Error.max=max(Std_Error), Std.Error.mean=mean(Std_Error)   )     
return( result_populationlevel[1]*141 )
}

MY_METHODOLOGY_for_patient_specific_value <- function(data){ x1=data$x1; x2=data$x2; x3=data$x3; x4=data$x4; x = data.frame(x1, x2, x3, x4); y=data$y; t=data$t
#quantile regression plots for quantile q
fit <- rqss(t~x1+x2+x3+x4, data = data.frame(x,t), tau =q)  #here, tau is the threshold. Its a bit confusing

# Estimation of theta
t_excess = c()
t_excess_index = c()
for (i in 1:length(y)) {
  if (t[i]> fitted(fit)[i]) {t_excess_index = c(t_excess_index, i); t_excess = c(t_excess, t[i]- fitted(fit)[i])}
} 

# GPD sigma and xi estimation. Note that A = T_e and B=X_e, that is, observations in the set S
A = t_excess
B=c()
for (i in 1:d) {
  B = cbind(B, x[t_excess_index, i] )
}
B = data.frame(B); names(B) = paste0("B", 1:d)

fit2 = evgam( list(A~s(B1)+s(B2)+s(B3)+s(B4), ~ 1), data =data.frame(A, B), family="gpd"  )


#Final regression
Y = y[t_excess_index]
Tr = t[t_excess_index]
Theta = exp(fitted(fit2)[,1])
threshold =  fitted(fit)[t_excess_index]  

# fit3 = gam(Y ~ s(threshold, Theta) +s(threshold, Theta, by=Tr)) #Use if we have many datapoints (more than 500 in the tail)
fit3 = gam(Y ~ s(threshold) +s(threshold, by=Tr))  #Use if we have intermediate number of datapoints (between 50-500 in the tail) 
# fit3 = gam(Y ~ threshold +threshold*Tr) #Use if we have little number of datapoints (less than 50 points in the tail)


tau = predict(fit, newdata = data.frame(x1 = 239, x2 = 0, x3 = 185, x4 = 0))
theta= predict(fit2, newdata = data.frame(B1 = 239, B2 = 0, B3 = 185, B4 = 0))[1]

result_for_x_star= predict(fit3,  type="terms" ,se.fit = TRUE, newdata = data.frame(threshold = tau, Theta = theta, Tr = 1)) 
result_for_x_star 


return( 41*result_for_x_star$fit[2] )
}

