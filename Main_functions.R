#HELLO AND WELOCME
#This R code for the methodology presented in our manuscirpt 'Extremal treatment effect'
#Here, you can find the following functions:
#                   'estimate_omega', where you can plug in your data and it estimates \hat{\omega} := \lim_{t\to\infty} \hat{\mu}(t+1) - \hat{\mu}(t) 
#                   'bootstrap', which computes the confidence intervals for estimate_omega. 
#                   'estimate_mu' where you can plug in your data and it estimates \hat{\mu}(value) 
#                   'bootstrap_mu', which computes the confidence intervals for estimate_mu. 
#For explanation of detailed steps, see the script 'Simple Example' where we explain step-by-step each part of the code, this code here serves just for you to use for your data without any nuisances

library(quantreg)
library(evgam)
library(mgcv)
library(Pareto)
library(evmix)
library(MASS)


#Example

#set.seed(1)
#n=2000 #sample size
#q = 0.9 #We estimate threshold tau as a q-quantile. This is the only hyperparameter that we have to choose. It represents the bias-variance trade off. You should always try a few different values to see how sensitive the result is on q: if by changing q by 0.01 the result changes a lot, something is wrong

###########################      Lets do it. First, generate data as follows:
#x = rbinom(prob = 0.75, size=1,  n); 
#t = x + rnorm(n)
#y = 1*t*(x==1)*(t>1) + 2*t*(x==0)*(t>1) -(1*t-1)*(t<1) +rnorm(n, 0, 1) #Note that omega = 1.25 since omega = 1*0.75 + 2*(1-0.75)

#x=data.frame(x); data = data.frame(t, y,x) ; x1=x

############################   plot them (this is the figure from Section 'Simple Example')
#plot_me <- function(){
#  fit <- rqss(t~x, data = data.frame(x,t), tau =q)
#  t_excess_index = c()
#  for (i in 1:(n-1)) {
#    if (t[i]> fitted(fit)[i]) {t_excess_index = c(t_excess_index, i)}
#  } 
#  plot(y~t, lwd=2)
#  points(y[x1==0]~t[x1==0], col = 'red', pch = 1, lwd = 2, cex=1)
##  points(y[t_excess_index]~t[t_excess_index], col = 'blue', pch = 22, lwd = 2, cex=1.5)
#  legend("bottomright", legend = c('x1=1', "x1=0"), col = c("black", "red"), pch = 16)
#}
#plot_me()

############################   Finally, estimate_omega should give you a number close to 1.25 +- 0.5
#estimate_omega(y, t, x, q, constant_sigma = TRUE, smooth_parameters = FALSE) #constant_sigma=TRUE since we have only one confounder and adding sigma would lead to not full rank

############################  Here is the function for bootstrap. Note that it works well only if n>=2000, since we need to bootstrap extremes and that require at least ~200 data above the threshold n*(1-q)>~200
#bootstrap(y,t,x, q=0.9, Bootstrap_size =500, constant_sigma=TRUE, smooth_parameters = FALSE, show_time=TRUE)


############################   Here is a function to estimate mu(5), should give you a number close to 5*1.25 +- 0.75
#estimate_mu(value = 5, y, t, x, q, constant_sigma = TRUE, smooth_parameters = FALSE) #constant_sigma=TRUE since we have only one confounder and adding sigma would lead to not full rank


############################   Finally, here is the function to compute the confidence intervals for mu(5). Note that it works well only if n>=2000, since we need to bootstrap extremes and that require at least ~200 data above the threshold n*(1-q)>~200
#bootstrap_for_mu(value =5, y, t, x, q,Bootstrap_size=100, constant_sigma = TRUE, smooth_parameters = FALSE)  #Increase Bootstrap_size for better results  













estimate_omega <- function(y,t,x, q='automatic', constant_sigma=FALSE, smooth_parameters = FALSE){ 
  
  n = length(y); d = ncol(x) ; if(is.null(d)){d=1}; 
  
  automatic_q <- function(n) {
    if (n <= 2000) return(0.9)
    if (n <= 10000) return(0.95)
    return(0.99)
  }
  if(q=='automatic') q=get_q(n)
  
  
  #tau estimation
  if(d==1){  
    if(smooth_parameters==FALSE){  fit <- rqss(t~x, data = data.frame(x,t), tau =q) }
    if(smooth_parameters==TRUE){   fit <- rqss(t~qss(x), data = data.frame(x,t), tau =q)}
  }
  if(d>1){   
    if(smooth_parameters==FALSE){formula <- as.formula(paste("t~", paste("x[,", 1:d, "]", collapse = " + ")))}
    if(smooth_parameters==TRUE){  formula <- as.formula(paste("t~", paste0("qss(X", 1:d, ")", collapse = " + ")))}
    fit <- rqss(formula, data = data.frame(x,t), tau =q) #fit = rqss(t~x1+x2+x3)
  }
  
  
  t_excess = c()
  t_excess_index = c()
  tau = fitted(fit) #tau is the estimated \hat{\tau}(x_i)
  for (i in 1:length(y)) {
    if (t[i]> tau[i]) {t_excess_index = c(t_excess_index, i); t_excess = c(t_excess, t[i]- tau[i])}
  } 
  
  Y_S = y[t_excess_index]
  T_S = t[t_excess_index]
  tau_S = tau[t_excess_index]  
  
  # GPD sigma and xi estimation
  if(d>1 & constant_sigma==FALSE){   
    A = t_excess
    B=c()
    for (i in 1:d) {
      B = cbind(B, x[t_excess_index, i] )
    }
    B = data.frame(B); names(B) = paste0("B", 1:d)
    
    if(smooth_parameters==FALSE){formula <- as.formula(paste("A~", paste0("B", 1:d, "", collapse = " + ")))}
    if(smooth_parameters==TRUE){formula <- as.formula(paste("A~", paste0("s(B", 1:d, ")", collapse = " + ")))}
    
    fit2 = evgam( list(formula, ~ 1), data =data.frame(A, B), family="gpd"  )
    sigma_S = exp(fitted(fit2)[,1])
  }
  
  #Step 2: Final regression
  
  if(constant_sigma==TRUE & smooth_parameters == FALSE)   fit3 = gam(Y_S ~ tau_S +tau_S*T_S)
  if(constant_sigma==TRUE & smooth_parameters == TRUE)    fit3 = gam(Y_S ~ s(tau_S, by=T_S))  
  if(constant_sigma==FALSE & smooth_parameters == FALSE)  fit3 = gam(Y_S ~ tau_S*sigma_S +tau_S*sigma_S*T_S)
  if(constant_sigma==FALSE & smooth_parameters == TRUE)   fit3 = gam(Y_S ~ s(tau_S, sigma_S) +s(tau_S, sigma_S, by=T_S))
  
  result = c()
  for (i in 1:length(Y_S)) {
    
    if(constant_sigma==TRUE & smooth_parameters == FALSE)   { patient_specific_result = fit3$coefficients[3] + fit3$coefficients[4]*tau_S[i] 
    result = c(result,  patient_specific_result)}
    if(constant_sigma==FALSE & smooth_parameters == FALSE)   { patient_specific_result = fit3$coefficients[4] + fit3$coefficients[6]*tau_S[i]  +  fit3$coefficients[7]*sigma_S[i] +  fit3$coefficients[8]*tau_S[i]*sigma_S[i]
    result = c(result,  patient_specific_result)}
    if(smooth_parameters == TRUE)   {predict_fit = predict(fit3,  type="terms" ,se.fit = FALSE, newdata = data.frame(tau_S = tau_S[i], sigma_S=sigma_S[i], T_S = 1000000000))
    result = c(result,  predict_fit[2]/1000000000)}
    
  }
  
  return(    mean(result)       )
}









#Main functions
estimate_mu <- function(value, y,t,x,q, constant_sigma=FALSE, smooth_parameters = FALSE){ 
  d = ncol(x) ; if(is.null(d)){d=1}; if(d==1){constant_sigma = TRUE}
  #tau estimation
  if(d==1){  
    if(smooth_parameters==FALSE){  fit <- rqss(t~x, data = data.frame(x,t), tau =q) }
    if(smooth_parameters==TRUE){   fit <- rqss(t~qss(x), data = data.frame(x,t), tau =q)}
  }
  if(d>1){   
    if(smooth_parameters==FALSE){formula <- as.formula(paste("t~", paste("x[,", 1:d, "]", collapse = " + ")))}
    if(smooth_parameters==TRUE){  formula <- as.formula(paste("t~", paste0("qss(X", 1:d, ")", collapse = " + ")))}
    fit <- rqss(formula, data = data.frame(x,t), tau =q) #fit = rqss(t~x1+x2+x3)
  }
  
  
  t_excess = c()
  t_excess_index = c()
  tau = fitted(fit) #tau is the estimated \hat{\tau}(x_i)
  for (i in 1:length(y)) {
    if (t[i]> tau[i]) {t_excess_index = c(t_excess_index, i); t_excess = c(t_excess, t[i]- tau[i])}
  } 
  
  Y_S = y[t_excess_index]
  T_S = t[t_excess_index]
  tau_S = tau[t_excess_index]  
  
  # GPD sigma and xi estimation
  if(d>1 & constant_sigma==FALSE){   
    A = t_excess
    B=c()
    for (i in 1:d) {
      B = cbind(B, x[t_excess_index, i] )
    }
    B = data.frame(B); names(B) = paste0("B", 1:d)
    
    if(smooth_parameters==FALSE){formula <- as.formula(paste("A~", paste0("B", 1:d, "", collapse = " + ")))}
    if(smooth_parameters==TRUE){formula <- as.formula(paste("A~", paste0("s(B", 1:d, ")", collapse = " + ")))}
    
    fit2 = evgam( list(formula, ~ 1), data =data.frame(A, B), family="gpd"  )
    sigma_S = exp(fitted(fit2)[,1])
  }
  
  #Step 2: Final regression
  
  if(constant_sigma==TRUE & smooth_parameters == FALSE)   fit3 = gam(Y_S ~ tau_S +tau_S*T_S)
  if(constant_sigma==TRUE & smooth_parameters == TRUE)    fit3 = gam(Y_S ~ s(tau_S) +s(tau_S, by=T_S))  
  if(constant_sigma==FALSE & smooth_parameters == FALSE)  fit3 = gam(Y_S ~ tau_S*sigma_S +tau_S*sigma_S*T_S)
  if(constant_sigma==FALSE & smooth_parameters == TRUE)   fit3 = gam(Y_S ~ s(tau_S, sigma_S) +s(tau_S, sigma_S, by=T_S))
  
  result = c()
  for (i in 1:length(Y_S)) {
    
    if(constant_sigma==TRUE)   { patient_specific_result = predict(fit3, se.fit = FALSE, newdata = data.frame(tau_S = tau_S[i], T_S = value))
    result = c(result,  patient_specific_result)}
    if(constant_sigma==FALSE)   { patient_specific_result =predict(fit3, se.fit = FALSE, newdata = data.frame(tau_S = tau_S[i], sigma_S=sigma_S[i], T_S = value))
    result = c(result,  patient_specific_result)}
  }
  
  return(    mean(result)       )
}









bootstrap <- function(y,t,x, q=0.9, Bootstrap_size=100, constant_sigma=FALSE, smooth_parameters = FALSE, show_time=TRUE){
  
  f <- function(data){ estimate_omega(data$y, data$t, data$x, q = q, constant_sigma, smooth_parameters)[1]  }
  data = data.frame(y,t,x)
  
  result = c()
  for (i in 1:Bootstrap_size) {
    z = sample_n(data, size = length(y), replace = TRUE) ;     z$x = z[, -c(1,2)]
    result = c(result, f(z))
    if(show_time==TRUE & i%%10 ==0)  cat('Step: ', "\r",Bootstrap_size-i)
  }
  result = as.numeric( result ); 
  
  result2 = result[abs(result)<median(abs(result)) + 5*quantile(abs(result), 0.75)]
  
  data = data.frame(y,t,x); data$x = data[, -c(1,2)]
  
  Conf_int_upper = quantile(result2, 0.95);
  Conf_int_lower = quantile(result2, 0.05);
  
  return( data.frame(Lower_conf_int=as.numeric(Conf_int_lower ),Upper_conf_int=as.numeric(Conf_int_upper)))
  
}






bootstrap_for_mu <- function(value, y,t,x, q=0.9, Bootstrap_size=100, constant_sigma=FALSE, smooth_parameters = FALSE,show_time=TRUE){
  
  f <- function(data){ estimate_mu(value, data$y, data$t, data$x, q = q, constant_sigma, smooth_parameters)[1]  }
  data = data.frame(y,t,x)
  
  result = c()
  for (i in 1:Bootstrap_size) {
    z = sample_n(data, size = length(y), replace = TRUE) ;     z$x = z[, -c(1,2)]
    result = c(result, f(z))
    if(show_time==TRUE & i%%10 ==0)  cat('Step: ', "\r",Bootstrap_size-i)
  }
  result = as.numeric( result ); 
  
  result2 = result[abs(result)<median(abs(result)) + 5*quantile(abs(result), 0.75)]
  
  data = data.frame(y,t,x); data$x = data[, -c(1,2)]
  
  Conf_int_upper = quantile(result2, 0.95);
  Conf_int_lower = quantile(result2, 0.05);
  
  return( data.frame(Lower_conf_int=as.numeric(Conf_int_lower ),Upper_conf_int=as.numeric(Conf_int_upper)))
  
}




