#HELLO AND WELOCME
#This R code explains in detail the code for the methodology presented in our manuscript 'Extremal treatment effect'
#We explain it on an example from section 'Simple Example'
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
library(causaldrf)

set.seed(0)
n=500 #sample size
q = 0.9 #We estimate threshold tau as a q-quantile. This is the only hyperparameter that we have to choose. It represents the bias-variance trade off. 

#Lets do it. First, generate data as follows:
x = rbinom(prob = 0.75, size=1,  n); 
t = x + rnorm(n)
y = 1*t*(x==1)*(t>1) + 2*t*(x==0)*(t>1) -(1*t-1)*(t<1) +rnorm(n, 0, 1)

#plot them (this is the figure from Section 'Simple Example')
plot_me <- function(){
  fit <- rqss(t~x, data = data.frame(x,t), tau =q)
  t_excess_index = c()
  for (i in 1:(n-1)) {
    if (t[i]> fitted(fit)[i]) {t_excess_index = c(t_excess_index, i)}
  } 
  plot(y~t, lwd=2)
  points(y[x==0]~t[x==0], col = 'red', pch = 1, lwd = 2, cex=1)
  points(y[t_excess_index]~t[t_excess_index], col = 'blue', pch = 22, lwd = 2, cex=1.5)
  legend("bottomright", legend = c('x1=1', "x1=0"), col = c("black", "red"), pch = 16)
}
plot_me()



#Here we compute other classical estimates from the literature
#We mostly just used library(causaldrf) and Kennedy et al code (see library(npcausal) )
#If you want to run Kennedy_estimate, you have to first run the code from file Kennedy_et_al_method located in the Simulations file


grid_val = seq(-1, 6, by = 0.5) #What values we want to estimate
x = jitter(x) #This is just because some estimates dont like binary con-founder

gam_list <- gam_est(Y = y,treat = t, data = data,grid_val = grid_val,treat_mod = "Normal",
                    treat_formula = t ~ x)

add_spl_list <- add_spl_est(Y = y,treat = t,data = data,grid_val = grid_val,knot_num = 3, treat_mod = "Normal",
                            treat_formula = t ~ x)

iptw_estimate <- iptw_est(Y = y,  treat = t,numerator_formula = t ~ 1,  data = data,degree = 2, treat_mod = "Normal", link_function = "inverse",
                          treat_formula = t ~ x)

#Kennedy_estimate = ctseff(y, a=t, x=matrix(c(x, rnorm(n)), ncol = 2),a.rng = c(min(t), 6), bw.seq = seq(0.1, 5, length.out = 100)) 



#Just plot the estimates
plot(t, y, xlab = "T", ylab = "Y", main = "DRF estimate", xlim = c(min(grid_val), max(grid_val-0.2)), ylim = c(min(y), max(y+3)), pch = 20, cex = 0.7)
for (i in seq(8.1, 10.1, by=0.1)) {  
  segments(4, 5, 7, i, col = "lightblue", lty = 1, lwd = 5)
}

segments(4, 5.03, 7, 9.1, col = "blue", lty = 1, lwd = 4) #estimate_mu(4) = 5.03  and omega = 1.26 +- 0.39
segments(-5, 4, 1, 1, col = "orange", lty = 1, lwd = 3)
segments( 1, 1,40,50 , col = "orange", lty = 1, lwd = 3)

lines(grid_val, gam_list$param, lty = 6,lwd = 3,col = "darkgreen")
lines(grid_val, add_spl_list$param,lty = 4,lwd = 3, col = "red")

wtf=iptw_estimate$param[1] + grid_val*iptw_estimate$param[2] + grid_val^2*iptw_estimate$param[3]
lines(grid_val, wtf, lty = 5,lwd = 3,col = "purple")

#lines(Kennedy_estimate$res$est~Kennedy_estimate$res$a.vals, lty = 64,lwd = 3,col = "grey", type = 'l')

legend("topleft",
       legend = c("True line", "Our estimate","CI our estimate", 'Kennedy et al', "Bia et al.", 'HI using GAM', 'IPTW'),
       col = c( "orange","blue","lightblue",'grey', "red","darkgreen", 'purple' ),
       lty = c(1, 1, 1, 64, 4, 6, 5),
       lwd = c(3,3,3,3,3,3),
       cex = 0.9)




#Now, we estimate the \mu(t) and \omega

# We compute omega by the following steps: 
#                 quantile regression
#                 estimation of theta
#                 GAM regression
#                 final result handling for 1) concrete with x* 2) population level result
#       then, we compute bootstrap confidence intervals

x=data.frame(x); data = data.frame(t, y,x) ; x1=x
#Estimation of the threshold tau as the 90% quantile
#fit <- rqss(t~qss(x), data = data.frame(x,t), tau =q) #Normally, we use smooth nonparametric regression. However, since we have just one confounder that is binary, we get Error 'A term has fewer unique covariate combinations than specified maximum degrees of freedom'. 
fit <- rqss(t~x, data = data.frame(x,t), tau =q) #here, tau and q are a bit confusing. Here, tau is just a name of parameter of the function rqss
estimated_threshold = fitted(fit) #THIS IS tau(x), 


#Computing the set S
t_excess = c()
t_excess_index = c()
for (i in 1:length(y)) {
  if (t[i]> estimated_threshold[i]) {t_excess_index = c(t_excess_index, i); t_excess = c(t_excess, t[i]- estimated_threshold[i])}
} 


#GPD sigma and xi estimation
A = x[t_excess_index,1]
#fit2 = evgam( list(t_excess~ s(A), ~ 1), family="gpd", data=data.frame(A, t_excess)  ) #Normally, we use smooth nonparametric regression. However, since we have just one confounder that is binary, we get Error 'A term has fewer unique covariate combinations than specified maximum degrees of freedom'. 
fit2 = evgam( list(t_excess~ A, ~ 1), family="gpd", data=data.frame(A, t_excess)  ) #Hence, we use only linear approximation since x is binary

#alpha and beta estimation
threshold =  fitted(fit)[t_excess_index]  
Y = y[t_excess_index]
Tr = t[t_excess_index]
sigma = exp(fitted(fit2)[,1])
# fit3 = gam(Y ~ s(threshold, sigma) +s(threshold, sigma, by=Tr)) #Normally we would use this. However, since we have just one confounder that is binary, we get Error 'A term has fewer unique covariate combinations than specified maximum degrees of freedom'. 
fit3 = gam(Y ~ threshold +threshold*Tr) #Hence, we use this. Note that if you compute sigma, you will end up with a vector of just two values sigma=(1.22, 3.44, 1.22,1.22,1.22, 3.44, 1.22, 3.44, 1.22, 3.44 ) while threshold=(0,1,0,0,0,1,0,1,0,1). Hence, in the regression it will be the same and we need just threshold 

#mu(t) expression
result = c(); 
for (i in 1:length(Y)) {
  patient_specific_result = fit3$coefficients[3] + fit3$coefficients[4]*threshold[i] 
  
  result = c(result,  patient_specific_result); 
}
mean(result) #you should get a number close to 1.25+-0.3



#Now we explain the confidence intervals computation
bootstrap_of_omega_for_binary_case <- function(data, bootstrap_size=500, print_progress=TRUE){
  
  #The following function is just everything explained previously in an automatic function
  estimate_omega_binary_case <- function(data){
    x=data$x; y=data$y; t=data$t
    
    fit <- rqss(t~x, data = data.frame(x,t), tau =q) #here, tau and q are a bit confusing. Here, tau is just a name of parameter of the function rqss
    estimated_threshold = fitted(fit) #THIS IS tau(x), 
    
    
    #Computing the set S
    t_excess = c()
    t_excess_index = c()
    for (i in 1:length(y)) {
      if (t[i]> estimated_threshold[i]) {t_excess_index = c(t_excess_index, i); t_excess = c(t_excess, t[i]- estimated_threshold[i])}
    } 
    
    
    A = x[t_excess_index]
    #GPD sigma and xi estimation
    fit2 = evgam( list(t_excess~ A, ~ 1), family="gpd", data=data.frame(A, t_excess)  )
    
    threshold =  fitted(fit)[t_excess_index]  
    Y = y[t_excess_index]
    Tr = t[t_excess_index]
    sigma = exp(fitted(fit2)[,1])
    
    fit3 = gam(Y ~ threshold +threshold*Tr) 
    
    result = c() 
    for (i in 1:length(Y)) {
      x = fit3$coefficients[3] + fit3$coefficients[4]*threshold[i] 
      
      result = c(result,  x); 
    }
    #mean(result)
    return(mean(result))} 
  
  #bootstrap
  result = c()
  for (i in 1:bootstrap_size) {
    z = sample_n(data, size = length(y), replace = TRUE) ;     z$x = z[, -c(1,2)]
    result = c(result, estimate_omega_binary_case(z))
    if(print_progress & i%%10 ==0) print(i)
  }
  result = as.numeric( result ); 
  
  result2 = result[abs(result)<median(abs(result)) + 5*quantile(abs(result), 0.75)] #Okay, this is a little bit ugly and dishonest. Sometimes, we obtain values such as follows: (1.69,  1.45,  8683196440,  1.60,  1.01,  1.38, 1.17). This is because numerical issues sometimes arise in the resampled data if the largest value is repeated multiple times 
  
  
  data = data.frame(y,t,x); data$x = data[, -c(1,2)]
  
  hat_omega = estimate_omega_binary_case(data)
  upper = as.numeric(quantile(result2, 0.9) )
  lower = as.numeric(quantile(result2, 0.1) )
  output = data.frame(hat_omega,  lower, upper)
  
  return(output)
  
}

bootstrap_of_omega_for_binary_case(data, 500)#you should get a number close to ~0.5


