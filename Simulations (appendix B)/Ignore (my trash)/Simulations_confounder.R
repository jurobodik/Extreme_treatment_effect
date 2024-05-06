#This R code contains two cases: Simulations from Section 'Simulations with a hidden confounder' and 'Simple Example'
#'Simple Example' is a special case when $delta = 0$ and $omega = 1.25$
#To reproduce Table from the Section 'Simulations with a hidden confounder', simply choose $n$, $delta$, $omega$ and run the code. 
final_result = c(); number_of_repetitions = 20 #In the end, we chose only 20, since the results were remarkably similar
set.seed(1)
for (i in 1:number_of_repetitions) {

  n=1000
  omega = 0
  delta = 0
  q = 0.9
  
  
  h = rnorm(n, 1, 1)
  x = rbinom(prob = 0.75, size=1,  n); 
  t = delta*h + x + rnorm(n)
  y = delta*h +omega*2/3*t*(x==1)*(t>1) + omega*2*t*(x==0)*(t>1) -1*t*(t<1) +1*(t<1) +rnorm(n, 0, 1)
  
  x=data.frame(x); data = data.frame(t, y,x) ; x1=x; x=data$x; y=data$y; t=data$t
 
  #Estimation of the threshold tau as the 90% quantile
  fit <- rqss(t~x, data = data.frame(x,t), tau =q) #here, tau is the threshold. Its a bit confusing
  
  #Computing the set S
  t_excess = c()
  t_excess_index = c()
  estimated_threshold = fitted(fit)
  for (i in 1:length(y)) {
    if (t[i]> estimated_threshold[i]) {t_excess_index = c(t_excess_index, i); t_excess = c(t_excess, t[i]- estimated_threshold[i])}
  } 
  
  
  #GPD sigma and xi estimation
  A = x[t_excess_index]
  fit2 = evgam( list(t_excess~ A, ~ 1), family="gpd", data=data.frame(A, t_excess)  )
  
  #alpha and beta estimation
  threshold =  fitted(fit)[t_excess_index]  
  Y = y[t_excess_index]
  Tr = t[t_excess_index]
  Theta = exp(fitted(fit2)[,1])
  # fit3 = gam(Y ~ s(threshold, Theta) +s(threshold, Theta, by=Tr)) #Normally we would use this. However, since we have just one confounder that is binary, we get Error 'A term has fewer unique covariate combinations than specified maximum degrees of freedom'. 
  fit3 = gam(Y ~ threshold +threshold*Tr) #Hence, we use this. Note that if you compute theta, you will end up with a vector of just two values theta=(1.22, 3.44, 1.22,1.22,1.22, 3.44, 1.22, 3.44, 1.22, 3.44 ) while threshold=(0,1,0,0,0,1,0,1,0,1). Hence, in the regression it will be the same and we need just threshold 

  #mu(t) expression
  result = c(); 
  for (i in 1:length(Y)) {
    patient_specific_result = fit3$coefficients[3] + fit3$coefficients[4]*threshold[i] 
    
    result = c(result,  patient_specific_result); 
  }
  final_result = c(final_result,  mean(result))
}

summary(final_result)
cat( 'hat(omega) = ',mean(final_result), ',     and 95% quantile = ', quantile(final_result, 0.95)-mean(final_result))






















