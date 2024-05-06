#Simulations with varying extremal region

#we denote by 'g' the patient-specific effect curve \mu_x
g = function(t, c, slope){  if(t>=c){return(5 - abs(slope)*(t-c))}
                            if(t<c){return(5 + abs(slope)*(t-c))}}

#First, we generate the plot of the function \mu_x
{n=10000

X = runif(n, 0, 4) #This is just for plot, not X as covariate
Y=c(); for (i in 1:n) {  Y=c(Y, g(X[i], 2, 3))}
plot(Y~X, ylab = expression(mu[x](t)), xlab = 't', cex = 0.5, cex.lab = 1.1, main = 'Effect of T on Y')
text(2, +0, 'c', pos = 1, cex = 2, col = 'red')
abline(v = 2, col = 'blue', lty = 2)
text(2.8, 3, '-Slope(x)', col = 'red', pos = 4)
text(0.2, 3, 'Slope(x)', col = 'red', pos = 4)
}

#This is just a nice visualisation of data. You can ignore it if you are just replicating data
plot_me <- function(){
  fit <- rqss(t~x, data = data.frame(x,t), tau =q)
  t_excess_index = c()
  for (i in 1:(n-1)) {
    if (t[i]> fitted(fit)[i]) {t_excess_index = c(t_excess_index, i)}
  } 
  plot(y~t, lwd=2)
  points(y[t_excess_index]~t[t_excess_index], col = 'blue', pch = 22, lwd = 2, cex=1.5)
  legend("bottomright", legend = c('x1=1', "x1=0"), col = c("black", "red"), pch = 16)
}

#Next, we define the automatic procedure for the esimation of omega
hat_omega <- function(y,t,x, tau=0.95, constant_theta=FALSE){ 
  
  fit <- rqss(t~x, data = data.frame(x,t), tau =q) #here, tau and q are a bit confusing. Here, tau is just a name of parameter of the function rqss
  estimated_threshold = fitted(fit) #THIS IS tau(x). We can check that it correctly estimates the coefficient beta=1
  
  
  #Computing the set S
  t_excess = c()
  t_excess_index = c()
  for (i in 1:length(y)) {
    if (t[i]> estimated_threshold[i]) {t_excess_index = c(t_excess_index, i); t_excess = c(t_excess, t[i]- estimated_threshold[i])}
  } 
  
  
  #GPD sigma and xi estimation
  A = x[t_excess_index]
  fit2 = evgam( list(t_excess~ s(A), ~ 1), family="gpd", data=data.frame(A, t_excess)  ) 
  
  #alpha and beta estimation
  threshold =  fitted(fit)[t_excess_index]  
  Y = y[t_excess_index]
  Tr = t[t_excess_index]
  Theta = exp(fitted(fit2)[,1])

  fit3 = gam(Y ~ threshold +threshold*Tr) #Again, we use this instead of fit3 = gam(Y ~ s(threshold, Theta) +s(threshold, Theta, by=Tr)). Why? In the case with one confounder,  
  
  #mu(t) expression
  result = c(); 
  for (i in 1:length(Y)) {
    patient_specific_result = fit3$coefficients[3] + fit3$coefficients[4]*threshold[i] 
    
    result = c(result,  patient_specific_result); 
  }
  return( mean(result) )
}

set.seed(1)
c=2
nu = 2
n=10000
q=0.9

result=c()
for (j in 1:50) {
x = rnorm(n, 0, 1)
t = x + rt(n,nu,1)
y=c(); for (i in 1:n) {y=c(y,        g(t[i], c, x[i])+rnorm(1)            )}

result=c(result, hat_omega(y,t,x,q))
if (j%%10==0) {print(j)}
  
}
#plot_me()
#lm(y~t+x)


#result
cat(mean(result), quantile(result, 0.95)[1] - mean(result)) #you should get \hat{\omega} = -0.76\pm 0.23
