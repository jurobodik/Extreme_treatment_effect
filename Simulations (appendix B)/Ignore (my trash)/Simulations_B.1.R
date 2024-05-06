# Simulations 1
library(raster)
library(ambient)
library(dplyr)
library(rgl)
library(quantreg)
library(raster)
library(ambient)
library(rgl)
library(evgam)
library(mgcv)
library(Pareto)
library(evmix)
library(boot)
library(copula)


random_function_3d<-function(){
  f <- long_grid(seq(-10, 10, length.out = 50),seq(-10, 10, length.out = 50), seq(-10, 10, length.out = 50))
  f$noise <- 1000*(gen_perlin(f$x, f$y, f$z, frequency = 0.01, fractal = 'fbm'))^2+ 5000*(gen_perlin(f$x, f$y, f$z, frequency = 0.01, fractal = 'fbm'))^2
  f$noise = f$noise +abs(f$x*f$y*f$z)^(0.5)+f$x*f$y*f$z #This is here to normalize gen_perlin a bit. 
  return(f) }


evaluation_of_f_3d<-function(f,X1,X2,X3){
  
  function_evaluated_in_x_y_z<-function(f,x,y,z){
    index <- which.min( sqrt((f$x - x)^2 + (f$y - y)^2 + (f$z - z)^2) )
    return(f$noise[index])
  }
  
  result=c()
  for (i in 1:length(X1)) {
    result = c(result, function_evaluated_in_x_y_z(f, X1[i], X2[i],X3[i])  )
  }
  return(result)
}                      

#EXPLANATION VIA EXAMPLE FOR GENERATING random function f(x1, x2, x3): 
#x1=rnorm(n); x2 = rnorm(n); x3=rnorm(n)
#f1 = random_function_3d()
#evaluation_of_f_3d(f1, x1, x2, x3) #THIS IS f(X1, X2, X3), WHERE f IS RANDOMLY GENERATED

set.seed(0)
omega=10
n=1000
alpha = 1
number_of_repetitions=10 #This is how much we want to wait
bootstrap_size = 20 #This is for bootstrap
q = 0.9

result = c()
result2 = c()
for (k in 1:number_of_repetitions) {
  
  gumbel.cop <- evCopula(family = 'gumbel', alpha, dim = 4)
  myMvd <- mvdc(copula=gumbel.cop, margins=c("exp", "norm", "norm", "norm"),
                paramMargins=list(list(1),
                                  list(mean = 0, sd = 1), list(mean = 0, sd = 1), list(mean = 0, sd = 1))
  )
  Z <- rMvdc(n, myMvd)
  t=Z[,1]; x1 = Z[,2]; x2 = Z[,3]; x3 = Z[,4]; 
  x= data.frame(x1, x2, x3)
  f1 = random_function_3d()
  y = 0.5*omega*t*(x1>0) + 1.5*omega*t*(x1<0) + evaluation_of_f_3d(f1, x1, x2, x3) +  rnorm(n)
  
  
  estimate_omega(y,t,x, q = q)
  
  result = c(result, estimate_omega(y,t,x, q = q, constant_sigma = FALSE, smooth_parameters = FALSE))
  result2 = c(result2, bootstrap(y,t,x,Bootstrap_size=bootstrap_size, q =q, show_time = FALSE)$Upper_conf_int)
  
  if(k%%2==0)cat("\r",k)
}
result = as.numeric(result)
result2 = as.numeric(result2)

mean(result)
quantile(result, 0.95) - mean(result)
mean(result2) -mean(result) 


