#####################################
# Assignment Multivariate Statistics
# Roel Veth 
# Carlos de Cloet


## Option 2

library('mvtnorm')

rho = 0.5
mu = c(0,0)
Sigma = matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)
n = 100
R = 500
eps = 0.2

mu_null = c(0,0)

set.seed(2)

results = replicate(R, {
  y = rmvnorm(n, mean = mu, sigma = Sigma) ## bivariate normal model
  contamination = rbinom(n, 1, eps)
  y = (1-contamination)*y + contamination*cbind(rep(50,n),rep(50,n)) ## Contamination model is here
  
  Test = n*colMeans((y-mu_null))%*%solve(var(y))%*%colMeans((y-mu_null)) ## gebruiken variantie van contaminated model
  signif = qf(0.95, 2, 98)*2*99/98 # critical value, where p=2, n=100
  c(Test,signif)
})


#Hier gaan we nog nuttige dingen typen