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

set.seed(2)

results = replicate(R, {
  y = rmvnorm(n, mean = mu, sigma = Sigma)
  contamination = rbinom(n, 1, eps)
  y = (1-contamination)*y + contamination*cbind(rep(100,n),rep(100,n)) ## Contamination model is here
  
})


#Hier gaan we nog nuttige dingen typen