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

set.seed(1)

results = replicate(R, {
  y <- rmvnorm(n, mean = mu, sigma = Sigma)
})


