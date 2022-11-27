#####################################
# Assignment Multivariate Statistics
# Option 2: Hotelling's T2 statistic under contaminated data
# Roel Veth 
# Carlos de Cloet


### Libraries
library('mvtnorm')
library('robustbase')


### Functions
HotellingsTestStat <- function(n, mean, sigmaInverse, hypothesis) {
  n*(mean-hypothesis)%*%sigmaInverse%*%(mean-hypothesis)
}




### Model Parameters
p <- 2 # Dimension of multivariate distrubution
rho <- 0.5 # Correlation between parameters
mu <- rep(0,p) # Location parameter
eps <- 0.5 # part of data which will be contaminated
significanceLevel = 0.01 # alpha = 0.05

# Two sigma options (comment one to choose the other)
# (1) rho on all off-diagonals
Sigma <- diag(1,p)+matrix(rho,nrow=p,ncol=p)-diag(rho,p)
# (2) rho^(number of steps from diagonal)
#Sigma <- rho^t(sapply(1:p, function(i, j) abs(i-j), 1:p))


### Contamination design: Give the parameters for the contamination distribution
muCont <- rep(50,p)
SigmaCont <- diag(1,p)


# #Probably remove this, not necessary for the mixedcontDirec <- 1 # 1 for all positive, -1 for all negative, 0 for random direction
# contMag <- 1 # 1 for always max, 0 for random magnitude. (random magnitude not yet implemented)
# contMagMax <- 500 # magnitude of contamination
# contFracPosi <- 0.5 # If contDirec = 0, the fraction of contaminations which are positive








### Simulation parameters
n <- 100 # Number of draws each simulation
R <- 500 # Number of simulations
mu_null <- rep(0,p) # Null hypothesis



### MCD can be used from the covMcd() command from robustbase


set.seed(1) # Set seed for replicability
results = replicate(R, {
  draws = rmvnorm(n, mean = mu, sigma = Sigma) ## bivariate normal model
  contaminatedDraws = rbern(n,eps) # Decide which draws are to be contaminated
  
  # Todo: Try to make contaminations more efficient, currently a contamination is made
  # for each draw, not just those who are randomly selected for contamination
  contamination = rmvnorm(n,mean = muCont, sigma = SigmaCont)
  
  # # Todo: Check if this can be removed, it is most likely not a good way to make the contaminations
  # if (contDirec == 0) { # contDirec == 0 is random direction
  #   contamination <- (2*rbern(n,contFracPosi)-1) * contMagMax 
  # } else { # contDirec = +1, -1, thus all contaminations are positive or negative
  #   contamination <- contMagMax
  # }
  
  data = (1-contaminatedDraws)*draws + contaminatedDraws*contamination # The contaminated data
  
  # First calculate statistics using regular estimator
  means = colMeans(data)
  sigmaInverse = solve(var(data))
  testStat = HotellingsTestStat(n,means,sigmaInverse,mu_null)
  critVal = qf(1-significanceLevel,p,n-p)*p*n-1/(n-p) #Volgens mij is de critical value niet goed, lijkt erg hoog
  c(testStat, critVal, (testStat>critVal))
  
  
  # Now calculate using MCD estimator
  
  
  
  # Test = n*colMeans((y-mu_null))%*%solve(var(y))%*%colMeans((y-mu_null)) ## gebruiken variantie van contaminated model
  # signif = qf(0.95, 2, 98)*2*99/98 # critical value, where p=2, n=100
  # c(Test,signif)
})

NumberOfRejections = sum(results[3,])

#Hier gaan we nog nuttige dingen typen