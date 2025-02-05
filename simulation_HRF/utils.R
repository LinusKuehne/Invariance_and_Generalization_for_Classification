# this script contains utility functions for the experiments on synthetic data


# generate uniformly dist. samples from [-max, -min] union [min, max]
# n: number of samples
# min: non-negative
# max: non-negative, satisfies max > min
runifstrong <- function(n, min, max){
  
  # sample the absolute value
  abs.samples <- runif(n, min, max)
  
  # randomly sample the sign
  sign <- sample(c(-1,1), size = n, replace = T)
  
  return(sign*abs.samples)
}




# computes weighted binary cross-entropy
# y: label in {0,1} 
# y.hat: probability in (0,1) (probability that Y=1 given some predictors)
BCE.weighted <- function(y, y.hat){
  
  # proportion of 1's
  ones <- mean(y) 
  
  # proportion of 0's
  zeros <- 1-ones 
  
  # compute weights
  w1 <- 1/(2*ones)
  w0 <- 1/(2*zeros)
  
  # cap extreme values such that the log is defined
  y.hat.norm <- y.hat
  y.hat.norm[y.hat > 0.99999999999] <- 0.99999999999
  y.hat.norm[y.hat < 0.00000000001] <- 0.00000000001
  
  nbce <- w1*y*log(y.hat.norm) + w0*(1-y)*log(1-y.hat.norm)
  
  return(-mean(nbce))
}

