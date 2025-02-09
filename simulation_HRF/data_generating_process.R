# This script contains functions to generate data from our 'standard' SCM,

#-------------------------------------------------------------------------------
# DGP for standard SCM
#-------------------------------------------------------------------------------


# returns a sample from standard SCM of one environment with different possible structural equations for Y
# n: sample size
# c: specific value of E / intervention
# env: name of the environment (string)
# mod: which model Y follows. Implemented: logistic regression ("logreg"), probit regression ("probit"), non-linear logistic regression ("nonlin"), bump model ("bump")
sim.SCM.mod <- function(n, c, env, mod = "logreg"){

  # parameters in SCM
  alpha <- -1  
  beta <- 1
  gamma <- 1.5
  
  # the intervention
  E <- rep(c,n)
  
  # we shouldn't use the value E directly, but we can know to which experiment data belongs
  Env <- as.factor(rep(env, n))
  
  X1 <- alpha*E + rnorm(n)
  
  Y <- numeric(n)
  
  # logistic regression model for Y
  if(mod == "logreg"){
    Ny <- rlogis(n)
    Y <- ifelse(Ny < beta*X1, 1, 0)
  }
  
  # probit regression model for Y
  if(mod == "probit"){
    Ny <- rnorm(n, mean = 0, sd = pi/sqrt(3))
    Y <- ifelse(Ny < beta*X1, 1, 0)
  }
  
  # non-linear logistic regression model for Y
  if(mod == "nonlin"){
    Ny <- rlogis(n)
    
    inpt <- beta*X1
    fx <- (1/20)*(0.75*inpt^3-5*inpt) + 2*sin(3*inpt) 
    Y <- ifelse(Ny < 3*fx + 1, 1, 0)
  }
  
  # bump model for Y
  if(mod == "bump"){
    Ny <- rlogis(n)
    
    inpt <- beta*X1
    
    fx <- abs(inpt + Ny)
    
    Y <- ifelse(fx < 2.5, 1, 0)
  }
  
  X2 <- Y + gamma*E + rnorm(n)
  
  X3 <- Y + rnorm(n)
  
  sample <- data.frame(X1 = X1, X2 = X2, X3 = X3, Y = Y, E = E, Env = Env)
  return(sample)
}


# returns five "training" samples (with weaker interventions) and ten "test" samples (with stronger interventions)
# n: number of observations per training environment => total 5*n data points in training sample
# n.test: number of observations per test environment => total 10*n.test data points in test sample
# int.strength.train: max magnitude of interventions in training data
# int.strength.test: max magnitude of interventions in test data
# mod: which model Y follows. Implemented: logistic regression ("logreg"), probit regression ("probit"), non-linear logistic regression ("nonlin"), bump model ("bump")
gen.sample.fixed <- function(n, n.test, int.strength.train = 1, int.strength.test = 5, mod = "logreg"){
  
  # sample intervention values for the five environments
  train.int <- runif(n = 10, min = -int.strength.train, max = int.strength.train)
  
  s.1 <- sim.SCM.mod(n = n, c = train.int[1], env = "train1", mod = mod)
  s.2 <- sim.SCM.mod(n = n, c = train.int[2], env = "train2", mod = mod)
  s.3 <- sim.SCM.mod(n = n, c = train.int[3], env = "train3", mod = mod)
  s.4 <- sim.SCM.mod(n = n, c = train.int[4], env = "train4", mod = mod)
  s.5 <- sim.SCM.mod(n = n, c = train.int[5], env = "train5", mod = mod)
  s.6 <- sim.SCM.mod(n = n, c = train.int[6], env = "train6", mod = mod)
  s.7 <- sim.SCM.mod(n = n, c = train.int[7], env = "train7", mod = mod)
  s.8 <- sim.SCM.mod(n = n, c = train.int[8], env = "train8", mod = mod)
  s.9 <- sim.SCM.mod(n = n, c = train.int[9], env = "train9", mod = mod)
  s.10 <- sim.SCM.mod(n = n, c = train.int[10], env = "train10", mod = mod)
  
  sample_train <- rbind(s.1, s.2, s.3, s.4, s.5, s.6, s.7, s.8, s.9, s.10)
  
  # interventions for testing sampled uniformly from [-int.strength.test, -int.strength.train] union [int.strength.train, int.strength.test]
  test.int <- runifstrong(n = 10, min = int.strength.train, max = int.strength.test)
  
  t.1 <- sim.SCM.mod(n = n.test, c = test.int[1], env = "test1", mod = mod)
  t.2 <- sim.SCM.mod(n = n.test, c = test.int[2], env = "test2", mod = mod)
  t.3 <- sim.SCM.mod(n = n.test, c = test.int[3], env = "test3", mod = mod)
  t.4 <- sim.SCM.mod(n = n.test, c = test.int[4], env = "test4", mod = mod)
  t.5 <- sim.SCM.mod(n = n.test, c = test.int[5], env = "test5", mod = mod)
  t.6 <- sim.SCM.mod(n = n.test, c = test.int[6], env = "test6", mod = mod)
  t.7 <- sim.SCM.mod(n = n.test, c = test.int[7], env = "test7", mod = mod)
  t.8 <- sim.SCM.mod(n = n.test, c = test.int[8], env = "test8", mod = mod)
  t.9 <- sim.SCM.mod(n = n.test, c = test.int[9], env = "test9", mod = mod)
  t.10 <- sim.SCM.mod(n = n.test, c = test.int[10], env = "test10", mod = mod)
  
  sample_test <- rbind(t.1, t.2, t.3, t.4, t.5, t.6, t.7, t.8, t.9, t.10)
  
  return(list("sample_train" = sample_train, "sample_test" = sample_test))
}

#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# simulate from a bigger SCM
#-------------------------------------------------------------------------------


sim.env.bigger <- function(n, d.1, d.2, env.name){
  
  sd <- 1
  a.1 <- 1.5
  a.2 <- 2
  
  X1 <- rnorm(n, sd = sd)
  X2 <- X1 + a.1*d.1 + rnorm(n, sd = sd)
  
  input <- X1 + 0.75*X2
  Ny <- rlogis(n)
  fx <- 3*((1/20)*(0.5*input^3-5*input) + 2*sin(input))
  Y <- ifelse(Ny < fx, 1, 0)
  
  X5 <- rnorm(n, sd = sd)
  X3 <- -Y + X5 + rnorm(n, sd = sd)
  
  X7 <- rnorm(n, sd = sd)
  X4 <- Y - 0.5*X3 + X7 + a.2*d.2 + rnorm(n, sd = sd/2)
  
  X6 <- Y - X4 + rnorm(n, sd = sd)
  
  X8 <- rnorm(n, sd = sd)
  X9 <- rnorm(n, sd = sd)
  X10 <- rnorm(n, sd = sd)
  
  Env <- as.factor(rep(env.name, n))
  
  sample <- data.frame(X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8, X9 = X9, X10 = X10, Y = Y, Env = Env)
  
  return(sample)
}



gen.sample.bigger <- function(n.train, n.test, int.strength.train = 1, int.strength.test = 1.5){
  
  # sample intervention values for the five environments
  d.1.train <- runif(n = 10, min = -int.strength.train, max = int.strength.train)
  d.2.train <- runif(n = 10, min = -int.strength.train, max = int.strength.train)
  
  s.1 <- sim.env.bigger(n = n.train, d.1 = d.1.train[1], d.2 = d.2.train[1], env.name = "train1")
  s.2 <- sim.env.bigger(n = n.train, d.1 = d.1.train[2], d.2 = d.2.train[2], env.name = "train2")
  s.3 <- sim.env.bigger(n = n.train, d.1 = d.1.train[3], d.2 = d.2.train[3], env.name = "train3")
  s.4 <- sim.env.bigger(n = n.train, d.1 = d.1.train[4], d.2 = d.2.train[4], env.name = "train4")
  s.5 <- sim.env.bigger(n = n.train, d.1 = d.1.train[5], d.2 = d.2.train[5], env.name = "train5")
  # s.6 <- sim.env.bigger(n = n.train, d.1 = d.1.train[6], d.2 = d.2.train[6], env.name = "train6")
  # s.7 <- sim.env.bigger(n = n.train, d.1 = d.1.train[7], d.2 = d.2.train[7], env.name = "train7")
  # s.8 <- sim.env.bigger(n = n.train, d.1 = d.1.train[8], d.2 = d.2.train[8], env.name = "train8")
  # s.9 <- sim.env.bigger(n = n.train, d.1 = d.1.train[9], d.2 = d.2.train[9], env.name = "train9")
  # s.10 <- sim.env.bigger(n = n.train, d.1 = d.1.train[10], d.2 = d.2.train[10], env.name = "train10")
  
  # sample_train <- rbind(s.1, s.2, s.3, s.4, s.5, s.6, s.7, s.8, s.9, s.10)
  sample_train <- rbind(s.1, s.2, s.3, s.4, s.5)
  
  # interventions for testing sampled uniformly from [-int.strength.test, -int.strength.train] union [int.strength.train, int.strength.test]
  d.1.test <- runifstrong(n = 10, min = int.strength.train, max = int.strength.test)
  d.2.test <- runifstrong(n = 10, min = int.strength.train, max = int.strength.test)
  
  t.1 <- sim.env.bigger(n = n.test, d.1 = d.1.test[1], d.2 = d.2.test[1], env.name = "test1")
  t.2 <- sim.env.bigger(n = n.test, d.1 = d.1.test[2], d.2 = d.2.test[2], env.name = "test2")
  t.3 <- sim.env.bigger(n = n.test, d.1 = d.1.test[3], d.2 = d.2.test[3], env.name = "test3")
  t.4 <- sim.env.bigger(n = n.test, d.1 = d.1.test[4], d.2 = d.2.test[4], env.name = "test4")
  t.5 <- sim.env.bigger(n = n.test, d.1 = d.1.test[5], d.2 = d.2.test[5], env.name = "test5")
  t.6 <- sim.env.bigger(n = n.test, d.1 = d.1.test[6], d.2 = d.2.test[6], env.name = "test6")
  t.7 <- sim.env.bigger(n = n.test, d.1 = d.1.test[7], d.2 = d.2.test[7], env.name = "test7")
  t.8 <- sim.env.bigger(n = n.test, d.1 = d.1.test[8], d.2 = d.2.test[8], env.name = "test8")
  t.9 <- sim.env.bigger(n = n.test, d.1 = d.1.test[9], d.2 = d.2.test[9], env.name = "test9")
  t.10 <- sim.env.bigger(n = n.test, d.1 = d.1.test[10], d.2 = d.2.test[10], env.name = "test10")
  
  sample_test <- rbind(t.1, t.2, t.3, t.4, t.5, t.6, t.7, t.8, t.9, t.10)
  
  return(list("sample_train" = sample_train, "sample_test" = sample_test))
}






