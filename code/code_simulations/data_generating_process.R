# This script contains functions to generate data from our 'standard' SCM,
# from the semirandom SCM with fixed DAG and random interventions and response variable,
# and from the fully random SCM



library(pcalg)
library(RBGL)
library(dagitty)
library(gtools)
library(rje)
library(igraph)






#-------------------------------------------------------------------------------
# DGP for standard SCM
#-------------------------------------------------------------------------------




# returns a sample from standard SCM of one environment
# n: sample size
# c: specific value of E / intervention
# env: name of the environment (string)
sim.SCM <- function(n, c, env){

  # parameters in SCM
  alpha <- -1
  beta <- 1
  gamma <- 1.5

  # the intervention
  E <- rep(c,n)

  # we shouldn't use the value E directly, but we can know to which experiment data belongs
  Env <- as.factor(rep(env, n))

  X1 <- alpha*E + rnorm(n)

  
  # logistic regression model for Y
  Ny <- rlogis(n)
  Y <- ifelse(Ny < beta*X1, 1, 0)

  X2 <- Y + gamma*E + rnorm(n)

  X3 <- Y + rnorm(n)

  sample <- data.frame(X1 = X1, X2 = X2, X3 = X3, Y = Y, E = E, Env = Env)
  return(sample)
}




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







# returns a sample from 5 environments from the standard SCM with different possible structural equations for Y
# n: number of observations per environment => total 5*n data points in sample
# t: scales the intervention strength (0 < t < 1)
# mod: which model Y follows. Implemented: logistic regression ("logreg"), probit regression ("probit"), non-linear logistic regression ("nonlin"), bump model ("bump")
gen.sample.standard <- function(n, t=1, mod = "logreg"){
  
  sample.0 <- sim.SCM.mod(n = n, c = t*0, env = "noInt", mod = mod)
  sample.1 <- sim.SCM.mod(n = n, c = t*1, env = "weakposInt", mod = mod)
  sample.2 <- sim.SCM.mod(n = n, c = t*2, env = "strongposInt", mod = mod)
  sample.3 <- sim.SCM.mod(n = n, c = t*(-1), env = "weaknegInt", mod = mod)
  sample.4 <- sim.SCM.mod(n = n, c = t*(-2), env = "strongnegInt", mod = mod)
  
  sample <- rbind(sample.0, sample.1, sample.2, sample.3, sample.4)
  return(sample)
}









# returns five "training" samples (with weaker interventions) and ten "test" samples (with stronger interventions)
# n: number of observations per training environment => total 5*n data points in training sample
# n.test: number of observations per test environment => total 10*n.test data points in test sample
# int.strength.train: max magnitude of interventions in training data
# int.strength.test: max magnitude of interventions in test data
gen.sample.fixed <- function(n, n.test, int.strength.train = 1, int.strength.test = 5){
  
  # sample intervention values for the five environments
  train.int <- runif(n = 5, min = -int.strength.train, max = int.strength.train)
  
  s.1 <- sim.SCM(n = n, c = train.int[1], env = "train1")
  s.2 <- sim.SCM(n = n, c = train.int[2], env = "train2")
  s.3 <- sim.SCM(n = n, c = train.int[3], env = "train3")
  s.4 <- sim.SCM(n = n, c = train.int[4], env = "train4")
  s.5 <- sim.SCM(n = n, c = train.int[5], env = "train5")
  
  sample_train <- rbind(s.1, s.2, s.3, s.4, s.5)
  
  # interventions for testing sampled uniformly from [-int.strength.test, -int.strength.train] union [int.strength.train, int.strength.test]
  test.int <- runifstrong(n = 10, min = int.strength.train, max = int.strength.test)
  
  t.1 <- sim.SCM(n = n.test, c = test.int[1], env = "test1")
  t.2 <- sim.SCM(n = n.test, c = test.int[2], env = "test2")
  t.3 <- sim.SCM(n = n.test, c = test.int[3], env = "test3")
  t.4 <- sim.SCM(n = n.test, c = test.int[4], env = "test4")
  t.5 <- sim.SCM(n = n.test, c = test.int[5], env = "test5")
  t.6 <- sim.SCM(n = n.test, c = test.int[6], env = "test6")
  t.7 <- sim.SCM(n = n.test, c = test.int[7], env = "test7")
  t.8 <- sim.SCM(n = n.test, c = test.int[8], env = "test8")
  t.9 <- sim.SCM(n = n.test, c = test.int[9], env = "test9")
  t.10 <- sim.SCM(n = n.test, c = test.int[10], env = "test10")
  
  sample_test <- rbind(t.1, t.2, t.3, t.4, t.5, t.6, t.7, t.8, t.9, t.10)
  
  return(list("sample_train" = sample_train, "sample_test" = sample_test))
}

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# DGP for random DAG
#-------------------------------------------------------------------------------






# generate a training samples from five environments and ten testing environments from the same random SCM
# n: number of observations per training environment => total 5*n data points in training sample
# n.test: number of observations per test environment => total 10*n.test data points in test sample
# d: number of covariates (including Y) for the randomly sampled DAG
# max.pa: maximum number of parents a node can have in the DAG
# num.int: number of interventions in the DAG
# t: scales the intervention values (0 < t < 1) 
# int.strength.train: max magnitude of interventions in training data
# int.strength.test: max magnitude of interventions in test data
# mod: which model Y follows. Implemented: logistic regression ("logreg"), probit regression ("probit"), non-linear logistic regression ("nonlin"), bump model ("bump")
generate.samples.random <- function(n, n.test, d, max.pa, num.int, t = 1, int.strength.train = 2, int.strength.test = 10, mod = "logreg"){
  
  # decide which node is response
  Y.ind <- sample(1:d, size = 1)
  
  # build adjacency matrix
  A.obj <- build.adj.mat(d = d, max.pa = max.pa, Y.ind = Y.ind, num.int = num.int)
  A <- A.obj$adj.mat

  
  # we now check if the resulting graph is disconnected
  adj.conn <- A + t(A)
  graph.conn <- graph.adjacency(adj.conn, mode = "undirected")
  
  # resample the adjacency matrix until associated DAG is connected
  while(!is.connected(graph.conn)){
    
    # build adjacency matrix
    A.obj <- build.adj.mat(d = d, max.pa = max.pa, Y.ind = Y.ind, num.int = num.int)
    A <- A.obj$adj.mat

    adj.conn <- A + t(A)
    graph.conn <- graph.adjacency(adj.conn, mode = "undirected")
  }
  

  
  # generate list with variable names in the correct order
  var.names <- generate.varnames(d = d, Y.ind = Y.ind, num.int = num.int)
  
  
  # turn A into DAG objects:
  
  # including I
  dag.full <- pcalg2dagitty(amat = t(A), labels = var.names, type = "dag")
  
  # excluding I
  dag.cov <- pcalg2dagitty(amat = t(A[1:d,1:d]), labels = var.names[1:d], type = "dag")
  
  
  # generate data from this DAG (with randomly sampled weights)
  samples <- dag2sample(d = d, n = n, n.test = n.test, num.int = num.int, dag.full = dag.full, dag.cov = dag.cov, var.names = var.names, t = t, int.strength.train = int.strength.train, int.strength.test = int.strength.test, mod = mod)
  
  
  returnlist <- list("sample_train" = samples$sample_train,
                     "sample_test" = samples$sample_test,
                     "dag.full" = dag.full,
                     "dag.cov"= dag.cov,
                     "num.int" = num.int)
  
  return(returnlist)
}







# given a DAG, generate training samples from five training environments and ten test environments (with stronger interventions)
# d: number of covariates (including Y) for the randomly sampled DAG
# n: number of observations per training environment => total 5*n data points in training sample
# n.test: number of observations per test environment => total 10*n.test data points in test sample
# num.int: number of interventions in the DAG
# dag.full: dagitty object of the DAG involving covariates, response, and environment variables
# dag.cov: dagitty object of the DAG involving just covariates and response
# var.names: vector of variable names in the correct order
# t: scales the intervention values (0 < t < 1) 
# int.strength.train: max magnitude of interventions in training data
# int.strength.test: max magnitude of interventions in test data
# mod: which model Y follows. Implemented: logistic regression ("logreg"), probit regression ("probit"), non-linear logistic regression ("nonlin"), bump model ("bump")
dag2sample <- function(d, n, n.test, num.int, dag.full, dag.cov, var.names, t, int.strength.train = 2, int.strength.test = 10, mod){
  
  # extract topological order
  top.order <- topologicalOrdering(dag.full) 
  
  # sample weights
  weights.SCM <- generate.weights(dag.full, top.order)
  
  
  # generate the training samples from five environments
  s1 <- generate.one.env(env.name = "env1", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
  s2 <- generate.one.env(env.name = "env2", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
  s3 <- generate.one.env(env.name = "env3", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
  s4 <- generate.one.env(env.name = "env4", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
  s5 <- generate.one.env(env.name = "env5", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
  
  # generate the test samples from ten environments
  t1 <- generate.one.env(env.name = "test1", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t2 <- generate.one.env(env.name = "test2", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t3 <- generate.one.env(env.name = "test3", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t4 <- generate.one.env(env.name = "test4", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t5 <- generate.one.env(env.name = "test5", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t6 <- generate.one.env(env.name = "test6", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t7 <- generate.one.env(env.name = "test7", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t8 <- generate.one.env(env.name = "test8", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t9 <- generate.one.env(env.name = "test9", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t10 <- generate.one.env(env.name = "test10", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  
  
  # if the randomly sampled DAG, weights, and interventions result in "degenerate" data 
  # (where at least 95% of response observations have the same value 0/1), we resample the weights and try again
  resamp.weights <- F
  
  # if the data from one environment is degenerate, generate.one.env is defined to be NULL
  if(is.null(s1) || is.null(s2) || is.null(s3) || is.null(s4) || is.null(s5) || is.null(t1) || is.null(t2) || is.null(t3) || is.null(t4) || is.null(t5) || is.null(t6) || is.null(t7) || is.null(t8) || is.null(t9) || is.null(t10)){
    resamp.weights <- T
  }
  
  # this counter avoids an infinite loop (but this won't happen in practice)
  num.resamp.weights <- 0
  
  # resample until we get non-degenerate samples
  while(resamp.weights){
    
    # increase counter
    num.resamp.weights <- num.resamp.weights + 1
    
    # sample weights
    weights.SCM <- generate.weights(dag.full, top.order)
    
    
    # generate the training samples from five environments
    s1 <- generate.one.env(env.name = "env1", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
    s2 <- generate.one.env(env.name = "env2", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
    s3 <- generate.one.env(env.name = "env3", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
    s4 <- generate.one.env(env.name = "env4", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
    s5 <- generate.one.env(env.name = "env5", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
    
    # generate the test samples from ten environments
    t1 <- generate.one.env(env.name = "test1", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t2 <- generate.one.env(env.name = "test2", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t3 <- generate.one.env(env.name = "test3", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t4 <- generate.one.env(env.name = "test4", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t5 <- generate.one.env(env.name = "test5", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t6 <- generate.one.env(env.name = "test6", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t7 <- generate.one.env(env.name = "test7", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t8 <- generate.one.env(env.name = "test8", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t9 <- generate.one.env(env.name = "test9", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t10 <- generate.one.env(env.name = "test10", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    
    
    resamp.weights <- F
    
    # avoid infinite loop (the probability of not finding a non-degenerate Y after 1000 resamplings is essentially 0)
    if(num.resamp.weights < 1000){
      
      if(is.null(s1) || is.null(s2) || is.null(s3) || is.null(s4) || is.null(s5) || is.null(t1) || is.null(t2) || is.null(t3) || is.null(t4) || is.null(t5) || is.null(t6) || is.null(t7) || is.null(t8) || is.null(t9) || is.null(t10)){
        resamp.weights <- T
      }
      
    }
    
  }
  
  
  
  samp <- rbind(s1, s2, s3, s4, s5)
  testsamp <- rbind(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10)
  
  ret.list <- list("sample_train" = samp, "sample_test" = testsamp)
  return(ret.list)
}









# this function generates a data set of size n from ONE environment for a given DAG
# env.name: name given to this environment to distinguish the environments later
# d: number of covariates (including Y) for the randomly sampled DAG
# n: number of observations per training environment => total 5*n data points in training sample
# num.int: number of interventions in the DAG
# dag.full: dagitty object of the DAG involving covariates, response, and environment variables
# dag.cov: dagitty object of the DAG involving just covariates and response
# top.order: topological order for the DAG equal to topologicalOrdering(dag.full) (dagitty)
# weights.SCM: list of edge weights (ordered along the topological order, for each variable 
# the corresponding entry contains the weights for the arrows pointing from the parents to that variable)
# var.names: vector of variable names in the correct order
# t: scales the intervention values (0 < t < 1) 
# usage: whether this sample is intended for testing (stronger interventions) or training (weaker interventions)
# int.strength.train: max magnitude of interventions in training data
# int.strength.test: max magnitude of interventions in test data
# mod: which model Y follows. Implemented: logistic regression ("logreg"), probit regression ("probit"), non-linear logistic regression ("nonlin"), bump model ("bump")
generate.one.env <- function(env.name, d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage, int.strength.train = 2, int.strength.test = 10, mod){
  


  # covariate names (just "X1", "X2", ...)
  cov.names <- rep("A", d-1)
  for(i in 1:(d-1)){
    number <- as.character(i)
    name <- paste("X", number, sep="") 
    cov.names[i] <- name  
  }

  # initialize data matrix
  # d-1 covariates + Y + Env + interventions
  data.mat <- matrix(NA, nrow = n, ncol = d+1+num.int)
  
  df.env <- as.data.frame(data.mat)
  
  names(df.env) <- append(cov.names, append(c("Y", "Env"), var.names[(d+1):length(var.names)]))
  
  
  # indicate from which environment this sample stems
  df.env$Env <- as.factor(rep(env.name, n))
  

  
  
  
  
  
  
  
  
  # sample the intervention values (differentiating between training (weak) and test (stronger))
  for(i in (d+2):ncol(df.env)){
    
    int.val <- 0
    
    if(usage == "train"){
      int.val <- t*runif(n=1, min = -int.strength.train, max = int.strength.train)
    }
    
    if(usage == "test"){
      int.val <- t*runifstrong(1, 1, int.strength.test)
    }
    
    df.env[,i] <- rep(int.val, n)
  }
  
  
  

  
  
  # sample from the SCM along the topological order
  for(i in 1:length(top.order)){
    
    # extract name of current variable
    name <- names(which(top.order==i))
    
    # see whether it is a covariate ("X"), the response ("Y"), or an intervention variable ("I")
    name.type <- substr(name, start = 1, stop=1)
    
    
    # do nothing if name.type == "I"
    
    
    
    if(name.type == "Y"){
      pa <- parents(x=dag.cov, v="Y")
      Ny <- rlogis(n)
      if(length(pa)==0){
        df.env$Y <- ifelse(Ny < 0, 1, 0)
      } else{
        
        # linear combination for the values of all parents with the corresponding weights
        comb <- lincomb(dat = df.env[, pa, drop = F], coef = weights.SCM[[i]])
        
        if(mod == "logreg"){
          df.env$Y <- ifelse(Ny < comb, 1, 0)
        }
        
        if(mod == "probit"){
          Ny <- rnorm(n, mean = 0, sd = pi/sqrt(3))
          df.env$Y <- ifelse(Ny < comb, 1, 0)
        }
        
        if(mod == "nonlin"){
          fx <- (1/20)*(0.75*comb^3 - 5*comb) + 2*sin(3*comb) 
          df.env$Y <- ifelse(Ny < 2*fx + 1, 1, 0)
        }
        
        
        if(mod == "bump"){
          fx <- abs(comb + Ny)
          df.env$Y <- ifelse(fx < 2.5, 1, 0)
        }
        
        
      }
    }
    
    
    if(name.type == "X"){
      pa <- parents(x=dag.full, v=name)
      error <- rnorm(n, mean = 0, sd = 0.5)
      if(length(pa)==0){
        df.env[,name] <- error
      } else{
        # linear combination for the values of all parents with the corresponding weights
        comb <- lincomb(dat = df.env[, pa, drop = F], coef = weights.SCM[[i]])
        df.env[,name] <- comb + error
      }
    }
    
  }
  
  
 
  # check whether the data is "degenerate" (> 95% of Y observations have the same value)
  sum.zero <- (sum(df.env$Y) <= n*0.05)
  sum.one <- (sum(df.env$Y) >= n*0.95)
  
  # if degenerate, we sample new interventions and try again. If still degenerate, 
  # we return NULL, which will then cause new weights to be sampled.
  resamp <- sum.zero || sum.one
  
  # intervention values are resampled at most max.resamp times
  max.resamp <- 5
  
  # count number of resampling
  num.resamp <- 0
  
  while(resamp){
    
    num.resamp <- num.resamp + 1
    
    # sample the intervention values (differentiating between training (weak) and test (stronger))
    for(i in (d+2):ncol(df.env)){
      
      int.val <- 0
      
      if(usage == "train"){
        int.val <- t*runif(n=1, min = -int.strength.train, max = int.strength.train)
      }
      
      if(usage == "test"){
        int.val <- t*runifstrong(1, 1, int.strength.test)
      }
      
      df.env[,i] <- rep(int.val, n)
    }
    
    
    
    
    
    
    # sample from the SCM along the topological order
    for(i in 1:length(top.order)){
      
      # extract name of current variable
      name <- names(which(top.order==i))
      
      # see whether it is a covariate ("X"), the response ("Y"), or an intervention variable ("I")
      name.type <- substr(name, start = 1, stop=1)
      
      
      # do nothing if(name.type == "I")
      
      
      
      if(name.type == "Y"){
        pa <- parents(x=dag.cov, v="Y")
        Ny <- rlogis(n)
        if(length(pa)==0){
          df.env$Y <- ifelse(Ny < 0, 1, 0)
        } else{
          
          # linear combination for the values of all parents with the corresponding weights
          comb <- lincomb(dat = df.env[, pa, drop = F], coef = weights.SCM[[i]])
          
          if(mod == "logreg"){
            df.env$Y <- ifelse(Ny < comb, 1, 0)
          }
          
          if(mod == "probit"){
            Ny <- rnorm(n, mean = 0, sd = pi/sqrt(3))
            df.env$Y <- ifelse(Ny < comb, 1, 0)
          }
          
          if(mod == "nonlin"){
            fx <- (1/20)*(0.75*comb^3 - 5*comb) + 2*sin(3*comb) 
            df.env$Y <- ifelse(Ny < 3*fx + 1, 1, 0)
          }
          
          if(mod == "bump"){
            fx <- abs(comb + Ny)
            df.env$Y <- ifelse(fx < 2.5, 1, 0)
          }
        }
      }
      
      
      if(name.type == "X"){
        pa <- parents(x=dag.full, v=name)
        error <- rnorm(n, mean = 0, sd = 0.5)
        if(length(pa)==0){
          df.env[,name] <- error
        } else{
          # linear combination for the values of all parents with the corresponding weights
          comb <- lincomb(dat = df.env[, pa, drop = F], coef = weights.SCM[[i]])
          df.env[,name] <- comb + error
        }
      }
      
    }
    
    
    # we only want to resample the interventions at mot max.resamp times (then, we resample the weights)
    resamp.more <- num.resamp < max.resamp
    
    
    sum.zero <- (sum(df.env$Y) <= n*0.05)
    sum.one <- (sum(df.env$Y) >= n*0.95)
    
    resamp <- sum.zero || sum.one
    
    
    
    resamp <- resamp && resamp.more
    
    
  }
  
  
  # here we give the order to resample all the weights if it is difficult to get a non-degenerate dist. for Y
  if(num.resamp >= max.resamp){
    return(NULL)
  }
  
  
  return(df.env)
}











# generate the weights for each variable used in the linear combinations for the SCM (no weights needed for interventions)
# for every node, we sample one weight per parent for the SCM
# dag.full: dagitty object of the full DAG including covariates, response Y and interventions I1, I2, ...
# top.order: topological order of the graph in dag.full
generate.weights <- function(dag.full, top.order){
  
  weights.SCM <- list()
  
  for(i in 1:length(top.order)){
    
    # extract name of current variable
    name <- names(which(top.order==i))
    
    # find its parents
    pa <- parents(x=dag.full, v=name)
    
    # no parents => no weights needed in this entry => insert NA
    if(length(pa)==0){
      weights.SCM[[i]] <- NA
    } else{
      weights.SCM[[i]] <- runifstrong(n=length(pa), min = 0.5, max = 1.5)
    }
    
  }
  return(weights.SCM)
}













# build the adjacency matrix A such that A(i,j) = 1 <=> node i -> j <=> j has parent i.
# Note that indices 1:d correspond to predictors, indices (d+1):(d+num.int) to intervention variables I
# d: number of covariates (including Y)
# max.pa: maximal number of parents of a variable
# Y.ind: number in 1:d which denotes the response Y
# num.int: number of interventions
build.adj.mat <- function(d, max.pa, Y.ind, num.int){
  
  # generate the parents of each variable
  parents <- generate.parents(d, max.pa)
  
  # can't intervene on response
  int <- sample((1:d)[-Y.ind], size = num.int)
  
  # on which nodes we intervene
  intervened.nodes <- int[order(int)]
  
  # initialize adjacency matrix
  A <- matrix(0, nrow = (d+num.int), ncol = (d+num.int))
  
  # entries for covariates including Y
  for(node in 1:d){
    if(sum(parents[[node]] != 0)){
      for(pa in parents[[node]]){
        A[pa,node] <- 1
      }
    }
  }
  
  
  # entries for interventions
  for(i in 1:num.int){
    A[(d+i), intervened.nodes[i]] <- 1
  }
  
  returnlist <- list("adj.mat" = A, "num.int" = num.int)
  return(returnlist)
}




# generate a list where each entry i is a vector of parents of variable i
# d: number of covariates including Y 
# max.pa: maximum number of parents a variable can have
generate.parents <- function(d, max.pa){
  
  parents <- list()
  
  for(node in 1:d){
    
    # the first node has no parents
    if(node == 1){
      parents[[1]] <- 0
    } 
    
    # the first few nodes in the top. ordering can't have max.pa parents
    else if(node-1 < max.pa){
      num.pa <- sample(x=0:(node-1), size = 1)
      
      if(num.pa==0){
        parents[[node]] <- 0
      } else{
        pa <- sample(x=1:(node-1), size = num.pa)
        parents[[node]] <- pa[order(pa)]
      }
    } 
    
    else{
      num.pa <- sample(x=0:max.pa, size = 1) 
      
      if(num.pa==0){
        parents[[node]] <- 0
      } else{
        pa <- sample(x=1:(node-1), size = num.pa)
        parents[[node]] <- pa[order(pa)]
      }
    }
    
  }
  
  return(parents)
}

#-------------------------------------------------------------------------------










#-------------------------------------------------------------------------------
# DGP for random SCM with fixed DAG ("semirandom SCM")
#-------------------------------------------------------------------------------



# generate a training samples from five environments and ten testing environments from the same random SCM with fixed DAG
# n: number of observations per training environment => total 5*n data points in training sample
# n.test: number of observations per test environment => total 10*n.test data points in test sample
# t: scales the intervention values (0 < t < 1) 
# int.strength.train: max magnitude of interventions in training data
# int.strength.test: max magnitude of interventions in test data
# mod: which model Y follows. Implemented: logistic regression ("logreg"), probit regression ("probit"), non-linear logistic regression ("nonlin"), bump model ("bump")
generate.samples.semirandom <- function(n, n.test, t = 1, int.strength.train = 2, int.strength.test = 8, mod = "logreg"){
  
  # decide which node is response
  Y.ind <- sample(1:6, size = 1)
  
  # sample two variables to intervene on
  int.vars <- sample((1:6)[-Y.ind], size = 2)
  
  # build adjacency matrix 
  A <- matrix(0, nrow = 8, ncol = 8)
  A[1,2] <- A[1,3] <- A[2,3] <- A[3,4] <- A[3,6] <- A[4,6] <- A[5,6] <- 1
  A[7, int.vars[1]] <- A[8, int.vars[2]] <- 1
  
  
  # list with variable names in correct order
  var.names <- generate.varnames(d = 6, Y.ind = Y.ind, num.int = 2)
  
  
  # turn A into DAG objects:
  
  # including I
  dag.full <- pcalg2dagitty(amat = t(A), labels = var.names, type = "dag")
  
  # excluding I
  dag.cov <- pcalg2dagitty(amat = t(A[1:6,1:6]), labels = var.names[1:6], type = "dag")
  
  
  # generate data from this DAG
  samples <- dag2sample.semirandom(d = 6, n = n, n.test = n.test, num.int = 2, dag.full = dag.full, dag.cov = dag.cov, var.names = var.names, t = t, int.strength.train = int.strength.train, int.strength.test = int.strength.test, mod = mod)
  
  
  returnlist <- list("sample_train" = samples$sample_train,
                     "sample_test" = samples$sample_test,
                     "dag.full" = dag.full,
                     "dag.cov"= dag.cov,
                     "num.int" = 2,
                     "response" = Y.ind)
  
  return(returnlist)
}






# generate the weights for each variable used in the linear combinations for the semirandom SCM (no weights needed for interventions)
# for every node, we sample one weight per parent for the SCM
# dag.full: dagitty object of the full DAG including covariates, response Y and interventions I1, I2, ...
# top.order: topological order of the graph in dag.full
generate.weights.semirandom <- function(dag.full, top.order){
  
  weights.SCM <- list()
  
  for(i in 1:length(top.order)){
    
    # extract name of current variable
    name <- names(which(top.order==i))
    
    # find its parents
    pa <- parents(x=dag.full, v=name)
    
    if(length(pa)==0){
      weights.SCM[[i]] <- NA
    } else{
      weights.SCM[[i]] <- sample(c(-1,1), size = length(pa), replace = T)
    }
    
  }
  return(weights.SCM)
}






# given a fixed DAG with interventions, generate training samples from five training environments and ten test environments (with stronger interventions)
# d: number of covariates (including Y) for the randomly sampled DAG
# n: number of observations per training environment => total 5*n data points in training sample
# n.test: number of observations per test environment => total 10*n.test data points in test sample
# num.int: number of interventions in the DAG
# dag.full: dagitty object of the DAG involving covariates, response, and environment variables
# dag.cov: dagitty object of the DAG involving just covariates and response
# var.names: vector of variable names in the correct order
# t: scales the intervention values (0 < t < 1) 
# int.strength.train: max magnitude of interventions in training data
# int.strength.test: max magnitude of interventions in test data
# mod: which model Y follows. Implemented: logistic regression ("logreg"), probit regression ("probit"), non-linear logistic regression ("nonlin"), bump model ("bump")
dag2sample.semirandom <- function(d, n, n.test, num.int, dag.full, dag.cov, var.names, t = 1, int.strength.train = 2, int.strength.test = 8, mod){
  
  # topological order
  top.order <- topologicalOrdering(dag.full) 
  
  
  # sample weights
  weights.SCM <- generate.weights.semirandom(dag.full, top.order)
  
  
  # generate the training samples from five environments
  s1 <- generate.one.env.semirandom(env.name = "env1", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
  s2 <- generate.one.env.semirandom(env.name = "env2", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
  s3 <- generate.one.env.semirandom(env.name = "env3", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
  s4 <- generate.one.env.semirandom(env.name = "env4", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
  s5 <- generate.one.env.semirandom(env.name = "env5", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
  
  # generate the test samples from ten environments
  t1 <- generate.one.env.semirandom(env.name = "test1", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t2 <- generate.one.env.semirandom(env.name = "test2", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t3 <- generate.one.env.semirandom(env.name = "test3", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t4 <- generate.one.env.semirandom(env.name = "test4", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t5 <- generate.one.env.semirandom(env.name = "test5", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t6 <- generate.one.env.semirandom(env.name = "test6", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t7 <- generate.one.env.semirandom(env.name = "test7", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t8 <- generate.one.env.semirandom(env.name = "test8", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t9 <- generate.one.env.semirandom(env.name = "test9", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  t10 <- generate.one.env.semirandom(env.name = "test10", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
  
  # if the data is degenerate (more than 95% or response observations have same value), we resample the weights and try again
  resamp.weights <- F
  
  # if a sample is degenerate, this is encoded by NULL
  if(is.null(s1) || is.null(s2) || is.null(s3) || is.null(s4) || is.null(s5) || is.null(t1) || is.null(t2) || is.null(t3) || is.null(t4) || is.null(t5) || is.null(t6) || is.null(t7) || is.null(t8) || is.null(t9) || is.null(t10)){
    resamp.weights <- T
  }
  
  # counter to avoid infinite loop (the probability of not finding a non-degenerate Y after 1000 resamplings is essentially 0)
  num.resamp.weights <- 0
  
  while(resamp.weights){
    
    # increase counter
    num.resamp.weights <- num.resamp.weights + 1
    
    # sample weights
    weights.SCM <- generate.weights.semirandom(dag.full, top.order)
    
    
    
    s1 <- generate.one.env.semirandom(env.name = "env1", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
    s2 <- generate.one.env.semirandom(env.name = "env2", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
    s3 <- generate.one.env.semirandom(env.name = "env3", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
    s4 <- generate.one.env.semirandom(env.name = "env4", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
    s5 <- generate.one.env.semirandom(env.name = "env5", d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "train", int.strength.train, int.strength.test, mod = mod)
    
    
    t1 <- generate.one.env.semirandom(env.name = "test1", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t2 <- generate.one.env.semirandom(env.name = "test2", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t3 <- generate.one.env.semirandom(env.name = "test3", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t4 <- generate.one.env.semirandom(env.name = "test4", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t5 <- generate.one.env.semirandom(env.name = "test5", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t6 <- generate.one.env.semirandom(env.name = "test6", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t7 <- generate.one.env.semirandom(env.name = "test7", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t8 <- generate.one.env.semirandom(env.name = "test8", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t9 <- generate.one.env.semirandom(env.name = "test9", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    t10 <- generate.one.env.semirandom(env.name = "test10", d, n = n.test, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage = "test", int.strength.train, int.strength.test, mod = mod)
    
    
    resamp.weights <- F
    
    # avoid infinite loop (the probability of not finding a non-degenerate Y after 1000 resamplings is essentially 0)
    if(num.resamp.weights < 1000){
      
      if(is.null(s1) || is.null(s2) || is.null(s3) || is.null(s4) || is.null(s5) || is.null(t1) || is.null(t2) || is.null(t3) || is.null(t4) || is.null(t5) || is.null(t6) || is.null(t7) || is.null(t8) || is.null(t9) || is.null(t10)){
        resamp.weights <- T
      }
      
    }
    
  }
  
  
  samp <- rbind(s1, s2, s3, s4, s5)
  testsamp <- rbind(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10)
  
  ret.list <- list("sample_train" = samp, "sample_test" = testsamp)
  return(ret.list)
}





# this function generates a data set of size n from ONE environment for a fixed DAG with some interventions
# env.name: name given to this environment to distinguish the environments later
# d: number of covariates (including Y) for the randomly sampled DAG
# n: number of observations per training environment => total 5*n data points in training sample
# num.int: number of interventions in the DAG
# dag.full: dagitty object of the DAG involving covariates, response, and environment variables
# dag.cov: dagitty object of the DAG involving just covariates and response
# top.order: topological order for the DAG equal to topologicalOrdering(dag.full) (dagitty)
# weights.SCM: list of edge weights (ordered along the topological order, for each variable 
# the corresponding entry contains the weights for the arrows pointing from the parents to that variable)
# var.names: vector of variable names in the correct order
# t: scales the intervention values (0 < t < 1) 
# usage: whether this sample is intended for testing (stronger interventions) or training (weaker interventions)
# int.strength.train: max magnitude of interventions in training data
# int.strength.test: max magnitude of interventions in test data
# mod: which model Y follows. Implemented: logistic regression ("logreg"), probit regression ("probit"), non-linear logistic regression ("nonlin"), bump model ("bump")
generate.one.env.semirandom <- function(env.name, d, n, num.int, dag.full, dag.cov, top.order, weights.SCM, var.names, t, usage, int.strength.train = 8, int.strength.test = 2, mod){

  
  # covariate names (just "X1", "X2", ...)
  cov.names <- rep("A", d-1)
  for(i in 1:(d-1)){
    number <- as.character(i)
    name <- paste("X", number, sep="") 
    cov.names[i] <- name  
  }

  
  
  # d-1 covariates + Y + Env + interventions
  data.mat <- matrix(NA, nrow = n, ncol = d+1+num.int)
  
  df.env <- as.data.frame(data.mat)
  
  names(df.env) <- append(cov.names, append(c("Y", "Env"), var.names[(d+1):length(var.names)]))
  
  df.env$Env <- as.factor(rep(env.name, n))
  
  
  
  
  # sample intervention value and distinguish between train/test
  for(i in (d+2):ncol(df.env)){
    
    int.val <- 0
    
    if(usage == "train"){
      int.val <- t*runif(n=1, min = -int.strength.train, max = int.strength.train)
    }
    
    if(usage == "test"){
      int.val <- t*runifstrong(1, 1, int.strength.test)
    }
    
    df.env[,i] <- rep(int.val, n)
  }
  
  
  
  for(i in 1:length(top.order)){
    
    # extract name of current variable
    name <- names(which(top.order==i))
    
    # see whether it is a covariate ("X"), the response ("Y"), or an intervention variable ("I")
    name.type <- substr(name, start = 1, stop=1)
    
    
    # do nothing if(name.type == "I")
    
    
    
    if(name.type == "Y"){
      pa <- parents(x=dag.cov, v="Y")
      Ny <- rlogis(n)
      if(length(pa)==0){
        df.env$Y <- ifelse(Ny < 0, 1, 0)
      } else{
        
        # linear combination of the values of the parents and the corresponding weights 
        comb <- lincomb(dat = df.env[, pa, drop = F], coef = weights.SCM[[i]])
        
        if(mod == "logreg"){
          df.env$Y <- ifelse(Ny < comb, 1, 0)
        }
        
        if(mod == "probit"){
          Ny <- rnorm(n, mean = 0, sd = pi/sqrt(3))
          df.env$Y <- ifelse(Ny < comb, 1, 0)
        }
        
        if(mod == "nonlin"){
          fx <- (1/20)*(0.75*comb^3 - 5*comb) + 2*sin(3*comb) 
          df.env$Y <- ifelse(Ny < 3*fx + 1, 1, 0)
        }
        
        if(mod == "bump"){
          fx <- abs(comb + Ny)
          df.env$Y <- ifelse(fx < 2.5, 1, 0)
        }

      }
    }
    
    
    if(name.type == "X"){
      pa <- parents(x=dag.full, v=name)
      error <- rnorm(n, mean = 0, sd = 0.5)
      if(length(pa)==0){
        df.env[,name] <- error
      } else{
        comb <- lincomb(dat = df.env[, pa, drop = F], coef = weights.SCM[[i]])
        df.env[,name] <- comb + error
      }
    }
    
  }
  
  
  # we resample the intervention values if the resulting sample is degenerate as most max.resample times
  max.resamp <- 5
  
  
  sum.zero <- (sum(df.env$Y) <= n*0.05)
  sum.one <- (sum(df.env$Y) >= n*0.95)
  
  
  resamp <- sum.zero || sum.one
  
  
  
  num.resamp <- 0
  
  while(resamp){
    
    num.resamp <- num.resamp + 1

    for(i in (d+2):ncol(df.env)){
      
      int.val <- 0
      
      if(usage == "train"){
        int.val <- t*runif(n=1, min = -int.strength.train, max = int.strength.train)
      }
      
      if(usage == "test"){
        int.val <- t*runifstrong(1, 1, int.strength.test)
      }
      
      df.env[,i] <- rep(int.val, n)
    }
    

    
    for(i in 1:length(top.order)){
      
      # extract name of current variable
      name <- names(which(top.order==i))
      
      # see whether it is a covariate ("X"), the response ("Y"), or an intervention variable ("I")
      name.type <- substr(name, start = 1, stop=1)
      
      
      # do nothing if(name.type == "I")
      
      if(name.type == "Y"){
        pa <- parents(x=dag.cov, v="Y")
        Ny <- rlogis(n)
        if(length(pa)==0){
          df.env$Y <- ifelse(Ny < 0, 1, 0)
        } else{
          
          comb <- lincomb(dat = df.env[, pa, drop = F], coef = weights.SCM[[i]])
          
          if(mod == "logreg"){
            df.env$Y <- ifelse(Ny < comb, 1, 0)
          }
          
          if(mod == "probit"){
            Ny <- rnorm(n, mean = 0, sd = pi/sqrt(3))
            df.env$Y <- ifelse(Ny < comb, 1, 0)
          }
          
          if(mod == "nonlin"){
            fx <- (1/20)*(0.75*comb^3 - 5*comb) + 2*sin(3*comb) 
            df.env$Y <- ifelse(Ny < 3*fx + 1, 1, 0)
          }
          
          if(mod == "bump"){
            fx <- abs(comb + Ny)
            df.env$Y <- ifelse(fx < 2.5, 1, 0)
          }

        }
      }
      
      
      if(name.type == "X"){
        pa <- parents(x=dag.full, v=name)
        error <- rnorm(n, mean = 0, sd = 0.5)
        if(length(pa)==0){
          df.env[,name] <- error
        } else{
          comb <- lincomb(dat = df.env[, pa, drop = F], coef = weights.SCM[[i]])
          df.env[,name] <- comb + error
        }
      }

    }
    
    
    # we only want to resample the interventions at mot max.resamp times (then, we resample the weights)
    resamp.more <- num.resamp < max.resamp
    
    sum.zero <- (sum(df.env$Y) <= n*0.05)
    sum.one <- (sum(df.env$Y) >= n*0.95)
    
    resamp <- sum.zero || sum.one
    
    resamp <- resamp && resamp.more
    
    
  }
  
  
  # here we give the order to resample all the weights if it is difficult to get a non-degenerate dist. for Y
  if(num.resamp >= max.resamp){
    return(NULL)
  }
  
  
  return(df.env)
}


#-------------------------------------------------------------------------------



