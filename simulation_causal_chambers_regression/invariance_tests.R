# this scripts contains all invariance tests used.
# the functions "pval.<test name>.set" computes the p-value for a specific set
# the functions "pval.<test name>" computes the p-value for all sets of predictors


# In all of these functions, set is a subset of 1:d, where d is the number of predictors,
# and sample is a dataframe where the first d columns correspond to the d predictors,
# the response is sample$Y, and sample$Env is a factor with the environment indices.

# There is a global variable called sets, which is a list of subsets of 1:d; 
# usually powerSet(1:d) (using the library rje)




library(pROC)
library(tramicp)
library(ranger)





#-------------------------------------------------------------------------------
# TRAM-GCM (RF)
#-------------------------------------------------------------------------------

# p-values for a single subset "set"
# icp is the output of rangerICP from library tramicp
pval.tram.rf.set <- function(set, icp){
  
  # the set is empty
  if(sum(set)<0.00001){
    return(pvalues(icp, "set")["Empty"])
  }
  
  # if the set is not empty, compute a vector of predictor names
  # corresponding to the output of rangerICP
  relevant.cov <- rep("A", length(set))
  
  if(sum(set)>0.00001){
    for(w in 1:length(set)){
      number <- as.character(set[w])
      name <- paste("X", number, sep="") 
      relevant.cov[w] <- name  
    }
    
    # compute the correct string to extract the output of rangerICP
    rhs <- paste(relevant.cov, collapse = "+")
    
    return(pvalues(icp, "set")[rhs])
  }
  
}


# compute the p-values for all subsets of predictors
pvalues.tram.rf <- function(sample){
  
  
  # we set up the correct formula and then apply rangerICP from tramicp
  
  num.Xs <- 0
  for(name in names(sample)){
    if(substr(name, start = 1, stop = 1) == "X"){
      num.Xs <- num.Xs+1
    }
  }
  
  cov.names <- rep("A", num.Xs)
  for(w in 1:num.Xs){
    number <- as.character(w)
    name <- paste("X", number, sep="") 
    cov.names[w] <- name  
  }
  
  # compute a string which is a "sum" of the predictors
  rhs <- paste(cov.names, collapse = " + ")
  
  # transform this into a formula
  form <- as.formula(paste("Y", "~", rhs))
  
  
  # compute the output with rangerICP
  icp <- rangerICP(formula = form, data = sample, env = ~ Env,
                   test = "gcm.test", verbose = FALSE)
  
  
  # initialize vector for the p-values for each subset
  pvals <- numeric(length(sets))
  
  # iterate over all subsets of predictors
  for(s in 1:length(sets)){
    
    # extract current set
    set <- sets[[s]]

    # compute p-value for set with the TRAM-GCM (RF) test
    pval <- pval.tram.rf.set(set, icp)
    
    # in rare cases, it can happen that a NA is returned. 
    # In this case, we randomly sample the p-value
    if(is.na(pval)){
      pvals[s] <- runif(1)
    } else{
      pvals[s] <- pval
    }
  }
  return(pvals)
  
}


#-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------
# Residual test
#-------------------------------------------------------------------------------

# p-values for a single subset "set"
pval.residual.set <- function(sample, set){
  
  # convert to factor
  y <- sample$Y
  
  # if set is empty
  res.y.rf <- y - mean(y)
  
  # if set is not empty
  if(sum(set)>0.00001){
    x <- sample[,set, drop = F]
    mod.y.rf <- ranger(y = y, x = x, num.trees = 1000, num.threads = 0)
    res.y.rf <- sample$Y - mod.y.rf$predictions
  }
  
  # put everything into a dataframe
  residual.dat <- data.frame(res = res.y.rf, group = sample$Env)
  
  # test for different means across environments
  test.means <- oneway.test(res ~ group, data = residual.dat)
  
  # extract p-value
  pval.mean <- test.means$p.value
  return(pval.mean)
  
}




# compute the p-values for all subsets of predictors
pvalues.residual <- function(sample){
  
  # initialize vector for the p-values for each subset
  pvals <- numeric(length(sets))
  
  # iterate over all subsets of predictors
  for(s in 1:length(sets)){
    
    if(s %% 10 == 0){
      print(paste0("Computing p-value for set ", s, " out of ", length(sets)))
    }
    
    # extract current set
    set <- sets[[s]]

    
    # compute p-value for set with the residual test
    pval <- pval.residual.set(sample, set)
    
    # in rare cases, it can happen that a NA is returned. 
    # In this case, we randomly sample the p-value
    if(is.na(pval)){
      pvals[s] <- runif(1)
    } else{
      pvals[s] <- pval
    }
  }
  return(pvals)
  
}


#-------------------------------------------------------------------------------












#-------------------------------------------------------------------------------
# combine all methods
#-------------------------------------------------------------------------------


pvalues.all <- function(sample){

  pvals <- data.frame(tram.rf = pvalues.tram.rf(sample),
                      residual = pvalues.residual(sample))

  return(pvals)
}


#-------------------------------------------------------------------------------










