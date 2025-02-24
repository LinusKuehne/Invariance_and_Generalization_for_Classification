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
# DeLong (RF)
#-------------------------------------------------------------------------------


# p-values for a single subset "set"
pval.delong.rf.set <- function(sample, set){
  
  # first, only extract the Env column in case set is empty
  Xmat <- sample[, "Env", drop = F]
  
  # if set is not empty, add the corresponding columns to Xmat
  if(sum(set)>0.0001){
    Xmat <- sample[, c(set,which(names(sample) == "Env")), drop = F]
  }
  
  # random forest using the environment
  RF.Env <- ranger(y = as.factor(sample$Y), x = Xmat, probability = T)
  
  # now we permute the environments
  Xmat$Env <- Xmat$Env[sample(1:nrow(Xmat), size = nrow(Xmat), replace = F)]
  
  # this random forest doesn't have access to environment information
  RF.noEnv <- ranger(y = as.factor(sample$Y), x = Xmat, probability = T)
  
  
  # we use OOB predictions to compute the ROC curves
  roc.noEnv <- roc(response = sample$Y, predictor = RF.noEnv$predictions[,"1"], quiet = T, direction="<")
  roc.Env <- roc(response = sample$Y, predictor = RF.Env$predictions[,"1"], quiet = T, direction="<")
  
  # run DeLong's test on the ROC curves
  test <- roc.test(roc1 = roc.noEnv, roc2 = roc.Env, method = "delong", alternative = "less")
  
  
  return(test$p.value)
}


# compute the p-values for all subsets of predictors
pvalues.delong.rf <- function(sample){
  
  # initialize vector for the p-values for each subset
  pvals <- numeric(length(sets))
  
  # iterate over all subsets of predictors
  for(s in 1:length(sets)){
    
    # extract current set
    set <- sets[[s]]
    
    # compute p-value for set with the DeLong (RF) test
    pval <- pval.delong.rf.set(sample, set)
    
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
# DeLong (GLM)
#-------------------------------------------------------------------------------


# p-values for a single subset "set"
pval.delong.glm.set <- function(sample, set){
  
  # first only include the response and Env column
  dat <- sample[, c("Y", "Env"), drop = F]
  
  # if set is not empty, we also add columns corresponding to the predictors
  if(sum(set)>0.0001){
    dat <- sample[, c(set, which(names(sample) == "Y"), which(names(sample) == "Env")), drop = F]
  }
  
  # initialize vectors for predictions WITH E and WITHOUT E
  preds.Env <- numeric(nrow(sample))
  preds.noEnv <- numeric(nrow(sample))
  
  # make predictions with K-fold CV (since we can't use OOB predictions as for random forests)
  K <- 10
  
  # compute folds
  folds <- base::sample(cut(1:nrow(sample), breaks = K, labels = F), replace = F)
  
  # cross-validation over the K folds for the model including E
  for(k in 1:K){
    test.ind <- which(folds == k)
    
    # fit and predict with glm model including E
    glm.Env.k <- glm(Y ~ ., data = dat[-test.ind, , drop = F], family = binomial(link = "logit"))
    preds.Env[test.ind] <- predict(glm.Env.k, newdata = dat[test.ind, , drop = F], type = "response")
  }
  


  # now we permute the environments
  dat$Env <- dat$Env[sample(1:nrow(dat), size = nrow(dat), replace = F)]
  
  # make predictions with K-fold CV (since we can't use OOB predictions as for random forests)
  for(k in 1:K){
    test.ind <- which(folds == k)
    
    # fit and predict with glm model with permuted E
    glm.noEnv.k <- glm(Y ~ ., data = dat[-test.ind, , drop = F], family = binomial(link = "logit"))
    preds.noEnv[test.ind] <- predict(glm.noEnv.k, newdata = dat[test.ind, , drop = F], type = "response")
  }
  
  # compute the ROC curves
  roc.noEnv <- roc(response = sample$Y, predictor = preds.noEnv, quiet = T, direction="<")
  roc.Env <- roc(response = sample$Y, predictor = preds.Env, quiet = T, direction="<")
  
  # compute the DeLong's test
  test <- roc.test(roc1 = roc.noEnv, roc2 = roc.Env, method = "delong", alternative = "less")
  
  return(test$p.value)
}


# compute the p-values for all subsets of predictors
pvalues.delong.glm <- function(sample){
  
  # initialize vector for the p-values for each subset
  pvals <- numeric(length(sets))
  
  # iterate over all subsets of predictors
  for(s in 1:length(sets)){
    
    # extract current set
    set <- sets[[s]]
    
    # compute p-value for set with the DeLong (GLM) test
    pval <- pval.delong.glm.set(sample, set)
    
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
# TRAM-GCM (GLM)
#-------------------------------------------------------------------------------

# p-values for a single subset "set"
# icp is the output of rangerICP from library tramicp
pval.tram.glm.set <- function(set, icp){
  
  # if the set is empty, return the corresponding output of glmICP
  if(sum(set)<0.00001){
    return(pvalues(icp, "set")["Empty"])
  }
  
  # if the set is not empty, we compute a string to extract the output of glmICP
  # for the corresponding subset
  relevant.cov <- rep("A", length(set))
  
  if(sum(set)>0.00001){
    for(w in 1:length(set)){
      number <- as.character(set[w])
      name <- paste("X", number, sep="") 
      relevant.cov[w] <- name  
    }
    
    rhs <- paste(relevant.cov, collapse = "+")
    
    return(pvalues(icp, "set")[rhs])
  }
  
}


# compute the p-values for all subsets of predictors
pvalues.tram.glm <- function(sample){
  
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
  
  rhs <- paste(cov.names, collapse = " + ")
  
  # transform into a formula
  form <- as.formula(paste("Y", "~", rhs))
  
  # run glmICP 
  icp <- glmICP(formula = form, data = sample, env = ~ Env,
                family = "binomial", verbose = FALSE)
  
  
  # initialize vector for the p-values for each subset
  pvals <- numeric(length(sets))
  
  # iterate over all subsets of predictors
  for(s in 1:length(sets)){
    
    # extract current set
    set <- sets[[s]]
    
    # compute p-value for set with the TRAM-GCM (GLM) test
    pval <- pval.tram.glm.set(set, icp)
    
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
# Correlation test
#-------------------------------------------------------------------------------

# p-values for a single subset "set"
pval.correlation.set <- function(sample, set){
  
  # convert factor to numeric
  y <- as.numeric(sample$Y)
  e <- sample$Env
  
  # if set is empty:
  res.y.rf <- y - mean(y)
  
  # if set is not empty
  if(sum(set)>0.00001){
    x <- sample[,set, drop = F]
    mod.y.rf <- ranger(y = y, x = x)
    res.y.rf <- y - mod.y.rf$predictions
  }
  
  # initialize one-hot encoding matrix to compute the environment residuals
  indicator.matrix <- matrix(0, nrow = nrow(sample), ncol = length(levels(sample$Env)))
  for(j in 1:nrow(sample)){
    indicator.matrix[j,] <- as.numeric(levels(e) == e[j])
  }
  
  # random guessing if set is empty (here we assume we have the same number of observations per environment)
  predictions.mat <- matrix((1/length(levels(e))), nrow = nrow(indicator.matrix), ncol = ncol(indicator.matrix))
  
  # if set is not empty
  if(sum(set)>0.000001){
    mod.e <- ranger(y = e, x = x, probability = T)
    predictions.mat <- mod.e$predictions
  }
  
  # compute environment residuals
  e.res.mat <- indicator.matrix - predictions.mat
  
  # initialize vector of p-values with respect to each column of e.res.mat
  pvals.vec <- numeric(ncol(e.res.mat))
  
  # iterate over the columns
  for(l in 1:length(pvals.vec)){
    
    # extract column 
    e.res <- e.res.mat[,l]
    
    # run correlation test
    test <- cor.test(res.y.rf, e.res)
    
    # extract p-value
    pvals.vec[l] <- test$p.value
  }
  
  # aggregate p-values with Bonferroni
  p.value <- min(length(pvals.vec)*min(pvals.vec),1)
  
  return(p.value)
}






# compute the p-values for all subsets of predictors
pvalues.correlation <- function(sample){
  
  # initialize vector for the p-values for each subset
  pvals <- numeric(length(sets))
  
  # iterate over all subsets of predictors
  for(s in 1:length(sets)){
    
    # extract current set
    set <- sets[[s]]

    # compute p-value for set with the correlation test
    pval <- pval.correlation.set(sample, set)
    
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
  y <- as.factor(sample$Y)
  
  # if set is empty
  res.y.rf <- sample$Y - mean(sample$Y)
  
  # if set is not empty
  if(sum(set)>0.00001){
    x <- sample[,set, drop = F]
    mod.y.rf <- ranger(y = y, x = x, probability = T)
    res.y.rf <- sample$Y - mod.y.rf$predictions[,"1"]
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
    
    if(s %% 100 == 0){
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

  pvals <- data.frame(delong.rf = pvalues.delong.rf(sample),
                      delong.glm = pvalues.delong.glm(sample),
                      tram.rf = pvalues.tram.rf(sample),
                      tram.glm = pvalues.tram.glm(sample),
                      correlation = pvalues.correlation(sample),
                      residual = pvalues.residual(sample))

  return(pvals)
}


#-------------------------------------------------------------------------------










