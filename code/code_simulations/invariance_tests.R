# this scripts contains all invariance tests used.
# the functions "pval.<test name>.set" computes the p-value for a specific set
# the functions "pval.<test name>" computes the p-value for all sets of predictors


# In all of these functions, set is a subset of 1:p, where p is the number of predictors,
# and sample is a dataframe where the first p columns correspond to the p covariates,
# the response is sample$Y, and sample$Env is a factor with the environment indeces.

# There is a global variable called sets, which is a list of subsets of 1:p; 
# usually powerSet(1:p) (using the library rje)




library(pROC)
library(tramicp)
library(ranger)






#-------------------------------------------------------------------------------
# DeLong (RF)
#-------------------------------------------------------------------------------



pval.delong.rf.set <- function(sample, set){
  Xmat <- sample[, "Env", drop = F]
  
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
  
  test <- roc.test(roc1 = roc.noEnv, roc2 = roc.Env, method = "delong", alternative = "less")
  
  
  return(test$p.value)
}


pvalues.delong.rf <- function(sample){
  pvals <- numeric(length(sets))
  
  for(s in 1:length(sets)){
    set <- sets[[s]]
    
    pval <- pval.delong.rf.set(sample, set)
    
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



pval.delong.glm.set <- function(sample, set){
  dat <- sample[, c("Y", "Env"), drop = F]
  
  if(sum(set)>0.0001){
    dat <- sample[, c(set, which(names(sample) == "Y"), which(names(sample) == "Env")), drop = F]
  }
  

  preds.Env <- numeric(nrow(sample))
  preds.noEnv <- numeric(nrow(sample))
  
  # make predictions with K-fold CV (since we can't use OOB predictions as for random forests)
  K <- 10
  folds <- base::sample(cut(1:nrow(sample), breaks = K, labels = F), replace = F)
  
  for(k in 1:K){
    test.ind <- which(folds == k)
    glm.Env.k <- glm(Y ~ ., data = dat[-test.ind, , drop = F], family = binomial(link = "logit"))
    preds.Env[test.ind] <- predict(glm.Env.k, newdata = dat[test.ind, , drop = F], type = "response")
  }
  


  # now we permute the environments
  dat$Env <- dat$Env[sample(1:nrow(dat), size = nrow(dat), replace = F)]
  
  # make predictions with K-fold CV (since we can't use OOB predictions as for random forests)
  for(k in 1:K){
    test.ind <- which(folds == k)
    
    glm.noEnv.k <- glm(Y ~ ., data = dat[-test.ind, , drop = F], family = binomial(link = "logit"))
    preds.noEnv[test.ind] <- predict(glm.noEnv.k, newdata = dat[test.ind, , drop = F], type = "response")
  }
  
  # compute the ROC curves
  roc.noEnv <- roc(response = sample$Y, predictor = preds.noEnv, quiet = T, direction="<")
  roc.Env <- roc(response = sample$Y, predictor = preds.Env, quiet = T, direction="<")
  
  test <- roc.test(roc1 = roc.noEnv, roc2 = roc.Env, method = "delong", alternative = "less")
  
  return(test$p.value)
}


pvalues.delong.glm <- function(sample){
  pvals <- numeric(length(sets))
  
  for(s in 1:length(sets)){
    set <- sets[[s]]
    
    pval <- pval.delong.glm.set(sample, set)
    
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

# icp is the output of rangerICP from library tramicp
pval.tram.rf.set <- function(set, icp){
  
  
  if(sum(set)<0.00001){
    return(pvalues(icp, "set")["Empty"])
  }
  
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
  
  rhs <- paste(cov.names, collapse = " + ")
  
  
  form <- as.formula(paste("Y", "~", rhs))
  
  
  icp <- rangerICP(formula = form, data = sample, env = ~ Env,
                   test = "gcm.test", verbose = FALSE)
  
  
  pvals <- numeric(length(sets))
  
  for(s in 1:length(sets)){
    set <- sets[[s]]

    pval <- pval.tram.rf.set(set, icp)
    
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


# icp is the output of rangerICP from library tramicp
pval.tram.glm.set <- function(set, icp){
  
  if(sum(set)<0.00001){
    return(pvalues(icp, "set")["Empty"])
  }
  
  
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
  
  
  form <- as.formula(paste("Y", "~", rhs))
  
  
  icp <- glmICP(formula = form, data = sample, env = ~ Env,
                family = "binomial", verbose = FALSE)
  
  
  pvals <- numeric(length(sets))
  
  for(s in 1:length(sets)){
    set <- sets[[s]]
    
    pval <- pval.tram.glm.set(set, icp)
    
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

pval.correlation.set <- function(sample, set){
  
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
  
  
  e.res.mat <- indicator.matrix - predictions.mat
  
  
  pvals.vec <- numeric(ncol(e.res.mat))
  for(l in 1:length(pvals.vec)){
    
    e.res <- e.res.mat[,l]
    test <- cor.test(res.y.rf, e.res)
    pvals.vec[l] <- test$p.value
    
  }
  
  p.value <- min(length(pvals.vec)*min(pvals.vec),1)
  return(p.value)
}







pvalues.correlation <- function(sample){
  pvals <- numeric(length(sets))
  
  for(s in 1:length(sets)){
    set <- sets[[s]]

    pval <- pval.correlation.set(sample, set)
    
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

pval.residual.set <- function(sample, set){
  
  y <- as.factor(sample$Y)
  
  # if set is empty
  res.y.rf <- sample$Y - mean(sample$Y)
  
  # if set is not empty
  if(sum(set)>0.00001){
    x <- sample[,set, drop = F]
    mod.y.rf <- ranger(y = y, x = x, probability = T)
    res.y.rf <- sample$Y - mod.y.rf$predictions[,"1"]
  }
  
  
  residual.dat <- data.frame(res = res.y.rf, group = sample$Env)
  test.means <- oneway.test(res ~ group, data = residual.dat)
  pval.mean <- test.means$p.value
  return(pval.mean)
  
}





pvalues.residual <- function(sample){
  pvals <- numeric(length(sets))
  
  for(s in 1:length(sets)){
    set <- sets[[s]]

    
    
    pval <- pval.residual.set(sample, set)
    
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










