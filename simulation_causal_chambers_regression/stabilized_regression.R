library(ranger)



# computes the prediciveness cutoff with a bootstrap procedure
# Smax: most predictive invariant subset
# sample: dataframe where the first d columns are the covariates named X1, ..., Xd. Contains also the response Y
# a.pred: predictiveness tuning parameter
# B: number of bootstrap iterations
c.pred <- function(Smax, sample, a.pred, B = 100){
  
  s.pred.vec <- numeric(B)
  
  for(b in 1:B){
    
    # create bootstrap sample
    boot.ind <- base::sample(x = 1:nrow(sample), size = nrow(sample), replace = T)
    boot.sample <- sample[boot.ind,]
    
    # if Smax is the empty set
    preds <- rep(mean(boot.sample$Y), nrow(sample[-boot.ind,]))
    
    # if Smax is NOT the empty set
    if(sum(Smax)>0.0001){
      
      # fit model on bootstrap sample
      boot.model <- ranger(y = boot.sample$Y, x = boot.sample[, Smax, drop = F], num.threads = 0)
      
      # predict for the observations NOT in the bootstrap sample used for training
      preds <- predict(boot.model, data = sample[-boot.ind, Smax, drop = F])$predictions
    }
    
    # compute the MAE for the current bootstrap iteration
    s.pred.vec[b] <- mean(abs(sample[-boot.ind, "Y"] - preds))
  }
  # return the quantile defining the cutoff parameter c(a.pred)
  return(quantile(x = s.pred.vec, probs = 1-a.pred))
}









# function to fit models on the invariant sets to evaluate predictiveness and later use these models for actual predictions
# sets.train: subset of the glob. variable "sets", denotes the subsets for which we should fit the models
# sample: dataframe where the first d columns are the covariates named X1, ..., Xd. Contains also the response Y
# usage: either "train" (compute out-of-sample MAE losses) or "predict" (don't compute out-of-sample MAE losses, saves time)
model.trainer <- function(sets.train, sample, usage){
  
  # store the fitted models in this list
  models.out <- list()
  MAE.losses <- numeric(length(sets.train))
  
  # train models.out
  for(s in 1:length(sets.train)){
    
    # extract the current set
    set <- sets.train[[s]]
    
    # distinuish whether set is empty or not
    if(sum(set) < 0.0001){
      
      # model trained on empty set is just encoded as the mean of the responses (the best we can do)
      models.out[[s]] <- mean(sample$Y)
      
      pred <- rep(mean(sample$Y), nrow(sample))
      
      # compute score
      MAE.losses[s] <- mean(abs(sample$Y - pred))
      
    } else{
      
      # fit random forest
      rf <- ranger(y = sample$Y, x = sample[, set, drop = F], num.threads = 0, num.trees = 1000)
      
      # add it to list
      models.out[[s]] <- rf
      
      # compute score on OOB samples
      preds <- rf$predictions
      MAE.losses[s] <- mean(abs(sample$Y - preds))
      
    }
  }
  
  list.out <- list("models" = models.out, "MAE.losses" = MAE.losses)
  
  return(list.out)
}





# main function for stabilized regression. Returns the invariant sets, the most predictive inv. sets, and models fitted on these sets
# sample: dataframe where the first d columns are the covariates named X1, ..., Xd. Contains also the response Y and the Env
# test: invariance test, one of ("delong.rf", "delong.glm", "tram.rf", "tram.glm", "corr", "residual")
# a.inv: tuning parameter regarding invariance (set is "invariant" if its p-value is >a.inv)
# a.pred: tuning parameter to compute the predictiveness cutoff parameter
# B: number of bootstrap iterations used to compute the predictiveness cutoff
# verbose: (boolean) whether the function should print out statements regarding the progress
stabilizedRegression <- function(sample, test, a.inv = 0.05, a.pred = 0.05, B = 100, verbose = T){
  
  
  # find invariant sets --------------------------------------------------------
  
  if(verbose){
    print(paste0("Find invariant sets: ", length(sets), " sets to test"))
  }
  
  
  # determine which invariance test we use 
  pval.test.set <- switch(as.character(test),
                          "tram.rf" = pvalues.tram.rf,
                          "residual" = pvalues.residual,
                          stop("Invalid input: no such test implemented"))
  
  
  # compute the p-values for this dataset for all subsets of predictors
  pvals <- pval.test.set(sample)
  
  # determine the empirically invariant subsets
  inv.sets <- sets[which(pvals > a.inv)]
  
  
  # if no set is classified as invariant, we take the "most invariant" one (with largest p-value)
  if(length(inv.sets)==0){
    inv.sets[[1]] <- sets[[which.max(pvals)]]
    
    if(verbose){
      print("no set S satisfies s.inv(S) >= a.inv. Hence, predictions are based on the subset of predictors with largest p-value for the chosen invariance test.")
    }
    
  }

  #-----------------------------------------------------------------------------
  
  
  
  # find most predictive invariant sets ----------------------------------------
  
  if(verbose){
    print("Find most predictive invariant sets")
  }
  
  
  # fit models on all invariant sets and compute s.pred
  m <- model.trainer(sets.train = inv.sets, sample = sample, usage = "train")
  vec.s.pred <- m$MAE.losses
  
  # find most predictive invariant model
  Smax.ind <- which.min(vec.s.pred)
  Smax <- inv.sets[[Smax.ind]]
  
  
  if(verbose){
    cat("The most predictive invariant set is", paste(Smax, collapse = " "), "\n")
    print("Compute cutoff c_pred(a_pred)")  
    }
  
  # compute the cutoff c(a.pred)
  c <- c.pred(Smax, sample, a.pred = a.pred, B = B)
  
  # extract the invariant subsets which are most predictive
  opt.sets <- inv.sets[which(vec.s.pred < c)]
  
  # if no set yields a loss lower than c, we take the best one (Smax)
  if(length(opt.sets) == 0){
    if(verbose){
      print("no set satisfies s.pred(S) >= c.pred")
    }
    opt.sets[[1]] <- Smax
  }
  
  # We have already fitted the models for the final predictions.
  
  # extract fitted models which are predictive enough
  opt.models <- m$models[which(vec.s.pred < c)]
  
  # if no model is chosen by this, we use the one corresponding to the set Smax
  if(length(opt.models) == 0){
    opt.models[[1]] <- m$models[[Smax.ind]]
  }
  
  result <- list("inv.sets" = inv.sets, "opt.sets" = opt.sets, "opt.models" = opt.models, "train.sample" = sample)
  
  return(result)
}







# used to make predictions for new data using a fitted stabilizedRegression object 
# stabClassMod: fitted stabilizedClassificaiton object 
# newsample: new data, same structure as sample used to fit stabClassMod (just needs the covariates)
predict.stabClass <- function(stabClassMod, newsample){
  
  # extract the fitted models and most predictive invariant subsets
  opt.models <- stabClassMod$opt.models
  opt.sets <- stabClassMod$opt.sets
  
  # initialize vector for predicted probabilities
  preds <- matrix(0, nrow = nrow(newsample), ncol = length(opt.models))
  
  # iterate over the most predictive invariant subsets
  for(m in 1:length(opt.models)){

    # if the corresp. set of predictors is empty, opt.models[[...]] is just the mean of the response over the sample
    if(!is.numeric(opt.models[[m]])){
      
      preds[, m] <- predict(opt.models[[m]], data = newsample[, opt.sets[[m]], drop = F])$predictions
      
    } else{
      preds[, m] <- rep(opt.models[[m]], length(preds[, m]))
    }
  }
  
  # aggregate the predicted probabilities of the models in the ensemble
  prediction <- rowMeans(preds)
  
  # return both results
  result.pred <- list("predictions" = prediction)
  
  return(result.pred)
}




# used to make predictions for new data using a fitted stabilizedRegression object 
# stabClassMod: fitted stabilizedClassificaiton object 
# newsample: new data, same structure as sample used to fit stabClassMod (just needs the covariates)
predict.stabClass.hrf <- function(stabClassMod, newsample){
  
  # extract the fitted models and most predictive invariant subsets
  opt.models <- stabClassMod$opt.models
  opt.sets <- stabClassMod$opt.sets
  train.sample <- stabClassMod$train.sample
  
  # initialize vector for predicted probabilities
  preds <- matrix(0, nrow = nrow(newsample), ncol = length(opt.models))
  
  # iterate over the most predictive invariant subsets
  for(m in 1:length(opt.models)){
    
    # if the corresp. set of predictors is empty, opt.models[[...]] is just the mean of the response over the sample
    if(!is.numeric(opt.models[[m]])){
      
      current.rf.fit <- opt.models[[m]]
      current.hrf.fit <- hrf.orig(y = train.sample$Y, x = train.sample[, opt.sets[[m]], drop = F], rf.fit = current.rf.fit)
      
      pred.current.hrf <- predict.hedgedrf(current.hrf.fit, data = newsample[, opt.sets[[m]], drop = F])
      
      preds[, m] <- pred.current.hrf
      
    } else{
      preds[, m] <- rep(opt.models[[m]], length(preds[, m]))
    }
  }
  
  # aggregate the predicted probabilities of the models in the ensemble
  prediction <- rowMeans(preds)
  
  # return both results
  result.pred <- list("predictions" = prediction)
  
  return(result.pred)
}


