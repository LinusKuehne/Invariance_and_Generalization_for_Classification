# in this script, the function "stabilizedClassification" and the corresponding
# prediction function are implemented, plus several auxiliary functions used 
# by stabilizedClassification


library(ranger)



# computes the prediciveness cutoff with a bootstrap procedure
# Smax: most predictive invariant subset
# sample: dataframe where the first d columns are the covariates named X1, ..., Xd. Contains also the response Y
# a.pred: predictiveness tuning parameter
# B: number of bootstrap iterations
# mod: model used for fitting (either "RF" or "GLM")
c.pred <- function(Smax, sample, a.pred, B = 100, mod){
  
  s.pred.vec <- numeric(B)
  
  for(b in 1:B){
    
    # create bootstrap sample
    boot.ind <- base::sample(x = 1:nrow(sample), size = nrow(sample), replace = T)
    boot.sample <- sample[boot.ind,]
    
    # if Smax is the empty set
    preds <- rep(mean(boot.sample$Y), nrow(sample[-boot.ind,]))
    
    # if Smax is NOT the empty set
    if(sum(Smax)>0.0001){
      
      if(mod == "RF"){
        boot.model <- ranger(y = as.factor(boot.sample$Y), x = boot.sample[, Smax, drop = F], probability = T)
        preds <- predict(boot.model, data = sample[-boot.ind, Smax, drop = F])$predictions[,"1"]
      }
      
      if(mod == "GLM"){
       
        # create the correct model formula for the glm interface
        x.strings <- paste0("X", Smax)
        rhs <- paste(x.strings, collapse = " + ")
        form <- as.formula(paste("Y", "~", rhs))
        
        boot.model <- glm(formula = form, family = binomial(link = "logit"), data = boot.sample)
        
        preds <- predict(boot.model, newdata = sample[-boot.ind, Smax, drop = F], type = "response")
        
      }
    }
    
    # compute the BCE score
    s.pred.vec[b] <- BCE.weighted(y = sample[-boot.ind, "Y"], y.hat = preds)
    
  }
  
  return(quantile(x = s.pred.vec, probs = 1-a.pred))
}









# function to fit models on the invariant sets to evaluate predictiveness and later use these models for actual predictions
# sets.train: subset of the glob. variable "sets", denotes the subsets for which we should fit the models
# sample: dataframe where the first d columns are the covariates named X1, ..., Xd. Contains also the response Y
# mod: model used for fitting (either "RF" or "GLM")
# usage: either "train" (compute out-of-sample wBCE scores) or "predict" (don't compute out-of-sample wBCE scores, saves time)
model.trainer <- function(sets.train, sample, mod, usage){
  
  # store the fitted models in this list
  models.out <- list()
  BCE.scores <- numeric(length(sets.train))
  
  # train models.out
  for(s in 1:length(sets.train)){
    set <- sets.train[[s]]
    
    # distinuish whether set is empty or not
    if(sum(set) < 0.0001){
      
      # model trained on empty set is just encoded as the mean of the responses (the best we can do)
      models.out[[s]] <- mean(sample$Y)
      
      pred <- rep(mean(sample$Y), nrow(sample))
      
      # compute score
      BCE.scores[s] <- BCE.weighted(y = sample$Y, y.hat = pred)
      
    } else{
      
      if(mod == "RF"){
        rf <- ranger(y = as.factor(sample$Y), x = sample[, set, drop = F], probability = T)
        models.out[[s]] <- rf
        
        
        # compute score on OOB samples
        preds <- rf$predictions[,"1"]
        BCE.scores[s] <- BCE.weighted(y = sample$Y, y.hat = preds)
      }
      
      if(mod == "GLM"){
        
        x.strings <- paste0("X", set)
        rhs <- paste(x.strings, collapse = " + ")
        form <- as.formula(paste("Y", "~", rhs))
        
        lr <- glm(formula = form, family = binomial(link = "logit"), data = sample)
        
        models.out[[s]] <- lr
        
        # we will compute the scores later if usage == "train", because we need CV 
      }
    }
  }
  
  
  
  if(mod == "GLM" && usage == "train"){
    
    # get K-fold CV folds
    K <- 10
    folds <- base::sample(cut(1:nrow(sample), breaks = K, labels = F), replace = F)
    
    score.mat <- matrix(-100, nrow = K, ncol = length(sets.train))
    
    for(k in 1:K){
      
      # separate train and test set in CV
      test.ind <- which(folds == k)
      samp.train <- sample[-test.ind, ]
      samp.test <- sample[test.ind, ]
      
      
      for(s in 1:length(sets.train)){
        set <- sets.train[[s]]
        
        # distinuish whether set is empty or not
        if(sum(set) < 0.0001){
          pred <- rep(mean(samp.train$Y), length(test.ind))
          
          # compute score for this set and this fold
          score.mat[k, s] <- BCE.weighted(y = samp.test$Y, y.hat = pred)
          
        } else{
          x.strings <- paste0("X", set)
          rhs <- paste(x.strings, collapse = " + ")
          form <- as.formula(paste("Y", "~", rhs))
          lr <- glm(formula = form, family = binomial(link = "logit"), data = samp.train)
            
          pred <- predict(lr, newdata = samp.test[, set, drop = F], type = "response")
            
          # compute score for this set and this fold
          score.mat[k, s] <- BCE.weighted(y = samp.test$Y, y.hat = pred)
        }
      }
    }
    
    # average the scores over all folds
    BCE.scores <- colMeans(score.mat)
  }
  
  
  list.out <- list("models" = models.out, "BCE.scores" = BCE.scores)
  
  return(list.out)
  
}





# main function for stabilized classification. Returns the invariant sets, the most predictive inv. sets, and models fitted on these sets
# sample: dataframe where the first d columns are the covariates named X1, ..., Xd. Contains also the response Y and the Env
# test: invariance test, one of ("delong.rf", "delong.glm", "tram.rf", "tram.glm", "corr", "residual")
# a.inv: tuning parameter regarding invariance (set is "invariant" if its p-value is >a.inv)
# a.pred: tuning parameter to compute the predictiveness cutoff parameter
# mod.internal: model (either "RF" or "GLM") used to rank the subsets in terms of predictiveness
# mod.output: model (either "RF" or "GLM") which are fitted on the most predictive invariant subsets and are returned
# B: number of bootstrap iterations used to compute the predictiveness cutoff
# verbose: (boolean) whether the function should print out statements regarding the progress
stabilizedClassification <- function(sample, test, a.inv = 0.05, a.pred = 0.05, mod.internal, mod.output, B = 100, verbose = T){
  
  
  # find invariant sets --------------------------------------------------------
  
  if(verbose){
    print(paste0("Find invariant sets: ", length(sets), " sets to test"))
  }
  
  
  
  pval.test.set <- switch(as.character(test),
                          "delong.rf" = pvalues.delong.rf,
                          "delong.glm" = pvalues.delong.glm,
                          "tram.rf" = pvalues.tram.rf,
                          "tram.glm" = pvalues.tram.glm,
                          "corr" = pvalues.correlation,
                          "residual" = pvalues.residual,
                          stop("Invalid input: no such test implemented"))
  
  
  pvals <- pval.test.set(sample)
  
  inv.sets <- sets[which(pvals > a.inv)]
  
  
  # if no set is classified as invariant, we take the "most invariant" one
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
  m <- model.trainer(sets.train = inv.sets, sample = sample, mod = mod.internal, usage = "train")
  vec.s.pred <- m$BCE.scores
  
  # find most predictive invariant model
  Smax.ind <- which.min(vec.s.pred)
  Smax <- inv.sets[[Smax.ind]]
  
  
  if(verbose){
    cat("The most predictive invariant set is", paste(Smax, collapse = " "), "\n")
    print("Compute cutoff c_pred(a_pred)")  
    }
  
  
  c <- c.pred(Smax, sample, a.pred = a.pred, B = B, mod = mod.internal)
  
  
  opt.sets <- inv.sets[which(vec.s.pred < c)]
  
  if(length(opt.sets) == 0){
    if(verbose){
      print("no set satisfies s.pred(S) >= c.pred")
    }
    opt.sets[[1]] <- Smax
  }
  
  
  if(mod.internal == mod.output){
    opt.models <- m$models[which(vec.s.pred < c)]
    
    if(length(opt.models) == 0){
      opt.models[[1]] <- m$models[[Smax.ind]]
    }
  } else{
    
    mm <- model.trainer(sets.train = opt.sets, sample = sample, mod = mod.output, usage = "test")
    opt.models <- mm$models
    
  }
  
  
  
  result <- list("inv.sets" = inv.sets, "opt.sets" = opt.sets, "opt.models" = opt.models, "mod.output" = mod.output)
  
  return(result)
}






# mod must be the same as mod.output for stabClassMod!!

# used to make predictions for new data using a fitted stabilizedClassificaiton object 
# stabClassMod: fitted stabilizedClassificaiton object 
# newsample: new data, same structure as sample used to fit stabClassMod (just needs the covariates)
predict.stabClass <- function(stabClassMod, newsample){
  
  mod <- stabClassMod$mod.output
  
  opt.models <- stabClassMod$opt.models
  opt.sets <- stabClassMod$opt.sets
  
  preds <- matrix(0, nrow = nrow(newsample), ncol = length(opt.models))
  
  for(m in 1:length(opt.models)){
    

    # if the corresp. set of predictors is empty, opt.models[[...]] is just the mean of the response over the sample
    if(!is.numeric(opt.models[[m]])){
      
      if(mod == "RF"){
        preds[, m] <- predict(opt.models[[m]], data = newsample[, opt.sets[[m]], drop = F])$predictions[,"1"]
      }
      
      if(mod == "GLM"){
        preds[, m] <- predict(opt.models[[m]], newdata = newsample[, opt.sets[[m]], drop = F], type = "response")
      }
      
    } else{
      preds[, m] <- rep(opt.models[[m]], length(preds[, m]))
    }
  }
  
  pred.prob <- rowMeans(preds)
  pred.class <- ifelse(pred.prob > 0.5, 1, 0)
  
  result.pred <- list("pred.probs" = pred.prob, "pred.class" = pred.class)
  
  return(result.pred)
}



