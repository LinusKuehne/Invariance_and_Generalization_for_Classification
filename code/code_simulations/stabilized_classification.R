




library(ranger)




c.pred <- function(Q, sample, a.pred, B = 100, mod){
  
  s.pred.vec <- numeric(B)
  
  for(b in 1:B){
    
    boot.ind <- base::sample(x = 1:nrow(sample), size = nrow(sample), replace = T)
    boot.sample <- sample[boot.ind,]
    
    preds <- rep(mean(boot.sample$Y), nrow(sample[-boot.ind,]))
    
    if(sum(Q)>0.0001){
      
      if(mod == "RF"){
        boot.model <- ranger(y = as.factor(boot.sample$Y), x = boot.sample[, Q, drop = F], probability = T)
        preds <- predict(boot.model, data = sample[-boot.ind, Q, drop = F])$predictions[,"1"]
      }
      
      if(mod == "GLM"){
       
        x.strings <- paste0("X", Q)
        rhs <- paste(x.strings, collapse = " + ")
        form <- as.formula(paste("Y", "~", rhs))
        boot.model <- glm(formula = form, family = binomial(link = "logit"), data = boot.sample)
        
        preds <- predict(boot.model, newdata = sample[-boot.ind, Q, drop = F], type = "response")
        
      }
    }
    
    s.pred.vec[b] <- nBCE(y = sample[-boot.ind, "Y"], y.hat = preds)
    
  }
  
  return(quantile(x = s.pred.vec, probs = 1-a.pred))
}










model.trainer <- function(sets.train, sample, mod, usage){
  
  models.out <- list()
  nBCE.scores <- numeric(length(sets.train))
  
  # train models.out
  for(s in 1:length(sets.train)){
    set <- sets.train[[s]]
    
    # empty set?
    if(sum(set) < 0.0001){
      
      models.out[[s]] <- mean(sample$Y)
      
      pred <- rep(mean(sample$Y), nrow(sample))
      
      nBCE.scores[s] <- nBCE(y = sample$Y, y.hat = pred)
      
    } else{
      
      if(mod == "RF"){
        rf <- ranger(y = as.factor(sample$Y), x = sample[, set, drop = F], probability = T)
        models.out[[s]] <- rf
        
        preds <- rf$predictions[,"1"]
        
        nBCE.scores[s] <- nBCE(y = sample$Y, y.hat = preds)
      }
      
      if(mod == "GLM"){
        
        x.strings <- paste0("X", set)
        rhs <- paste(x.strings, collapse = " + ")
        form <- as.formula(paste("Y", "~", rhs))
        
        lr <- glm(formula = form, family = binomial(link = "logit"), data = sample)
        
        models.out[[s]] <- lr
        
      }
    }
  }
  
  
  
  if(mod == "GLM" && usage == "train"){
    
    K <- 10
    folds <- base::sample(cut(1:nrow(sample), breaks = K, labels = F), replace = F)
    
    score.mat <- matrix(-100, nrow = K, ncol = length(sets.train))
    
    for(k in 1:K){
      
      test.ind <- which(folds == k)
      samp.train <- sample[-test.ind, ]
      samp.test <- sample[test.ind, ]
      
      
      for(s in 1:length(sets.train)){
        set <- sets.train[[s]]
        
        # empty set?
        if(sum(set) < 0.0001){
          pred <- rep(mean(samp.train$Y), length(test.ind))
          
          score.mat[k, s] <- nBCE(y = samp.test$Y, y.hat = pred)
          
        } else{
          x.strings <- paste0("X", set)
          rhs <- paste(x.strings, collapse = " + ")
          form <- as.formula(paste("Y", "~", rhs))
          lr <- glm(formula = form, family = binomial(link = "logit"), data = samp.train)
            
          pred <- predict(lr, newdata = samp.test[, set, drop = F], type = "response")
            
          score.mat[k, s] <- nBCE(y = samp.test$Y, y.hat = pred)
          
        }
      }
    }
    
    nBCE.scores <- colMeans(score.mat)
  }
  
  
  list.out <- list("models" = models.out, "nBCE.scores" = nBCE.scores)
  
  return(list.out)
  
}






# either sets (list of sets in [d] which should be tested for stability) or d
# (number of covariates) must be given. If only d is given, sets is set to 
# powerSet(1:d). If sets is given, d is ignored
# The first d columns of sample are the covariates named as X1, ..., Xd

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
  vec.s.pred <- m$nBCE.scores
  
  # find most predictive invariant model
  Q.ind <- which.min(vec.s.pred)
  Q <- inv.sets[[Q.ind]]
  
  
  if(verbose){
    cat("The most predictive invariant set is", paste(Q, collapse = " "), "\n")
    print("Compute cutoff c_pred(a_pred)")  
    }
  
  
  c <- c.pred(Q, sample, a.pred = a.pred, B = B, mod = mod.internal)
  
  
  opt.sets <- inv.sets[which(vec.s.pred < c)]
  
  if(length(opt.sets) == 0){
    if(verbose){
      print("no set satisfies s.pred(S) >= c.pred")
    }
    opt.sets[[1]] <- Q
  }
  
  
  if(mod.internal == mod.output){
    opt.models <- m$models[which(vec.s.pred < c)]
    
    if(length(opt.sets) == 0){
      opt.models[[1]] <- m$models[[Q.ind]]
    }
  } else{
    
    mm <- model.trainer(sets.train = opt.sets, sample = sample, mod = mod.output, usage = "test")
    opt.models <- mm$models
    
  }
  
  
  
  result <- list("inv.sets" = inv.sets, "opt.sets" = opt.sets, "opt.models" = opt.models)
  
  return(result)
}






# mod must be the same as mod.output for stabClassMod!!

predict.stabClass <- function(stabClassMod, newsample, mod){
  
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






# I think RFs are bad under strong interventions, because they'd have to extrapolate. Data in these 
# regions never observed have no splits -> bad performance because outside of the observed training 
# data, the RF will be essentially constant.

# If Y were continuous, then linear regression would also work outside of observed data on stable sets
# because E[Y|S] would be linear in x^S. In analogy, I think the GLM should also still work (better at least)




