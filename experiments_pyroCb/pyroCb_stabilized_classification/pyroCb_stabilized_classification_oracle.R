# in this script, we run stabilized classification on the pyroCb dataset using oracle
# invariance information from the DeLong (RF) test with continuous E


library(ranger)
library(pROC)
library(rje)



# get the path of this script
script_dir <- getwd()


# load in the dataset and the grouping into different environments
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))


# load the stabilized classification utils for the pyroCb data
source("../../code/code_pyroCb/pyroCb_stabilized_classification_utils.R")



# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# sets to check invariance
sets <- powerSet(varincl)
sets[[1]] <- c(0)



# tuning parameters
a.inv <- 0.05
a.pred <- 0.05

# number of bootstrap samples
B <- 100

# convert factor to numeric vector
y.num <- as.numeric(labels)-1


# we use five environments
envs <- env5


# load in the oracle p-values from the continuous DeLong's test
load(file.path(script_dir, "../saved_data/pyroCb_ICP_delong_cont.rdata"))


# use oracle p-values computed with data from all environments
pvals <- pvals.delong.cont.1tail




set.seed(1)


# initialize the losses
wbce.per.env.sc <- data.frame(internal_mean = numeric(length(levels(envs))),
                              internal_worst = numeric(length(levels(envs))))






# worst case loss to rank the predictive subsets s_pred,2(S) -------------------

# iterate over the environments with LOEO CV
for(e in 1:length(levels(envs))){

  # find indices for held-out environments
  i.test <- which(envs == levels(envs)[e])
  i.train <- -i.test

  # split data into train/test according to the held-out environment e
  X.train <- cube[i.train, ]
  X.val <- cube[i.test, ]
  
  labels.train <- labels[i.train]
  labels.test <- labels[i.test]
  
  y.num.train <- y.num[i.train]
  y.num.test <- y.num[i.test]
  
  # exclude the held-out environment
  train.env <- droplevels(envs[i.train])
  
  

  # invariant sets -------------------------------------------------------------
  inv.sets <- sets[pvals > a.inv]
  
  # if no set is invariant, use the one with maximal p-value
  if(length(inv.sets) < 0.0001){
    inv.sets[[1]] <- sets[[which.max(pvals)]]
  }
  
  
  # invariant and class. optimal -----------------------------------------------
  
  # initialize vector for losses corresponding to all empirically invariant subsets
  wbce.worst.vec <- numeric(length(inv.sets))
  
  # iterate over all invariant subsets to compute s_pred
  for(s in 1:length(inv.sets)){
    
    print(paste0("find class opt. sets for env ",e, ": set ",   s, " out of ", length(inv.sets), " for worst-case error"))
    
    # extract current invariant subset
    set <- inv.sets[[s]]
    
    if(sum(set) < 0.01){
      stop("empty set is invariant")
    }
    
    # get corresponding column indices for the predictor matrix
    set.indx <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
    
    # get training data for current invariant subset 
    dat.train <- X.train[, set.indx]
    
    
    # conduct LOEO CV to find worst-case error
    wbce.train.envs <- numeric(length(levels(train.env)))
    for(ee in 1:length(levels(train.env))){
    
      # hold out a further environment
      test.indx <- which(train.env == levels(train.env)[ee])
      train.indx <- -test.indx

      # make RF predictions
      rf <- ranger(y = labels.train[train.indx], x = dat.train[train.indx, ], probability = T)
      probs <- predict(rf, data = dat.train[test.indx, ])$predictions[,"1"]
      
      # compute weighted BCE loss of these predictions
      wbce.train.envs[ee] <- BCE.weighted(y = y.num.train[test.indx], y.hat = probs)
    }
    
    # extract worst-case error
    wbce.worst.vec[s] <- max(wbce.train.envs)
  }
  

  
  # determine c.pred with bootstrap --------------------------------------------
  
  # find best invariant subset and corresponding columns in predictor matrix
  Q.set <- inv.sets[[which.min(wbce.worst.vec)]]
  Q.indx <- as.vector(unlist(sapply(X = Q.set, function(i) posts[i]:(posts[i+1]-1))))
  

  # compute bootstrap errors of Q
  boot.vec <- numeric(B)
  for(b in 1:B){

    # generate bootstrap indices
    boot.ind <- base::sample(x = 1:nrow(X.train), size = nrow(X.train), replace = T)
    
    # get bootstrap sample
    dat.boot <- X.train[boot.ind, Q.indx]
    boot.env.train <- droplevels(train.env[boot.ind])
    
    # initialize loss vector for current iteration
    wbce.boot <- numeric(length(levels(boot.env.train)))
    
    # run LOEO CV for the bootstrap sample (same procedure as above)
    for(eee in 1:length(levels(boot.env.train))){
      test.indx <- which(boot.env.train == levels(boot.env.train)[eee])
      train.indx <- -test.indx
      rf <- ranger(y = (labels.train[boot.ind])[train.indx], x = dat.boot[train.indx, ], probability = T)
      probs <- predict(rf, data = dat.boot[test.indx, ])$predictions[,"1"]
      wbce.boot[eee] <- BCE.weighted(y = (y.num.train[boot.ind])[test.indx], y.hat = probs)
    }
    
    boot.vec[b] <- max(wbce.boot)
  }

  # compute cutoff c(a.pred)
  c.pred <- quantile(x = boot.vec, probs = 1-a.pred)

  # extract the subsets with lower error than c.pred
  inv.class.opt.sets <- inv.sets[which(wbce.worst.vec < c.pred)]
  
  
  # if no set is selected, take the best one
  if(length(inv.class.opt.sets) < 0.0001){
    inv.class.opt.sets[[1]] <- Q.set
  }
  
  
  
  
  # now we make predictions using inv.class.opt.sets ---------------------------

  # initialize vector for predicted probabilities
  prob.preds <- numeric(length(labels.test))
  
  # iterate over all subsets in the ensemble
  for(set.indx in 1:length(inv.class.opt.sets)){
    print(paste0("make predictions for inv.class.opt sets for env ", e, ": set ", set.indx, " out of ", length(inv.class.opt.sets)))
    
    # extract subset and obtain corresponding columns in the predictor matrix
    set <- inv.class.opt.sets[[set.indx]]
    ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
    
    # take care of unbalanced labels
    table_y <- table(labels.train)  
    weights <- length(labels.train)/(2*table_y) 
    
    # make predictions
    rf.mod <- ranger(y = labels.train, x = X.train[, ind.set], probability = T, class.weights = weights)
    prob.preds <- prob.preds + predict(rf.mod, data = X.val[, ind.set])$predictions[,"1"]
  }
  
  # average over the used subsets
  prob.preds <- prob.preds/length(inv.class.opt.sets)
  
  # compute weighted BCE loss
  wbce.per.env.sc$internal_worst[e] <- BCE.weighted(y = y.num.test, y.hat = prob.preds)
}








# mean LOEO CV loss to rank the predictive subsets with s_pred,1(S) ------------

# iterate over the environments with LOEO CV
for(e in 1:length(levels(envs))){
  
  # find indices of held-out environments
  i.test <- which(envs == levels(envs)[e])
  i.train <- -i.test
  
  # split data into train/test according to held-out environment
  X.train <- cube[i.train, ]
  X.val <- cube[i.test, ]
  
  labels.train <- labels[i.train]
  labels.test <- labels[i.test]
  
  y.num.train <- y.num[i.train]
  y.num.test <- y.num[i.test]
  
  
  # generate 5 folds for CV while making sure all observations for the same wildfires
  # are placed in the same fold
  event_df_train <- unlist((event_df$wildfire_id)[i.train])
  event_df_train <- data.frame("wildfire_id" = event_df_train)
  seg_base <- event_df_train[!duplicated(event_df_train), , drop = F]
  folds <- sample(cut(1:nrow(seg_base), breaks = 5, labels = F), replace = F)
  seg_base$clusters_env <- folds
  event_df_train <- merge(x = event_df_train, y = seg_base, all.x = T)
  
  
  # exclude the held-out environment
  train.env <- droplevels(envs[i.train])
  
  
  # extract the empirically invariant subsets ----------------------------------
  inv.sets <- sets[pvals > a.inv]
  
  # if no set is invariant, use the one with the largest p-value
  if(length(inv.sets) < 0.0001){
    inv.sets[[1]] <- sets[[which.max(pvals)]]
  }
  

  
  # find invariant and class. optimal subsets ----------------------------------
  
  # initialize score vector for the subsets
  wbce.vec <- numeric(length(inv.sets))
  
  # compute CV with these clusters/folds
  clusters <- as.factor(event_df_train$clusters_env)
  
  # iterate over the invariant sets and compute predictiveness score
  for(s in 1:length(inv.sets)){
    
    print(paste0("find class opt. sets for env ",e, ": set ",   s, " out of ", length(inv.sets), " for mean error"))
    
    # extract current set
    set <- inv.sets[[s]]
    
    # compute CV losses
    probs <- get.probs.sc(set, X.train, labels.train, envs = clusters, posts)
    
    # compute weighted BCE loss
    wbce.vec[s] <- BCE.weighted(y = y.num.train, y.hat = probs)
  }
  

  
  # determine the cutoff c.pred(a.pred) with bootstrap -------------------------
  
  # get best performing set
  Q.set <- inv.sets[[which.min(wbce.vec)]]

  
  # find bootstrap performance scores of Q
  boot.vec <- numeric(B)
  for(b in 1:B){
    
    # get bootstrap indices
    boot.ind <- base::sample(x = 1:nrow(X.train), size = nrow(X.train), replace = T)
    
    # define bootstrap sample
    dat.boot <- X.train[boot.ind, ]

    # compute prediction score as before on the bootstrap sample
    probs <- get.probs.sc(set = Q.set, cube = dat.boot, labels = labels.train[boot.ind], envs = clusters, posts)
    boot.vec[b] <- BCE.weighted(y = y.num.train[boot.ind], y.hat = probs)
  }
  
  # compute c(a_pred) as the quantile
  c.pred <- quantile(x = boot.vec, probs = 1-a.pred)
  
  # extract the subsets used for prediction
  inv.class.opt.sets <- inv.sets[which(wbce.vec < c.pred)]
  
  
  # if no set is selected, take the best one
  if(length(inv.class.opt.sets) < 0.0001){
    inv.class.opt.sets[[1]] <- Q.set
  }
  
  
  # initialize the final prediction probabilities
  prob.preds <- numeric(length(labels.test))
  
  # iterate over the sets in inv.class.opt.sets and make predictions
  for(set.indx in 1:length(inv.class.opt.sets)){
    print(paste0("make predictions for inv.class.opt sets for env ", e, ": set ", set.indx, " out of ", length(inv.class.opt.sets)))
    
    # extract subset for prediction and get the corresponding columns in the predictor matrix
    set <- inv.class.opt.sets[[set.indx]]
    ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
    
    # take care of unbalanced labels
    table_y <- table(labels.train) 
    weights <- length(labels.train)/(2*table_y) 
    
    # compute predicted probabilities
    rf.mod <- ranger(y = labels.train, x = X.train[, ind.set], probability = T, class.weights = weights)
    prob.preds <- prob.preds + predict(rf.mod, data = X.val[, ind.set])$predictions[,"1"]
  }
  
  
  # average over all subsets in the ensemble
  prob.preds <- prob.preds/length(inv.class.opt.sets)
  
  # compute loss
  wbce.per.env.sc$internal_mean[e] <- BCE.weighted(y = y.num.test, y.hat = prob.preds)
}





print("save data")

save(wbce.per.env.sc, file = "../saved_data/pyroCb_stabilized_classification_oracle.rdata")





# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/pyroCb_stabilized_classification_oracle.txt"))





