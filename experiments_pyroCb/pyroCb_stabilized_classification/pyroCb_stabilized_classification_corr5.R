# in this script, we run stabilized classification on the pyroCb dataset with the correlation test


library(ranger)
library(pROC)
library(rje)



# get the path of this script
script_dir <- getwd()


# load in the dataset and the grouping into different environments
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))

# load the invariance tests and stabilized classification utils for the pyroCb data
source(file.path(script_dir, "../../code/code_pyroCb/pyroCb_invariance_tests.R"))
source(file.path(script_dir, "../../code/code_pyroCb/pyroCb_stabilized_classification_utils.R"))



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

# divide the data into five environments
envs <- env5

# initialize
pvals <- matrix(NA, nrow = length(sets), ncol = length(levels(envs)))

# find invariant subsets when excluding each environment once
for(e in 1:length(levels(envs))){

  # set the seed
  set.seed(1)

  # determine which observations are held out
  i.test <- which(envs == levels(envs)[e])
  i.train <- -i.test

  # find remaining environments
  train.env <- droplevels(envs[i.train])

  # split data into train and test according to held-out environment
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


  # initialize vector for p-values
  pvals.e <- numeric(length(sets))


  # iterate over all subsets of predictors
  for(s in 1:length(sets)){
    
    #extract current set
    set <- sets[[s]]
    
    print(paste0("Compute set ", s, " for env ", e))
    
    # run invariance test
    test.corr <- corr.single(set = set, cube = X.train, labels = labels.train, y.num = y.num.train, env_test = train.env, cluster.assoc = event_df_train$clusters_env, posts)

    
    # if the environments are predicted perfectly, there is no variations in the 
    # residual matrix corresponding to the environments. This yields an NA for the
    # p-value of the correlation test. In this case, since there is no variation, 
    # we cannot conduct the test. Therefore, we also cannot reject the invariance
    # and we therefore set the p-value to 1
    if(is.na(test.corr)){
      test.corr <- 1
    }
    
    # store p-value for this subset
    pvals.e[s] <- test.corr
  }

  # store p-values vector corresponding to the current environment e being held out
  pvals[, e] <- pvals.e

}



# compute how many subsets are invariant at the significance level a.inv for each environment being held out once
# If there are no invariant sets for one environment, we will take the one with highest p-value
num.inv.sets <- numeric(length(levels(envs)))

for(e in 1:length(levels(envs))){
  pvals.e <- pvals[,e]
  num.inv.sets[e] <- sum(pvals.e>a.inv)
}




set.seed(1)

# initialize losses
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
  inv.sets <- sets[pvals[,e] > a.inv]
  
  # if no set is invariant, use the one with maximal p-value
  if(length(inv.sets) < 0.0001){
    inv.sets[[1]] <- sets[[which.max(pvals[,e])]]
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
    # extract worst-case loss
    wbce.worst.vec[s] <- max(wbce.train.envs)
  }
  
  
  
  
  
  # determine c.pred with bootstrap --------------------------------------------
  
  # find best invariant subset and corresponding columns in predictor matrix
  Q.set <- inv.sets[[which.min(wbce.worst.vec)]]
  Q.indx <- as.vector(unlist(sapply(X = Q.set, function(i) posts[i]:(posts[i+1]-1))))
  

  # compute bootstrap errors of Q
  boot.vec <- numeric(B)
  for(b in 1:B){

    # get bootstrap indices
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

  # most predictive invariant subsets
  inv.class.opt.sets <- inv.sets[which(wbce.worst.vec < c.pred)]
  
  
  # if no set is selected, take the best one
  if(length(inv.class.opt.sets) < 0.0001){
    inv.class.opt.sets[[1]] <- Q.set
  }
  
  
  # now we make predictions using inv.class.opt.sets ---------------------------
  
  # initialize vector for predicted probabilities
  prob.preds <- numeric(length(labels.test))
  
  # iterate over all subsets in the ensemlbe
  for(set.indx in 1:length(inv.class.opt.sets)){
    print(paste0("make predictions for inv.class.opt sets for env ", e, ": set ", set.indx, " out of ", length(inv.class.opt.sets)))
    
    # extract subset and obtain corresponding columns in the predictor matrix
    set <- inv.class.opt.sets[[set.indx]]
    ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
    
    # take care of unbalanced labels
    table_y <- table(labels.train)  # frequency of each class
    weights <- length(labels.train)/(2*table_y)  # one half times inverse of class frequency 
    
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
  
  # only consider remaining environments
  train.env <- droplevels(envs[i.train])

  
  # extract the empirically invariant subsets ----------------------------------
  inv.sets <- sets[pvals[,e] > a.inv]
  
  # if no set is invariant, use the one with the largest p-value
  if(length(inv.sets) < 0.0001){
    inv.sets[[1]] <- sets[[which.max(pvals[,e])]]
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
    
    # compute CV predictions
    probs <- get.probs.sc(set, X.train, labels.train, envs = clusters, posts)
    
    # compute weighted BCE loss
    wbce.vec[s] <- BCE.weighted(y = y.num.train, y.hat = probs)
  }
  

  
  # determine the cutoff c.pred(a.pred) with bootstrap -------------------------
  
  # get best performing set
  Q.set <- inv.sets[[which.min(wbce.vec)]]

  
  # compute bootstrap prediction scores of the best set
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
  
  # compute cutoff c(a.pred)
  c.pred <- quantile(x = boot.vec, probs = 1-a.pred)
  
  # define the most predictive invariant subsets
  inv.class.opt.sets <- inv.sets[which(wbce.vec < c.pred)]
  
  
  # if no set is selected, take the best one
  if(length(inv.class.opt.sets) < 0.0001){
    inv.class.opt.sets[[1]] <- Q.set
  }
  
  
  # initialize vector for predicted probabilities
  prob.preds <- numeric(length(labels.test))
  
  
  
  # iterate over all subsets in the ensemble
  for(set.indx in 1:length(inv.class.opt.sets)){
    print(paste0("make predictions for inv.class.opt sets for env ", e, ": set ", set.indx, " out of ", length(inv.class.opt.sets)))
    
    # extract subset for prediction and get the corresponding columns in the predictor matrix
    set <- inv.class.opt.sets[[set.indx]]
    ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
    
    # take care of unbalanced lables
    table_y <- table(labels.train)  # frequency of each class
    weights <- length(labels.train)/(2*table_y)  # one half times inverse of class frequency 
    
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

save(wbce.per.env.sc,
     num.inv.sets,
     file = file.path(script_dir, "../saved_data/pyroCb_stabilized_classification_corr5.rdata"))




# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/pyroCb_stabilized_classification_corr5.txt"))





