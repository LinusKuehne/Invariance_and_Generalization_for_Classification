# in this script, we compute the mean and worst-case generalization error for 
# all subsets of 13 variables with random forests

library(ranger)
library(rje)


# get the path of this script
script_dir <- getwd()


# load in the dataset and the division into different environments
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))

# load in utilities
source("../../code/code_pyroCb/pyroCb_stabilized_classification_utils.R")


# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# sets to check invariance
sets <- powerSet(varincl)
sets[[1]] <- c(0)


# convert factor to numeric
y.num <- as.numeric(labels)-1


# initialize vectors to store results for each set
wbce.avg <- wbce.worst <- numeric(length(sets))


# use grouping into five environments
env <- env5


set.seed(1)



# empty set: -------------------------------------------------------------------


wbce.vec.empty <- numeric(length(levels(env)))

# iterate over the environments
for(e in 1:length(levels(env))){
  
  # compute indices for held-out environments
  test.indx <- which(env == levels(env)[e])
  train.indx <- -test.indx
  
  # compute predictions using the empty set of predictors
  probs <- rep(mean(y.num[train.indx]), length(test.indx))
  
  # compute weighted BCE loss for this environment
  wbce.vec.empty[e] <- BCE.weighted(y = y.num[test.indx], y.hat = probs)
}

# aggregate the losses for the empty set
wbce.avg[1] <- mean(wbce.vec.empty)
wbce.worst[1] <- max(wbce.vec.empty)



#-------------------------------------------------------------------------------









# iterate over all non-empty subsets 
for(s in 2:length(sets)){
  
  print(paste0("working on set ", s))
  
  # extract current subset
  set <- sets[[s]]
  
  # obtain columns corresponding to the set
  set.indx <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
  
  # extract the corresponding columns from the predictor matrix
  dat.set <- cube[, set.indx]
  
  
  # initialize a vector for the per-environment losses
  wbce.vec <- numeric(length(levels(env)))
  
  
  # iterate over the different environments
  for(e in 1:length(levels(env))){
    
    # compute the indices for the held-out environment
    test.indx <- which(env == levels(env)[e])
    train.indx <- -test.indx
    
    # compute predictions
    rf <- ranger(y = labels[train.indx], x = dat.set[train.indx, ], probability = T)
    probs <- predict(rf, data = dat.set[test.indx, ])$predictions[,"1"]
    
    # compute weighted BCE loss
    wbce.vec[e] <- BCE.weighted(y = y.num[test.indx], y.hat = probs)
  }
  
  # L^mean
  wbce.avg[s] <- mean(wbce.vec)
  
  # L^worst
  wbce.worst[s] <- max(wbce.vec)
  
}





# save the results
save(wbce.avg, wbce.worst, file = "../saved_data/pyroCb_subsets_generalization.rdata")






# store the sessionInfo:
writeLines(capture.output(sessionInfo()), "../sessionInfo/pyroCb_subsets_generalization.txt")











