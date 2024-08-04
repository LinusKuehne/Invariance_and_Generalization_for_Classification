# in this script, we use simple random forests instead of stabilized classification. 
# These are the baseline methods.


library(ranger)
library(rje)



# get the path of this script
script_dir <- getwd()


# load in the dataset and the division into different environments
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))


# load in the utils
source("../../code/code_pyroCb/pyroCb_stabilized_classification_utils.R")



# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# we consider all subsets of the 13 variables
sets <- powerSet(varincl)
sets[[1]] <- c(0)



# convert factor to numeric vector
y.num <- as.numeric(labels)-1


# group into five environments
envs <- env5


# initialize losses for each environment for the LOEO CV
wbce.per.env.rf.selected.vars <- wbce.per.env.rf.all.vars <- numeric(length(levels(envs)))

set.seed(1)


# iterate over all environments in LOEO CV
for(e in 1:length(levels(envs))){
  
  # get indices corresponding to the held-out environment
  i.test <- which(envs == levels(envs)[e])
  i.train <- -i.test
  
  # split data into train/validation
  X.train <- cube[i.train, ]
  X.val <- cube[i.test, ]
  
  labels.train <- labels[i.train]
  labels.test <- labels[i.test]
  
  y.num.train <- y.num[i.train]
  y.num.test <- y.num[i.test]
  
  
  # first, use all screened variables 
  set <- sets[[length(sets)]]
  
  # get the corresponding columns in the predictor matrix
  ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
  
  # train model on these variables
  rf.mod.sel <- ranger(y = labels.train, x = X.train[, ind.set], probability = T)
  
  # compute predictions of probabilities
  pred.probs.sel <- predict(rf.mod.sel, data = X.val[, ind.set])$predictions[,"1"]
  
  # compute weighted BCE losses
  wbce.per.env.rf.selected.vars[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.sel)
  
  
  # use all available predictors for the second model
  rf.mod.all <- ranger(y = labels.train, x = X.train, probability = T)
  pred.probs.all <- predict(rf.mod.all, data = X.val)$predictions[,"1"]
  
  # compute losses for this model as well
  wbce.per.env.rf.all.vars[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.all)
}





# now we repeat this B times to find the effect of randomness in random forest 
# training on the loss

set.seed(1)

# number of iterations
B <- 200


# initialize loss vectors
wbce.mean.sel <- wbce.mean.all <- numeric(B)

wbce.worst.sel <- wbce.worst.all <- numeric(B)


# iterate over the simulations
for(b in 1:B){
  
  # initialize the loss vectors for this iteration
  wbce.per.env.rf.sel.b <- wbce.per.env.rf.all.b <- numeric(length(levels(envs)))
  
  # iterate over all environments in LOEO CV
  for(e in 1:length(levels(envs))){
    
    # get indices corresponding to the held-out environment
    i.test <- which(envs == levels(envs)[e])
    i.train <- -i.test
    
    # split data into train/validation
    X.train <- cube[i.train, ]
    X.val <- cube[i.test, ]
    
    labels.train <- labels[i.train]
    labels.test <- labels[i.test]
    
    y.num.train <- y.num[i.train]
    y.num.test <- y.num[i.test]
    
    
    # first, use all screened variables
    set <- sets[[length(sets)]]
    
    # get the corresponding columns in the predictor matrix
    ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
    
    # train model on these variables
    rf.mod.sel <- ranger(y = labels.train, x = X.train[, ind.set], probability = T)
    
    # compute predictions of probabilities
    pred.probs.sel <- predict(rf.mod.sel, data = X.val[, ind.set])$predictions[,"1"]
    
    # compute weighted BCE losses
    wbce.per.env.rf.sel.b[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.sel)
    
    
    # use all available predictors for the second model
    rf.mod.all <- ranger(y = labels.train, x = X.train, probability = T)
    pred.probs.all <- predict(rf.mod.all, data = X.val)$predictions[,"1"]
    
    # compute losses for this model as well
    wbce.per.env.rf.all.b[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.all)
  }
  
  # compute results for this simulation iteration
  wbce.mean.sel[b] <- mean(wbce.per.env.rf.sel.b)
  wbce.mean.all[b] <- mean(wbce.per.env.rf.all.b)
  
  wbce.worst.sel[b] <- max(wbce.per.env.rf.sel.b)
  wbce.worst.all[b] <- max(wbce.per.env.rf.all.b)

  
}


# compute the standard deviations 
sd(wbce.mean.sel)
sd(wbce.mean.all)

sd(wbce.worst.sel)
sd(wbce.worst.all)




# save the results
save(wbce.per.env.rf.selected.vars,
     wbce.per.env.rf.all.vars,
     wbce.worst.all,
     wbce.worst.sel,
     wbce.mean.all,
     wbce.mean.sel,
     file = "../saved_data/pyroCb_random_forest.rdata")





# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/pyroCb_random_forest.txt"))





