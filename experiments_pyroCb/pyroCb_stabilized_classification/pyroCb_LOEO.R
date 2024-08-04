# in this script, we run LOEO on the pyroCb dataset. For the thesis, we consider
# top.n = 5 and top.n = 25


# tuning parameter
top.n <- 5




set.seed(1)

library(ranger)
library(pROC)
library(rje)


# get the path of this script
script_dir <- getwd()


# load in the dataset and the division into the different environments
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))

# load in the utils
source("../../code/code_pyroCb/pyroCb_stabilized_classification_utils.R")


# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# list of all subsets
sets <- powerSet(varincl)
sets[[1]] <- c(0)


# convert factor to numeric vector
y.num <- as.numeric(labels)-1


# use five environments
envs <- env5


# store results
wbce.per.env.loeo <- numeric(length(levels(envs)))
names(wbce.per.env.loeo) <- levels(envs)




# iterate over all environments in LOEO CV
for(e in 1:length(levels(envs))){
  
  # get the indices corresponding to environment e
  i.test <- which(envs == levels(envs)[e])
  i.train <- -i.test
  
  # split data into train/test according to the held-out environment
  X.train <- cube[i.train, ]
  X.val <- cube[i.test, ]
  
  labels.train <- labels[i.train]
  labels.test <- labels[i.test]
  
  y.num.train <- y.num[i.train]
  y.num.test <- y.num[i.test]
  
  # only consider the remaining environments
  train.env <- droplevels(envs[i.train])
  
  # initialize vector to store the out-of-distribution losses corresponding to each
  # subset for the current held-out environment e
  err.e <- numeric(length(sets))
  
  
  # iterate over all subsets
  for(s in 1:length(sets)){
    
    print(paste0("working on set ", s, " for environment ",   e))
    
    # extract current set
    set <- sets[[s]]
    
    # if set is empty, permute the rows 
    dat.train <- X.train[sample(1:nrow(X.train), size = nrow(X.train), replace = F), ]
    
    
    # if set is not empty, use the unpermuted data
    if(sum(set)>0.001){
      
      # get columns corresponding to the current subset
      set.indx <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
      
      # get the columns for this subset
      dat.train <- X.train[, set.indx]
    }
    
    # initialize vector to store the losses for each held out environment from the 
    # remaining environments for this subset
    wbce.train.envs <- numeric(length(levels(train.env)))
    
    # find subset with best out-of-distribution prediction performance on current training folds
    for(ee in 1:length(levels(train.env))){
      
      # get row indices for the held-out environment
      test.indx <- which(train.env == levels(train.env)[ee])
      train.indx <- -test.indx
      
      # compute out-of-distribution predictions
      rf <- ranger(y = labels.train[train.indx], x = dat.train[train.indx, ], probability = T)
      probs <- predict(rf, data = dat.train[test.indx, ])$predictions[,"1"]
      
      # compute the loss
      wbce.train.envs[ee] <- BCE.weighted(y = y.num.train[test.indx], y.hat = probs)
    }
    # aggregate the LOEO CV losses over the environments with the worst-case out-of-distribution loss
    err.e[s] <- max(wbce.train.envs)
  }
  

  # determine the top.n best sets
  opt.sets.e <- (sort(err.e, index.return=TRUE, decreasing=FALSE)$ix)[1:top.n]
  
  # store probability predictions
  prob.preds <- numeric(length(labels.test))
  
  # iterate over the top.n best sets
  for(set.indx in 1:length(opt.sets.e)){
    print(paste0("make predictions for env ", e))
    
    # extract subset
    set <- sets[[opt.sets.e[set.indx]]]
    
    # use all columns if the subset is empty (and later randomize)
    ind.set <- 1:ncol(X.train)
    
    # if the subset is not empty, get the columns 
    if(sum(set)>0.001){
      ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
    }
    
    # if subset is empty, randomize rows
    randomizer <- 1:nrow(X.train)
    if(sum(set)< 0.001){
      randomizer <- sample(1:nrow(X.train), size = nrow(X.train), replace = F)
    }

    # take care of unbalanced classes in current training data
    table_y <- table(labels.train)  # frequency of each class
    weights <- length(labels.train)/(2*table_y)  # one half times inverse of class frequency 
    
    # compute predictions
    rf.mod <- ranger(y = labels.train, x = X.train[randomizer, ind.set], probability = T, class.weights = weights)
    prob.preds <- prob.preds + predict(rf.mod, data = X.val[, ind.set])$predictions[,"1"]
  }
  
  # averate the predictions over the top.n subsets
  prob.preds <- prob.preds/length(opt.sets.e)
  
  # compute the loss
  wbce.per.env.loeo[e] <- BCE.weighted(y = y.num.test, y.hat = prob.preds)
}



# get the file name depending on the value of top.n
file.name <- paste0("pyroCb_LOEO_top_n_", top.n)

print("save data")

# save data
save(wbce.per.env.loeo, file = paste0("../saved_data/", file.name, ".rdata"))


# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, paste0("../sessionInfo/", file.name, ".txt")))













