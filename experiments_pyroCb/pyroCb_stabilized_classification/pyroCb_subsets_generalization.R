# in this script, we compute the mean and worst-case generalization error for 
# all subsets of 13 variables




library(ranger)
library(rje)



# get the path of this script
script_dir <- getwd()


# load in the dataset
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))


source("../../code/code_pyroCb/pyroCb_stabilized_classification_utils.R")



# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# sets to check stability
sets <- powerSet(varincl)
sets[[1]] <- c(0)






y.num <- as.numeric(labels)-1


wbce.avg <- wbce.worst <- numeric(length(sets))



env <- env5


set.seed(1)



# empty set: -------------------------------------------------------------------


wbce.vec.empty <- numeric(length(levels(env)))

for(e in 1:length(levels(env))){
  
  test.indx <- which(env == levels(env)[e])
  train.indx <- -test.indx
  
  probs <- rep(mean(y.num[train.indx]), length(test.indx))
  
  wbce.vec.empty[e] <- BCE.weighted(y = y.num[test.indx], y.hat = probs)
  
}

wbce.avg[1] <- mean(wbce.vec.empty)
wbce.worst[1] <- max(wbce.vec.empty)



#-------------------------------------------------------------------------------











for(s in 2:length(sets)){
  
  print(paste0("working on set ", s))
  
  set <- sets[[s]]
  
  set.indx <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
  
  dat.set <- cube[, set.indx]
  
  
  wbce.vec <- numeric(length(levels(env)))
  
  for(e in 1:length(levels(env))){
    
    
    test.indx <- which(env == levels(env)[e])
    train.indx <- -test.indx
    
    rf <- ranger(y = labels[train.indx], x = dat.set[train.indx, ], probability = T)
    probs <- predict(rf, data = dat.set[test.indx, ])$predictions[,"1"]
    
    wbce.vec[e] <- BCE.weighted(y = y.num[test.indx], y.hat = probs)
    
  }
  
  wbce.avg[s] <- mean(wbce.vec)
  wbce.worst[s] <- max(wbce.vec)
  
}






save(wbce.avg, wbce.worst, file = "../saved_data/pyroCb_subsets_generalization.rdata")






# store the sessionInfo:
writeLines(capture.output(sessionInfo()), "../sessionInfo/pyroCb_subsets_generalization.txt")











