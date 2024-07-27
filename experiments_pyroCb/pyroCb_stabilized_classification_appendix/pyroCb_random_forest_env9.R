# in this script, we use simple random forests instead of stabilized classification


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





envs <- env9



wbce.per.env.rf.selected.vars <- wbce.per.env.rf.all.vars <- numeric(length(levels(envs)))

set.seed(1)



for(e in 1:length(levels(envs))){
  
  i.test <- which(envs == levels(envs)[e])
  i.train <- -i.test
  
  X.train <- cube[i.train, ]
  X.val <- cube[i.test, ]
  
  labels.train <- labels[i.train]
  labels.test <- labels[i.test]
  
  y.num.train <- y.num[i.train]
  y.num.test <- y.num[i.test]
  
  
  # use all variables this time
  set <- sets[[length(sets)]]
  
  ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
  
  
  rf.mod.sel <- ranger(y = labels.train, x = X.train[, ind.set], probability = T)
  
  pred.probs.sel <- predict(rf.mod.sel, data = X.val[, ind.set])$predictions[,"1"]
  
  wbce.per.env.rf.selected.vars[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.sel)
  
  
  
  rf.mod.all <- ranger(y = labels.train, x = X.train, probability = T)
  pred.probs.all <- predict(rf.mod.all, data = X.val)$predictions[,"1"]
  
  wbce.per.env.rf.all.vars[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.all)

}






save(wbce.per.env.rf.selected.vars,
     wbce.per.env.rf.all.vars,
     file = file.path(script_dir, "saved_data/pyroCb_random_forest_env9.rdata"))





# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/pyroCb_random_forest_env9.txt"))





