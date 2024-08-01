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


# group into five environments
envs <- env5



wbce.per.env.rf.selected.vars <- wbce.per.env.rf.all.vars <- numeric(length(levels(envs)))

set.seed(1)


# LOEO CV
for(e in 1:length(levels(envs))){
  
  i.test <- which(envs == levels(envs)[e])
  i.train <- -i.test
  
  X.train <- cube[i.train, ]
  X.val <- cube[i.test, ]
  
  labels.train <- labels[i.train]
  labels.test <- labels[i.test]
  
  y.num.train <- y.num[i.train]
  y.num.test <- y.num[i.test]
  
  
  # use all screened variables 
  set <- sets[[length(sets)]]
  
  ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
  
  rf.mod.sel <- ranger(y = labels.train, x = X.train[, ind.set], probability = T)
  
  pred.probs.sel <- predict(rf.mod.sel, data = X.val[, ind.set])$predictions[,"1"]
  
  wbce.per.env.rf.selected.vars[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.sel)
  
  
  # use all available predictors
  rf.mod.all <- ranger(y = labels.train, x = X.train, probability = T)
  pred.probs.all <- predict(rf.mod.all, data = X.val)$predictions[,"1"]
  
  wbce.per.env.rf.all.vars[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.all)

}





# now we repeat this B times to find the effect of randomness in random forest 
# training on the loss

set.seed(1)


B <- 200

wbce.mean.sel <- wbce.mean.all <- numeric(B)

wbce.worst.sel <- wbce.worst.all <- numeric(B)


for(b in 1:B){
  
  
  wbce.per.env.rf.sel.b <- wbce.per.env.rf.all.b <- numeric(length(levels(envs)))
  
  # LOEO CV 
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
    
    wbce.per.env.rf.sel.b[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.sel)
    
    
    
    rf.mod.all <- ranger(y = labels.train, x = X.train, probability = T)
    pred.probs.all <- predict(rf.mod.all, data = X.val)$predictions[,"1"]
    
    wbce.per.env.rf.all.b[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.all)
    
  }
  
  wbce.mean.sel[b] <- mean(wbce.per.env.rf.sel.b)
  wbce.mean.all[b] <- mean(wbce.per.env.rf.all.b)
  
  wbce.worst.sel[b] <- max(wbce.per.env.rf.sel.b)
  wbce.worst.all[b] <- max(wbce.per.env.rf.all.b)

  
}


sd(wbce.mean.sel)
sd(wbce.mean.all)

sd(wbce.worst.sel)
sd(wbce.worst.all)





save(wbce.per.env.rf.selected.vars,
     wbce.per.env.rf.all.vars,
     wbce.worst.all,
     wbce.worst.sel,
     wbce.mean.all,
     wbce.mean.sel,
     file = "../saved_data/pyroCb_random_forest.rdata")





# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/pyroCb_random_forest.txt"))





