# in this script, we run LOEO on the pyroCb dataset


# tuning parameter
top.n <- 25








library(ranger)
library(pROC)
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



envs <- env5

wbce.per.env.loeo <- numeric(length(levels(envs)))
names(wbce.per.env.loeo) <- levels(envs)











for(e in 1:length(levels(envs))){
  
  i.test <- which(envs == levels(envs)[e])
  i.train <- -i.test
  
  X.train <- cube[i.train, ]
  X.val <- cube[i.test, ]
  
  labels.train <- labels[i.train]
  labels.test <- labels[i.test]
  
  y.num.train <- y.num[i.train]
  y.num.test <- y.num[i.test]
  
  train.env <- droplevels(envs[i.train])
  
  
  
  
  
  err.e <- numeric(length(sets))
  
  for(s in 1:length(sets)){
    
    print(paste0("working on set ", s, " for environment ",   e))
    
    set <- sets[[s]]
    
    # if set is empty
    dat.train <- X.train[sample(1:nrow(X.train), size = nrow(X.train), replace = F), ]
    
    
    # if set is not empty
    if(sum(set)>0.001){
      set.indx <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
      dat.train <- X.train[, set.indx]
    }
    

    wbce.train.envs <- numeric(length(levels(train.env)))
    
    for(ee in 1:length(levels(train.env))){
      
      
      test.indx <- which(train.env == levels(train.env)[ee])
      train.indx <- -test.indx
      
      
      
      rf <- ranger(y = labels.train[train.indx], x = dat.train[train.indx, ], probability = T)
      probs <- predict(rf, data = dat.train[test.indx, ])$predictions[,"1"]
      
      wbce.train.envs[ee] <- BCE.weighted(y = y.num.train[test.indx], y.hat = probs)
      
    }
    
    err.e[s] <- max(wbce.train.envs)
    
  }
  
  
  
  
  # top n best sets
  opt.sets.e <- (sort(err.e, index.return=TRUE, decreasing=FALSE)$ix)[1:top.n]
  
  
  
  
  prob.preds <- numeric(length(labels.test))
  
  for(set.indx in 1:length(opt.sets.e)){
    print(paste0("make predictions for env ", e))
    
    set <- sets[[opt.sets.e[set.indx]]]
    
    ind.set <- 1:ncol(X.train)
    
    if(sum(set)>0.001){
      ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
    }
    
    randomizer <- 1:nrow(X.train)
    if(sum(set)< 0.001){
      randomizer <- sample(1:nrow(X.train), size = nrow(X.train), replace = F)
    }
    
  
    
    
    
    
    
    
    table_y <- table(labels.train)  # frequency of each class
    
    weights <- length(labels.train)/(2*table_y)  # one half times inverse of class frequency 
    
    
    rf.mod <- ranger(y = labels.train, x = X.train[randomizer, ind.set], probability = T, class.weights = weights)
    
    prob.preds <- prob.preds + predict(rf.mod, data = X.val[, ind.set])$predictions[,"1"]
    
  }
  
  prob.preds <- prob.preds/length(opt.sets.e)
  
  
  wbce.per.env.loeo[e] <- BCE.weighted(y = y.num.test, y.hat = prob.preds)
  
}




file.name <- paste0("pyroCb_LOEO_top_n_", top.n)




print("save data")

save(wbce.per.env.loeo, file = paste0("../saved_data/", file.name, ".rdata"))





# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, paste0("../sessionInfo/", file.name, ".txt")))













