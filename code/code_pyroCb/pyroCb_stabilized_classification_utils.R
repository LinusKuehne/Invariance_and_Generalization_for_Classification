# this script contains utility functions for stabilized classification on the 
# pyroCb dataset




#-------------------------------------------------------------------------------
# compute out-of-sample predictions
#-------------------------------------------------------------------------------


# calculate out-of-distribution probability of pyrocb occurrence for given cluster using a 
# random forest trained on all other clusters
# data: a dataframe containing the observations to use
# labels: factor with pyroCb occurrence (1) or non-occurrence (0)
# test.ind: row indices for observations which should be used for testing
cluster.probs.sc <- function(data, labels, test.ind){
  
  X_train <- data[-test.ind, ]
  X_val <- data[test.ind, ]
  
  y_train <- labels[-test.ind]
  y_val <- labels[test.ind]
  
  table_y <- table(y_train)  # frequency of each class
  
  weights <- length(y_train)/(2*table_y)  # one half times inverse of class frequency 
  
  rf <- ranger(y = y_train, x = X_train, probability = T)
  
  y_pred <- predict(rf, data = X_val)$predictions[,"1"]
  
  return(y_pred)
}






# wrapper function which calculates the out-of-distribution probability of pyrocb occurrence 
# by calling cluster.probs.sc for each cluster
# set: subset of indxIncl, encodes the variables used for predicition
# cube: observation dataframe for the predictors
# labels: factor with pyroCb occurrence (1) or non-occurrence (0)
# envs: environment vector used to compute the LOEO folds
# posts: vector containing the start indices of every variable in terms of the columns of cube
get.probs.sc <- function(set, cube, labels, envs, posts){
  
  if(sum(set) < 0.0001){
    stop("Is the empty set really invariant?")
  }
  
  ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
  
  dat <- cube[, ind.set]
  
  y_pred <- rep(-1, nrow(cube))
  
  for(e in levels(envs)){
    
    test.ind <- which(envs == e)
    
    y_pred_clust <- cluster.probs.sc(data = dat, labels = labels, test.ind = test.ind)
    
    if(length(y_pred_clust) != length(test.ind)){
      stop("error")
    }
    
    y_pred[test.ind] <- y_pred_clust
  }
  
  return(y_pred)
}

#-------------------------------------------------------------------------------










#-------------------------------------------------------------------------------
# loss functions
#-------------------------------------------------------------------------------


# computes the binary cross-entropy loss
# y in {0,1} (label)
# y.hat in (0,1) (probability that Y=1 given some predictors)
BCE <- function(y, y.hat){
  
  y.hat.norm <- y.hat
  
  # cap extreme values such that the log is defined
  y.hat.norm[y.hat == 1] <- 0.99999999999
  y.hat.norm[y.hat == 0] <- 0.00000000001
  
  bce <- y*log(y.hat.norm) + (1-y)*log(1-y.hat.norm)
  
  return(-mean(bce))
}


# computes the weighted binary cross-entropy loss
# y in {0,1} (label)
# y.hat in (0,1) (probability that Y=1 given some predictors)
BCE.weighted <- function(y, y.hat){

  ones <- mean(y) # proportion of 1's
  zeros <- 1-ones # proportion of 0's
  w1 <- 1/(2*ones)
  w0 <- 1/(2*zeros)
  
  y.hat.norm <- y.hat
  
  # cap extreme values such that the log is defined
  y.hat.norm[y.hat > 0.99999999999] <- 0.99999999999
  y.hat.norm[y.hat < 0.00000000001] <- 0.00000000001
  
  bce <- w1*y*log(y.hat.norm) + w0*(1-y)*log(1-y.hat.norm)
  
  return(-mean(bce))
}


#-------------------------------------------------------------------------------








