

# calculate out-of-sample prob. of pyrocb occurrence for given cluster using RF trained on all other clusters
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












# wrapper function which calculates out-of-sample prob. of pyrocb occurrence by calling cluster.probs for each cluster
get.probs.sc <- function(set, cube, labels, envs, posts){
  
  if(sum(set) < 0.0001){
    stop("Is the empty set really invariant?")
  }
  
  
  ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
  
  dat <- cube[, ind.set]
  
  
  
  y_pred <- rep(-1, nrow(cube))
  

  
  for(e in levels(envs)){
    
    test.ind <- which(envs == e)
    
    #print(paste0("Cluster ", clust))
    
    y_pred_clust <- cluster.probs.sc(data = dat, labels = labels, test.ind = test.ind)
    
    #indclust <- which(cluster.assoc == clust)
    
    if(length(y_pred_clust) != length(test.ind)){
      stop("error")
    }
    
    y_pred[test.ind] <- y_pred_clust
  }
  
  return(y_pred)
  
}






BCE <- function(y, y.hat){
  # computes binary cross-entropy
  # y in {0,1} (label)
  # y.hat in (0,1) (probability that Y=1 given some predictors)
  
  y.hat.norm <- y.hat
  y.hat.norm[y.hat == 1] <- 0.99999999999
  y.hat.norm[y.hat == 0] <- 0.00000000001
  
  bce <- y*log(y.hat.norm) + (1-y)*log(1-y.hat.norm)
  
  return(-mean(bce))
}



BCE.weighted <- function(y, y.hat){
  # computes weighted binary cross-entropy
  # y in {0,1} (label)
  # y.hat in (0,1) (probability that Y=1 given some predictors)
  
  
  ones <- mean(y) # proportion of 1's
  zeros <- 1-ones # proportion of 0's
  w1 <- 1/(2*ones)
  w0 <- 1/(2*zeros)
  
  
  y.hat.norm <- y.hat
  y.hat.norm[y.hat > 0.99999999999] <- 0.99999999999
  y.hat.norm[y.hat < 0.00000000001] <- 0.00000000001
  
  bce <- w1*y*log(y.hat.norm) + w0*(1-y)*log(1-y.hat.norm)
  
  return(-mean(bce))
}










