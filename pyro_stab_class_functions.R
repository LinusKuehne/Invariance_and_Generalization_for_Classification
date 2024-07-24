

# calculate out-of-sample prob. of pyrocb occurrence for given cluster using RF trained on all other clusters
cluster.probs <- function(data, labels, test.ind){
  
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
get.probs <- function(set, cube, labels, envs, posts){
  
  if(sum(set) < 0.0001){
    stop("Is the empty set really invariant?")
  }
  
  
  ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
  
  dat <- cube[, ind.set]
  
  
  
  y_pred <- rep(-1, nrow(cube))
  

  
  for(e in levels(envs)){
    
    test.ind <- which(envs == e)
    
    #print(paste0("Cluster ", clust))
    
    y_pred_clust <- cluster.probs(data = dat, labels = labels, test.ind = test.ind)
    
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





# calculate out-of-sample prob. of pyrocb occurrence for given cluster using RF trained on all other clusters
cluster.probs.delong <- function(data, labels, envVar, cluster.assoc, cluster){
  
  indx_val <- which(cluster.assoc == cluster)
  indx_train <- -indx_val
  
  X_train <- data[indx_train, ]
  X_val <- data[indx_val, ]
  y_train <- labels[indx_train]
  y_val <- labels[indx_val]
  
  
  if(!is.null(envVar)){
    envVar_train <- envVar[indx_train, , drop = F]
    envVar_val <- envVar[indx_val, , drop = F]
    
    X_train <- cbind(X_train, envVar_train)
    X_val <- cbind(X_val, envVar_val)
    
    if(ncol(envVar_train) == 1){
      names(X_train)[ncol(X_train)] <- "Env"
      names(X_val)[ncol(X_val)] <- "Env"
    }
    
  }
  
  
  table_y <- table(y_train)  # frequency of each class
  
  weights <- length(y_train)/(2*table_y)  # one half times inverse of class frequency 
  
  
  
  rf <- ranger(y = y_train, x = X_train, probability = T, class.weights = weights)
  
  y_pred <- predict(rf, data = X_val)$predictions[,"1"]
  
  return(y_pred)
}










# wrapper function which calculates out-of-sample prob. of pyrocb occurrence by calling cluster.probs for each cluster
get.probs.delong <- function(set, cube, labels, envVar = NULL, cluster.assoc, posts){
  
  dat <- NA
  
  if(sum(set)<0.00001){
    permute <- sample(1:nrow(cube), size = nrow(cube), replace = F)
    dat <- cube[permute, ]
  } else{
    ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
    dat <- cube[, ind.set]
  }
  
  
  clusts <- unique(cluster.assoc)[order(unique(cluster.assoc))]
  
  
  
  y_pred <- rep(-1, nrow(cube))
  
  
  
  
  
  for(clust in clusts){
    #print(paste0("Cluster ", clust))
    y_pred_clust <- cluster.probs.delong(data = dat, 
                                         labels = labels, 
                                         envVar = envVar,
                                         cluster.assoc = cluster.assoc, 
                                         cluster = clust)
    
    indclust <- which(cluster.assoc == clust)
    
    if(length(y_pred_clust) != length(indclust)){
      stop("error")
    }
    
    y_pred[indclust] <- y_pred_clust
  }
  
  
  return(y_pred)
  
}




# performs DeLong test between the ROC curves of RF including environment and RF excluding environment
delong.test <- function(y_pred_noE, y_pred_E, labels){
  
  roc_noE <- roc(response = labels, predictor = y_pred_noE, quiet = T, direction="<")
  roc_E <- roc(response = labels, predictor = y_pred_E, quiet = T, direction="<")
  
  
  test.less <- roc.test(roc1 = roc_noE, roc2 = roc_E, method = "delong", alternative = "less")
  test.two.sided <- roc.test(roc1 = roc_noE, roc2 = roc_E, method = "delong", alternative = "two.sided")
  
  ret.list <- list("stat" = test.less$statistic, 
                   "pval_1tail" = test.less$p.value, 
                   "pval_2tail" = test.two.sided$p.value,
                   "auc_E" = test.less$estimate[2],
                   "auc_noE" = test.less$estimate[1])
  
  return(ret.list)
}







# wrapper which performs delong test given raw data (first, calculates probabilities and then calls delong.test)
delong <- function(set, cube, labels, envVar, cluster.assoc, posts){
  
  # permuted indeces for y_pred_noE
  perm <- sample(1:nrow(envVar), size = nrow(envVar), replace = F)
  
  y_pred_noE <- get.probs.delong(set, cube, labels, envVar = envVar[perm, , drop = F], cluster.assoc, posts)
  y_pred_E <- get.probs.delong(set, cube, labels, envVar = envVar, cluster.assoc, posts)
  
  res <- delong.test(y_pred_noE, y_pred_E, labels)
  
  return(res)
}


















