
# calculate out-of-sample prob. of pyrocb occurrence for given cluster using RF trained on all other clusters
cluster.probs.rf <- function(data, labels, envVar, cluster.assoc, cluster){
  
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











# calculate out-of-sample prob. of environment for given cluster using RF trained on all other clusters
cluster.probs.env <- function(data, envs, cluster.assoc, cluster){
  
  indx_val <- which(cluster.assoc == cluster)
  indx_train <- -indx_val
  
  X_train <- data[indx_train, ]
  X_val <- data[indx_val, ]
  env_train <- envs[indx_train]
  env_val <- envs[indx_val]
  
  
  table_env <- table(env_train)  # frequency of each class
  
  weights <- length(env_train)/(length(levels(envs))*table_env)  # 1/length(levels(envs)) times inverse of class frequency 
  
  rf <- ranger(y = env_train, x = X_train, probability = T, class.weights = weights)
  
  env_pred <- predict(rf, data = X_val)$predictions
  
  return(env_pred)
}








# wrapper function which calculates out-of-sample prob. of environment by calling cluster.probs for each cluster
get.probs.env <- function(set, cube, envs, cluster.assoc, posts){
  
  dat <- NA
  
  if(sum(set)<0.00001){
    permute <- sample(1:nrow(cube), size = nrow(cube), replace = F)
    dat <- cube[permute, ]
  } else{
    ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
    dat <- cube[, ind.set]
  }
  
  
  clusts <- unique(cluster.assoc)[order(unique(cluster.assoc))]
  
  
  
  env_pred <- matrix(-1, nrow = nrow(cube), ncol = length(levels(envs)))
  
  for(clust in clusts){
    #print(paste0("Cluster ", clust))
    env_pred_clust <- cluster.probs.env(data = dat, 
                                        envs = envs, 
                                        cluster.assoc = cluster.assoc, 
                                        cluster = clust)
    
    indclust <- which(cluster.assoc == clust)
    
    if(nrow(env_pred_clust) != length(indclust)){
      stop("error")
    }
    
    env_pred[indclust, ] <- env_pred_clust
  }
  
  
  return(env_pred)
}












# wrapper function which calculates out-of-sample prob. of pyrocb occurrence by calling cluster.probs for each cluster
get.probs <- function(set, cube, labels, envVar = NULL, cluster.assoc, posts, model = "RF"){
  
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
  
  
  
  cluster.probs.mod <- cluster.probs.rf
  
  if(model == "GLM"){
    cluster.probs.mod <- cluster.probs.glm
  }
  
  for(clust in clusts){
    #print(paste0("Cluster ", clust))
    y_pred_clust <- cluster.probs.mod(data = dat, 
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

















################################################################################
# hypothesis tests:
################################################################################





# wrapper which performs delong test given raw data (first, calculates probabilities and then calls delong.test)
delong <- function(set, cube, labels, envVar, cluster.assoc, posts){
  
  # permuted indeces for y_pred_noE
  perm <- sample(1:nrow(envVar), size = nrow(envVar), replace = F)
  
  y_pred_noE <- get.probs(set, cube, labels, envVar = envVar[perm, , drop = F], cluster.assoc, posts, model = "RF")
  y_pred_E <- get.probs(set, cube, labels, envVar = envVar, cluster.assoc, posts, model = "RF")
  
  res <- delong.test(y_pred_noE, y_pred_E, labels)
  
  return(res)
}








inv.res.dist <- function(set, cube, labels, y.num, group5, group9, cluster.assoc, posts){
  

  probs <- get.probs(set, cube, labels, envVar = NULL, cluster.assoc, posts, model = "RF")
  
  res.y.rf <- y.num - probs
  
  residual.dat <- data.frame(res = res.y.rf, group5 = group5, group9 = group9)
  
  

  test.means.5 <- oneway.test(res ~ group5, data = residual.dat)
  test.means.9 <- oneway.test(res ~ group9, data = residual.dat)
  
  
  ret.list <- list("p.val_fiveEnv" = test.means.5$p.value, "p.val_nineEnv" = test.means.9$p.value)
  
  return(ret.list)
  
}











rangerICP.paper <- function(set, cube, labels, y.num, envs, cluster.assoc, posts){
  
  probs.y <- get.probs(set, cube, labels, envVar = NULL, cluster.assoc, posts, model = "RF")
  
  r <- matrix(y.num - probs.y, ncol = 1)
  
  indicator.matrix <- matrix(0, nrow = nrow(cube), ncol = length(levels(envs)))
  
  for(j in 1:nrow(cube)){
    indicator.matrix[j,] <- as.numeric(levels(envs) == envs[j])
  }
  
  indicator.matrix <- indicator.matrix[,-1]
  
  
  probs.env <- get.probs.env(set, cube, envs, cluster.assoc, posts)
  predictions.mat <- probs.env[,-1]
  
  
  e.res <- indicator.matrix - predictions.mat
  
  
  L <- matrix(r, nrow = nrow(cube), ncol = ncol(e.res)) * e.res
  
  
  
  sigma <- crossprod(L)/nrow(cube) - tcrossprod(colMeans(L))
  
  eig <- eigen(sigma)
  
  siginvhalf <- eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
  
  tstat <- siginvhalf %*% colSums(L)/sqrt(nrow(cube))
  
  p.value <- pchisq(sum(tstat^2), df = ncol(L), lower.tail = FALSE)
  
  return(p.value)
}










