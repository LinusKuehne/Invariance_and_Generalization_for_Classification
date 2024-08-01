# this script implements the invariance tests specifically for the pyroCb dataset






#-------------------------------------------------------------------------------
# out-of-distribution prediction functions used by several tests
#-------------------------------------------------------------------------------



# calculate out-of-distribution probability of pyrocb occurrence for given cluster/LOEO CV fold using a 
# random forest trained on all other clusters
# data: a dataframe containing the observations to use
# labels: factor with pyroCb occurrence (1) or non-occurrence (0)
# envVar: environment variable (can be NULL if no environment variable should be used for prediction)
# cluster.assoc: association of every observation to a particular cluster/LOEO CV fold used for the out-of-distribution predictions
# cluster: for which cluster we should make predictions
cluster.probs <- function(data, labels, envVar, cluster.assoc, cluster){
  
  indx_val <- which(cluster.assoc == cluster)
  indx_train <- -indx_val
  
  X_train <- data[indx_train, ]
  X_val <- data[indx_val, ]
  y_train <- labels[indx_train]
  y_val <- labels[indx_val]
  
  # if envVar is not null, we include it in predictions
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
  
  rf <- ranger(y = y_train, x = X_train, probability = T, class.weights = weights, num.trees = 200, max.depth = 20)
  
  y_pred <- predict(rf, data = X_val)$predictions[,"1"]
  
  return(y_pred)
}














# calculate out-of-distribution predictions for the environment variable for given cluster using a 
# random forest trained on all other clusters
# data: a dataframe containing the observations to use
# envs: environment variable 
# cluster.assoc: association of every observation to a particular cluster/LOEO CV fold used for the out-of-distribution predictions
# cluster: for which cluster we should make predictions
cluster.probs.env <- function(data, envs, cluster.assoc, cluster){
  
  indx_val <- which(cluster.assoc == cluster)
  indx_train <- -indx_val
  
  X_train <- data[indx_train, ]
  X_val <- data[indx_val, ]
  env_train <- envs[indx_train]
  env_val <- envs[indx_val]
  
  table_env <- table(env_train)  # frequency of each class
  
  weights <- length(env_train)/(length(levels(envs))*table_env)  # 1/length(levels(envs)) times inverse of class frequency 
  
  rf <- ranger(y = env_train, x = X_train, probability = T, class.weights = weights, num.trees = 200, max.depth = 20)
  
  env_pred <- predict(rf, data = X_val)$predictions
  
  return(env_pred)
}








# wrapper function which calculates out-of-distribution predictions of environment by calling cluster.probs.env for each cluster
# set: subset of predictors which should be used for prediction
# cube: a dataframe containing all observations
# envs: environment variable 
# cluster.assoc: association of every observation to a particular cluster/LOEO CV fold used for the out-of-distribution predictions
# cluster: for which cluster we should make predictions
# posts: vector containing the start indices of every variable in terms of the columns of cube
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





# wrapper function which calculates out-of-distribution predictions of pyrocb occurrence by calling cluster.probs for each cluster
# set: subset of predictors which should be used for prediction
# cube: a dataframe containing all observations
# labels: observations of pyroCb occurrence (0/1)
# envVar: environment variable observations. If it is not NULL, the environment is used for prediction, too
# cluster.assoc: association of every observation to a particular cluster/LOEO CV fold used for the out-of-distribution predictions
# cluster: for which cluster we should make predictions
# posts: vector containing the start indices of every variable in terms of the columns of cube
get.probs <- function(set, cube, labels, envVar = NULL, cluster.assoc, posts){
  
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
    y_pred_clust <- cluster.probs(data = dat, 
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



#-------------------------------------------------------------------------------















#-------------------------------------------------------------------------------
# invariance tests
#-------------------------------------------------------------------------------




# performs DeLong test between the ROC curves of a random forest including environment and 
# a random forest excluding environment
# y_pred_noE: probability prediction of Y without using E by a random forest
# y_pred_E: probability prediction of Y using E by a random forest
# labels: observations of pyroCb occurrence (0/1)
delong.test <- function(y_pred_noE, y_pred_E, labels){
  
  # compute ROC curves
  roc_noE <- roc(response = labels, predictor = y_pred_noE, quiet = T, direction="<")
  roc_E <- roc(response = labels, predictor = y_pred_E, quiet = T, direction="<")
  
  # delong test with pROC
  test.less <- roc.test(roc1 = roc_noE, roc2 = roc_E, method = "delong", alternative = "less")
  test.two.sided <- roc.test(roc1 = roc_noE, roc2 = roc_E, method = "delong", alternative = "two.sided")
  
  ret.list <- list("stat" = test.less$statistic, 
                   "pval_1tail" = test.less$p.value, 
                   "pval_2tail" = test.two.sided$p.value,
                   "auc_E" = test.less$estimate[2],
                   "auc_noE" = test.less$estimate[1])
  
  return(ret.list)
}





# wrapper function which performs the complete delong test given raw data (first, calculates probabilities and then calls delong.test)
# set: subset of predictors which should be used for prediction
# cube: a dataframe containing all observations
# labels: observations of pyroCb occurrence (0/1)
# envVar: environment variable observations. If it is not NULL, the environment is used for prediction, too
# cluster.assoc: association of every observation to a particular cluster/LOEO CV fold used for the out-of-distribution predictions
# posts: vector containing the start indices of every variable in terms of the columns of cube
delong <- function(set, cube, labels, envVar, cluster.assoc, posts){
  
  # permute indices for y_pred_noE
  perm <- sample(1:nrow(envVar), size = nrow(envVar), replace = F)
  
  # compute the out-of-dist. predictions
  y_pred_noE <- get.probs(set, cube, labels, envVar = envVar[perm, , drop = F], cluster.assoc, posts)
  y_pred_E <- get.probs(set, cube, labels, envVar = envVar, cluster.assoc, posts)
  
  # run the DeLong's test
  res <- delong.test(y_pred_noE, y_pred_E, labels)
  
  return(res)
}











# invariant residual distribution test (residual test)
# set: subset of predictors which should be used for prediction
# cube: a dataframe containing all observations
# labels: observations of pyroCb occurrence (0/1)
# y.num: labels converted to numeric vector instead of factor
# group5: grouping into 5 distinct environments
# group9: grouping into 9 distinct environments
# cluster.assoc: association of every observation to a particular cluster/LOEO CV fold used for the out-of-distribution predictions
# posts: vector containing the start indices of every variable in terms of the columns of cube
residual <- function(set, cube, labels, y.num, group5, group9, cluster.assoc, posts){
  
  # compute predictions
  probs <- get.probs(set, cube, labels, envVar = NULL, cluster.assoc, posts)
  
  # compute residuals
  res.y.rf <- y.num - probs
  
  # store in dataframe
  residual.dat <- data.frame(res = res.y.rf, group5 = group5, group9 = group9)
  
  # ANOVA test
  test.means.5 <- oneway.test(res ~ group5, data = residual.dat)
  test.means.9 <- oneway.test(res ~ group9, data = residual.dat)
  
  ret.list <- list("p.val_fiveEnv" = test.means.5$p.value, "p.val_nineEnv" = test.means.9$p.value)
  
  return(ret.list)
}






# correlation test combined for 5 and 9 environments
# set: subset of predictors which should be used for prediction
# cube: a dataframe containing all observations
# labels: observations of pyroCb occurrence (0/1)
# y.num: labels converted to numeric vector instead of factor
# group5: grouping into 5 distinct environments
# group9: grouping into 9 distinct environments
# cluster.assoc: association of every observation to a particular cluster/LOEO CV fold used for the out-of-distribution predictions
# posts: vector containing the start indices of every variable in terms of the columns of cube
corr <- function(set, cube, labels, y.num, group5, group9, cluster.assoc, posts){
  
  # compute predictions
  probs <- get.probs(set, cube, labels, envVar = NULL, cluster.assoc, posts)
  
  # residuals
  res.y.rf <- y.num - probs
  
  # initialize "ground truth" matrix for environment predictions
  indicator.matrix.5 <- matrix(0, nrow = nrow(cube), ncol = length(levels(group5)))
  for(j in 1:nrow(cube)){
    indicator.matrix.5[j,] <- as.numeric(levels(group5) == group5[j])
  }
  indicator.matrix.9 <- matrix(0, nrow = nrow(cube), ncol = length(levels(group9)))
  for(j in 1:nrow(cube)){
    indicator.matrix.9[j,] <- as.numeric(levels(group9) == group9[j])
  }
  
  # compute predictions for environments
  env.res.5 <- get.probs.env(set, cube, group5, cluster.assoc, posts)
  env.res.9 <- get.probs.env(set, cube, group9, cluster.assoc, posts)
  
  # compute environment residuals
  e.res.mat.5 <- indicator.matrix.5 - env.res.5
  e.res.mat.9 <- indicator.matrix.9 - env.res.9
  
  
  pvals.vec.5 <- numeric(ncol(e.res.mat.5))
  
  for(l in 1:length(pvals.vec.5)){
    e.res.5 <- e.res.mat.5[,l]
    test.5 <- cor.test(res.y.rf, e.res.5)
    pvals.vec.5[l] <- test.5$p.value
  }
  
  # bonferroni correction
  p.value.5 <- min(length(pvals.vec.5)*min(pvals.vec.5),1)
  

  pvals.vec.9 <- numeric(ncol(e.res.mat.9))
  
  for(l in 1:length(pvals.vec.9)){
    
    e.res.9 <- e.res.mat.9[,l]
    test.9 <- cor.test(res.y.rf, e.res.9)
    pvals.vec.9[l] <- test.9$p.value
  }
  
  # bonferroni correction
  p.value.9 <- min(length(pvals.vec.9)*min(pvals.vec.9),1)
  

  ret.list <- list("p.val_fiveEnv" = p.value.5, "p.val_nineEnv" = p.value.9)
  
  return(ret.list)
}





# correlation test combined for a single user-specified grouping
# set: subset of predictors which should be used for prediction
# cube: a dataframe containing all observations
# labels: observations of pyroCb occurrence (0/1)
# y.num: labels converted to numeric vector instead of factor
# env_test: user-specified grouping into environments
# cluster.assoc: association of every observation to a particular cluster/LOEO CV fold used for the out-of-distribution predictions
# posts: vector containing the start indices of every variable in terms of the columns of cube
corr.single <- function(set, cube, labels, y.num, env_test, cluster.assoc, posts){
  
  # this code is analogous to the code for corr()
  
  probs <- get.probs(set, cube, labels, envVar = NULL, cluster.assoc, posts)
  
  res.y.rf <- y.num - probs
  
  indicator.matrix <- matrix(0, nrow = nrow(cube), ncol = length(levels(env_test)))
  
  for(j in 1:nrow(cube)){
    indicator.matrix[j,] <- as.numeric(levels(env_test) == env_test[j])
  }
  
  env.res <- get.probs.env(set, cube, env_test, cluster.assoc, posts)
  
  e.res.mat <- indicator.matrix - env.res

  pvals.vec <- numeric(ncol(e.res.mat))
  
  for(l in 1:length(pvals.vec)){
    
    e.res <- e.res.mat[,l]
    test <- cor.test(res.y.rf, e.res)
    pvals.vec[l] <- test$p.value
  }
  
  p.value <- min(length(pvals.vec)*min(pvals.vec),1)

  return(p.value)
}








# TRAM-GCM (RF) test as implemented in arXiv:2309.12833v3
# set: subset of predictors which should be used for prediction
# cube: a dataframe containing all observations
# labels: observations of pyroCb occurrence (0/1)
# y.num: labels converted to numeric vector instead of factor
# env: grouping into distinct environments
# cluster.assoc: association of every observation to a particular cluster/LOEO CV fold used for the out-of-distribution predictions
# posts: vector containing the start indices of every variable in terms of the columns of cube
rangerICP.paper <- function(set, cube, labels, y.num, envs, cluster.assoc, posts){
  
  probs.y <- get.probs(set, cube, labels, envVar = NULL, cluster.assoc, posts)
  
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




#-------------------------------------------------------------------------------



