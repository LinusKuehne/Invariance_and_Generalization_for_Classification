# In this script, we compare the in-distribution and out-of-distribution losses
# of small invariant subsets and of much larger subsets.


library(ranger)
library(rje)



# get the path of this script
script_dir <- getwd()


# load in the dataset
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))


# get utilities
source("../../code/code_pyroCb/pyroCb_stabilized_classification_utils.R")



# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# sets to check stability
sets <- powerSet(varincl)
sets[[1]] <- c(0)




y.num <- as.numeric(labels)-1



# factor describing the association to one of two continents
env_continent <- as.factor(ifelse(event_df$lon > 0, "aus", "am"))


# choose which environment grouping is to be used
#envs <- env5
envs <- env_continent


# define subsets of predictors to use
set.1 <- c(12, 28, 29)  # corr test intersection and invariant by corr test with largest pval
set.2 <- c(12, 13, 28, 29, 30) # invariant by corr test with second largest pval
set.3 <- sets[[length(sets)]] # 13 screened variables
set.4 <- indxIncl # all available variables






# store in-dist. and out-of-dist. errors
err.out.1 <- err.out.2 <- err.out.3 <- err.out.4 <- numeric(length(levels(envs)))
err.in.1 <- err.in.2 <- err.in.3 <- err.in.4 <- numeric(length(levels(envs)))



set.seed(1)



for(e in 1:length(levels(envs))){
  
  print(e)
  
  i.test <- which(envs == levels(envs)[e])
  i.train <- -i.test
  
  
  # generate 5 folds for CV while making sure all observations for the same wildfires
  # are placed in the same fold
  event_df_train <- unlist((event_df$wildfire_id)[i.train])
  event_df_train <- data.frame("wildfire_id" = event_df_train)
  seg_base <- event_df_train[!duplicated(event_df_train), , drop = F]
  folds <- sample(cut(1:nrow(seg_base), breaks = 5, labels = F), replace = F)
  seg_base$clusters_env <- folds
  event_df_train <- merge(x = event_df_train, y = seg_base, all.x = T)
  clusters <- as.factor(event_df_train$clusters_env)
  
  
  X.train <- cube[i.train, ]
  X.val <- cube[i.test, ]
  
  labels.train <- labels[i.train]
  labels.test <- labels[i.test]
  
  y.num.train <- y.num[i.train]
  y.num.test <- y.num[i.test]
  
  
  # get indeces corresponding to columns in "cube"
  ind.set.1 <- as.vector(unlist(sapply(X = set.1, function(i) posts[i]:(posts[i+1]-1))))
  ind.set.2 <- as.vector(unlist(sapply(X = set.2, function(i) posts[i]:(posts[i+1]-1))))
  ind.set.3 <- as.vector(unlist(sapply(X = set.3, function(i) posts[i]:(posts[i+1]-1))))
  ind.set.4 <- as.vector(unlist(sapply(X = set.4, function(i) posts[i]:(posts[i+1]-1))))
  
  
  # get in-dist. CV predictions
  in.probs.1 <- get.probs.sc(set.1, X.train, labels.train, envs = clusters, posts)
  in.probs.2 <- get.probs.sc(set.2, X.train, labels.train, envs = clusters, posts)
  in.probs.3 <- get.probs.sc(set.3, X.train, labels.train, envs = clusters, posts)
  in.probs.4 <- get.probs.sc(set.4, X.train, labels.train, envs = clusters, posts)
  
  # get in-dist. errors
  err.in.1[e] <- BCE.weighted(y = y.num.train, y.hat = in.probs.1)
  err.in.2[e] <- BCE.weighted(y = y.num.train, y.hat = in.probs.2)
  err.in.3[e] <- BCE.weighted(y = y.num.train, y.hat = in.probs.3)
  err.in.4[e] <- BCE.weighted(y = y.num.train, y.hat = in.probs.4)
  
  # get out-of-dist. predictions
  rf.1 <- ranger(y = labels.train, x = X.train[, ind.set.1], probability = T)
  rf.2 <- ranger(y = labels.train, x = X.train[, ind.set.2], probability = T)
  rf.3 <- ranger(y = labels.train, x = X.train[, ind.set.3], probability = T)
  rf.4 <- ranger(y = labels.train, x = X.train[, ind.set.4], probability = T)
  out.probs.1 <- predict(rf.1, data = X.val[, ind.set.1])$predictions[,"1"]
  out.probs.2 <- predict(rf.2, data = X.val[, ind.set.2])$predictions[,"1"]
  out.probs.3 <- predict(rf.3, data = X.val[, ind.set.3])$predictions[,"1"]
  out.probs.4 <- predict(rf.4, data = X.val[, ind.set.4])$predictions[,"1"]
  
  # get out-of-dist. errors
  err.out.1[e] <- BCE.weighted(y = y.num.test, y.hat = out.probs.1)
  err.out.2[e] <- BCE.weighted(y = y.num.test, y.hat = out.probs.2)
  err.out.3[e] <- BCE.weighted(y = y.num.test, y.hat = out.probs.3)
  err.out.4[e] <- BCE.weighted(y = y.num.test, y.hat = out.probs.4)
}

# mean in-dist. error
e.in.1 <- mean(err.in.1)
e.in.2 <- mean(err.in.2)
e.in.3 <- mean(err.in.3)
e.in.4 <- mean(err.in.4)
e.in <- c(e.in.1, e.in.2, e.in.3, e.in.4)

# mean out-of-dist. error
e.out.mean.1 <- mean(err.out.1)
e.out.mean.2 <- mean(err.out.2)
e.out.mean.3 <- mean(err.out.3)
e.out.mean.4 <- mean(err.out.4)
e.out.mean <- c(e.out.mean.1, e.out.mean.2, e.out.mean.3, e.out.mean.4)

# increase when going from mean IN-dist. error to MEAN OUT-of-dist. error
increase.mean.1 <- (e.out.mean.1-e.in.1)/e.in.1
increase.mean.2 <- (e.out.mean.2-e.in.2)/e.in.2
increase.mean.3 <- (e.out.mean.3-e.in.3)/e.in.3
increase.mean.4 <- (e.out.mean.4-e.in.4)/e.in.4
increase.mean <- c(increase.mean.1, increase.mean.2, increase.mean.3, increase.mean.4)

# worst-case out-of-dist. error
e.out.max.1 <- max(err.out.1)
e.out.max.2 <- max(err.out.2)
e.out.max.3 <- max(err.out.3)
e.out.max.4 <- max(err.out.4)
e.out.max <- c(e.out.max.1, e.out.max.2, e.out.max.3, e.out.max.4)

# increase when going from mean IN-dist. error to WORST-CASE OUT-of-dist. error
increase.max.1 <- (e.out.max.1-e.in.1)/e.in.1
increase.max.2 <- (e.out.max.2-e.in.2)/e.in.2
increase.max.3 <- (e.out.max.3-e.in.3)/e.in.3
increase.max.4 <- (e.out.max.4-e.in.4)/e.in.4
increase.max <- c(increase.max.1, increase.max.2, increase.max.3, increase.max.4)


# collect everything 
df <- data.frame("L_in" = e.in,
                 "L_out_mean" = e.out.mean,
                 "increase_mean" = increase.mean,
                 "L_out_max" = e.out.max,
                 "increase_max" = increase.max)













# test what happens to the increase for randomly generated "small sets"

set.seed(1)

# number of iterations
B <- 200

vec.increase.mean <- vec.increase.max <- numeric(B)


for(b in 1:B){
  
  set <- sample(x = varincl, size = 3)
  ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))
  
  
  err.out <- numeric(length(levels(envs)))
  err.in <- numeric(length(levels(envs)))
  
  for(e in 1:length(levels(envs))){
    
    i.test <- which(envs == levels(envs)[e])
    i.train <- -i.test
    
    
    # generate 5 folds for CV while making sure all observations for the same wildfires
    # are placed in the same fold
    event_df_train <- unlist((event_df$wildfire_id)[i.train])
    event_df_train <- data.frame("wildfire_id" = event_df_train)
    seg_base <- event_df_train[!duplicated(event_df_train), , drop = F]
    folds <- sample(cut(1:nrow(seg_base), breaks = 5, labels = F), replace = F)
    seg_base$clusters_env <- folds
    event_df_train <- merge(x = event_df_train, y = seg_base, all.x = T)
    clusters <- as.factor(event_df_train$clusters_env)
    
    
    X.train <- cube[i.train, ]
    X.val <- cube[i.test, ]
    
    labels.train <- labels[i.train]
    labels.test <- labels[i.test]
    
    y.num.train <- y.num[i.train]
    y.num.test <- y.num[i.test]
    
    
    # in-dist. error
    in.probs <- get.probs.sc(set, X.train, labels.train, envs = clusters, posts)
    err.in[e] <- BCE.weighted(y = y.num.train, y.hat = in.probs)
    
    
    # out-of-dist. error
    rf <- ranger(y = labels.train, x = X.train[, ind.set], probability = T)
    out.probs <- predict(rf, data = X.val[, ind.set])$predictions[,"1"]
    err.out[e] <- BCE.weighted(y = y.num.test, y.hat = out.probs)
  }

  vec.increase.mean[b] <- (mean(err.out) - mean(err.in))/mean(err.in)
  vec.increase.max[b] <- (max(err.out) - mean(err.in))/mean(err.in)
}



file.name <- paste0("pyroCb_small_invariant_subsets_E_", length(levels(envs)))


save(df, vec.increase.mean, vec.increase.max, file = file.path(script_dir, paste0("saved_data/", file.name, ".rdata")))


# get confidence intervals
t.test(vec.increase.mean)
t.test(vec.increase.max)




# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, paste0("sessionInfo/", file.name, ".txt")))




