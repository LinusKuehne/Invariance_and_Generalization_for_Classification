# in this script, we show that making predictions using an invariant set can be 
# much better than using all available variables in case of (other) interventions 
# in the test distribution


library(ranger)


# get the path of this script
script_dir <- getwd()




# computes binary cross-entropy
# y in {0,1} (label)
# y.hat in (0,1) (probability that Y=1 given some predictors)
BCE <- function(y, y.hat){
  
  y.hat.norm <- y.hat
  y.hat.norm[y.hat == 1] <- 0.99999999999
  y.hat.norm[y.hat == 0] <- 0.00000000001
  
  nbce <- y*log(y.hat.norm) + (1-y)*log(1-y.hat.norm)
  
  return(-mean(nbce))
}





set.seed(1)

# number of training samples
n.train <- 100

# number of test samples
n.test <- 100

# sum of samples
n <- n.train + n.test

# number of experiment repetitions
B <- 200

# store the scores
BCEs.rf <- ACCs.rf <- BCEs.glm <- ACCs.glm <- matrix(0, nrow = B, ncol = 3)
names(BCEs.rf) <- names(ACCs.rf) <- names(BCEs.glm) <- names(ACCs.glm) <- c("{1,2}", "{1}", "{2}")


for(b in 1:B){
  
  print(b)
  
  # generate a data set
  X.1 <- rnorm(n)
  e.y <- rlogis(n)
  Y <- ifelse(e.y < 2*X.1, 1, 0)
  e.2 <- rnorm(n)
  X.2 <- 2*Y + e.2 + c(rep(0, n.train), rep(3, n.test))
  
  # put data into a train and test dataframe
  d.train <- data.frame(Y = as.factor(Y[1:n.train]), X.1 = X.1[1:n.train], X.2 = X.2[1:n.train])
  d.test <- data.frame(Y = as.factor(Y[(1+n.train):n]), X.1 = X.1[(1+n.train):n], X.2 = X.2[(1+n.train):n])
  
  # fit random forest on variables {1,2}
  rf.12 <- ranger(Y ~ ., data = d.train, probability = T)
  preds.rf.12 <- predict(rf.12, data = d.test)$predictions[, "1"]
  BCEs.rf[b, 1] <- BCE(y = Y[(1+n.train):n], y.hat = preds.rf.12)
  ACCs.rf[b, 1] <- mean(d.test$Y == ifelse(preds.rf.12 > 0.5, 1,0))
  
  # fit random forest on variable {1}
  rf.1 <- ranger(Y ~ X.1, data = d.train, probability = T)
  preds.rf.1 <- predict(rf.1, data = d.test)$predictions[, "1"]
  BCEs.rf[b, 2] <- BCE(y = Y[(1+n.train):n], y.hat = preds.rf.1)
  ACCs.rf[b, 2] <- mean(d.test$Y == ifelse(preds.rf.1 > 0.5, 1,0))
  
  # fit random forest on variable {2}
  rf.2 <- ranger(Y ~ X.2, data = d.train, probability = T)
  preds.rf.2 <- predict(rf.2, data = d.test)$predictions[, "1"]
  BCEs.rf[b, 3] <- BCE(y = Y[(1+n.train):n], y.hat = preds.rf.2)
  ACCs.rf[b, 3] <- mean(d.test$Y == ifelse(preds.rf.2 > 0.5, 1,0))
  
  # fit logistic regression on variables {1,2}
  glm.12 <- glm(Y ~ ., data = d.train, family = binomial(link = "logit"))
  preds.glm.12 <- predict(glm.12, newdata = d.test, type = "response")
  BCEs.glm[b, 1] <- BCE(y = Y[(1+n.train):n], y.hat = preds.glm.12)
  ACCs.glm[b, 1] <- mean(d.test$Y == ifelse(preds.glm.12 > 0.5, 1,0))
  
  # fit logistic regression on variable {1}
  glm.1 <- glm(Y ~ X.1, data = d.train, family = binomial(link = "logit"))
  preds.glm.1 <- predict(glm.1, newdata = d.test, type = "response")
  BCEs.glm[b, 2] <- BCE(y = Y[(1+n.train):n], y.hat = preds.glm.1)
  ACCs.glm[b, 2] <- mean(d.test$Y == ifelse(preds.glm.1 > 0.5, 1,0))
  
  # fit logistic regression on variable {2}
  glm.2 <- glm(Y ~ X.2, data = d.train, family = binomial(link = "logit"))
  preds.glm.2 <- predict(glm.2, newdata = d.test, type = "response")
  BCEs.glm[b, 3] <- BCE(y = Y[(1+n.train):n], y.hat = preds.glm.2)
  ACCs.glm[b, 3] <- mean(d.test$Y == ifelse(preds.glm.2 > 0.5, 1,0))
}

# average over all simulation repetitions
colMeans(BCEs.glm)
colMeans(ACCs.glm)

colMeans(BCEs.rf)
colMeans(ACCs.rf)




# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/example_predicting_on_inv_sets.txt"))


