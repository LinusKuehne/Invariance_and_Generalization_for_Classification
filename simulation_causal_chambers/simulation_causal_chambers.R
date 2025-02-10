library(ggplot2)
library(rje)
library(viridis)
library(ranger)

# load in the functions needed
source("simulation_causal_chambers/invariance_tests.R")
source("simulation_causal_chambers/utils.R")
source("simulation_causal_chambers/hrf.R")
source("simulation_causal_chambers/stabilized_classification.R")


#-------------------------------------------------------------------------------
# Data preparation
#-------------------------------------------------------------------------------

df <- read.table("simulation_causal_chambers/dataframe_cont.csv", stringsAsFactors = TRUE, sep=";", header = TRUE, row.names = 1)

boxplot(ir_1 ~ intervention, data = df)

df <- df[df$intervention != "uniform_red_mid",]
df$intervention <- droplevels(df$intervention)

thresh <- 12500
df$Y <- ifelse(df$ir_1 > thresh, 1, 0)

df <- df[,c("red", "green", "blue", "vis_1", "vis_2", "vis_3", "ir_2", "ir_3", "Y", "intervention")]

colnames(df) <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "Y", "Env")


#-------------------------------------------------------------------------------
# Parameters for the experiment
#-------------------------------------------------------------------------------

# set the seed
set.seed(1)

# number of predictors
d <- 8
env.col.idx <- which(names(df) == "Env")

# sets to check stability
sets <- powerSet(1:d)
sets[[1]] <- c(0)

# number of bootstrap samples to compute c.pred
B <- 50

# the invariance test to be used
inv_test = "delong.glm"

# number of environments
n.env <- length(levels(df$Env))


#-------------------------------------------------------------------------------
# Run the experiment and evaluate on data from held-out environments/domains
#-------------------------------------------------------------------------------

accuracies <- data.frame(oracle = numeric(n.env),
                         rf = numeric(n.env),
                         hrf.orig = numeric(n.env),
                         hrf.ood = numeric(n.env), 
                         sc = numeric(n.env),
                         sc.hrf = numeric(n.env))

waccuracies <- accuracies
wBCEscores <- accuracies

for(e in 1:n.env){
  
  print(paste0("------ Environment ", e, " ------"))
  
  i.test = which(df$Env == levels(df$Env)[e])
  i.train = which(df$Env != levels(df$Env)[e])
  
  df.test = df[i.test,]
  df.train = df[i.train,]

  
  # ----------------------------------------------
  # standard random forest
  # ----------------------------------------------
  
  rf.fit <- ranger(y = as.factor(df.train$Y), x = df.train[, 1:d], probability = T, num.threads = 0, num.trees = 1000)
  pred.rf <- predict(rf.fit, data = df.test[,1:d])$predictions[,"1"]
  
  mean(df.test$Y == ifelse(pred.rf>0.5, 1, 0))
  weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.rf>0.5, 1, 0))
  BCE.weighted(y = df.test$Y, y.hat = pred.rf)
  
  accuracies$rf[e] <- mean(df.test$Y == ifelse(pred.rf>0.5, 1, 0))
  waccuracies$rf[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.rf>0.5, 1, 0))
  wBCEscores$rf[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.rf)
  
  
  # ----------------------------------------------
  # oracle RF fitted on largest invariant subset
  # ----------------------------------------------
  
  oracle.fit <- ranger(y = as.factor(df.train$Y), x = df.train[, c(1,2,3,5,6,7,8)], probability = T, num.threads = 0, num.trees = 1000)
  pred.oracle <- predict(oracle.fit, data = df.test[, c(1,2,3,5,6,7,8)])$predictions[,"1"]

  mean(df.test$Y == ifelse(pred.oracle>0.5, 1, 0))
  weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.oracle>0.5, 1, 0))
  BCE.weighted(y = df.test$Y, y.hat = pred.oracle)
  
  accuracies$oracle[e] <- mean(df.test$Y == ifelse(pred.oracle>0.5, 1, 0))
  waccuracies$oracle[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.oracle>0.5, 1, 0))
  wBCEscores$oracle[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.oracle)
  
  
  # ----------------------------------------------
  # original HRF
  # ----------------------------------------------
  
  hrf.orig.fit <- hrf.orig(y = as.factor(df.train$Y), x = df.train[, 1:d], rf.fit = rf.fit)
  pred.hrf.orig <- predict.hedgedrf(hrf.orig.fit, data = df.test[,1:d])
  
  mean(df.test$Y == ifelse(pred.hrf.orig>0.5, 1, 0))
  weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.hrf.orig>0.5, 1, 0))
  BCE.weighted(y = df.test$Y, y.hat = pred.hrf.orig)
  
  accuracies$hrf.orig[e] <- mean(df.test$Y == ifelse(pred.hrf.orig>0.5, 1, 0))
  waccuracies$hrf.orig[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.hrf.orig>0.5, 1, 0))
  wBCEscores$hrf.orig[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.hrf.orig)
  
  
  # ----------------------------------------------
  # ood HRF
  # ----------------------------------------------
  
  hrf.ood.fit <- hrf.ood(y = as.factor(df.train$Y), x = df.train[, c(1:d,env.col.idx)])
  pred.hrf.ood <- predict.hedgedrf(hrf.ood.fit, data = df.test[,1:d])
  
  mean(df.test$Y == ifelse(pred.hrf.ood>0.5, 1, 0))
  weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.hrf.ood>0.5, 1, 0))
  BCE.weighted(y = df.test$Y, y.hat = pred.hrf.ood)
  
  accuracies$hrf.ood[e] <- mean(df.test$Y == ifelse(pred.hrf.ood>0.5, 1, 0))
  waccuracies$hrf.ood[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.hrf.ood>0.5, 1, 0))
  wBCEscores$hrf.ood[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.hrf.ood)
  
  
  # ----------------------------------------------
  # stabilized classification
  # ----------------------------------------------
  
  sc.fit <- stabilizedClassification(sample = df.train, test = inv_test, B = B, verbose = F)
  pred.sc <- predict.stabClass(sc.fit, newsample = df.test[,1:d])

  mean(df.test$Y == pred.sc$pred.class)
  weighted_accuracy(labels = df.test$Y, predictions = pred.sc$pred.class)
  BCE.weighted(y = df.test$Y, y.hat = pred.sc$pred.probs)

  accuracies$sc[e] <- mean(df.test$Y == pred.sc$pred.class)
  waccuracies$sc[e] <- weighted_accuracy(labels = df.test$Y, predictions = pred.sc$pred.class)
  wBCEscores$sc[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.sc$pred.probs)


  # ----------------------------------------------
  # stabilized classification with HRF classifier
  # ----------------------------------------------

  pred.sc.hrf <- predict.stabClass.hrf(sc.fit, newsample = df.test[,1:d])

  mean(df.test$Y == pred.sc.hrf$pred.class)
  weighted_accuracy(labels = df.test$Y, predictions = pred.sc.hrf$pred.class)
  BCE.weighted(y = df.test$Y, y.hat = pred.sc.hrf$pred.probs)

  accuracies$sc.hrf[e] <- mean(df.test$Y == pred.sc.hrf$pred.class)
  waccuracies$sc.hrf[e] <- weighted_accuracy(labels = df.test$Y, predictions = pred.sc.hrf$pred.class)
  wBCEscores$sc.hrf[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.sc.hrf$pred.probs)
}


save(accuracies, waccuracies, wBCEscores, file = "simulation_causal_chambers/results.rdata")
