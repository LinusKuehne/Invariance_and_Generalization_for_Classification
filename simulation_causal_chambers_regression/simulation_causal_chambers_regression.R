library(ggplot2)
library(rje)
library(viridis)
library(ranger)

# load in the functions needed
source("simulation_causal_chambers_regression/invariance_tests.R")
source("simulation_causal_chambers_regression/hrf.R")
source("simulation_causal_chambers_regression/stabilized_regression.R")


#-------------------------------------------------------------------------------
# Data preparation (maybe try with dataframe_cont_ext.csv !!!!!!!!!!!!!!!!!!!!)
#-------------------------------------------------------------------------------

df <- read.table("simulation_causal_chambers_regression/dataframe_cont.csv", stringsAsFactors = TRUE, sep=";", header = TRUE, row.names = 1)

boxplot(ir_1 ~ intervention, data = df)

df <- df[,c("red", "green", "blue", "vis_1", "vis_2", "vis_3", "ir_2", "ir_3", "ir_1", "intervention")]

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
inv_test = "residual"

# number of environments
n.env <- length(levels(df$Env))


#-------------------------------------------------------------------------------
# Run the experiment and evaluate on data from held-out environments/domains
#-------------------------------------------------------------------------------

mae <- data.frame(oracle = numeric(n.env),
                  rf = numeric(n.env),
                  hrf.orig = numeric(n.env),
                  hrf.ood = numeric(n.env), 
                  sr = numeric(n.env),
                  sr.hrf = numeric(n.env))

mse <- mae


for(e in 1:n.env){
  
  print(paste0("------ Environment ", e, " ------"))
  
  i.test = which(df$Env == levels(df$Env)[e])
  i.train = which(df$Env != levels(df$Env)[e])
  
  df.test = df[i.test,]
  df.test$Env = droplevels(df.test$Env)
  
  df.train = df[i.train,]
  df.train$Env = droplevels(df.train$Env)


  
  # ----------------------------------------------
  # standard random forest
  # ----------------------------------------------
  
  rf.fit <- ranger(y = df.train$Y, x = df.train[, 1:d], num.threads = 0, num.trees = 1000)
  pred.rf <- predict(rf.fit, data = df.test[,1:d])$predictions
  
  mse$rf[e] <- mean((df.test$Y - pred.rf)^2)
  mae$rf[e] <- mean(abs(df.test$Y - pred.rf))
  
  
  # ----------------------------------------------
  # oracle RF fitted on largest invariant subset
  # ----------------------------------------------
  
  oracle.fit <- ranger(y = df.train$Y, x = df.train[, c(1,2,3,5,6,7,8)], num.threads = 0, num.trees = 1000)
  pred.oracle <- predict(oracle.fit, data = df.test[, c(1,2,3,5,6,7,8)])$predictions

  mse$oracle[e] <- mean((df.test$Y - pred.oracle)^2)
  mae$oracle[e] <- mean(abs(df.test$Y - pred.oracle))
  
  
  # ----------------------------------------------
  # original HRF
  # ----------------------------------------------
  
  hrf.orig.fit <- hrf.orig(y = df.train$Y, x = df.train[, 1:d], rf.fit = rf.fit)
  pred.hrf.orig <- predict.hedgedrf(hrf.orig.fit, data = df.test[,1:d])
  
  mse$hrf.orig[e] <- mean((df.test$Y - pred.hrf.orig)^2)
  mae$hrf.orig[e] <- mean(abs(df.test$Y - pred.hrf.orig))
  
  
  # ----------------------------------------------
  # ood HRF
  # ----------------------------------------------
  
  hrf.ood.fit <- hrf.ood(y = df.train$Y, x = df.train[, c(1:d,env.col.idx)])
  pred.hrf.ood <- predict.hedgedrf(hrf.ood.fit, data = df.test[,1:d])
  
  mse$hrf.ood[e] <- mean((df.test$Y - pred.hrf.ood)^2)
  mae$hrf.ood[e] <- mean(abs(df.test$Y - pred.hrf.ood))
  

  # ----------------------------------------------
  # stabilized regression
  # ----------------------------------------------
  
  sr.fit <- stabilizedRegression(sample = df.train, test = inv_test, a.inv = 0.0001, B = B, verbose = F)
  pred.sr <- predict.stabClass(sr.fit, newsample = df.test[,1:d])$predictions

  mse$sr[e] <- mean((df.test$Y - pred.sr)^2)
  mae$sr[e] <- mean(abs(df.test$Y - pred.sr))

  
  # ----------------------------------------------
  # stabilized regression with HRF predictor
  # ----------------------------------------------

  pred.sr.hrf <- predict.stabClass.hrf(sr.fit, newsample = df.test[,1:d])$predictions

  mse$sr.hrf[e] <- mean((df.test$Y - pred.sr.hrf)^2)
  mae$sr.hrf[e] <- mean(abs(df.test$Y - pred.sr.hrf))
}


save(mae, mse, file = "simulation_causal_chambers_regression/results.rdata")
