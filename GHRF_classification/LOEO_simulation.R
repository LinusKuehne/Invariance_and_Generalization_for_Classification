library(rje)
library(ranger)
library(ggplot2)
library(tidyr)
library(dplyr)




# load in the functions needed
source("GHRF_classification/invariance_tests.R")
source("GHRF_classification/utils.R")
source("GHRF_classification/hrf.R")
source("GHRF_classification/stabilized_classification.R")




#-------------------------------------------------------------------------------
# Data preparation
#-------------------------------------------------------------------------------

df_obs <- read.table("data_genome/df_obs.csv", stringsAsFactors = TRUE, sep=";", header = TRUE, row.names = 1)
df_inter_obs <- read.table("data_genome/df_inter_obs.csv", stringsAsFactors = TRUE, sep=";", header = TRUE, row.names = 1)
df_total <- rbind(df_obs, df_inter_obs)


boxplot(Y ~ Env, data = df_total)
abline(h=0, col="red")
df_total$Y <- ifelse(df_total$Y > 0, 1, 0)

for(e in levels(df_total$Env)){
  y.vec <- df_total[df_total$Env == e, "Y"]
  print(e)
  print(mean(y.vec))
  print(length(y.vec))
}


#-------------------------------------------------------------------------------
# Parameters for the experiment
#-------------------------------------------------------------------------------

# set the seed
set.seed(1)

# number of predictors
d <- 9
env.col.idx <- which(names(df_total) == "Env")

# sets to check stability
sets <- powerSet(1:d)
sets[[1]] <- c(0)

# number of bootstrap samples to compute c.pred
B <- 50

# the invariance test to be used
inv_test = "residual"

# number of environments
n.env <- length(levels(df_total$Env))

n.sim <- 10


#-------------------------------------------------------------------------------
# Run the experiment and evaluate on data from held-out environments/domains
#-------------------------------------------------------------------------------


# Initialize empty lists to store results per simulation
accuracies_list <- list()
waccuracies_list <- list()
wBCEscores_list <- list()


for(sim in 1:n.sim){
  print(sim)
  
  idx <- sample(1:nrow(df_obs), size = 750, replace = FALSE)
  df <- rbind(df_obs[idx,], df_inter_obs)
  df$Y <- ifelse(df$Y > 0, 1, 0)
  
  
  accuracies.sim <- waccuracies.sim <- wBCEscores.sim <- data.frame(
    sim = sim,  # Store the simulation number
    env = 1:n.env,  # Store the environment index
    rf = numeric(n.env),
    hrf.orig = numeric(n.env),
    ghrf.1 = numeric(n.env),
    ghrf.2 = numeric(n.env),
    ghrf.3 = numeric(n.env),
    ghrf.4 = numeric(n.env)
    # sc = numeric(n.env),
    # sc.hrf = numeric(n.env)
  )
  
  
  
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
    
    rf.fit <- ranger(y = as.factor(df.train$Y), x = df.train[, 1:d], probability = T, num.threads = 0)
    pred.rf <- predict(rf.fit, data = df.test[,1:d])$predictions[,"1"]
    
    accuracies.sim$rf[e] <- mean(df.test$Y == ifelse(pred.rf>0.5, 1, 0))
    waccuracies.sim$rf[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.rf>0.5, 1, 0))
    wBCEscores.sim$rf[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.rf)
    
    
    # ----------------------------------------------
    # original HRF
    # ----------------------------------------------
    
    hrf.orig.fit <- hrf.orig(y = as.factor(df.train$Y), x = df.train[, 1:d], rf.fit = rf.fit)
    pred.hrf.orig <- predict.hedgedrf(hrf.orig.fit, data = df.test[,1:d])

    accuracies.sim$hrf.orig[e] <- mean(df.test$Y == ifelse(pred.hrf.orig>0.5, 1, 0))
    waccuracies.sim$hrf.orig[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.hrf.orig>0.5, 1, 0))
    wBCEscores.sim$hrf.orig[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.hrf.orig)

    
    # ----------------------------------------------
    # stabilized classification
    # ----------------------------------------------
    
    # fitting and predicting once takes about 4min on my MB
    # sc.fit <- stabilizedClassification(sample = df.train, test = inv_test, B = B, verbose = F)
    # pred.sc <- predict.stabClass(sc.fit, newsample = df.test[,1:d])
    # 
    # accuracies.sim$sc[e] <- mean(df.test$Y == pred.sc$pred.class)
    # waccuracies.sim$sc[e] <- weighted_accuracy(labels = df.test$Y, predictions = pred.sc$pred.class)
    # wBCEscores.sim$sc[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.sc$pred.probs)


    # ----------------------------------------------
    # stabilized classification with HRF classifier
    # ----------------------------------------------

    # pred.sc.hrf <- predict.stabClass.hrf(sc.fit, newsample = df.test[,1:d])
    # 
    # accuracies.sim$sc.hrf[e] <- mean(df.test$Y == pred.sc.hrf$pred.class)
    # waccuracies.sim$sc.hrf[e] <- weighted_accuracy(labels = df.test$Y, predictions = pred.sc.hrf$pred.class)
    # wBCEscores.sim$sc.hrf[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.sc.hrf$pred.probs)


    # ----------------------------------------------
    # GHRF 1
    # ----------------------------------------------

    ghrf.1.fit <- ghrf.1(y = as.factor(df.train$Y), x = df.train[, c(1:d,env.col.idx)])
    pred.ghrf.1 <- predict.ghrf(object = ghrf.1.fit, data = df.test[,1:d])

    accuracies.sim$ghrf.1[e] <- mean(df.test$Y == ifelse(pred.ghrf.1>0.5, 1, 0))
    waccuracies.sim$ghrf.1[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.ghrf.1>0.5, 1, 0))
    wBCEscores.sim$ghrf.1[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.ghrf.1)


    # ----------------------------------------------
    # GHRF 2
    # ----------------------------------------------

    ghrf.2.fit <- ghrf.2(y = as.factor(df.train$Y), x = df.train[, c(1:d,env.col.idx)])
    pred.ghrf.2 <- predict.ghrf(object = ghrf.2.fit, data = df.test[,1:d])

    accuracies.sim$ghrf.2[e] <- mean(df.test$Y == ifelse(pred.ghrf.2>0.5, 1, 0))
    waccuracies.sim$ghrf.2[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.ghrf.2>0.5, 1, 0))
    wBCEscores.sim$ghrf.2[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.ghrf.2)



    # ----------------------------------------------
    # GHRF 3
    # ----------------------------------------------

    ghrf.3.fit <- ghrf.3(y = as.factor(df.train$Y), x = df.train[, c(1:d,env.col.idx)])
    pred.ghrf.3 <- predict.ghrf(object = ghrf.3.fit, data = df.test[,1:d])

    accuracies.sim$ghrf.3[e] <- mean(df.test$Y == ifelse(pred.ghrf.3>0.5, 1, 0))
    waccuracies.sim$ghrf.3[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.ghrf.3>0.5, 1, 0))
    wBCEscores.sim$ghrf.3[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.ghrf.3)


    # ----------------------------------------------
    # GHRF 4
    # ----------------------------------------------

    ghrf.4.fit <- ghrf.4(y = as.factor(df.train$Y), x = df.train[, c(1:d,env.col.idx)])
    pred.ghrf.4 <- predict.ghrf(object = ghrf.4.fit, data = df.test[,1:d])

    accuracies.sim$ghrf.4[e] <- mean(df.test$Y == ifelse(pred.ghrf.4>0.5, 1, 0))
    waccuracies.sim$ghrf.4[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.ghrf.4>0.5, 1, 0))
    wBCEscores.sim$ghrf.4[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.ghrf.4)
    
  }
  
  # Store simulation results
  accuracies_list[[sim]] <- accuracies.sim
  waccuracies_list[[sim]] <- waccuracies.sim
  wBCEscores_list[[sim]] <- wBCEscores.sim
}

# Combine results into long format
accuracies_df <- do.call(rbind, accuracies_list)
waccuracies_df <- do.call(rbind, waccuracies_list)
wBCEscores_df <- do.call(rbind, wBCEscores_list)


save(accuracies_df, waccuracies_df, wBCEscores_df, file = "GHRF_classification/LOEO_classification_results.rdata")




#-------------------------------------------------------------------------------
# analyze the results
#-------------------------------------------------------------------------------

accuracies_avg <- accuracies_df %>%
  group_by(env) %>%
  summarise(across(-sim, mean))

waccuracies_avg <- waccuracies_df %>%
  group_by(env) %>%
  summarise(across(-sim, mean))

wBCEscores_avg <- wBCEscores_df %>%
  group_by(env) %>%
  summarise(across(-sim, mean))



# Convert to long format
waccuracies_long <- pivot_longer(waccuracies_df, cols = -c(sim, env), names_to = "method", values_to = "waccuracy")


# Create boxplot
ggplot(waccuracies_long, aes(x = method, y = waccuracy, fill = method)) +
  geom_boxplot() +
  facet_wrap(~env) +
  theme_bw() +
  labs(title = "Weighted Accuracy Distribution Across Simulations", x = "Method", y = "Weighted Accuracy") +
  theme(axis.text.x = element_text(hjust = 1, angle=45))

ggsave("GHRF_classification/boxplot_LOEO_classification.png")

# create bar plot
ggplot(waccuracies_long, aes(x = factor(env), y = waccuracy, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Mean Accuracy Across Simulations",
       x = "Environment",
       y = "Mean Weighted Accuracy",
       fill = "Method") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1))  # Rotate x-axis labels if needed

ggsave("GHRF_classification/barplot_LOEO_classification.png")

