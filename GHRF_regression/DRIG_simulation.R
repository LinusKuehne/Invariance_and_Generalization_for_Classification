library(rje)
library(ranger)
library(ggplot2)
library(tidyr)
library(dplyr)



# load in the functions needed
source("GHRF_regression/hrf.R")




#-------------------------------------------------------------------------------
# Data preparation
#-------------------------------------------------------------------------------

df_obs <- read.table("data_genome/df_obs.csv", stringsAsFactors = TRUE, sep=";", header = TRUE, row.names = 1)
df_inter_obs <- read.table("data_genome/df_inter_obs.csv", stringsAsFactors = TRUE, sep=";", header = TRUE, row.names = 1)
df_inter_hidden <- read.table("data_genome/df_inter_hidden.csv", stringsAsFactors = TRUE, sep=";", header = TRUE, row.names = 1)

idx <- sample(1:nrow(df_obs), size = 750, replace = FALSE)
df <- rbind(df_obs[idx,], df_inter_obs)

boxplot(Y ~ Env, data = df_inter_hidden)
abline(h=0, col="red")

df$Y <- ifelse(df$Y > 0, 1, 0)
df_inter_hidden$Y <- ifelse(df_inter_hidden$Y > 0, 1, 0)


for(e in levels(df_inter_hidden$Env)){
  y.vec <- df_inter_hidden[df_inter_hidden$Env == e, "Y"]
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
env.col.idx <- which(names(df) == "Env")

# sets to check stability
sets <- powerSet(1:d)
sets[[1]] <- c(0)

# number of bootstrap samples to compute c.pred
B <- 50

# the invariance test to be used
inv_test = "residual"

# number of environments
n.env <- length(levels(df_inter_hidden$Env))


#-------------------------------------------------------------------------------
# train models
#-------------------------------------------------------------------------------

# standard random forest
rf.fit <- ranger(y = as.factor(df$Y), x = df[, 1:d], probability = T, num.threads = 0)

# original HRF
hrf.orig.fit <- hrf.orig(y = as.factor(df$Y), x = df[, 1:d], rf.fit = rf.fit)

# stabilized classification
sc.fit <- stabilizedClassification(sample = df, test = inv_test, B = B, verbose = F)

# GHRF 1
ghrf.1.fit <- ghrf.1(y = as.factor(df$Y), x = df[, c(1:d,env.col.idx)])

# GHRF 2
ghrf.2.fit <- ghrf.2(y = as.factor(df$Y), x = df[, c(1:d,env.col.idx)])

# GHRF 3
ghrf.3.fit <- ghrf.3(y = as.factor(df$Y), x = df[, c(1:d,env.col.idx)])

# GHRF 4
ghrf.4.fit <- ghrf.4(y = as.factor(df$Y), x = df[, c(1:d,env.col.idx)])



#-------------------------------------------------------------------------------
# make predictions
#-------------------------------------------------------------------------------

accuracies <- data.frame(rf = numeric(n.env),
                         hrf.orig = numeric(n.env),
                         ghrf.1 = numeric(n.env),
                         ghrf.2 = numeric(n.env),
                         ghrf.3 = numeric(n.env),
                         ghrf.4 = numeric(n.env),
                         sc = numeric(n.env),
                         sc.hrf = numeric(n.env)
)

waccuracies <- accuracies
wBCEscores <- accuracies



for(e in 1:n.env){
  
  print(paste0("------ Environment ", e, " ------"))
  
  i.test = which(df_inter_hidden$Env == levels(df_inter_hidden$Env)[e])

  df.test = df_inter_hidden[i.test,]
  df.test$Env = droplevels(df.test$Env)
  
  
  # ----------------------------------------------
  # standard random forest
  # ----------------------------------------------
  
  pred.rf <- predict(rf.fit, data = df.test[,1:d])$predictions[,"1"]
  
  accuracies$rf[e] <- mean(df.test$Y == ifelse(pred.rf>0.5, 1, 0))
  waccuracies$rf[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.rf>0.5, 1, 0))
  wBCEscores$rf[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.rf)
  
  
  # ----------------------------------------------
  # original HRF
  # ----------------------------------------------
  
  pred.hrf.orig <- predict.hedgedrf(hrf.orig.fit, data = df.test[,1:d])
  
  accuracies$hrf.orig[e] <- mean(df.test$Y == ifelse(pred.hrf.orig>0.5, 1, 0))
  waccuracies$hrf.orig[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.hrf.orig>0.5, 1, 0))
  wBCEscores$hrf.orig[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.hrf.orig)
  
  
  # ----------------------------------------------
  # stabilized classification
  # ----------------------------------------------

  pred.sc <- predict.stabClass(sc.fit, newsample = df.test[,1:d])

  accuracies$sc[e] <- mean(df.test$Y == pred.sc$pred.class)
  waccuracies$sc[e] <- weighted_accuracy(labels = df.test$Y, predictions = pred.sc$pred.class)
  wBCEscores$sc[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.sc$pred.probs)


  # ----------------------------------------------
  # stabilized classification with HRF classifier
  # ----------------------------------------------

  pred.sc.hrf <- predict.stabClass.hrf(sc.fit, newsample = df.test[,1:d])

  accuracies$sc.hrf[e] <- mean(df.test$Y == pred.sc.hrf$pred.class)
  waccuracies$sc.hrf[e] <- weighted_accuracy(labels = df.test$Y, predictions = pred.sc.hrf$pred.class)
  wBCEscores$sc.hrf[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.sc.hrf$pred.probs)


  # ----------------------------------------------
  # GHRF 1
  # ----------------------------------------------

  pred.ghrf.1 <- predict.ghrf(object = ghrf.1.fit, data = df.test[,1:d])

  accuracies$ghrf.1[e] <- mean(df.test$Y == ifelse(pred.ghrf.1>0.5, 1, 0))
  waccuracies$ghrf.1[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.ghrf.1>0.5, 1, 0))
  wBCEscores$ghrf.1[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.ghrf.1)


  # ----------------------------------------------
  # GHRF 2
  # ----------------------------------------------

  pred.ghrf.2 <- predict.ghrf(object = ghrf.2.fit, data = df.test[,1:d])

  accuracies$ghrf.2[e] <- mean(df.test$Y == ifelse(pred.ghrf.2>0.5, 1, 0))
  waccuracies$ghrf.2[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.ghrf.2>0.5, 1, 0))
  wBCEscores$ghrf.2[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.ghrf.2)

  
  
  # ----------------------------------------------
  # GHRF 3
  # ----------------------------------------------
  
  pred.ghrf.3 <- predict.ghrf(object = ghrf.3.fit, data = df.test[,1:d])

  accuracies$ghrf.3[e] <- mean(df.test$Y == ifelse(pred.ghrf.3>0.5, 1, 0))
  waccuracies$ghrf.3[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.ghrf.3>0.5, 1, 0))
  wBCEscores$ghrf.3[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.ghrf.3)
  
  
  # ----------------------------------------------
  # GHRF 4
  # ----------------------------------------------

  pred.ghrf.4 <- predict.ghrf(object = ghrf.4.fit, data = df.test[,1:d])

  accuracies$ghrf.4[e] <- mean(df.test$Y == ifelse(pred.ghrf.4>0.5, 1, 0))
  waccuracies$ghrf.4[e] <- weighted_accuracy(labels = df.test$Y, predictions = ifelse(pred.ghrf.4>0.5, 1, 0))
  wBCEscores$ghrf.4[e] <- BCE.weighted(y = df.test$Y, y.hat = pred.ghrf.4)
  
}


df <- waccuracies


df_long <- df %>%
  mutate(ID = row_number()) %>%  # Create an identifier for each row
  pivot_longer(cols = -ID, names_to = "Model", values_to = "Accuracy") %>%
  select(-ID)  # Remove ID if it's no longer needed


# Compute summary statistics
summary_df <- df_long %>%
  group_by(Model) %>%
  summarise(
    min = min(Accuracy),
    q25 = quantile(Accuracy, 0.25),
    median = median(Accuracy),
    q75 = quantile(Accuracy, 0.75),
    max = max(Accuracy)
  ) %>%
  pivot_longer(cols = -Model, names_to = "Quantile", values_to = "Value")



# Define color mapping for quantiles
quantile_levels <- c("min", "q25", "median", "q75", "max")
summary_df$Quantile <- factor(summary_df$Quantile, levels = quantile_levels)
colors <- scales::viridis_pal()(length(quantile_levels))

# Plot using ggplot2
ggplot(summary_df, aes(x = Model, y = Value, color = Quantile, group = Quantile)) +
  geom_point(size = 3) +
  geom_line(aes(group = Quantile), size = 1) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(
    title = "Comparison of Accuracy Distributions Across Models",
    y = "Weighted Accuracy",
    color = "Quantile"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave("GHRF_regression/quantiles_DRIG.png")




save(accuracies, waccuracies, wBCEscores, file = "simulation_causal_chambers_GHRF/results.rdata")






metric <- waccuracies

plot(metric$rf, type = "l", col="black")
# points(metric$oracle, type ="l", col="red")
points(metric$hrf.orig, type ="l", col="red")
# points(metric$hrf.ood, type ="l", col="blue")
points(metric$ghrf.1, type ="l", col="green")
points(metric$ghrf.2, type ="l", col="blue")
points(metric$ghrf.3, type ="l", col="orange")
points(metric$ghrf.4, type ="l", col="darkgray")



legend("bottomleft", legend=c("RF", "HRF", "GHRF 1", "GHRF 2", "GHRF 3", "GHRF 4"), col = c("black", "red", "green", "blue", "orange", "darkgray"), lty = 1)
