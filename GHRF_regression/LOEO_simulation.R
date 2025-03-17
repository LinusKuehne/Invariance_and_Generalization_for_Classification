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

df_inter_obs <- df_inter_obs[df_inter_obs$Env != "Y",]
df_inter_obs$Env <- droplevels(df_inter_obs$Env)

df_total <- rbind(df_obs, df_inter_obs)


boxplot(Y ~ Env, data = df_total)


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

# number of environments
n.env <- length(levels(df_total$Env))

n.sim <- 20


#-------------------------------------------------------------------------------
# Run the experiment and evaluate on data from held-out environments/domains
#-------------------------------------------------------------------------------


# Initialize empty lists to store results per simulation
mse_list <- list()
mae_list <- list()


for(sim in 1:n.sim){
  print(sim)
  
  idx <- sample(1:nrow(df_obs), size = 750, replace = FALSE)
  df <- rbind(df_obs[idx,], df_inter_obs)

  
  mse.sim <- mae.sim <- data.frame(
    sim = sim,  # Store the simulation number
    env = 1:n.env,  # Store the environment index
    rf = numeric(n.env),
    hrf.orig = numeric(n.env),
    ghrf.1 = numeric(n.env),
    ghrf.2 = numeric(n.env),
    ghrf.3 = numeric(n.env),
    ghrf.4 = numeric(n.env)
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
    
    rf.fit <- ranger(y = df.train$Y, x = df.train[, 1:d], num.threads = 0)
    pred.rf <- predict(rf.fit, data = df.test[,1:d])$predictions
    
    mse.sim$rf[e] <- mean((df.test$Y - pred.rf)^2)
    mae.sim$rf[e] <- mean(abs(df.test$Y - pred.rf))
    
    
    # ----------------------------------------------
    # original HRF
    # ----------------------------------------------
    
    hrf.orig.fit <- hrf.orig(y = df.train$Y, x = df.train[, 1:d], rf.fit = rf.fit)
    pred.hrf.orig <- predict.hedgedrf(hrf.orig.fit, data = df.test[,1:d])

    mse.sim$hrf.orig[e] <- mean((df.test$Y - pred.hrf.orig)^2)
    mae.sim$hrf.orig[e] <- mean(abs(df.test$Y - pred.hrf.orig))
    

    # ----------------------------------------------
    # GHRF 1
    # ----------------------------------------------

    ghrf.1.fit <- ghrf.1(y = df.train$Y, x = df.train[, c(1:d,env.col.idx)])
    pred.ghrf.1 <- predict.ghrf(object = ghrf.1.fit, data = df.test[,1:d])

    mse.sim$ghrf.1[e] <- mean((df.test$Y - pred.ghrf.1)^2)
    mae.sim$ghrf.1[e] <- mean(abs(df.test$Y - pred.ghrf.1))


    # ----------------------------------------------
    # GHRF 2
    # ----------------------------------------------

    ghrf.2.fit <- ghrf.2(y = df.train$Y, x = df.train[, c(1:d,env.col.idx)])
    pred.ghrf.2 <- predict.ghrf(object = ghrf.2.fit, data = df.test[,1:d])

    mse.sim$ghrf.2[e] <- mean((df.test$Y - pred.ghrf.2)^2)
    mae.sim$ghrf.2[e] <- mean(abs(df.test$Y - pred.ghrf.2))


    # ----------------------------------------------
    # GHRF 3
    # ----------------------------------------------

    ghrf.3.fit <- ghrf.3(y = df.train$Y, x = df.train[, c(1:d,env.col.idx)])
    pred.ghrf.3 <- predict.ghrf(object = ghrf.3.fit, data = df.test[,1:d])

    mse.sim$ghrf.3[e] <- mean((df.test$Y - pred.ghrf.3)^2)
    mae.sim$ghrf.3[e] <- mean(abs(df.test$Y - pred.ghrf.3))


    # ----------------------------------------------
    # GHRF 4
    # ----------------------------------------------

    ghrf.4.fit <- ghrf.4(y = df.train$Y, x = df.train[, c(1:d,env.col.idx)])
    pred.ghrf.4 <- predict.ghrf(object = ghrf.4.fit, data = df.test[,1:d])

    mse.sim$ghrf.4[e] <- mean((df.test$Y - pred.ghrf.4)^2)
    mae.sim$ghrf.4[e] <- mean(abs(df.test$Y - pred.ghrf.4))
    
  }
  
  # Store simulation results
  mse_list[[sim]] <- mse.sim
  mae_list[[sim]] <- mae.sim
}

# Combine results into long format
mae_df <- do.call(rbind, mae_list)
mse_df <- do.call(rbind, mse_list)


save(mae_df, mse_df, file = "GHRF_regression/LOEO_classification_results.rdata")




#-------------------------------------------------------------------------------
# analyze the results
#-------------------------------------------------------------------------------

mae_avg <- mae_df %>%
  group_by(env) %>%
  summarise(across(-sim, mean))

mse_avg <- mse_df %>%
  group_by(env) %>%
  summarise(across(-sim, mean))


# Convert to long format
mse_long <- pivot_longer(mse_df, cols = -c(sim, env), names_to = "method", values_to = "mse")


# Create boxplot
ggplot(mse_long, aes(x = method, y = mse, fill = method)) +
  geom_boxplot() +
  facet_wrap(~env) +
  theme_bw() +
  labs(title = "MSE Distribution Across Simulations", x = "Method", y = "MSE") +
  theme(axis.text.x = element_text(hjust = 1, angle=45))

ggsave("GHRF_regression/boxplot_LOEO_classification.png")

# create bar plot
ggplot(mse_long, aes(x = factor(env), y = mse, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "MSE Across Simulations",
       x = "Environment",
       y = "MSE",
       fill = "Method") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1))  # Rotate x-axis labels if needed

ggsave("GHRF_regression/barplot_LOEO_classification.png")

