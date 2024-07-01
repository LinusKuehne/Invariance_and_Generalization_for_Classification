# In this script, we test the performance of Stabilized Classification on data from 
# the standard SCM


library(ggplot2)
library(patchwork)
library(rje)
library(viridis)


# get the path of this script
script_dir <- getwd()

# load in the functions needed
source(file.path(script_dir, "../code/code_simulations/invariance_tests.R"))
source(file.path(script_dir, "../code/code_simulations/data_generating_process.R"))
source(file.path(script_dir, "../code/code_simulations/utils.R"))
source(file.path(script_dir, "../code/code_simulations/stabilized_classification.R"))







#-------------------------------------------------------------------------------
# Parameters for the simulation
#-------------------------------------------------------------------------------




# sets to check stability
sets <- powerSet(1:3)
sets[[1]] <- c(0)

# number of samples per environment (training: n, testing: n.test)
n.test <- 250
n <- 250

# number of simulation runs 
# n.sim = 500 should take around 1/2 day on my MB
n.sim <- 500

# number of bootstrap samples to compute c.pred
B <- 100

# the tests to be used
test.1 = "tram.glm"
test.2 = "residual"


name <- factor(c(rep("SC (a)",n.sim), rep("SC (b)",n.sim), rep("SC (c)",n.sim), rep("SC (d)",n.sim), rep("Log. Reg.", n.sim), rep("RF", n.sim)), levels = c("SC (a)", "SC (b)", "SC (c)", "SC (d)", "Log. Reg.", "RF"))


#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# run the simulation for the logistic regression model
#-------------------------------------------------------------------------------


# set the seed
set.seed(1)

# accuracies of the models
accuracies.mod1 <- data.frame(sc.1.rfglm = numeric(n.sim),
                         sc.1.rfrf = numeric(n.sim),
                         sc.2.rfglm = numeric(n.sim),
                         sc.2.rfrf = numeric(n.sim),
                         glm = numeric(n.sim),
                         rf = numeric(n.sim)
                         )


negwBCEscores.mod1 <- accuracies.mod1



for(sim in 1:n.sim){
  print(paste0("Simulation iteration ",sim, " out of ", n.sim))
  
  # generate a sample of the standard SCM
  s <- gen.sample.fixed(n = n, n.test = n.test, int.strength.train = 1/2, int.strength.test = 2.5, mod = "logreg")
  
  # extract generated datasets
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # stabilized classification with test 1 (RF, GLM)
  output.sc.1.rfglm <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.1.rfglm <- predict.stabClass(output.sc.1.rfglm, newsample = sample_test[,1:3])
  accuracies.mod1$sc.1.rfglm[sim] <- mean(sample_test$Y == pred.sc.1.rfglm$pred.class)
  negwBCEscores.mod1$sc.1.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfglm$pred.probs)
  
  # stabilized classification with test 1 (RF, RF)
  output.sc.1.rfrf <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.1.rfrf <- predict.stabClass(output.sc.1.rfrf, newsample = sample_test[,1:3])
  accuracies.mod1$sc.1.rfrf[sim] <- mean(sample_test$Y == pred.sc.1.rfrf$pred.class)
  negwBCEscores.mod1$sc.1.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfrf$pred.probs)
  
  # stabilized classification with test 2 (RF, GLM)
  output.sc.2.rfglm <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.2.rfglm <- predict.stabClass(output.sc.2.rfglm, newsample = sample_test[,1:3])
  accuracies.mod1$sc.2.rfglm[sim] <- mean(sample_test$Y == pred.sc.2.rfglm$pred.class)
  negwBCEscores.mod1$sc.2.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfglm$pred.probs)
  
  # stabilized classification with test 2 (RF, RF)
  output.sc.2.rfrf <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.2.rfrf <- predict.stabClass(output.sc.2.rfrf, newsample = sample_test[,1:3])
  accuracies.mod1$sc.2.rfrf[sim] <- mean(sample_test$Y == pred.sc.2.rfrf$pred.class)
  negwBCEscores.mod1$sc.2.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfrf$pred.probs)
  
  # standard logistic regression
  output.glm <- glm(Y ~ ., data = sample[, 1:4], family = binomial(link = "logit"))
  pred.glm <- predict(output.glm, newdata = sample_test[,1:3], type = "response")
  accuracies.mod1$glm[sim] <- mean(sample_test$Y == ifelse(pred.glm>0.5, 1, 0))
  negwBCEscores.mod1$glm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.glm)
  
  # standard random forest
  output.rf <- ranger(y = as.factor(sample$Y), x = sample[, 1:3], probability = T)
  pred.rf <- predict(output.rf, data = sample_test[,1:3])$predictions[,"1"]
  accuracies.mod1$rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))
  negwBCEscores.mod1$rf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.rf)
}





# create a dataframe
data.acc.mod1 <- data.frame(
  name=name,
  value=c(accuracies.mod1$sc.1.rfglm, accuracies.mod1$sc.1.rfrf, accuracies.mod1$sc.2.rfglm, accuracies.mod1$sc.2.rfrf, accuracies.mod1$glm, accuracies.mod1$rf)
)

data.wbce.mod1 <- data.frame(
  name=name,
  value=c(negwBCEscores.mod1$sc.1.rfglm, negwBCEscores.mod1$sc.1.rfrf, negwBCEscores.mod1$sc.2.rfglm, negwBCEscores.mod1$sc.2.rfrf, negwBCEscores.mod1$glm, negwBCEscores.mod1$rf)
)



#-------------------------------------------------------------------------------













#-------------------------------------------------------------------------------
# run the simulation for the probit regression model
#-------------------------------------------------------------------------------


# set the seed
set.seed(1)

# accuracies of the models
accuracies.mod2 <- data.frame(sc.1.rfglm = numeric(n.sim),
                              sc.1.rfrf = numeric(n.sim),
                              sc.2.rfglm = numeric(n.sim),
                              sc.2.rfrf = numeric(n.sim),
                              glm = numeric(n.sim),
                              rf = numeric(n.sim)
)


negwBCEscores.mod2 <- accuracies.mod2



for(sim in 1:n.sim){
  print(paste0("Simulation iteration ",sim, " out of ", n.sim))
  
  # generate a sample of the standard SCM
  s <- gen.sample.fixed(n = n, n.test = n.test, int.strength.train = 1/2, int.strength.test = 4, mod = "probit")
  
  # extract generated datasets
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # stabilized classification with test 1 (RF, GLM)
  output.sc.1.rfglm <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.1.rfglm <- predict.stabClass(output.sc.1.rfglm, newsample = sample_test[,1:3])
  accuracies.mod2$sc.1.rfglm[sim] <- mean(sample_test$Y == pred.sc.1.rfglm$pred.class)
  negwBCEscores.mod2$sc.1.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfglm$pred.probs)
  
  # stabilized classification with test 1 (RF, RF)
  output.sc.1.rfrf <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.1.rfrf <- predict.stabClass(output.sc.1.rfrf, newsample = sample_test[,1:3])
  accuracies.mod2$sc.1.rfrf[sim] <- mean(sample_test$Y == pred.sc.1.rfrf$pred.class)
  negwBCEscores.mod2$sc.1.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfrf$pred.probs)
  
  # stabilized classification with test 2 (RF, GLM)
  output.sc.2.rfglm <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.2.rfglm <- predict.stabClass(output.sc.2.rfglm, newsample = sample_test[,1:3])
  accuracies.mod2$sc.2.rfglm[sim] <- mean(sample_test$Y == pred.sc.2.rfglm$pred.class)
  negwBCEscores.mod2$sc.2.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfglm$pred.probs)
  
  # stabilized classification with test 2 (RF, RF)
  output.sc.2.rfrf <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.2.rfrf <- predict.stabClass(output.sc.2.rfrf, newsample = sample_test[,1:3])
  accuracies.mod2$sc.2.rfrf[sim] <- mean(sample_test$Y == pred.sc.2.rfrf$pred.class)
  negwBCEscores.mod2$sc.2.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfrf$pred.probs)
  
  # standard logistic regression
  output.glm <- glm(Y ~ ., data = sample[, 1:4], family = binomial(link = "logit"))
  pred.glm <- predict(output.glm, newdata = sample_test[,1:3], type = "response")
  accuracies.mod2$glm[sim] <- mean(sample_test$Y == ifelse(pred.glm>0.5, 1, 0))
  negwBCEscores.mod2$glm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.glm)
  
  # standard random forest
  output.rf <- ranger(y = as.factor(sample$Y), x = sample[, 1:3], probability = T)
  pred.rf <- predict(output.rf, data = sample_test[,1:3])$predictions[,"1"]
  accuracies.mod2$rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))
  negwBCEscores.mod2$rf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.rf)
}





# create a dataframe
data.acc.mod2 <- data.frame(
  name=name,
  value=c(accuracies.mod2$sc.1.rfglm, accuracies.mod2$sc.1.rfrf, accuracies.mod2$sc.2.rfglm, accuracies.mod2$sc.2.rfrf, accuracies.mod2$glm, accuracies.mod2$rf)
)

data.wbce.mod2 <- data.frame(
  name=name,
  value=c(negwBCEscores.mod2$sc.1.rfglm, negwBCEscores.mod2$sc.1.rfrf, negwBCEscores.mod2$sc.2.rfglm, negwBCEscores.mod2$sc.2.rfrf, negwBCEscores.mod2$glm, negwBCEscores.mod2$rf)
)



#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# run the simulation for the non-linear logistic regression model
#-------------------------------------------------------------------------------


# set the seed
set.seed(1)

# accuracies of the models
accuracies.mod3 <- data.frame(sc.1.rfglm = numeric(n.sim),
                              sc.1.rfrf = numeric(n.sim),
                              sc.2.rfglm = numeric(n.sim),
                              sc.2.rfrf = numeric(n.sim),
                              glm = numeric(n.sim),
                              rf = numeric(n.sim)
)


negwBCEscores.mod3 <- accuracies.mod3



for(sim in 1:n.sim){
  print(paste0("Simulation iteration ",sim, " out of ", n.sim))
  
  # generate a sample of the standard SCM
  s <- gen.sample.fixed(n = n, n.test = n.test, int.strength.train = 1/2, int.strength.test = 4, mod = "nonlin")
  
  # extract generated datasets
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # stabilized classification with test 1 (RF, GLM)
  output.sc.1.rfglm <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.1.rfglm <- predict.stabClass(output.sc.1.rfglm, newsample = sample_test[,1:3])
  accuracies.mod3$sc.1.rfglm[sim] <- mean(sample_test$Y == pred.sc.1.rfglm$pred.class)
  negwBCEscores.mod3$sc.1.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfglm$pred.probs)
  
  # stabilized classification with test 1 (RF, RF)
  output.sc.1.rfrf <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.1.rfrf <- predict.stabClass(output.sc.1.rfrf, newsample = sample_test[,1:3])
  accuracies.mod3$sc.1.rfrf[sim] <- mean(sample_test$Y == pred.sc.1.rfrf$pred.class)
  negwBCEscores.mod3$sc.1.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfrf$pred.probs)
  
  # stabilized classification with test 2 (RF, GLM)
  output.sc.2.rfglm <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.2.rfglm <- predict.stabClass(output.sc.2.rfglm, newsample = sample_test[,1:3])
  accuracies.mod3$sc.2.rfglm[sim] <- mean(sample_test$Y == pred.sc.2.rfglm$pred.class)
  negwBCEscores.mod3$sc.2.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfglm$pred.probs)
  
  # stabilized classification with test 2 (RF, RF)
  output.sc.2.rfrf <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.2.rfrf <- predict.stabClass(output.sc.2.rfrf, newsample = sample_test[,1:3])
  accuracies.mod3$sc.2.rfrf[sim] <- mean(sample_test$Y == pred.sc.2.rfrf$pred.class)
  negwBCEscores.mod3$sc.2.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfrf$pred.probs)
  
  # standard logistic regression
  output.glm <- glm(Y ~ ., data = sample[, 1:4], family = binomial(link = "logit"))
  pred.glm <- predict(output.glm, newdata = sample_test[,1:3], type = "response")
  accuracies.mod3$glm[sim] <- mean(sample_test$Y == ifelse(pred.glm>0.5, 1, 0))
  negwBCEscores.mod3$glm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.glm)
  
  # standard random forest
  output.rf <- ranger(y = as.factor(sample$Y), x = sample[, 1:3], probability = T)
  pred.rf <- predict(output.rf, data = sample_test[,1:3])$predictions[,"1"]
  accuracies.mod3$rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))
  negwBCEscores.mod3$rf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.rf)
}





# create a dataframe
data.acc.mod3 <- data.frame(
  name=name,
  value=c(accuracies.mod3$sc.1.rfglm, accuracies.mod3$sc.1.rfrf, accuracies.mod3$sc.2.rfglm, accuracies.mod3$sc.2.rfrf, accuracies.mod3$glm, accuracies.mod3$rf)
)

data.wbce.mod3 <- data.frame(
  name=name,
  value=c(negwBCEscores.mod3$sc.1.rfglm, negwBCEscores.mod3$sc.1.rfrf, negwBCEscores.mod3$sc.2.rfglm, negwBCEscores.mod3$sc.2.rfrf, negwBCEscores.mod3$glm, negwBCEscores.mod3$rf)
)




#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# run the simulation for the bump model
#-------------------------------------------------------------------------------


# set the seed
set.seed(1)

# accuracies of the models
accuracies.mod4 <- data.frame(sc.1.rfglm = numeric(n.sim),
                              sc.1.rfrf = numeric(n.sim),
                              sc.2.rfglm = numeric(n.sim),
                              sc.2.rfrf = numeric(n.sim),
                              glm = numeric(n.sim),
                              rf = numeric(n.sim)
)


negwBCEscores.mod4 <- accuracies.mod4



for(sim in 1:n.sim){
  print(paste0("Simulation iteration ",sim, " out of ", n.sim))
  
  # generate a sample of the standard SCM
  s <- gen.sample.fixed(n = n, n.test = n.test, int.strength.train = 1/2, int.strength.test = 4, mod = "bump")
  
  # extract generated datasets
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # stabilized classification with test 1 (RF, GLM)
  output.sc.1.rfglm <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.1.rfglm <- predict.stabClass(output.sc.1.rfglm, newsample = sample_test[,1:3])
  accuracies.mod4$sc.1.rfglm[sim] <- mean(sample_test$Y == pred.sc.1.rfglm$pred.class)
  negwBCEscores.mod4$sc.1.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfglm$pred.probs)
  
  # stabilized classification with test 1 (RF, RF)
  output.sc.1.rfrf <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.1.rfrf <- predict.stabClass(output.sc.1.rfrf, newsample = sample_test[,1:3])
  accuracies.mod4$sc.1.rfrf[sim] <- mean(sample_test$Y == pred.sc.1.rfrf$pred.class)
  negwBCEscores.mod4$sc.1.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfrf$pred.probs)
  
  # stabilized classification with test 2 (RF, GLM)
  output.sc.2.rfglm <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.2.rfglm <- predict.stabClass(output.sc.2.rfglm, newsample = sample_test[,1:3])
  accuracies.mod4$sc.2.rfglm[sim] <- mean(sample_test$Y == pred.sc.2.rfglm$pred.class)
  negwBCEscores.mod4$sc.2.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfglm$pred.probs)
  
  # stabilized classification with test 2 (RF, RF)
  output.sc.2.rfrf <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.2.rfrf <- predict.stabClass(output.sc.2.rfrf, newsample = sample_test[,1:3])
  accuracies.mod4$sc.2.rfrf[sim] <- mean(sample_test$Y == pred.sc.2.rfrf$pred.class)
  negwBCEscores.mod4$sc.2.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfrf$pred.probs)
  
  # standard logistic regression
  output.glm <- glm(Y ~ ., data = sample[, 1:4], family = binomial(link = "logit"))
  pred.glm <- predict(output.glm, newdata = sample_test[,1:3], type = "response")
  accuracies.mod4$glm[sim] <- mean(sample_test$Y == ifelse(pred.glm>0.5, 1, 0))
  negwBCEscores.mod4$glm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.glm)
  
  # standard random forest
  output.rf <- ranger(y = as.factor(sample$Y), x = sample[, 1:3], probability = T)
  pred.rf <- predict(output.rf, data = sample_test[,1:3])$predictions[,"1"]
  accuracies.mod4$rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))
  negwBCEscores.mod4$rf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.rf)
}





# create a dataframe
data.acc.mod4 <- data.frame(
  name=name,
  value=c(accuracies.mod4$sc.1.rfglm, accuracies.mod4$sc.1.rfrf, accuracies.mod4$sc.2.rfglm, accuracies.mod4$sc.2.rfrf, accuracies.mod4$glm, accuracies.mod4$rf)
)

data.wbce.mod4 <- data.frame(
  name=name,
  value=c(negwBCEscores.mod4$sc.1.rfglm, negwBCEscores.mod4$sc.1.rfrf, negwBCEscores.mod4$sc.2.rfglm, negwBCEscores.mod4$sc.2.rfrf, negwBCEscores.mod4$glm, negwBCEscores.mod4$rf)
)




#-------------------------------------------------------------------------------





save(data.acc.mod1, 
     data.acc.mod2, 
     data.acc.mod3, 
     data.acc.mod4, 
     data.wbce.mod1,
     data.wbce.mod2,
     data.wbce.mod3,
     data.wbce.mod4,
     file = file.path(script_dir, "saved_data/stabclass_standard.rdata"))







#-------------------------------------------------------------------------------
# generate plots for all models
#-------------------------------------------------------------------------------

min.wbce <- min(min(data.wbce.mod1$value), min(data.wbce.mod2$value), min(data.wbce.mod3$value), min(data.wbce.mod4$value))
max.wbce <- max(max(data.wbce.mod1$value), max(data.wbce.mod2$value), max(data.wbce.mod3$value), max(data.wbce.mod4$value))


size <- 10


# logistic regression

plt.acc.mod1 <- ggplot(data.acc.mod1, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  ggtitle("Logistic Regression Model") +
  xlab("") +
  ylim(0,1)+
  ylab("Test accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.acc.mod1


plt.wbce.mod1 <- ggplot(data.wbce.mod1, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  ggtitle("Logistic Regression Model") +
  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Test negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.wbce.mod1






# probit regression

plt.acc.mod2 <- ggplot(data.acc.mod2, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  ggtitle("Probit Regression Model") +
  xlab("") +
  ylim(0,1)+
  ylab("Test accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
plt.acc.mod2


plt.wbce.mod2 <- ggplot(data.wbce.mod2, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  ggtitle("Probit Regression Model") +
  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Test negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
plt.wbce.mod2






# non-linear logistic regression

plt.acc.mod3 <- ggplot(data.acc.mod3, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  ggtitle("Non-Linear Logistic Regression Model") +
  xlab("") +
  ylim(0,1)+
  ylab("Test accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
plt.acc.mod3


plt.wbce.mod3 <- ggplot(data.wbce.mod3, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  ggtitle("Non-Linear Logistic Regression Model") +
  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Test negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
plt.wbce.mod3







# Bump model

plt.acc.mod4 <- ggplot(data.acc.mod4, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  ggtitle("Bump Model") +
  xlab("") +
  ylim(0,1)+
  ylab("Test accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
plt.acc.mod4


plt.wbce.mod4 <- ggplot(data.wbce.mod4, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  ggtitle("Bump Model") +
  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Test negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
plt.wbce.mod4








# Combine the plots 
combined.acc <- plt.acc.mod1 + plt.acc.mod2 + plt.acc.mod3 + plt.acc.mod4 
combined.acc

ggsave(filename = file.path(script_dir, "saved_plots/stabclass_standard_acc.pdf"), width = 7, height = 7.5)




combined.wbce <- plt.wbce.mod1 + plt.wbce.mod2 + plt.wbce.mod3 + plt.wbce.mod4 
combined.wbce
ggsave(filename = file.path(script_dir, "saved_plots/stabclass_standard_wbce.pdf"), width = 7, height = 7.5)








#-------------------------------------------------------------------------------

# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/stabclass_standard.txt"))



