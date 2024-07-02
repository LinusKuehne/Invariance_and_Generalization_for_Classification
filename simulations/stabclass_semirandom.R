# In this script, we test the performance of Stabilized Classification on data from 
# the semi-random SCM


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
sets <- powerSet(1:5)
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

# for plotting
size <- 10

# max intervention strengths
strength.train <- 2
strength.test <- 8

# number of covariates (including Y)
d <- 6


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


# record whether the simulated DAG satisfies MB(Y) == SB_I(Y)
MB.equal.SB.mod1 <- rep(T, n.sim)




for(sim in 1:n.sim){
  print(paste0("Simulation iteration ",sim, " out of ", n.sim))
  
  # generate a sample of the semirandom SCM
  s <- generate.samples.semirandom(n = n, n.test = n.test, int.strength.train = strength.train, int.strength.test = strength.test, mod = "logreg")
  
  # extract generated datasets
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # extract used DAGs
  dag.full <- s$dag.full
  dag.cov <- s$dag.cov
  
  # compute stable blanket and markov blanket
  SB.Y <- stableBlanket(dag.full, dag.cov, num.int = 2, d = d)
  MB.Y <- markovBlanket(dag.cov)
  
  # are they the same?
  MB.equal.SB.mod1[sim] <- identical(SB.Y, MB.Y)
  
  # stabilized classification with test 1 (RF, GLM)
  output.sc.1.rfglm <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.1.rfglm <- predict.stabClass(output.sc.1.rfglm, newsample = sample_test[,1:(d-1)])
  accuracies.mod1$sc.1.rfglm[sim] <- mean(sample_test$Y == pred.sc.1.rfglm$pred.class)
  negwBCEscores.mod1$sc.1.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfglm$pred.probs)
  
  # stabilized classification with test 1 (RF, RF)
  output.sc.1.rfrf <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.1.rfrf <- predict.stabClass(output.sc.1.rfrf, newsample = sample_test[,1:(d-1)])
  accuracies.mod1$sc.1.rfrf[sim] <- mean(sample_test$Y == pred.sc.1.rfrf$pred.class)
  negwBCEscores.mod1$sc.1.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfrf$pred.probs)
  
  # stabilized classification with test 2 (RF, GLM)
  output.sc.2.rfglm <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.2.rfglm <- predict.stabClass(output.sc.2.rfglm, newsample = sample_test[,1:(d-1)])
  accuracies.mod1$sc.2.rfglm[sim] <- mean(sample_test$Y == pred.sc.2.rfglm$pred.class)
  negwBCEscores.mod1$sc.2.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfglm$pred.probs)
  
  # stabilized classification with test 2 (RF, RF)
  output.sc.2.rfrf <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.2.rfrf <- predict.stabClass(output.sc.2.rfrf, newsample = sample_test[,1:(d-1)])
  accuracies.mod1$sc.2.rfrf[sim] <- mean(sample_test$Y == pred.sc.2.rfrf$pred.class)
  negwBCEscores.mod1$sc.2.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfrf$pred.probs)
  
  # standard logistic regression
  output.glm <- glm(Y ~ ., data = sample[, 1:d], family = binomial(link = "logit"))
  pred.glm <- predict(output.glm, newdata = sample_test[,1:(d-1)], type = "response")
  accuracies.mod1$glm[sim] <- mean(sample_test$Y == ifelse(pred.glm>0.5, 1, 0))
  negwBCEscores.mod1$glm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.glm)
  
  # standard random forest
  output.rf <- ranger(y = as.factor(sample$Y), x = sample[, 1:(d-1)], probability = T)
  pred.rf <- predict(output.rf, data = sample_test[,1:(d-1)])$predictions[,"1"]
  accuracies.mod1$rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))
  negwBCEscores.mod1$rf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.rf)
}


# separate data depending on whether is was generated for a DAG with SB(Y) == MB(Y) or SB(Y) != MB(Y)
acc.eq.mod1 <- accuracies.mod1[MB.equal.SB.mod1, ]
acc.neq.mod1 <- accuracies.mod1[!MB.equal.SB.mod1, ]

negwBCE.eq.mod1 <- negwBCEscores.mod1[MB.equal.SB.mod1, ]
negwBCE.neq.mod1 <- negwBCEscores.mod1[!MB.equal.SB.mod1, ]

n.eq.mod1 <- sum(MB.equal.SB.mod1)
n.neq.mod1 <- n.sim-n.eq.mod1

name.eq.mod1 <- factor(c(rep("SC (a)",n.eq.mod1), rep("SC (b)",n.eq.mod1), rep("SC (c)",n.eq.mod1), rep("SC (d)",n.eq.mod1), rep("Log. Reg.",n.eq.mod1), rep("RF", n.eq.mod1)), levels = c("SC (a)", "SC (b)", "SC (c)", "SC (d)", "Log. Reg.", "RF"))
name.neq.mod1 <- factor(c(rep("SC (a)",n.neq.mod1), rep("SC (b)",n.neq.mod1), rep("SC (c)",n.neq.mod1), rep("SC (d)",n.neq.mod1), rep("Log. Reg.",n.neq.mod1), rep("RF", n.neq.mod1)), levels = c("SC (a)", "SC (b)", "SC (c)", "SC (d)", "Log. Reg.", "RF"))





# create dataframes for plotting
data.eq.acc.mod1 <- data.frame(
  name=name.eq.mod1,
  value=c(acc.eq.mod1$sc.1.rfglm, acc.eq.mod1$sc.1.rfrf, acc.eq.mod1$sc.2.rfglm, acc.eq.mod1$sc.2.rfrf, acc.eq.mod1$glm, acc.eq.mod1$rf)
)

data.neq.acc.mod1 <- data.frame(
  name=name.neq.mod1,
  value=c(acc.neq.mod1$sc.1.rfglm, acc.neq.mod1$sc.1.rfrf, acc.neq.mod1$sc.2.rfglm, acc.neq.mod1$sc.2.rfrf, acc.neq.mod1$glm, acc.neq.mod1$rf)
)


data.eq.wbce.mod1 <- data.frame(
  name=name.eq.mod1,
  value=c(negwBCE.eq.mod1$sc.1.rfglm, negwBCE.eq.mod1$sc.1.rfrf, negwBCE.eq.mod1$sc.2.rfglm, negwBCE.eq.mod1$sc.2.rfrf, negwBCE.eq.mod1$glm, negwBCE.eq.mod1$rf)
)

data.neq.wbce.mod1 <- data.frame(
  name=name.neq.mod1,
  value=c(negwBCE.neq.mod1$sc.1.rfglm, negwBCE.neq.mod1$sc.1.rfrf, negwBCE.neq.mod1$sc.2.rfglm, negwBCE.neq.mod1$sc.2.rfrf, negwBCE.neq.mod1$glm, negwBCE.neq.mod1$rf)
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

# record whether the simulated DAG satisfies MB(Y) == SB_I(Y)
MB.equal.SB.mod2 <- rep(T, n.sim)




for(sim in 1:n.sim){
  print(paste0("Simulation iteration ",sim, " out of ", n.sim))
  
  # generate a sample of the semirandom SCM
  s <- generate.samples.semirandom(n = n, n.test = n.test, int.strength.train = strength.train, int.strength.test = strength.test, mod = "probit")
  
  # extract generated datasets
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # extract used DAGs
  dag.full <- s$dag.full
  dag.cov <- s$dag.cov
  
  # compute stable blanket and markov blanket
  SB.Y <- stableBlanket(dag.full, dag.cov, num.int = 2, d = d)
  MB.Y <- markovBlanket(dag.cov)
  
  # are they the same?
  MB.equal.SB.mod2[sim] <- identical(SB.Y, MB.Y)
  
  # stabilized classification with test 1 (RF, GLM)
  output.sc.1.rfglm <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.1.rfglm <- predict.stabClass(output.sc.1.rfglm, newsample = sample_test[,1:(d-1)])
  accuracies.mod2$sc.1.rfglm[sim] <- mean(sample_test$Y == pred.sc.1.rfglm$pred.class)
  negwBCEscores.mod2$sc.1.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfglm$pred.probs)
  
  # stabilized classification with test 1 (RF, RF)
  output.sc.1.rfrf <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.1.rfrf <- predict.stabClass(output.sc.1.rfrf, newsample = sample_test[,1:(d-1)])
  accuracies.mod2$sc.1.rfrf[sim] <- mean(sample_test$Y == pred.sc.1.rfrf$pred.class)
  negwBCEscores.mod2$sc.1.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfrf$pred.probs)
  
  # stabilized classification with test 2 (RF, GLM)
  output.sc.2.rfglm <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.2.rfglm <- predict.stabClass(output.sc.2.rfglm, newsample = sample_test[,1:(d-1)])
  accuracies.mod2$sc.2.rfglm[sim] <- mean(sample_test$Y == pred.sc.2.rfglm$pred.class)
  negwBCEscores.mod2$sc.2.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfglm$pred.probs)
  
  # stabilized classification with test 2 (RF, RF)
  output.sc.2.rfrf <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.2.rfrf <- predict.stabClass(output.sc.2.rfrf, newsample = sample_test[,1:(d-1)])
  accuracies.mod2$sc.2.rfrf[sim] <- mean(sample_test$Y == pred.sc.2.rfrf$pred.class)
  negwBCEscores.mod2$sc.2.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfrf$pred.probs)
  
  # standard logistic regression
  output.glm <- glm(Y ~ ., data = sample[, 1:d], family = binomial(link = "logit"))
  pred.glm <- predict(output.glm, newdata = sample_test[,1:(d-1)], type = "response")
  accuracies.mod2$glm[sim] <- mean(sample_test$Y == ifelse(pred.glm>0.5, 1, 0))
  negwBCEscores.mod2$glm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.glm)
  
  # standard random forest
  output.rf <- ranger(y = as.factor(sample$Y), x = sample[, 1:(d-1)], probability = T)
  pred.rf <- predict(output.rf, data = sample_test[,1:(d-1)])$predictions[,"1"]
  accuracies.mod2$rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))
  negwBCEscores.mod2$rf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.rf)
}


# separate data depending on whether is was generated for a DAG with SB(Y) == MB(Y) or SB(Y) != MB(Y)
acc.eq.mod2 <- accuracies.mod2[MB.equal.SB.mod2, ]
acc.neq.mod2 <- accuracies.mod2[!MB.equal.SB.mod2, ]

negwBCE.eq.mod2 <- negwBCEscores.mod2[MB.equal.SB.mod2, ]
negwBCE.neq.mod2 <- negwBCEscores.mod2[!MB.equal.SB.mod2, ]

n.eq.mod2 <- sum(MB.equal.SB.mod2)
n.neq.mod2 <- n.sim-n.eq.mod2

name.eq.mod2 <- factor(c(rep("SC (a)",n.eq.mod2), rep("SC (b)",n.eq.mod2), rep("SC (c)",n.eq.mod2), rep("SC (d)",n.eq.mod2), rep("Log. Reg.",n.eq.mod2), rep("RF", n.eq.mod2)), levels = c("SC (a)", "SC (b)", "SC (c)", "SC (d)", "Log. Reg.", "RF"))
name.neq.mod2 <- factor(c(rep("SC (a)",n.neq.mod2), rep("SC (b)",n.neq.mod2), rep("SC (c)",n.neq.mod2), rep("SC (d)",n.neq.mod2), rep("Log. Reg.",n.neq.mod2), rep("RF", n.neq.mod2)), levels = c("SC (a)", "SC (b)", "SC (c)", "SC (d)", "Log. Reg.", "RF"))






# create dataframes for plotting
data.eq.acc.mod2 <- data.frame(
  name=name.eq.mod2,
  value=c(acc.eq.mod2$sc.1.rfglm, acc.eq.mod2$sc.1.rfrf, acc.eq.mod2$sc.2.rfglm, acc.eq.mod2$sc.2.rfrf, acc.eq.mod2$glm, acc.eq.mod2$rf)
)

data.neq.acc.mod2 <- data.frame(
  name=name.neq.mod2,
  value=c(acc.neq.mod2$sc.1.rfglm, acc.neq.mod2$sc.1.rfrf, acc.neq.mod2$sc.2.rfglm, acc.neq.mod2$sc.2.rfrf, acc.neq.mod2$glm, acc.neq.mod2$rf)
)


data.eq.wbce.mod2 <- data.frame(
  name=name.eq.mod2,
  value=c(negwBCE.eq.mod2$sc.1.rfglm, negwBCE.eq.mod2$sc.1.rfrf, negwBCE.eq.mod2$sc.2.rfglm, negwBCE.eq.mod2$sc.2.rfrf, negwBCE.eq.mod2$glm, negwBCE.eq.mod2$rf)
)

data.neq.wbce.mod2 <- data.frame(
  name=name.neq.mod2,
  value=c(negwBCE.neq.mod2$sc.1.rfglm, negwBCE.neq.mod2$sc.1.rfrf, negwBCE.neq.mod2$sc.2.rfglm, negwBCE.neq.mod2$sc.2.rfrf, negwBCE.neq.mod2$glm, negwBCE.neq.mod2$rf)
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


# record whether the simulated DAG satisfies MB(Y) == SB_I(Y)
MB.equal.SB.mod3 <- rep(T, n.sim)



for(sim in 1:n.sim){
  print(paste0("Simulation iteration ",sim, " out of ", n.sim))
  
  # generate a sample of the semirandom SCM
  s <- generate.samples.semirandom(n = n, n.test = n.test, int.strength.train = strength.train, int.strength.test = strength.test, mod = "nonlin")
  
  # extract generated datasets
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # extract used DAGs
  dag.full <- s$dag.full
  dag.cov <- s$dag.cov
  
  # compute stable blanket and markov blanket
  SB.Y <- stableBlanket(dag.full, dag.cov, num.int = 2, d = d)
  MB.Y <- markovBlanket(dag.cov)
  
  # are they the same?
  MB.equal.SB.mod3[sim] <- identical(SB.Y, MB.Y)
  
  # stabilized classification with test 1 (RF, GLM)
  output.sc.1.rfglm <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.1.rfglm <- predict.stabClass(output.sc.1.rfglm, newsample = sample_test[,1:(d-1)])
  accuracies.mod3$sc.1.rfglm[sim] <- mean(sample_test$Y == pred.sc.1.rfglm$pred.class)
  negwBCEscores.mod3$sc.1.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfglm$pred.probs)
  
  # stabilized classification with test 1 (RF, RF)
  output.sc.1.rfrf <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.1.rfrf <- predict.stabClass(output.sc.1.rfrf, newsample = sample_test[,1:(d-1)])
  accuracies.mod3$sc.1.rfrf[sim] <- mean(sample_test$Y == pred.sc.1.rfrf$pred.class)
  negwBCEscores.mod3$sc.1.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfrf$pred.probs)
  
  # stabilized classification with test 2 (RF, GLM)
  output.sc.2.rfglm <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.2.rfglm <- predict.stabClass(output.sc.2.rfglm, newsample = sample_test[,1:(d-1)])
  accuracies.mod3$sc.2.rfglm[sim] <- mean(sample_test$Y == pred.sc.2.rfglm$pred.class)
  negwBCEscores.mod3$sc.2.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfglm$pred.probs)
  
  # stabilized classification with test 2 (RF, RF)
  output.sc.2.rfrf <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.2.rfrf <- predict.stabClass(output.sc.2.rfrf, newsample = sample_test[,1:(d-1)])
  accuracies.mod3$sc.2.rfrf[sim] <- mean(sample_test$Y == pred.sc.2.rfrf$pred.class)
  negwBCEscores.mod3$sc.2.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfrf$pred.probs)
  
  # standard logistic regression
  output.glm <- glm(Y ~ ., data = sample[, 1:d], family = binomial(link = "logit"))
  pred.glm <- predict(output.glm, newdata = sample_test[,1:(d-1)], type = "response")
  accuracies.mod3$glm[sim] <- mean(sample_test$Y == ifelse(pred.glm>0.5, 1, 0))
  negwBCEscores.mod3$glm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.glm)
  
  # standard random forest
  output.rf <- ranger(y = as.factor(sample$Y), x = sample[, 1:(d-1)], probability = T)
  pred.rf <- predict(output.rf, data = sample_test[,1:(d-1)])$predictions[,"1"]
  accuracies.mod3$rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))
  negwBCEscores.mod3$rf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.rf)
}


# separate data depending on whether is was generated for a DAG with SB(Y) == MB(Y) or SB(Y) != MB(Y)
acc.eq.mod3 <- accuracies.mod3[MB.equal.SB.mod3, ]
acc.neq.mod3 <- accuracies.mod3[!MB.equal.SB.mod3, ]

negwBCE.eq.mod3 <- negwBCEscores.mod3[MB.equal.SB.mod3, ]
negwBCE.neq.mod3 <- negwBCEscores.mod3[!MB.equal.SB.mod3, ]

n.eq.mod3 <- sum(MB.equal.SB.mod3)
n.neq.mod3 <- n.sim-n.eq.mod3

name.eq.mod3 <- factor(c(rep("SC (a)",n.eq.mod3), rep("SC (b)",n.eq.mod3), rep("SC (c)",n.eq.mod3), rep("SC (d)",n.eq.mod3), rep("Log. Reg.",n.eq.mod3), rep("RF", n.eq.mod3)), levels = c("SC (a)", "SC (b)", "SC (c)", "SC (d)", "Log. Reg.", "RF"))
name.neq.mod3 <- factor(c(rep("SC (a)",n.neq.mod3), rep("SC (b)",n.neq.mod3), rep("SC (c)",n.neq.mod3), rep("SC (d)",n.neq.mod3), rep("Log. Reg.",n.neq.mod3), rep("RF", n.neq.mod3)), levels = c("SC (a)", "SC (b)", "SC (c)", "SC (d)", "Log. Reg.", "RF"))






# create dataframes for plotting
data.eq.acc.mod3 <- data.frame(
  name=name.eq.mod3,
  value=c(acc.eq.mod3$sc.1.rfglm, acc.eq.mod3$sc.1.rfrf, acc.eq.mod3$sc.2.rfglm, acc.eq.mod3$sc.2.rfrf, acc.eq.mod3$glm, acc.eq.mod3$rf)
)

data.neq.acc.mod3 <- data.frame(
  name=name.neq.mod3,
  value=c(acc.neq.mod3$sc.1.rfglm, acc.neq.mod3$sc.1.rfrf, acc.neq.mod3$sc.2.rfglm, acc.neq.mod3$sc.2.rfrf, acc.neq.mod3$glm, acc.neq.mod3$rf)
)


data.eq.wbce.mod3 <- data.frame(
  name=name.eq.mod3,
  value=c(negwBCE.eq.mod3$sc.1.rfglm, negwBCE.eq.mod3$sc.1.rfrf, negwBCE.eq.mod3$sc.2.rfglm, negwBCE.eq.mod3$sc.2.rfrf, negwBCE.eq.mod3$glm, negwBCE.eq.mod3$rf)
)

data.neq.wbce.mod3 <- data.frame(
  name=name.neq.mod3,
  value=c(negwBCE.neq.mod3$sc.1.rfglm, negwBCE.neq.mod3$sc.1.rfrf, negwBCE.neq.mod3$sc.2.rfglm, negwBCE.neq.mod3$sc.2.rfrf, negwBCE.neq.mod3$glm, negwBCE.neq.mod3$rf)
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



# record whether the simulated DAG satisfies MB(Y) == SB_I(Y)
MB.equal.SB.mod4 <- rep(T, n.sim)




for(sim in 1:n.sim){
  print(paste0("Simulation iteration ",sim, " out of ", n.sim))
  
  # generate a sample of the semirandom SCM
  s <- generate.samples.semirandom(n = n, n.test = n.test, int.strength.train = strength.train, int.strength.test = strength.test, mod = "bump")
  
  # extract generated datasets
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # extract used DAGs
  dag.full <- s$dag.full
  dag.cov <- s$dag.cov
  
  # compute stable blanket and markov blanket
  SB.Y <- stableBlanket(dag.full, dag.cov, num.int = 2, d = d)
  MB.Y <- markovBlanket(dag.cov)
  
  # are they the same?
  MB.equal.SB.mod4[sim] <- identical(SB.Y, MB.Y)
  
  # stabilized classification with test 1 (RF, GLM)
  output.sc.1.rfglm <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.1.rfglm <- predict.stabClass(output.sc.1.rfglm, newsample = sample_test[,1:(d-1)])
  accuracies.mod4$sc.1.rfglm[sim] <- mean(sample_test$Y == pred.sc.1.rfglm$pred.class)
  negwBCEscores.mod4$sc.1.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfglm$pred.probs)
  
  # stabilized classification with test 1 (RF, RF)
  output.sc.1.rfrf <- stabilizedClassification(sample = sample, test = test.1, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.1.rfrf <- predict.stabClass(output.sc.1.rfrf, newsample = sample_test[,1:(d-1)])
  accuracies.mod4$sc.1.rfrf[sim] <- mean(sample_test$Y == pred.sc.1.rfrf$pred.class)
  negwBCEscores.mod4$sc.1.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.1.rfrf$pred.probs)
  
  # stabilized classification with test 2 (RF, GLM)
  output.sc.2.rfglm <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "GLM", B = B, verbose = F)
  pred.sc.2.rfglm <- predict.stabClass(output.sc.2.rfglm, newsample = sample_test[,1:(d-1)])
  accuracies.mod4$sc.2.rfglm[sim] <- mean(sample_test$Y == pred.sc.2.rfglm$pred.class)
  negwBCEscores.mod4$sc.2.rfglm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfglm$pred.probs)
  
  # stabilized classification with test 2 (RF, RF)
  output.sc.2.rfrf <- stabilizedClassification(sample = sample, test = test.2, mod.internal = "RF", mod.output = "RF", B = B, verbose = F)
  pred.sc.2.rfrf <- predict.stabClass(output.sc.2.rfrf, newsample = sample_test[,1:(d-1)])
  accuracies.mod4$sc.2.rfrf[sim] <- mean(sample_test$Y == pred.sc.2.rfrf$pred.class)
  negwBCEscores.mod4$sc.2.rfrf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.sc.2.rfrf$pred.probs)
  
  # standard logistic regression
  output.glm <- glm(Y ~ ., data = sample[, 1:d], family = binomial(link = "logit"))
  pred.glm <- predict(output.glm, newdata = sample_test[,1:(d-1)], type = "response")
  accuracies.mod4$glm[sim] <- mean(sample_test$Y == ifelse(pred.glm>0.5, 1, 0))
  negwBCEscores.mod4$glm[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.glm)
  
  # standard random forest
  output.rf <- ranger(y = as.factor(sample$Y), x = sample[, 1:(d-1)], probability = T)
  pred.rf <- predict(output.rf, data = sample_test[,1:(d-1)])$predictions[,"1"]
  accuracies.mod4$rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))
  negwBCEscores.mod4$rf[sim] <- nBCE.weighted(y = sample_test$Y, y.hat = pred.rf)
}



# separate data depending on whether is was generated for a DAG with SB(Y) == MB(Y) or SB(Y) != MB(Y)
acc.eq.mod4 <- accuracies.mod4[MB.equal.SB.mod4, ]
acc.neq.mod4 <- accuracies.mod4[!MB.equal.SB.mod4, ]

negwBCE.eq.mod4 <- negwBCEscores.mod4[MB.equal.SB.mod4, ]
negwBCE.neq.mod4 <- negwBCEscores.mod4[!MB.equal.SB.mod4, ]

n.eq.mod4 <- sum(MB.equal.SB.mod4)
n.neq.mod4 <- n.sim-n.eq.mod4

name.eq.mod4 <- factor(c(rep("SC (a)",n.eq.mod4), rep("SC (b)",n.eq.mod4), rep("SC (c)",n.eq.mod4), rep("SC (d)",n.eq.mod4), rep("Log. Reg.",n.eq.mod4), rep("RF", n.eq.mod4)), levels = c("SC (a)", "SC (b)", "SC (c)", "SC (d)", "Log. Reg.", "RF"))
name.neq.mod4 <- factor(c(rep("SC (a)",n.neq.mod4), rep("SC (b)",n.neq.mod4), rep("SC (c)",n.neq.mod4), rep("SC (d)",n.neq.mod4), rep("Log. Reg.",n.neq.mod4), rep("RF", n.neq.mod4)), levels = c("SC (a)", "SC (b)", "SC (c)", "SC (d)", "Log. Reg.", "RF"))






# create dataframes for plotting
data.eq.acc.mod4 <- data.frame(
  name=name.eq.mod4,
  value=c(acc.eq.mod4$sc.1.rfglm, acc.eq.mod4$sc.1.rfrf, acc.eq.mod4$sc.2.rfglm, acc.eq.mod4$sc.2.rfrf, acc.eq.mod4$glm, acc.eq.mod4$rf)
)

data.neq.acc.mod4 <- data.frame(
  name=name.neq.mod4,
  value=c(acc.neq.mod4$sc.1.rfglm, acc.neq.mod4$sc.1.rfrf, acc.neq.mod4$sc.2.rfglm, acc.neq.mod4$sc.2.rfrf, acc.neq.mod4$glm, acc.neq.mod4$rf)
)


data.eq.wbce.mod4 <- data.frame(
  name=name.eq.mod4,
  value=c(negwBCE.eq.mod4$sc.1.rfglm, negwBCE.eq.mod4$sc.1.rfrf, negwBCE.eq.mod4$sc.2.rfglm, negwBCE.eq.mod4$sc.2.rfrf, negwBCE.eq.mod4$glm, negwBCE.eq.mod4$rf)
)

data.neq.wbce.mod4 <- data.frame(
  name=name.neq.mod4,
  value=c(negwBCE.neq.mod4$sc.1.rfglm, negwBCE.neq.mod4$sc.1.rfrf, negwBCE.neq.mod4$sc.2.rfglm, negwBCE.neq.mod4$sc.2.rfrf, negwBCE.neq.mod4$glm, negwBCE.neq.mod4$rf)
)





#-------------------------------------------------------------------------------


save(accuracies.mod1,
     accuracies.mod2,
     accuracies.mod3,
     accuracies.mod4,
     negwBCEscores.mod1,
     negwBCEscores.mod2,
     negwBCEscores.mod3,
     negwBCEscores.mod4,
     MB.equal.SB.mod1,
     MB.equal.SB.mod2,
     MB.equal.SB.mod3,
     MB.equal.SB.mod4,
     file = file.path(script_dir, "saved_data/stabclass_semirandom.rdata"))










#-------------------------------------------------------------------------------
# generate plots for all models
#-------------------------------------------------------------------------------


# plotting window
min.wbce <- -0.1 + min(min(data.eq.wbce.mod1$value), min(data.eq.wbce.mod2$value), min(data.eq.wbce.mod3$value), min(data.eq.wbce.mod4$value), min(data.neq.wbce.mod1$value), min(data.neq.wbce.mod2$value), min(data.neq.wbce.mod3$value), min(data.neq.wbce.mod4$value))
max.wbce <- 0.1 + max(max(data.eq.wbce.mod1$value), max(data.eq.wbce.mod2$value), max(data.eq.wbce.mod3$value), max(data.eq.wbce.mod4$value), max(data.neq.wbce.mod1$value), max(data.neq.wbce.mod2$value), max(data.neq.wbce.mod3$value), max(data.neq.wbce.mod4$value))



# logistic regression ----------------------------------------------------------

plt.eq.acc.mod1 <- ggplot(data.eq.acc.mod1, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Logistic Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" == "SB(Y)" * " (" * v * " repetitions)"), list(v = n.eq.mod1)))
  ) +
  xlab("") +
  ylim(0,1)+
  ylab("Accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.eq.acc.mod1



plt.neq.acc.mod1 <- ggplot(data.neq.acc.mod1, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Logistic Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" != "SB(Y)" * " (" * v * " repetitions)"), list(v = n.neq.mod1)))
  ) +  xlab("") +
  ylim(0,1)+
  ylab("Accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.neq.acc.mod1




plt.eq.wbce.mod1 <- ggplot(data.eq.wbce.mod1, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Logistic Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" == "SB(Y)" * " (" * v * " repetitions)"), list(v = n.eq.mod1)))
  ) +
  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.eq.wbce.mod1



plt.neq.wbce.mod1 <- ggplot(data.neq.wbce.mod1, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Logistic Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" != "SB(Y)" * " (" * v * " repetitions)"), list(v = n.neq.mod1)))
  ) +  
  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.neq.wbce.mod1







# probit regression ------------------------------------------------------------

plt.eq.acc.mod2 <- ggplot(data.eq.acc.mod2, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Probit Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" == "SB(Y)" * " (" * v * " repetitions)"), list(v = n.eq.mod2)))
  ) +
  xlab("") +
  ylim(0,1)+
  ylab("Accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.eq.acc.mod2



plt.neq.acc.mod2 <- ggplot(data.neq.acc.mod2, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Probit Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" != "SB(Y)" * " (" * v * " repetitions)"), list(v = n.neq.mod2)))
  ) +  xlab("") +
  ylim(0,1)+
  ylab("Accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.neq.acc.mod2




plt.eq.wbce.mod2 <- ggplot(data.eq.wbce.mod2, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Probit Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" == "SB(Y)" * " (" * v * " repetitions)"), list(v = n.eq.mod2)))
  ) +
  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.eq.wbce.mod2



plt.neq.wbce.mod2 <- ggplot(data.neq.wbce.mod2, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Probit Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" != "SB(Y)" * " (" * v * " repetitions)"), list(v = n.neq.mod2)))
  ) +  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.neq.wbce.mod2






# non-linear logistic regression -----------------------------------------------

plt.eq.acc.mod3 <- ggplot(data.eq.acc.mod3, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Non-Linear Logistic Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" == "SB(Y)" * " (" * v * " repetitions)"), list(v = n.eq.mod3)))
  ) +
  xlab("") +
  ylim(0,1)+
  ylab("Accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.eq.acc.mod3



plt.neq.acc.mod3 <- ggplot(data.neq.acc.mod3, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Non-Linear Logistic Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" != "SB(Y)" * " (" * v * " repetitions)"), list(v = n.neq.mod3)))
  ) +  xlab("") +
  ylim(0,1)+
  ylab("Accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.neq.acc.mod3




plt.eq.wbce.mod3 <- ggplot(data.eq.wbce.mod3, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Non-Linear Logistic Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" == "SB(Y)" * " (" * v * " repetitions)"), list(v = n.eq.mod3)))
  ) +
  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.eq.wbce.mod3



plt.neq.wbce.mod3 <- ggplot(data.neq.wbce.mod3, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Non-Linear Logistic Regression Model"),
    subtitle = eval(substitute(expression("MB(Y)" != "SB(Y)" * " (" * v * " repetitions)"), list(v = n.neq.mod3)))
  ) +  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.neq.wbce.mod3





# bump model -------------------------------------------------------------------

plt.eq.acc.mod4 <- ggplot(data.eq.acc.mod4, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Bump Model"),
    subtitle = eval(substitute(expression("MB(Y)" == "SB(Y)" * " (" * v * " repetitions)"), list(v = n.eq.mod4)))
  ) +
  xlab("") +
  ylim(0,1)+
  ylab("Accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.eq.acc.mod4



plt.neq.acc.mod4 <- ggplot(data.neq.acc.mod4, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Bump Model"),
    subtitle = eval(substitute(expression("MB(Y)" != "SB(Y)" * " (" * v * " repetitions)"), list(v = n.neq.mod4)))
  ) +  xlab("") +
  ylim(0,1)+
  ylab("Accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.neq.acc.mod4




plt.eq.wbce.mod4 <- ggplot(data.eq.wbce.mod4, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Bump Model"),
    subtitle = eval(substitute(expression("MB(Y)" == "SB(Y)" * " (" * v * " repetitions)"), list(v = n.eq.mod4)))
  ) +
  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.eq.wbce.mod4



plt.neq.wbce.mod4 <- ggplot(data.neq.wbce.mod4, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  labs(
    title = expression("Bump Model"),
    subtitle = eval(substitute(expression("MB(Y)" != "SB(Y)" * " (" * v * " repetitions)"), list(v = n.neq.mod4)))
  ) +  xlab("") +
  ylim(min.wbce,max.wbce)+
  ylab("Negative weighted BCE") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
#plt.neq.wbce.mod4














# combine plots ----------------------------------------------------------------


combined.acc <- plt.eq.acc.mod1 + plt.neq.acc.mod1 + plt.eq.acc.mod2 + plt.neq.acc.mod2 + plt.eq.acc.mod3 + plt.neq.acc.mod3 + plt.eq.acc.mod4 + plt.neq.acc.mod4 + plot_layout(ncol = 2, nrow = 4)
combined.acc
ggsave(filename = file.path(script_dir, "saved_plots/stabclass_semirandom_acc.pdf"), width = 8, height = 10)


combined.wbce <- plt.eq.wbce.mod1 + plt.neq.wbce.mod1 + plt.eq.wbce.mod2 + plt.neq.wbce.mod2 + plt.eq.wbce.mod3 + plt.neq.wbce.mod3 + plt.eq.wbce.mod4 + plt.neq.wbce.mod4 + plot_layout(ncol = 2, nrow = 4)
combined.wbce
ggsave(filename = file.path(script_dir, "saved_plots/stabclass_semirandom_wbce.pdf"), width = 8, height = 10)








#-------------------------------------------------------------------------------

# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/stabclass_semirandom.txt"))



