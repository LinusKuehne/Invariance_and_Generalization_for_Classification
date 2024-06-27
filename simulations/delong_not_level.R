# In this script, we show that DeLong (RF) where E is included/excluded
# (not just permuted) is NOT level, i.e. P_null(reject) > 0.05. We can 
# show this by running the test for a subset for which cond. independence should be
# satisfied, computing the pvalues for increasing sample size, and then computing 
# the rejection rate (sample size increases incrementally).
# We compare this to permutation of E (instead of exclusion) since this is 
# supposedly better for small S (which we have here) according to page 32 
# of (Heinze-Deml, Peters, Meinshausen, 2017) "Invariant Causal Prediction for 
# Nonlinear Models" (arXiv:1706.08576v2).

# We see from the plot that limsup_n sup_null P_null (reject) > 0.05. Hence, 
# the test doesn't have uniform asymptotic level.


# However, this method works much better if we permute E in RF.noEnv instead of 
# leaving E away completely. This was suggested on p. 32 of (Heinze-Deml, Peters, Meinshausen, 2017).
# Running the "level-test" below for set <- c(1) and set <- c(1,3) also shows that 
# the test is **approximately** uniformly asymptotically level. 



# get the path of this script
script_dir <- getwd()


library(ranger)
library(pROC)
library(ggplot2)
library(patchwork)







#-------------------------------------------------------------------------------
# data generating process
#-------------------------------------------------------------------------------


# returns a sample from standard SCM of one environment
# n: sample size
# c: specific value of E / intervention
# env: name of the environment (string)
sim.SCM1 <- function(n, c, env){

  # parameters in SCM
  alpha <- -1 
  beta <- 1
  gamma <- 1.5
  
  # the intervention
  E <- rep(c,n)
  
  # we shouldn't use the value E directly, but we can know to which experiment data belongs
  Env <- as.factor(rep(env, n))
  
  X1 <- alpha*E + rnorm(n)
  
  # logistic regression model for Y
  Ny <- rlogis(n)
  Y <- ifelse(Ny < beta*X1, 1, 0)
  
  X2 <- Y + gamma*E + rnorm(n)
  
  X3 <- Y + rnorm(n)
  
  sample <- data.frame(X1 = X1, X2 = X2, X3 = X3, Y = Y, E = E, Env = Env)
  return(sample)
}



# returns a sample from 5 environments from the standard SCM 
# n: number of observations per environment => total 5*n data points in sample
gen.sample <- function(n){
  sample.0 <- sim.SCM1(n,0,"noInt")
  sample.1 <- sim.SCM1(n,1,"weakposInt")
  sample.2 <- sim.SCM1(n,2,"strongposInt")
  sample.3 <- sim.SCM1(n,-1,"weaknegInt")
  sample.4 <- sim.SCM1(n,-2,"strongnegInt")
  sample <- rbind(sample.0, sample.1, sample.2, sample.3, sample.4)
  return(sample)
}



# these are the possible subsets of predictors 1:3
sets <- list(c(1,2,3), c(1,3), c(1,2), c(2,3), c(1), c(2), c(3), c(0))
#-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------
# invariance tests
#-------------------------------------------------------------------------------


# set is a subset of 1:3, where 3 is the number of predictors,
# and sample is a dataframe where the first 3 columns correspond to the 3 covariates,
# the response is sample$Y, and sample$Env is a factor with the environment indeces.


# DeLong (RF) test where the environment observations are permuted instead of omitted
# for the model which doesn't have access to the environment
pval.permute <- function(sample, set){
  Xmat <- sample[, "Env", drop = F]
  
  if(sum(set)>0.0001){
    Xmat <- sample[, c(set,which(names(sample) == "Env")), drop = F]
  }
  
  # random forest using the environment
  RF.Env <- ranger(y = as.factor(sample$Y), x = Xmat, probability = T)
  
  # now we permute the environments
  Xmat$Env <- Xmat$Env[sample(1:nrow(Xmat), size = nrow(Xmat), replace = F)]
  
  # this random forest doesn't have access to environment information
  RF.noEnv <- ranger(y = as.factor(sample$Y), x = Xmat, probability = T)
  
  # we use OOB predictions to compute the ROC curves
  roc.noEnv <- roc(response = sample$Y, predictor = RF.noEnv$predictions[,"1"], quiet = T, direction="<")
  roc.Env <- roc(response = sample$Y, predictor = RF.Env$predictions[,"1"], quiet = T, direction="<")
  
  test <- roc.test(roc1 = roc.noEnv, roc2 = roc.Env, method = "delong", alternative = "less")
  
  return(test$p.value)
}






# DeLong (RF) test where the environment observations are excluded
# for the model which doesn't have access to the environment
pval.exclude <- function(sample, set){
  
  # random forest using the environment
  RF.Env <- ranger(y = as.factor(sample$Y), 
                   x = sample[, c(set, which(names(sample) == "Env")), drop = F],
                   probability = T)
  
  # this random forest doesn't have access to environment information
  RF.noEnv <- ranger(y = as.factor(sample$Y), x = sample[, set, drop = F], probability = T)
  
  
  # we use OOB predictions to compute the ROC curves
  roc.noEnv <- roc(response = sample$Y, predictor = RF.noEnv$predictions[,2], quiet = T, direction="<")
  roc.Env <- roc(response = sample$Y, predictor = RF.Env$predictions[,2], quiet = T, direction="<")
  
  test <- roc.test(roc1 = roc.noEnv, roc2 = roc.Env, method = "delong", alternative = "less")
  
  return(test$p.value)
}


#-------------------------------------------------------------------------------









#-------------------------------------------------------------------------------
# function for running the simulation 
#-------------------------------------------------------------------------------


# set: subset of 1:3 (we will use invariant sets {1} and {1,3} here)
# n: number of observations per environment for the maximal sample size considered in the experiment
# nreps: number of simulation repetitions
# samp.sizes: vector of different sample sizes per environment to try (with values between 0 and n)
# pval.function: which invariance test to apply (either pval.permute or pval.exclude)
sim.run <- function(set, n, nreps, samp.sizes, pval.function){
  
  pvalues.mat <- matrix(0, nrow = nreps, ncol = length(samp.sizes))
  
  for(b in 1:nreps){
    print(paste0("Simulation iteration ", b, " out of ", nreps))
    
    # generate a sample of the maximal size considered (here, the rows are ordered by environment)
    sample.ordered <- gen.sample(n)
    
    # permute the rows such that the rows are not ordered by environment
    sample.total <- sample.ordered[sample(1:nrow(sample.ordered), size = nrow(sample.ordered), replace = F), ]
    
    
    # now we consider iteratively larger proportions of sample.total for the experiment to see how level changes
    for(m in 1:length(samp.sizes)){
      sample <- sample.total[1:samp.sizes[m],]
      pvalues.mat[b,m] <- pval.function(sample, set)
    }
  }
  
  # see whether test would accept or reject
  reject.mat <- ifelse(pvalues.mat < 0.05, 1, 0)
  
  
  # calculate the rejection rate over the simulation iterations
  rejection.rate <- colMeans(reject.mat)
  
  # calculate standard deviation
  sd.rejection.rate <- apply(X = reject.mat, MARGIN = 2, FUN = sd)
  
  
  df.output <- data.frame(sample.size = samp.sizes, rejection.rate = rejection.rate, sd.rej.rate = sd.rejection.rate)
  return(df.output)
}


#-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------
# running the simulation
#-------------------------------------------------------------------------------


# n observations per environment for maximally large sample => 5*n observations in total
n <- 200

# repeat experiment nreps times
# 1000 nreps should take around one hour on my MB
nreps <- 500


# samp.sizes are the iteratively larger number of observations from the total sample used
# We shouldn't start at a too low number such as 10, because then it sometimes happens
# than one label (0/1) for Y is not in the dataset
samp.sizes <- seq(from=25, to=5*n, by = 50)




set.seed(1)

# simulations for the invariant set {1,3}
df.exclude13 <- sim.run(set = sets[[2]], n=n, nreps = nreps, samp.sizes = samp.sizes, pval.function = pval.exclude)
df.permute13 <- sim.run(set = sets[[2]], n=n, nreps = nreps, samp.sizes = samp.sizes, pval.function = pval.permute)

# simulations for the invariant set {1}
df.exclude1 <- sim.run(set = sets[[5]], n=n, nreps = nreps, samp.sizes = samp.sizes, pval.function = pval.exclude)
df.permute1 <- sim.run(set = sets[[5]], n=n, nreps = nreps, samp.sizes = samp.sizes, pval.function = pval.permute)


save(df.exclude13, df.permute13, df.exclude1, df.permute1, file = file.path(script_dir, "saved_data/delong_not_level.rdata"))
#-------------------------------------------------------------------------------









#-------------------------------------------------------------------------------
# plotting
#-------------------------------------------------------------------------------



# prepare data sets for plotting
df.all13 <- df.exclude13[, c("sample.size", "rejection.rate")]
df.all13$rej.rate.permute <- df.permute13$rejection.rate
names(df.all13) <- c("sample.size", "rej.exclude", "rej.permute")
df.all13$level <- rep(0.05, nrow(df.all13))

# compute 95% t-test confidence interval
df.all13$ci.perm <- qt(p = 0.975, df = nreps-1)*df.permute13$sd.rej.rate / sqrt(nreps)
df.all13$ci.excl <- qt(p = 0.975, df = nreps-1)*df.exclude13$sd.rej.rate / sqrt(nreps)





df.all1 <- df.exclude1[, c("sample.size", "rejection.rate")]
df.all1$rej.rate.permute <- df.permute1$rejection.rate
names(df.all1) <- c("sample.size", "rej.exclude", "rej.permute")
df.all1$level <- rep(0.05, nrow(df.all1))

# compute 95% t-test confidence interval
df.all1$ci.perm <- qt(p = 0.975, df = nreps-1)*df.permute1$sd.rej.rate / sqrt(nreps)
df.all1$ci.excl <- qt(p = 0.975, df = nreps-1)*df.exclude1$sd.rej.rate / sqrt(nreps)



color.palette <- scales::hue_pal()(3)
size <- 10


# plot for set {1,3}
p13 <- ggplot(df.all13, aes(sample.size)) + 
  geom_line(aes(y = rej.exclude, colour = "exclude environment", linetype = "exclude environment")) + 
  geom_line(aes(y = rej.permute, colour = "permute environment", linetype = "permute environment")) +
  geom_line(aes(y = level, colour = "desired level (0.05)", linetype = "desired level (0.05)")) +
  geom_ribbon(aes(ymin = rej.exclude - ci.excl, ymax = rej.exclude + ci.excl), fill = color.palette[1], alpha = 0.2) +
  geom_ribbon(aes(ymin = rej.permute - ci.perm, ymax = rej.permute + ci.perm), fill = color.palette[3], alpha = 0.2) +
  coord_cartesian(ylim=c(0,1)) +
  scale_color_manual(values = c("desired level (0.05)" = "black", "exclude environment" = color.palette[1], "permute environment" = color.palette[3])) +
  scale_linetype_manual(values = c("desired level (0.05)" = "dotted", "exclude environment" = "solid", "permute environment" = "solid")) +
  labs(color = 'DeLong (RF)', linetype = 'DeLong (RF)') +
  xlab("Sample size") +
  ylab("Rejection rate") +
  theme_bw(base_size = size) +
  ggtitle("Subset S = {1,3}")
#p13



# plot for set {1}
p1 <- ggplot(df.all1, aes(sample.size)) + 
  geom_line(aes(y = rej.exclude, colour = "exclude environment", linetype = "exclude environment")) + 
  geom_line(aes(y = rej.permute, colour = "permute environment", linetype = "permute environment")) +
  geom_line(aes(y = level, colour = "desired level (0.05)", linetype = "desired level (0.05)")) +
  geom_ribbon(aes(ymin = rej.exclude - ci.excl, ymax = rej.exclude + ci.excl), fill = color.palette[1], alpha = 0.2) +
  geom_ribbon(aes(ymin = rej.permute - ci.perm, ymax = rej.permute + ci.perm), fill = color.palette[3], alpha = 0.2) +
  coord_cartesian(ylim=c(0,1)) +
  scale_color_manual(values = c("desired level (0.05)" = "black", "exclude environment" = color.palette[1], "permute environment" = color.palette[3])) +
  scale_linetype_manual(values = c("desired level (0.05)" = "dotted", "exclude environment" = "solid", "permute environment" = "solid")) +
  labs(color = 'DeLong (RF)', linetype = 'DeLong (RF)') +
  xlab("Sample size") +
  ylab("Rejection rate") +
  theme_bw(base_size = size) +
  ggtitle("Subset S = {1}")
#p1




combined <- p1 + p13 & theme(legend.position = "bottom", legend.title=element_blank(), legend.text = element_text(size=size)) 

combined + plot_layout(guides = "collect")


ggsave(filename = file.path(script_dir, "saved_plots/delong_not_level.pdf"), width = 6, height = 4)




#-------------------------------------------------------------------------------


writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/delong_not_level.txt"))








