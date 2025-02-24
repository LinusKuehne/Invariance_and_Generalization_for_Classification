library(ggplot2)
library(rje)
library(viridis)


# load in the functions needed
source("simulation_GHRF/invariance_tests.R")
source("simulation_GHRF/utils.R")
source("simulation_GHRF/data_generating_process.R")
source("simulation_GHRF/hrf.R")
#source("simulation_GHRF/stabilized_classification.R")


#-------------------------------------------------------------------------------
# Parameters for the simulation
#-------------------------------------------------------------------------------

# number of samples per domain (training: n, testing: n.test)
n.test <- 500
n <- 500

# find number of predictors
sample <- gen.sample.bigger(n.train = n, n.test = n.test)$sample_train
d <- length(grep("^X", names(sample)))
env.col.idx <- which(names(sample) == "Env")

# stable blanket
stable.blanket <- c(1,2,3,5)

# sets to check stability
sets <- powerSet(1:d)
sets[[1]] <- c(0)



# number of simulation runs 
n.sim <- 50

# number of bootstrap samples to compute c.pred
B <- 50

# the invariance test to be used
inv_test = "residual"

#-------------------------------------------------------------------------------

# set the seed
set.seed(1)



accuracies <- data.frame(sb = numeric(n.sim),
                         rf = numeric(n.sim),
                         hrf.orig = numeric(n.sim),
                         hrf.ood = numeric(n.sim),
                         ghrf.1 = numeric(n.sim),
                         ghrf.2 = numeric(n.sim),
                         ghrf.3 = numeric(n.sim),
                         ghrf.4 = numeric(n.sim),
                         sc = numeric(n.sim),
                         sc.hrf = numeric(n.sim))

sb.first.split <- data.frame(rf = numeric(n.sim),
                             hrf.orig = numeric(n.sim),
                             hrf.ood = numeric(n.sim),
                             ghrf.1 = numeric(n.sim),
                             ghrf.2 = numeric(n.sim))


for(sim in 1:n.sim){
  
  print(sim)
  s <- gen.sample.bigger(n.train = n, n.test = n.test)
  
  # extract generated datasets
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # # stabilized classification
  # sc.fit <- stabilizedClassification(sample = sample, test = inv_test, B = B, verbose = F)
  # pred.sc <- predict.stabClass(sc.fit, newsample = sample_test[,1:d])
  # accuracies$sc[sim] <- mean(sample_test$Y == pred.sc$pred.class)
  # 
  # # stabilized classification with HRF classifier
  # pred.sc.hrf <- predict.stabClass.hrf(sc.fit, newsample = sample_test[,1:d])
  # accuracies$sc.hrf[sim] <- mean(sample_test$Y == pred.sc.hrf$pred.class)
  
  # oracle RF fitted on the stable blanket
  sb.fit <- ranger(y = as.factor(sample$Y), x = sample[, stable.blanket], probability = T, num.threads = 0)
  pred.sb <- predict(sb.fit, data = sample_test[, stable.blanket])$predictions[,"1"]
  accuracies$sb[sim] <- mean(sample_test$Y == ifelse(pred.sb>0.5, 1, 0))

  # standard random forest
  rf.fit <- ranger(y = as.factor(sample$Y), x = sample[, 1:d], probability = T, num.threads = 0)
  pred.rf <- predict(rf.fit, data = sample_test[,1:d])$predictions[,"1"]
  accuracies$rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))

  first.splits.rf <- sapply(rf.fit$forest$split.varIDs, function(x) x[1])
  sb.first.split$rf[sim] <- mean(first.splits.rf %in% (stable.blanket-1))

  # original hedged RF
  hrf.orig.fit <- hrf.orig(y = as.factor(sample$Y), x = sample[, 1:d], rf.fit = rf.fit)
  pred.hrf.orig <- predict.hedgedrf(hrf.orig.fit, data = sample_test[,1:d])
  accuracies$hrf.orig[sim] <- mean(sample_test$Y == ifelse(pred.hrf.orig>0.5, 1, 0))

  hrf.orig.splitIDs <- hrf.orig.fit$rf.fit$forest$split.varIDs
  first.splits.hrf.orig <- sapply(hrf.orig.splitIDs, function(x) x[1])
  sb.first.split$hrf.orig[sim] <- sum(hrf.orig.fit$tree.weights * (first.splits.hrf.orig %in% (stable.blanket-1)))

  # ood hedged RF
  hrf.ood.fit <- hrf.ood(y = as.factor(sample$Y), x = sample[, c(1:d,env.col.idx)])
  pred.hrf.ood <- predict.hedgedrf(hrf.ood.fit, data = sample_test[,1:d])
  accuracies$hrf.ood[sim] <- mean(sample_test$Y == ifelse(pred.hrf.ood>0.5, 1, 0))

  hrf.ood.splitIDs <- hrf.ood.fit$rf.fit$forest$split.varIDs
  first.splits.hrf.ood <- sapply(hrf.ood.splitIDs, function(x) x[1])
  sb.first.split$hrf.ood[sim] <- t(hrf.ood.fit$tree.weights) %*% as.matrix(as.numeric(first.splits.hrf.ood %in% (stable.blanket-1)), ncol = 1)
  
  # GHRF 1
  ghrf.1.fit <- ghrf.1(y = as.factor(sample$Y), x = sample[, c(1:d,env.col.idx)])
  pred.ghrf.1 <- predict.ghrf(object = ghrf.1.fit, data = sample_test[,1:d])
  accuracies$ghrf.1[sim] <- mean(sample_test$Y == ifelse(pred.ghrf.1>0.5, 1, 0))
  
  #ghrf.1.splitIDs <- ghrf.1.fit$rf.fit$forest$split.varIDs
  #first.splits.ghrf.1 <- sapply(ghrf.1.splitIDs, function(x) x[1])
  #sb.first.split$ghrf.1[sim] <- t(ghrf.1.fit$tree.weights) %*% as.matrix(as.numeric(first.splits.ghrf.1 %in% (stable.blanket-1)), ncol = 1)
  
  # GHRF 2
  ghrf.2.fit <- ghrf.2(y = as.factor(sample$Y), x = sample[, c(1:d,env.col.idx)])
  pred.ghrf.2 <- predict.ghrf(object = ghrf.2.fit, data = sample_test[,1:d])
  accuracies$ghrf.2[sim] <- mean(sample_test$Y == ifelse(pred.ghrf.2>0.5, 1, 0))
  
  # GHRF 3
  ghrf.3.fit <- ghrf.3(y = as.factor(sample$Y), x = sample[, c(1:d,env.col.idx)])
  pred.ghrf.3 <- predict.ghrf(object = ghrf.3.fit, data = sample_test[,1:d])
  accuracies$ghrf.3[sim] <- mean(sample_test$Y == ifelse(pred.ghrf.3>0.5, 1, 0))
  
  # GHRF 4
  ghrf.4.fit <- ghrf.4(y = as.factor(sample$Y), x = sample[, c(1:d,env.col.idx)])
  pred.ghrf.4 <- predict.ghrf(object = ghrf.4.fit, data = sample_test[,1:d])
  accuracies$ghrf.4[sim] <- mean(sample_test$Y == ifelse(pred.ghrf.4>0.5, 1, 0))
}



#-------------------------------------------------------------------------------
# generate plot
#-------------------------------------------------------------------------------


name <- factor(c(rep("SB (oracle)",n.sim), rep("RF",n.sim), rep("HRF",n.sim), rep("HRF ood", n.sim), rep("GHRF 1", n.sim), rep("GHRF 2", n.sim), rep("GHRF 3", n.sim), rep("GHRF 4", n.sim)), levels = c("SB (oracle)", "RF", "HRF", "HRF ood", "GHRF 1", "GHRF 2", "GHRF 3", "GHRF 4"))

# create a dataframe
data.acc <- data.frame(
  name=name,
  value=c(accuracies$sb, accuracies$rf, accuracies$hrf.orig, accuracies$hrf.ood, accuracies$ghrf.1, accuracies$ghrf.2, accuracies$ghrf.3, accuracies$ghrf.4)
)





size <- 11.5


# 

plt.acc <- ggplot(data.acc, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  ggtitle("Prediction accuracy on 10 test domains\nSCM with 7 informative & 12 noise predictors") +
  xlab("") +
  ylim(0.75,0.91)+
  ylab("Accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
plt.acc

ggsave(filename = "simulation_GHRF/ghrf_scm_12_noise.pdf", width = 7.5, height = 7.5)


sb.ratio <- colMeans(sb.first.split)

save(sb.ratio, file = "simulation_GHRF/sb_ratio.rdata")


# TODO: weigh the domains by a difficulty measure




