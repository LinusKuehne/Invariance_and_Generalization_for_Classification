library(ggplot2)
library(rje)
library(viridis)


# load in the functions needed
source("invariance_tests.R")
source("utils.R")
source("data_generating_process.R")
source("hrf.R")
source("stabilized_classification.R")


#-------------------------------------------------------------------------------
# Parameters for the simulation
#-------------------------------------------------------------------------------

# sets to check stability
sets <- powerSet(1:3)
sets[[1]] <- c(0)

# number of samples per domain (training: n, testing: n.test)
n.test <- 500
n <- 500

# number of simulation runs 
n.sim <- 10

# number of bootstrap samples to compute c.pred
B <- 100

# the invariance test to be used
inv_test = "residual"

#-------------------------------------------------------------------------------

# set the seed
set.seed(1)



accuracies <- data.frame(rf = numeric(n.sim),
                         hrf.orig = numeric(n.sim),
                         hrf.ood = numeric(n.sim), 
                         sc = numeric(n.sim),
                         sc.hrf = numeric(n.sim))

sb.first.split <- data.frame(rf = numeric(n.sim),
                             hrf.orig = numeric(n.sim),
                             hrf.ood = numeric(n.sim))


for(sim in 1:n.sim){
  
  print(sim)
  s <- gen.sample.fixed(n = n, n.test = n.test, int.strength.train = 1/2, int.strength.test = 3, mod = "nonlin")
  
  # extract generated datasets
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # stabilized classification
  sc.fit <- stabilizedClassification(sample = sample, test = inv_test, B = B, verbose = F)
  pred.sc <- predict.stabClass(sc.fit, newsample = sample_test[,1:3])
  accuracies$sc[sim] <- mean(sample_test$Y == pred.sc$pred.class)
  
  # stabilized classification with HRF classifier
  pred.sc.hrf <- predict.stabClass.hrf(sc.fit, newsample = sample_test[,1:3])
  accuracies$sc.hrf[sim] <- mean(sample_test$Y == pred.sc.hrf$pred.class)
  
  
  # 
  # # standard random forest
  # rf.fit <- ranger(y = as.factor(sample$Y), x = sample[, 1:3], probability = T, num.threads = 0)
  # pred.rf <- predict(rf.fit, data = sample_test[,1:3])$predictions[,"1"]
  # accuracies$rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))
  # 
  # first.splits.rf <- sapply(rf.fit$forest$split.varIDs, function(x) x[1])
  # sb.first.split$rf[sim] <- mean(first.splits.rf %in% c(0,2))
  # 
  # # original hedged RF
  # hrf.orig.fit <- hrf.orig(y = as.factor(sample$Y), x = sample[, 1:3], rf.fit = rf.fit)
  # pred.hrf.orig <- predict.hedgedrf(hrf.orig.fit, data = sample_test[,1:3])
  # accuracies$hrf.orig[sim] <- mean(sample_test$Y == ifelse(pred.hrf.orig>0.5, 1, 0))
  # 
  # hrf.orig.splitIDs <- hrf.orig.fit$rf.fit$forest$split.varIDs
  # first.splits.hrf.orig <- sapply(hrf.orig.splitIDs, function(x) x[1])
  # sb.first.split$hrf.orig[sim] <- sum(hrf.orig.fit$tree.weights * (first.splits.hrf.orig %in% c(0,2)))
  # 
  # # ood hedged RF
  # hrf.ood.fit <- hrf.ood(y = as.factor(sample$Y), x = sample[, c(1,2,3,6)])
  # pred.hrf.ood <- predict.hedgedrf(hrf.ood.fit, data = sample_test[,1:3])
  # accuracies$hrf.ood[sim] <- mean(sample_test$Y == ifelse(pred.hrf.ood>0.5, 1, 0))
  # 
  # hrf.ood.splitIDs <- hrf.ood.fit$rf.fit$forest$split.varIDs
  # first.splits.hrf.ood <- sapply(hrf.ood.splitIDs, function(x) x[1])
  # sb.first.split$hrf.ood[sim] <- t(hrf.ood.fit$tree.weights) %*% as.matrix(as.numeric(first.splits.hrf.ood %in% c(0,2)), ncol = 1)
}



#-------------------------------------------------------------------------------
# generate plot
#-------------------------------------------------------------------------------


name <- factor(c(rep("RF",n.sim), rep("HRF orig",n.sim), rep("HRF ood", n.sim), rep("SC (RF)", n.sim), rep("SC (HRF)", n.sim)), levels = c("RF", "HRF orig", "HRF ood", "SC (RF)", "SC (HRF)"))

# create a dataframe
data.acc <- data.frame(
  name=name,
  value=c(accuracies$rf, accuracies$hrf.orig, accuracies$hrf.ood, accuracies$sc, accuracies$sc.hrf)
)





size <- 11.5


# non-linear logistic regression

plt.acc <- ggplot(data.acc, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(color="black", size=0.025, alpha=0.2) +
  ggtitle("Prediction Performance on Test Domains") +
  xlab("") +
  ylim(0.65,0.95)+
  ylab("Accuracy") +
  theme_bw(base_size = size) +
  theme(legend.position="none") 
plt.acc




colMeans(sb.first.split)




# next steps:
# try with more noise, more difficult DGP (maybe larger graph, maybe random SCM)
# in this context: combine SC with HRF
# use kappa = 1 => everything non-negative => can use weighted mean





