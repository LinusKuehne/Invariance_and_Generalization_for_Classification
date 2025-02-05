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
sets <- powerSet(1:7)
sets[[1]] <- c(0)

# number of samples per domain (training: n, testing: n.test)
n.test <- 500
n <- 500

# number of simulation runs 
n.sim <- 150

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
                         sc = numeric(n.sim),
                         sc.hrf = numeric(n.sim))

sb.first.split <- data.frame(rf = numeric(n.sim),
                             hrf.orig = numeric(n.sim),
                             hrf.ood = numeric(n.sim))


for(sim in 1:n.sim){
  
  print(sim)
  s <- gen.sample.bigger(n.train = n, n.test = n.test)
  
  # extract generated datasets
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # stabilized classification
  sc.fit <- stabilizedClassification(sample = sample, test = inv_test, B = B, verbose = F)
  pred.sc <- predict.stabClass(sc.fit, newsample = sample_test[,1:7])
  accuracies$sc[sim] <- mean(sample_test$Y == pred.sc$pred.class)

  # stabilized classification with HRF classifier
  pred.sc.hrf <- predict.stabClass.hrf(sc.fit, newsample = sample_test[,1:7])
  accuracies$sc.hrf[sim] <- mean(sample_test$Y == pred.sc.hrf$pred.class)
  
  # oracle RF fitted on the stable blanket
  sb.fit <- ranger(y = as.factor(sample$Y), x = sample[, c(1,2,3,5)], probability = T, num.threads = 0)
  pred.sb <- predict(sb.fit, data = sample_test[, c(1,2,3,5)])$predictions[,"1"]
  accuracies$sb[sim] <- mean(sample_test$Y == ifelse(pred.sb>0.5, 1, 0))

  # standard random forest
  rf.fit <- ranger(y = as.factor(sample$Y), x = sample[, 1:7], probability = T, num.threads = 0)
  pred.rf <- predict(rf.fit, data = sample_test[,1:7])$predictions[,"1"]
  accuracies$rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))

  first.splits.rf <- sapply(rf.fit$forest$split.varIDs, function(x) x[1])
  sb.first.split$rf[sim] <- mean(first.splits.rf %in% c(0,1,2,4))

  # original hedged RF
  hrf.orig.fit <- hrf.orig(y = as.factor(sample$Y), x = sample[, 1:7], rf.fit = rf.fit)
  pred.hrf.orig <- predict.hedgedrf(hrf.orig.fit, data = sample_test[,1:7])
  accuracies$hrf.orig[sim] <- mean(sample_test$Y == ifelse(pred.hrf.orig>0.5, 1, 0))

  hrf.orig.splitIDs <- hrf.orig.fit$rf.fit$forest$split.varIDs
  first.splits.hrf.orig <- sapply(hrf.orig.splitIDs, function(x) x[1])
  sb.first.split$hrf.orig[sim] <- sum(hrf.orig.fit$tree.weights * (first.splits.hrf.orig %in% c(0,1,2,4)))

  # ood hedged RF
  hrf.ood.fit <- hrf.ood(y = as.factor(sample$Y), x = sample[, c(1,2,3,4,5,6,7,9)])
  pred.hrf.ood <- predict.hedgedrf(hrf.ood.fit, data = sample_test[,1:7])
  accuracies$hrf.ood[sim] <- mean(sample_test$Y == ifelse(pred.hrf.ood>0.5, 1, 0))

  hrf.ood.splitIDs <- hrf.ood.fit$rf.fit$forest$split.varIDs
  first.splits.hrf.ood <- sapply(hrf.ood.splitIDs, function(x) x[1])
  sb.first.split$hrf.ood[sim] <- t(hrf.ood.fit$tree.weights) %*% as.matrix(as.numeric(first.splits.hrf.ood %in% c(0,1,2,4)), ncol = 1)
}



#-------------------------------------------------------------------------------
# generate plot
#-------------------------------------------------------------------------------


name <- factor(c(rep("SB (oracle)",n.sim), rep("RF",n.sim), rep("HRF orig",n.sim), rep("HRF ood", n.sim), rep("SC (RF)", n.sim), rep("SC (HRF)", n.sim)), levels = c("SB (oracle)", "RF", "HRF orig", "HRF ood", "SC (RF)", "SC (HRF)"))

# create a dataframe
data.acc <- data.frame(
  name=name,
  value=c(accuracies$sb, accuracies$rf, accuracies$hrf.orig, accuracies$hrf.ood, accuracies$sc, accuracies$sc.hrf)
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

ggsave(filename = "hrf_sc_acc.pdf", width = 7.5, height = 7.5)


sb.ratio <- colMeans(sb.first.split)

save(sb.ratio, file = "sb_ratio.rdata")
load("sb_ratio.rdata")






