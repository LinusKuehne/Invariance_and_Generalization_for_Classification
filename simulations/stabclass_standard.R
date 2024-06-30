# In this script, we test the performance of Stabilized Classification on data from 
# the standard SCM


library(ggplot2)
library(patchwork)
library(rje)
#library(tidyverse)
#library(hrbrthemes)
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
n.sim <- 20

# number of bootstrap samples to compute c.pred
B <- 20

# the test to be used
test = "tram.glm"

#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# run the simulation
#-------------------------------------------------------------------------------


# set the seed
set.seed(1)

# accuracies of the models
accuracies <- data.frame(acc.sc.glmglm = numeric(n.sim),
                         acc.sc.rfglm = numeric(n.sim),
                         acc.sc.pred.anova = numeric(n.sim),
                         acc.glm = numeric(n.sim),
                         acc.rf = numeric(n.sim))



for(sim in 1:n.sim){
  print(paste0("Simulation iteration ",sim, " out of ", n.sim))
  
  # generate a sample the fixed SCM
  s <- gen.sample.fixed(n = n, n.test = n.test, int.strength.train = 1/2, int.strength.test = 4)
  
  # extract generated objects
  sample <- s$sample_train
  sample_test <- s$sample_test
  
  # stabilized classification
  output.sc <- stabilizedClassification(sample = sample, test = test, mod.internal = "GLM", mod.output = "GLM", B = B, verbose = F)
  pred.sc <- predict.stabClass(output.sc, newsample = sample_test[,1:3], mod = "GLM")
  accuracies.fixed$acc.sc.anova[sim] <- mean(sample_test$Y == pred.sc$pred.class)
  
  
  
  # standard logistic regression
  output.glm <- glm(Y ~ ., data = sample[, 1:4], family = binomial(link = "logit"))
  pred.glm <- predict(output.glm, newdata = sample_test[,1:3], type = "response")
  accuracies.fixed$acc.glm[sim] <- mean(sample_test$Y == ifelse(pred.glm>0.5, 1, 0))
  
  # standard random forest
  output.rf <- ranger(y = as.factor(sample$Y), x = sample[, 1:3], probability = T)
  pred.rf <- predict(output.rf, data = sample_test[,1:3])$predictions[,"1"]
  accuracies.fixed$acc.rf[sim] <- mean(sample_test$Y == ifelse(pred.rf>0.5, 1, 0))
}





name <- factor(c(rep("SC",n.sim), rep("SCpred",n.sim), rep("Log. Reg.", n.sim), rep("RF", n.sim)), levels = c("SC", "SCpred", "Log. Reg.", "RF"))

# create a dataset
data.standard <- data.frame(
  name=name,
  value=c(accuracies.fixed$acc.sc.anova, accuracies.fixed$acc.sc.pred.anova, accuracies.fixed$acc.glm, accuracies.fixed$acc.rf )
)

# Plot

plt.standard <- ggplot(data.standard, aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.2, alpha=0.9) +
  ggtitle("Accuracy Standard SCM") +
  xlab("") +
  ylim(0,1)+
  ylab("Test accuracy") +
  theme_bw() +
  theme(legend.position="none") 

plt.standard


























