# in this script, we run ICP with the correlation test


library(ranger)
library(pROC)
library(rje)

# get the path of this script
script_dir <- getwd()


# load in the dataset
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))

source("../../code/code_pyroCb/pyroCb_invariance_tests.R")


# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# sets to check stability
sets <- powerSet(varincl)
sets[[1]] <- c(0)



# tuning parameter
a.inv <- 0.05

y.num <- as.numeric(labels)-1


pvals.corr.5 <- pvals.corr.9 <- numeric(length(sets))


set.seed(1)


for(s in 1:length(sets)){
  print(paste0("Working on set number ", s, " out of ", length(sets)))
  
  set <- sets[[s]]
  
  test.corr <- corr(set, cube, labels, y.num, group5 = env5, group9 = env9, cluster.assoc = event_df$cluster_random, posts)
  
  pvals.corr.5[s] <- test.corr$p.val_fiveEnv
  pvals.corr.9[s] <- test.corr$p.val_nineEnv
}


inv.sets.corr.5 <- sets[pvals.corr.5>a.inv]
inv.sets.corr.9 <- sets[pvals.corr.9>a.inv]


# Compute the intersection of all vectors in the list
intersection.corr.5 <- Reduce(intersect, inv.sets.corr.5)
intersection.corr.9 <- Reduce(intersect, inv.sets.corr.9)


save(pvals.corr.5, 
     pvals.corr.9, 
     inv.sets.corr.5,
     inv.sets.corr.9,
     intersection.corr.5,
     intersection.corr.9,
     file = "../saved_data/pyroCb_ICP_corr.rdata")




# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/pyroCb_ICP_corr.txt"))




























