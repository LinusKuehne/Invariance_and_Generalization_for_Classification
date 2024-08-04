# in this script, we run ICP with the continuous DeLong test on the pyroCb dataset


library(ranger)
library(pROC)
library(rje)

# get the path of this script
script_dir <- getwd()


# load in the dataset
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))


# load the invariance tests implemented for the pyroCb data
source("../../code/code_pyroCb/pyroCb_invariance_tests.R")


names(envVars) <- c("lat", "long", "date")


# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# sets to check invariance
sets <- powerSet(varincl)
sets[[1]] <- c(0)


# tuning parameter
a.inv <- 0.05

# convert labels to numeric vector
y.num <- as.numeric(labels)-1

# initialize
pvals.delong.cont.1tail <- numeric(length(sets))


set.seed(1)

# iterate over all subsets of predictors to compute the p-values with the invariance test
for(s in 1:length(sets)){
  print(paste0("Working on set number ", s, " out of ", length(sets)))
  
  # extract current set
  set <- sets[[s]]
  
  # run the invariance test
  test.delong.cont <- delong(set, cube, labels, envVar = envVars, cluster.assoc = event_df$cluster_random, posts)
  
  # store p-value
  pvals.delong.cont.1tail[s] <- test.delong.cont$pval_1tail
}

# find the subsets which are invariant at the significance level a.inv
inv.sets.delong.cont <- sets[pvals.delong.cont.1tail>a.inv]


# Compute the intersection of all vectors in the list
intersection.delong.cont <- Reduce(intersect, inv.sets.delong.cont)


# save the results
save(pvals.delong.cont.1tail, 
     inv.sets.delong.cont,
     intersection.delong.cont,
     file = "../saved_data/pyroCb_ICP_delong_cont.rdata")




# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/pyroCb_ICP_delong_cont.txt"))




























