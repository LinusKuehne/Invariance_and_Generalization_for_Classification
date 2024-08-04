# in this script, we run ICP with the tram-GCM (RF) test with 9 environments on the pyroCb dataset


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
pvals.tram.9 <- numeric(length(sets))


set.seed(1)



# iterate over all subsets of predictors to compute the p-values with the invariance test
for(s in 1:length(sets)){
  print(paste0("Working on set number ", s, " out of ", length(sets)))
  
  # extract current set
  set <- sets[[s]]
  
  # run the invariance test
  test.tram.9 <- rangerICP.paper(set, cube, labels, y.num, envs = env9, cluster.assoc = event_df$cluster_random, posts)  
  
  # store p-value
  pvals.tram.9[s] <- test.tram.9
}



# find the subsets which are invariant at the significance level a.inv
inv.sets.tram.9 <- sets[pvals.tram.9>a.inv]


# Compute the intersection of all vectors in the list
intersection.tram.9 <- Reduce(intersect, inv.sets.tram.9)


# save the results
save(pvals.tram.9, 
     inv.sets.tram.9,
     intersection.tram.9,
     file = "../saved_data/pyroCb_ICP_tram_9.rdata")




# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/pyroCb_ICP_tram_9.txt"))




























