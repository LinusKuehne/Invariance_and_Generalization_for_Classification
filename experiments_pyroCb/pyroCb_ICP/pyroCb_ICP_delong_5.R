# in this script, we run ICP with the DeLong test on five environments




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



pvals.delong.5.1tail <- numeric(length(sets))





set.seed(1)


for(s in 1:length(sets)){
  print(paste0("Working on set number ", s, " out of ", length(sets)))
  
  set <- sets[[s]]
  
  test.delong.5 <- delong(set, cube, labels, envVar = matrix(env5, ncol = 1), cluster.assoc = event_df$cluster_random, posts)
  
  
  pvals.delong.5.1tail[s] <- test.delong.5$pval_1tail

  
}




inv.sets.delong.5 <- sets[pvals.delong.5.1tail>a.inv]


# Compute the intersection of all vectors in the list
intersection.delong.5 <- Reduce(intersect, inv.sets.delong.5)





save(pvals.delong.5.1tail, 
     inv.sets.delong.5,
     intersection.delong.5,
     file = "../saved_data/pyroCb_ICP_delong_5.rdata")




# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/pyroCb_ICP_delong_5.txt"))




























