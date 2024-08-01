# in this script, we compute the prediction performance of random guessing 


library(rje)


# get the path of this script
script_dir <- getwd()


# load in the dataset
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))



source("../../code/code_pyroCb/pyroCb_stabilized_classification_utils.R")



# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# sets to check stability
sets <- powerSet(varincl)
sets[[1]] <- c(0)




y.num <- as.numeric(labels)-1

# group into five environments
envs <- env5



wbce.per.env.random <- wbce.per.env.empty <- numeric(length(levels(envs)))

set.seed(1)


# LOEO CV
for(e in 1:length(levels(envs))){
  
  i.test <- which(envs == levels(envs)[e])
  i.train <- -i.test
  
  labels.train <- labels[i.train]
  labels.test <- labels[i.test]
  
  y.num.train <- y.num[i.train]
  y.num.test <- y.num[i.test]
  
  pred.probs.random <- runif(n = length(y.num.test), min = 0, max = 1)
  
  pred.probs.empty <- rep(mean(y.num.train), length(y.num.test))
  
  wbce.per.env.random[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.random)
  wbce.per.env.empty[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.empty)

}

max(wbce.per.env.empty)
mean(wbce.per.env.empty)

max(wbce.per.env.random)
mean(wbce.per.env.random)







# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/pyroCb_random_guessing.txt"))





