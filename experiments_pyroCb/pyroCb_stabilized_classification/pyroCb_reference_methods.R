# in this script, we compute the prediction performance of the reference models, i.e.,
# random guessing and using the training labels (i.e. the empty set of predictors)


library(rje)


# get the path of this script
script_dir <- getwd()


# load in the dataset and the grouping into different environments
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))


# load in utilities
source("../../code/code_pyroCb/pyroCb_stabilized_classification_utils.R")



# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# sets to check invarinace
sets <- powerSet(varincl)
sets[[1]] <- c(0)



# convert factor into numeric vector
y.num <- as.numeric(labels)-1

# group into five environments
envs <- env5


# initialize vectors for losses
wbce.per.env.random <- wbce.per.env.empty <- numeric(length(levels(envs)))

set.seed(1)


# iterate over the environments with LOEO CV 
for(e in 1:length(levels(envs))){
  
  # get indices for currently held-out environment
  i.test <- which(envs == levels(envs)[e])
  i.train <- -i.test
  
  # get the train/test labels
  labels.train <- labels[i.train]
  labels.test <- labels[i.test]
  y.num.train <- y.num[i.train]
  y.num.test <- y.num[i.test]
  
  # sample the "predicted probabilities"
  pred.probs.random <- runif(n = length(y.num.test), min = 0, max = 1)
  
  # use the mean of the labels on the training data as the predicted probability
  pred.probs.empty <- rep(mean(y.num.train), length(y.num.test))
  
  # compute losses of both approaches
  wbce.per.env.random[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.random)
  wbce.per.env.empty[e] <- BCE.weighted(y = y.num.test, y.hat = pred.probs.empty)
}


# compute worst-case and mean out-of-distribution LOEO CV losses
max(wbce.per.env.empty)
mean(wbce.per.env.empty)

max(wbce.per.env.random)
mean(wbce.per.env.random)







# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/pyroCb_reference_methods.txt"))





