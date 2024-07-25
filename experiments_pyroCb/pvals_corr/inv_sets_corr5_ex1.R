library(ranger)
library(pROC)
library(rje)



# get the path of this script
script_dir <- getwd()


# load in the dataset
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))
load(file.path(script_dir, "../saved_data/discrete_envs.rdata"))

source(file.path(script_dir, "../../code/code_pyroCb/pyroCb_invariance_tests.R"))
source(file.path(script_dir, "../../code/code_pyroCb/pyroCb_stabilized_classification_utils.R"))






# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# sets to check stability
sets <- powerSet(varincl)
sets[[1]] <- c(0)



# tuning parameters
a.inv <- 0.05
a.pred <- 0.05

# number of bootstrap samples
B <- 100


y.num <- as.numeric(labels)-1





envs <- env5




# exclude1 ---------------------------------------------------------------------
e <- 1
# exclude1 ---------------------------------------------------------------------







pvals <- matrix(NA, nrow = length(sets), ncol = length(levels(envs)))

set.seed(1)

i.test <- which(envs == levels(envs)[e])
i.train <- -i.test

train.env <- droplevels(envs[i.train])

X.train <- cube[i.train, ]
X.val <- cube[i.test, ]

labels.train <- labels[i.train]
labels.test <- labels[i.test]

y.num.train <- y.num[i.train]
y.num.test <- y.num[i.test]






event_df_train <- unlist((event_df$wildfire_id)[i.train])
event_df_train <- data.frame("wildfire_id" = event_df_train)
seg_base <- event_df_train[!duplicated(event_df_train), , drop = F]
folds <- sample(cut(1:nrow(seg_base), breaks = 5, labels = F), replace = F)
seg_base$clusters_env <- folds
event_df_train <- merge(x = event_df_train, y = seg_base, all.x = T)



pvals.e <- numeric(length(sets))



for(s in 1:length(sets)){
  set <- sets[[s]]
  print(paste0("Compute set ", s, " for env ", e))
  test.corr <- corr.single(set = set, cube = X.train, labels = labels.train, y.num = y.num.train, env_test = train.env, cluster.assoc = event_df_train$clusters_env, posts)
  
  # if the environments are predicted perfectly, there is no variations in the 
  # residual matrix corresponding to the environments. This yields an NA for the
  # p-value of the correlation test. In this case, since there is no variation, 
  # we cannot conduct the test. Therefore, we also cannot reject the invariance
  # and we therefore set the p-value to 1
  if(is.na(test.corr)){
    test.corr <- 1
  }
  
  pvals.e[s] <- test.corr
}

pvals[, e] <- pvals.e






pvals_ex1 <- pvals



save(pvals_ex1, file = file.path(script_dir, "inv_sets_corr5_ex1.rdata"))



# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "inv_sets_corr5_ex1.txt"))




