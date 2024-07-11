# In this script, we perform variable screening for the pyroCb data. We compute the order
# of exclusion of all variables using the group lasso for logistic regression


library(grpreg)


# get the path of this script
script_dir <- getwd()


# load in the dataset
load(file.path(script_dir, "../../data/exported_pyrocb.rdata"))


# obtain data for the 28 variables we consider
indxincl <- as.vector(unlist(sapply(X = indxIncl, function(i) posts[i]:(posts[i+1]-1))))
dat <- cube[, indxincl]


# rename the data 
X <- dat
y <- as.numeric(labels)-1


# define the structure of the grouping (the last two variables don't have 11 summary statistics)
groups <- numeric(318)
first.part <- rep(1:28, each = 11)
groups[1:length(first.part)] <- first.part
groups[(length(first.part)+1):length(groups)] <- rep(29,10)
groups[(length(groups)-5):length(groups)] <- rep(30,6)
ind <- which(groups %in% indxIncl)
groups <- groups[ind]


# fit the group lasso for logistic regression
fit <- grpreg(X = X, y = y, group = groups, penalty="grLasso", family = "binomial")



# plot lasso traces if interested
#plot(fit, label = T)




# num.var: how many variables we want to include
for(num.var in 1:28){
  
  
  lambda.ind <- -1
  
  # find the lambda yielding the desired number of included variables
  for(i in 1:length(fit$lambda)){
    n.incl <- predict(fit, type="ngroups", lambda=fit$lambda[i])
    
    if(n.incl >= num.var){
      lambda.ind <- i
      break
    }
  }
  
  
  # determine the identity of nonzero groups
  id.nonzero <- predict(fit, type="groups", lambda=fit$lambda[lambda.ind])   
  
  # convert to correct variable number
  indx.nonzero <- as.numeric(levels(id.nonzero))[id.nonzero]
  
  # order them
  indx.nonzero <- indx.nonzero[order(indx.nonzero)]
  
  print(paste0("When screening to ", length(indx.nonzero), " variables, include the following:"))
  print(varss[indx.nonzero])
  
}



