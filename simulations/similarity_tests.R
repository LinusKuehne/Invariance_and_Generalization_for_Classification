# in this script, we evaluate to which degree the ranking of subsets according to 
# the p-values by the different tests is similar


library(rje)
library(ggplot2)
library(reshape2)



# get the path of this script
script_dir <- getwd()

# load in the functions needed
source(file.path(script_dir, "../code/code_simulations/invariance_tests.R"))
source(file.path(script_dir, "../code/code_simulations/data_generating_process.R"))
source(file.path(script_dir, "../code/code_simulations/utils.R"))



#-------------------------------------------------------------------------------
# parameters
#-------------------------------------------------------------------------------

# sets to check stability
sets <- powerSet(1:5)
sets[[1]] <- c(0)

# number of samples per environment 
n <- 250

# number of simulation runs 
n.sim <- 300

# for plotting
size <- 10

# max intervention strengths
strength.train <- 5
strength.test <- 8

# number of covariates (including Y)
d <- 6

# max number of parents
max.pa <- 3 

# number of interventions
num.int <- 2 


#-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------
# run simulation
#-------------------------------------------------------------------------------

set.seed(1)

# store results 
mat <- matrix(0, nrow = 7, ncol = 7)



for(sim in 1:n.sim){
  
  print(paste0("Simulation iteration ", sim, " out of ", n.sim))
  
  # generate a sample of the random SCM
  s <- generate.samples.random(n = n, n.test = 1000, d = d, max.pa = max.pa, num.int = num.int, int.strength.train = strength.train, int.strength.test = strength.test, mod = "logreg")
  
  # extract generated dataset
  sample <- s$sample_train
  
  # compute p-values for all subsets
  pvals <- pvalues.all(sample)
  
  
  
  # find the worst-case error for subsets --------------------------------------
  
  wbce.worst <- numeric(length(sets))
  
  env <- sample$Env
  
  # separately deal with empty set
  wbce.vec.empty <- numeric(length(levels(env)))
  for(e in 1:length(levels(env))){
    
    test.indx <- which(env == levels(env)[e])
    train.indx <- -test.indx
    
    probs <- rep(mean(sample$Y[train.indx]), length(test.indx))
    wbce.vec.empty[e] <- BCE.weighted(y = sample$Y[test.indx], y.hat = probs)
  }
  
  wbce.worst[1] <- max(wbce.vec.empty)
  

  # now for the non-empty sets
  for(s in 2:length(sets)){
    
    set <- sets[[s]]
    
    dat.set <- sample[, set, drop = F]
    
    wbce.vec <- numeric(length(levels(env)))
    
    # iterate over the environments
    for(e in 1:length(levels(env))){
      
      test.indx <- which(env == levels(env)[e])
      train.indx <- -test.indx
      
      
      table_y <- table(sample$Y[train.indx])  # frequency of each class
      
      weights <- length(sample$Y[train.indx])/(2*table_y)  # one half times inverse of class frequency 
      
       
      
      rf <- ranger(y = as.factor(sample$Y[train.indx]), x = dat.set[train.indx, , drop = F], class.weights = weights, probability = T)
      probs <- predict(rf, data = dat.set[test.indx, , drop = F])$predictions[,"1"]
      
      wbce.vec[e] <- BCE.weighted(y = sample$Y[test.indx], y.hat = probs)
      
    }
    
    wbce.worst[s] <- max(wbce.vec)
    
  }
  
  
  # we want a positive Kendall's tau if the ranking is similar:
  neg.wbce.worst <- -wbce.worst
  
  # ----------------------------------------------------------------------------
  
  
  # add the generalization scores to the dataframe
  pvals$neg.wbce.worst <- neg.wbce.worst
  

  # compute Kendall's tau between all vector pairs -----------------------------
  
  tau.mat <- matrix(0, nrow = ncol(pvals), ncol = ncol(pvals))
  
  
  
  for(j in 1:(ncol(tau.mat)-1)){
    for(i in (j+1):ncol(tau.mat)){
      tau.mat[i,j] <- cor(x = pvals[,i], y = pvals[,j], method = "kendall")
    }
  }

  # add this result to mat
  mat <- mat + tau.mat
  # ----------------------------------------------------------------------------
  
}


# divide by the number of simulations as we have always added mat in each iteration
mat <- mat/n.sim

#-------------------------------------------------------------------------------


save(mat, file = file.path(script_dir, "saved_data/similarity_tests.rdata"))





#-------------------------------------------------------------------------------
# plotting
#-------------------------------------------------------------------------------


# convert to a symmetric matrix
mat <- mat + t(mat)

diag(mat) <- NA

df <- melt(mat)


# generate grid
x <- c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual", "Neg. LOEO CV loss")
y <- c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual", "Neg. LOEO CV loss")
data <- expand.grid(X=x, Y=y)
data$Z <- df$value


# assign colors to kendall tau values and plot
plt <- ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", na.value = "black")  +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Consistency between different rankings", fill = "Kendall's Tau") +
  theme_bw(base_size = size) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_blank())

plt



ggsave(file.path(script_dir, "saved_plots/similarity_tests.pdf"), width = 6, height = 5.07)






#-------------------------------------------------------------------------------

# store sessionInfo()
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/similarity_tests.txt"))


























