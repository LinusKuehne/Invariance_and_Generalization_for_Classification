# In this script, we generate data from the random SCM, and from the semi-random SCM with 
# fixed DAG and random response and edge weights, and compute the p-values
# using all of our invariance tests for the invariant sets. 
# We then plot the empirical CDF (ECDF) of these p-values (p-values of all invariant
# sets pooled together). If the tests are level, the ECDFs should lie on or below the diagonal.




library(ggplot2)
library(patchwork)
library(rje)

source("invariance_tests.R")
source("DGP.R")
source("utils.R")




#-------------------------------------------------------------------------------
# parameters
#-------------------------------------------------------------------------------

# number of simulation repetitions
# 250 should take around 2*3 h on my MB
nreps <- 250

# n observations per environment => 5*n observations in total
n <- 200

# DAG has 6 variables (including Y)
d <- 6

# DAG interventions on two variables
num.int <- 2

# a variable should have at most 3 parents
max.pa <- 3

# predictor names
Xnames <- rep("A", d-1)
for(w in 1:(d-1)){
  number <- as.character(w)
  name <- paste("X", number, sep="") 
  Xnames[w] <- name  
}

# power set of all predictors (names)
ps <- powerSet(Xnames)

# all possible subsets of predictors 1:(d-1) 
sets <- powerSet(1:(d-1))
sets[[1]] <- 0


#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# run the simulation for the random SCM
#-------------------------------------------------------------------------------


set.seed(1)


# in these vectors, we pool together the pvalues for all invariant sets across all simulations
null.pvals.delong.rf.randDAG <- numeric(0)
null.pvals.delong.glm.randDAG <- numeric(0)
null.pvals.tram.rf.randDAG <- numeric(0)
null.pvals.tram.glm.randDAG <- numeric(0)
null.pvals.corr.randDAG <- numeric(0)
null.pvals.residual.randDAG <- numeric(0)


for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # generate a sample (n.test is arbitrary, we just use "training" data)
  s <- generate.samples.random(n = n, n.test = 1000, d = d, max.pa = max.pa, num.int = num.int, t = 1)
  sample <- s$sample_train
  
  # compute p-values
  pvals <- pvalues.all(sample)
  
  # find invariant sets from DAG 
  inv.sets <- int.stable.sets(s$dag.full, d = d, num.int = num.int)
  
  # find the indeces of the invariant sets in the list of all subsets
  ind.inv.sets <- (1:length(ps))[ps %in% inv.sets]
  
  # find the rows corresponding to the invariant sets
  pvals.inv <- pvals[ind.inv.sets, ]
  
  # append the pvalues for the invariant sets to the vectors
  null.pvals.delong.rf.randDAG <- c(null.pvals.delong.rf.randDAG, pvals.inv$delong.rf)
  null.pvals.delong.glm.randDAG <- c(null.pvals.delong.glm.randDAG, pvals.inv$delong.glm)
  null.pvals.tram.rf.randDAG <- c(null.pvals.tram.rf.randDAG, pvals.inv$tram.rf)
  null.pvals.tram.glm.randDAG <- c(null.pvals.tram.glm.randDAG, pvals.inv$tram.glm)
  null.pvals.corr.randDAG <- c(null.pvals.corr.randDAG, pvals.inv$correlation)
  null.pvals.residual.randDAG <- c(null.pvals.residual.randDAG, pvals.inv$residual)
  
}




#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# run the simulation for the semirandom SCM
#-------------------------------------------------------------------------------


set.seed(1)


# in these vectors, we pool together the pvalues for all invariant sets across all simulations
null.pvals.delong.rf.semirand <- numeric(0)
null.pvals.delong.glm.semirand <- numeric(0)
null.pvals.tram.rf.semirand <- numeric(0)
null.pvals.tram.glm.semirand <- numeric(0)
null.pvals.corr.semirand <- numeric(0)
null.pvals.residual.semirand <- numeric(0)



for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # generate a sample (n.test is arbitrary, we just use "training" data)
  s <- generate.samples.semirandom(n = n, n.test = 1000, t = 1)
  sample <- s$sample_train
  
  # compute p-values
  pvals <- pvalues.all(sample)
  
  # find invariant sets from DAG 
  inv.sets <- int.stable.sets(s$dag.full, d = d, num.int = num.int)
  
  # find the indeces of the invariant sets in the list of all subsets
  ind.inv.sets <- (1:length(ps))[ps %in% inv.sets]
  
  # find the rows corresponding to the invariant sets
  pvals.inv <- pvals[ind.inv.sets, ]
  
  # append the pvalues for the invariant sets to the vectors
  null.pvals.delong.rf.semirand <- c(null.pvals.delong.rf.semirand, pvals.inv$delong.rf)
  null.pvals.delong.glm.semirand <- c(null.pvals.delong.glm.semirand, pvals.inv$delong.glm)
  null.pvals.tram.rf.semirand <- c(null.pvals.tram.rf.semirand, pvals.inv$tram.rf)
  null.pvals.tram.glm.semirand <- c(null.pvals.tram.glm.semirand, pvals.inv$tram.glm)
  null.pvals.corr.semirand <- c(null.pvals.corr.semirand, pvals.inv$correlation)
  null.pvals.residual.semirand <- c(null.pvals.residual.semirand, pvals.inv$residual)
  
}


save(null.pvals.delong.rf.randDAG, 
     null.pvals.delong.glm.randDAG, 
     null.pvals.tram.rf.randDAG,
     null.pvals.tram.glm.randDAG,
     null.pvals.corr.randDAG,
     null.pvals.residual.randDAG,
     null.pvals.delong.rf.semirand, 
     null.pvals.delong.glm.semirand, 
     null.pvals.tram.rf.semirand,
     null.pvals.tram.glm.semirand,
     null.pvals.corr.semirand,
     null.pvals.residual.semirand,
     file = "ecdf_both_randomSCMs.rdata")

#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# plotting
#-------------------------------------------------------------------------------


# prepare data for plotting
df.invariant.ecdf.randDAG <- data.frame(
  x = c(null.pvals.delong.rf.randDAG, null.pvals.delong.glm.randDAG, null.pvals.tram.rf.randDAG, null.pvals.tram.glm.randDAG, null.pvals.corr.randDAG, null.pvals.residual.randDAG),
  g = gl(n = 6, k = length(null.pvals.delong.rf.randDAG))
)


df.invariant.ecdf.semirand <- data.frame(
  x = c(null.pvals.delong.rf.semirand, null.pvals.delong.glm.semirand, null.pvals.tram.rf.semirand, null.pvals.tram.glm.semirand, null.pvals.corr.semirand, null.pvals.residual.semirand),
  g = gl(n = 6, k = length(null.pvals.delong.rf.semirand))
)


size <- 10

p.randDAG <- ggplot(df.invariant.ecdf.randDAG, aes(x, colour = g)) +
  stat_ecdf() +
  geom_abline() +
  scale_color_hue(labels=c('DeLong (RF)', 'DeLong (GLM)', 'TRAM-GCM (RF)', 'TRAM-GCM (GLM)', 'Correlation', 'Residual')) +
  labs(color='') +
  xlab("p-value") +
  ylab("Empirical CDF") +
  ggtitle("Random SCM") +
  theme_bw(base_size = size) 
#p.randDAG

p.semirand <- ggplot(df.invariant.ecdf.semirand, aes(x, colour = g)) +
  stat_ecdf() +
  geom_abline() +
  scale_color_hue(labels=c('DeLong (RF)', 'DeLong (GLM)', 'TRAM-GCM (RF)', 'TRAM-GCM (GLM)', 'Correlation', 'Residual')) +
  labs(color='') +
  xlab("p-value") +
  ylab("Empirical CDF") +
  ggtitle("Semi-Random SCM") +
  theme_bw(base_size = size) 
#p.semirand




combined <- p.semirand + p.randDAG & theme(legend.position = "bottom", legend.title=element_blank(), legend.text = element_text(size=size)) 

combined + plot_layout(guides = "collect")

ggsave(filename = "ecdf_both_randomSCMs.pdf", width = 6, height = 4)




#-------------------------------------------------------------------------------








