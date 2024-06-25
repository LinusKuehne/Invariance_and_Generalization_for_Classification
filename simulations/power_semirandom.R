# In this script, we conduct two experiments with data from the semi-random SCM:
# First, for a fixed sample size, we investigate the effect of the intervention strength 
# (by scaling the intervention values by a factor t in (0,1)) on the surrogate measure Jaccard index
# between the invariant sets found by the invariance tests and the ground truth invariant sets. 
# Second, in a very similar setup, we fix the intervention strength and vary the sample size.




library(ggplot2)
library(reshape2)
library(patchwork)
library(pROC)
library(rje)


source("DGP.R")
source("invariance_tests.R")
source("utils.R")





#-------------------------------------------------------------------------------
# parameters shared by both experiments
#-------------------------------------------------------------------------------

# number of variables (including Y, excluding I)
d <- 6

# number of interventions
num.int <- 2 

# 0 encodes the empty set
sets <- powerSet(1:(d-1))
sets[[1]] <- 0


# number of simulation repetitions (nreps = 300 would take roughly 30 h for both experiments combined)
nreps <- 300


#-------------------------------------------------------------------------------






################################################################################
# Experiment 1: Investigate effect of intervention strength on power
################################################################################




#-------------------------------------------------------------------------------
# parameters for experiment 1
#-------------------------------------------------------------------------------


# fixed number of n observations per environment => 5*n observations in total
n <- 100

# number of different intervention strengths to consider
nt <- 10

# vector of intervention strength scalings between 0 and 1
int.strength <- seq(from=0, to=1, length.out = nt)


#-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------
# run experiment 1
#-------------------------------------------------------------------------------

set.seed(1)

# matrix used to initialize list for the scores
mat.init <- as.data.frame(matrix(0, nrow = nreps, ncol = nt))



# in this list we store the scores for each invariance test and each intervention strength
list.jaccard.int <- list(mat.delong.rf = mat.init,
                     mat.delong.glm = mat.init,
                     mat.tram.rf = mat.init,
                     mat.tram.glm = mat.init,
                     mat.correlation = mat.init,
                     mat.residual = mat.init)
                     



for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  
  for(i in 1:nt){
    
    # intervention strength scaling
    t <- int.strength[i]
    
    # generate data with this intervention strength
    s <- generate.samples.semirandom(n = n, n.test = 1000, t = t)
    
    # extract training sample
    sample <- s$sample_train
    
    # compute p-values for all sets using all invariance tests
    pvals.sample <- pvalues.all(sample)
    
    # find sets which satisfy the d-separation associated with invariance from the DAG
    stable.sets <- int.stable.sets(dag.full = s$dag.full, d = d, num.int = num.int)
    
    # compute the Jaccard index
    score.j <- apply(X = pvals.sample, MARGIN = 2, FUN = jaccard, thresh = 0.05, stable.sets = stable.sets, num.covariates = (d-1))

    # store these appropriately
    for(c in 1:ncol(pvals.sample)){
      list.jaccard.int[[c]][b,i] <- score.j[c]
    }
    
  }
}



save(list.jaccard.int, file = "power_semirandom_int.rdata")

#-------------------------------------------------------------------------------








#-------------------------------------------------------------------------------
# Plotting
#-------------------------------------------------------------------------------



# prepare data for plotting
df.jaccard.int <- data.frame(delong.rf = numeric(nt),
                             delong.glm = numeric(nt),
                             tram.rf = numeric(nt),
                             tram.glm = numeric(nt),
                             correlation = numeric(nt),
                             residual = numeric(nt),
                             intervention.strength = int.strength)


df.jaccard.sd.int <- df.jaccard.int



# compute mean and standard deviations of Jaccard index across simulations
for(u in 1:length(list.jaccard.int)){
  
  df.jaccard.int[,u] <- colMeans(list.jaccard.int[[u]])

  df.jaccard.sd.int[,u] <- apply(X = list.jaccard.int[[u]], MARGIN = 2, FUN = sd)

}




# bring data into a better format for plotting
df.plot.jaccard.int <- melt(df.jaccard.int,  id.vars = "intervention.strength", variable.name = 'Candidate')
sd.plot.jaccard.int <- melt(df.jaccard.sd.int,  id.vars = "intervention.strength", variable.name = 'Candidate')


# compute t-test confidence intervals
df.plot.jaccard.int$ci <- sd.plot.jaccard.int$value * qt(0.975, df = nreps - 1) / sqrt(nreps)







size <- 10
alpha <- 0.3





p.jaccard.int <- ggplot(df.plot.jaccard.int, aes(intervention.strength, value)) +
  scale_color_hue(labels=c('DeLong (RF)', 'DeLong (GLM)', 'TRAM-GCM (RF)', 'TRAM-GCM (GLM)', 'Correlation', 'Residual')) +
  geom_line(aes(colour = Candidate)) +
  geom_ribbon(aes(ymin = value - ci, ymax = value + ci, fill = Candidate), alpha = alpha, show.legend = F) +
  coord_cartesian(ylim=c(0.2,1)) +
  xlab("Intervention strength") +
  ylab("Average Jaccard index") +
  ggtitle("Fixed Sample Size") +
  theme_bw(base_size = size) 

#p.jaccard.int


#-------------------------------------------------------------------------------










################################################################################
# Experiment 2: Investigate effect of sample size on power
################################################################################




#-------------------------------------------------------------------------------
# parameters for experiment 2
#-------------------------------------------------------------------------------

# how many different sample sizes we want to consider
nt <- 10

# this is the vector of different samples per environment => in total, we have 5 times that number
sample.size <- floor(seq(from = 20, to = 200, length = nt))


# fixed intervention strength:
t <- 0.5

#-------------------------------------------------------------------------------








#-------------------------------------------------------------------------------
# run experiment 2
#-------------------------------------------------------------------------------

set.seed(1)


# matrix used to initialize lists for the scores
mat.init <- as.data.frame(matrix(0, nrow = nreps, ncol = nt))


# in this list we store the scores for each invariance test and each sample size
list.jaccard.samp <- list(mat.delong.rf = mat.init,
                          mat.delong.glm = mat.init,
                          mat.tram.rf = mat.init,
                          mat.tram.glm = mat.init,
                          mat.correlation = mat.init,
                          mat.residual = mat.init)




for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  
  for(i in 1:nt){
    
    # determine the sample size
    n <- sample.size[i]
    
    # generate data with this sample size
    s <- generate.samples.semirandom(n = n, n.test = 1000, t = t)
    
    # extract training sample
    sample <- s$sample_train
    
    # compute p-values for all sets using all invariance tests
    pvals.sample <- pvalues.all(sample)
    
    # find sets which satisfy the d-separation associated with invariance from the DAG
    stable.sets <- int.stable.sets(dag.full = s$dag.full, d = d, num.int = num.int)
    
    # compute Jaccard index
    score.j <- apply(X = pvals.sample, MARGIN = 2, FUN = jaccard, thresh = 0.05, stable.sets = stable.sets, num.covariates = (d-1))

    # save these values appropriately
    for(c in 1:ncol(pvals.sample)){
      list.jaccard.samp[[c]][b,i] <- score.j[c]
    }
    
  }
}



save(list.jaccard.samp, file = "power_semirandom_samp.rdata")
#-------------------------------------------------------------------------------









#-------------------------------------------------------------------------------
# Plotting
#-------------------------------------------------------------------------------




# prepare data for plotting
df.jaccard.samp <- data.frame(delong.rf = numeric(nt),
                              delong.glm = numeric(nt),
                              tram.rf = numeric(nt),
                              tram.glm = numeric(nt),
                              correlation = numeric(nt),
                              residual = numeric(nt),
                              sample.size = sample.size)


df.jaccard.sd.samp <- df.jaccard.samp


# compute mean and standard deviations of Jaccard index across simulations 
for(u in 1:length(list.jaccard.samp)){
  
  df.jaccard.samp[,u] <- colMeans(list.jaccard.samp[[u]])

  df.jaccard.sd.samp[,u] <- apply(X = list.jaccard.samp[[u]], MARGIN = 2, FUN = sd)

}




# bring data into a better format for plotting
df.plot.jaccard.samp <- melt(df.jaccard.samp,  id.vars = "sample.size", variable.name = 'Candidate')
sd.plot.jaccard.samp <- melt(df.jaccard.sd.samp,  id.vars = "sample.size", variable.name = 'Candidate')

# compute t-test confidence intervals
df.plot.jaccard.samp$ci <- sd.plot.jaccard.samp$value * qt(0.975, df = nreps - 1) / sqrt(nreps)







size <- 10
alpha <- 0.3



# sample.size is the sample size PER ENVIRONMENT
# For the total sample size: multiply by 5
df.plot.jaccard.samp$sample.size <- 5*df.plot.jaccard.samp$sample.size


p.jaccard.samp <- ggplot(df.plot.jaccard.samp, aes(sample.size, value)) +
  scale_color_hue(labels=c('DeLong (RF)', 'DeLong (GLM)', 'TRAM-GCM (RF)', 'TRAM-GCM (GLM)', 'Correlation', 'Residual')) +
  geom_line(aes(colour = Candidate)) +
  geom_ribbon(aes(ymin = value - ci, ymax = value + ci, fill = Candidate), alpha = alpha, show.legend = F) +
  coord_cartesian(ylim=c(0.3,1)) +
  xlab("Sample size") +
  ylab("Average Jaccard index") +
  ggtitle("Fixed Intervention Strength") +
  theme_bw(base_size = size)

#p.jaccard.samp


#-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------
# combine the plots from both experiments 
#-------------------------------------------------------------------------------

combined <- p.jaccard.int + p.jaccard.samp & theme(legend.position = "bottom", legend.title=element_blank(), legend.text = element_text(size=size))

combined + plot_layout(guides = "collect")

ggsave(filename = "power_semirandom.pdf", width = 6, height = 4)

#-------------------------------------------------------------------------------



