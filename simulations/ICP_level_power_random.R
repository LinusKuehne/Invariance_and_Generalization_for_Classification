# In this script, we investigate the effect of intervention strength and sample size 
# on the performance of ICP with regard to FWER and Jaccard index for data from the
# random SCM



library(ggplot2)
library(patchwork)
library(rje)
library(reshape2)


# get the path of this script
script_dir <- getwd()

# load in the functions needed
source(file.path(script_dir, "../code/code_simulations/invariance_tests.R"))
source(file.path(script_dir, "../code/code_simulations/data_generating_process.R"))
source(file.path(script_dir, "../code/code_simulations/utils.R"))






#-------------------------------------------------------------------------------
# parameters shared by both experiments
#-------------------------------------------------------------------------------

# number of variables (including Y, excluding I)
d <- 6

# number of interventions
num.int <- 2 

# max number of parents
max.pa <- 3 

# 0 encodes the empty set
sets <- powerSet(1:(d-1))
sets[[1]] <- 0


# number of simulation repetitions (nreps = 300 would take roughly 30 h for both experiments combined)
nreps <- 500


icp.methods <- c("delong.rf", "delong.glm", "tram.rf", "tram.glm", "correlation", "residual", "oracle")


#-------------------------------------------------------------------------------




################################################################################
# Experiment 1: Investigate effect of intervention strength on ICP
################################################################################




#-------------------------------------------------------------------------------
# parameters for experiment 1
#-------------------------------------------------------------------------------


# fixed number of n observations per environment => 5*n observations in total
n <- 250

# number of different intervention strengths to consider
nt <- 10

# vector of intervention strength scalings between 0 and 1
int.strength <- seq(from=0, to=1, length.out = nt)


#-------------------------------------------------------------------------------








#-------------------------------------------------------------------------------
# run experiment 1
#-------------------------------------------------------------------------------


set.seed(1)



# store the jaccard scores and family-wise errors for each intervention strength
df.init <- as.data.frame(matrix(0, nrow = nreps, ncol = length(icp.methods)))
names(df.init) <- icp.methods
list.jaccard.int <- list.FWE.int <- replicate(nt, df.init, simplify=FALSE)




for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  for(i in 1:nt){
    
    # set intervention strength
    t <- int.strength[i]
    
    # generate data from the random SCM
    s <- generate.samples.random(n = n, n.test = 1000, d = d, max.pa = max.pa, num.int = num.int, t = t)
    
    # extract training data
    sample <- s$sample_train
    
    # compute the pvalues for all sets using all invariance tests
    pvals.sample <- pvalues.all(sample)
    
    # find variable names of parents
    pars <- parents(x = s$dag.cov, v = "Y")
    
    # convert to a vector of variable indeces
    parents <- numeric(0)
    if(length(pars)>0){
      parents <- unname(sapply(X = pars, function(x) { as.numeric(substr(x, start = 2, stop = 2))}))
    }
    
    # compute the ICP oracle output in terms of covariate names for comparison
    s.ICP <- ICP.oracle(s$dag.cov, s$dag.full, num.int)
    
    # convert to a vector of variable indeces
    S.ICP <- numeric(0)
    if(length(s.ICP) > 0){
      S.ICP <- unname(sapply(X = s.ICP, function(x) { as.numeric(substr(x, start = 2, stop = 2))}))
    }
    
    # compute the ICP output from the p-values (i.e., intersection of all invariant sets)
    icp.output <- apply(X = pvals.sample, MARGIN = 2, FUN = calc.ICP.output, simplify = F)
    
    # append the oracle output
    icp.output$oracle <- S.ICP
    
    # compute family-wise errors
    FWE.i <- sapply(X = icp.output, FUN = FWE, parents = parents)
    list.FWE.int[[i]][b, ] <- FWE.i
    
    # compute jaccard index
    jaccard.i <- sapply(X = icp.output, FUN = jaccard.ICP, parents = parents)
    list.jaccard.int[[i]][b, ] <- jaccard.i
  }
}

save(list.FWE.int, list.jaccard.int, file = file.path(script_dir, "saved_data/ICP_level_power_random_int.rdata"))


#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# plotting
#-------------------------------------------------------------------------------


# compute the family-wise error RATE and standard deviations
FWER.list.int <- lapply(X = list.FWE.int, FUN = colMeans)
FWER.sd.list.int <- lapply(X = list.FWE.int, function(x) { apply(X = x, MARGIN = 2, FUN = sd) })

# convert to dataframe
df.FWER.int <- data.frame(matrix(unlist(FWER.list.int), nrow = nt, byrow = T))
names(df.FWER.int) <- icp.methods
df.FWER.int$intervention.strength <- int.strength

# convert to dataframe
df.FWER.sd.int <- data.frame(matrix(unlist(FWER.sd.list.int), nrow = nt, byrow = T))
names(df.FWER.sd.int) <- icp.methods
df.FWER.sd.int$intervention.strength <- int.strength



# compute the mean and standard devitation of Jaccard index
jaccard.means.list.int <- lapply(X = list.jaccard.int, FUN = colMeans)
jaccard.sd.list.int <- lapply(X = list.jaccard.int, function(x) { apply(X = x, MARGIN = 2, FUN = sd) })

# convert to dataframe
df.jaccard.int <- data.frame(matrix(unlist(jaccard.means.list.int), nrow = nt, byrow = T))
names(df.jaccard.int) <- icp.methods
df.jaccard.int$intervention.strength <- int.strength

# convert to dataframe
df.jaccard.sd.int <- data.frame(matrix(unlist(jaccard.sd.list.int), nrow = nt, byrow = T))
names(df.jaccard.sd.int) <- icp.methods
df.jaccard.sd.int$intervention.strength <- int.strength




# prepare data for plotting
df.plot.jaccard.int <- melt(df.jaccard.int,  id.vars = "intervention.strength", variable.name = 'Candidate')
df.plot.FWER.int <- melt(df.FWER.int,  id.vars = "intervention.strength", variable.name = 'Candidate')

sd.plot.jaccard.int <- melt(df.jaccard.sd.int, id.vars = "intervention.strength", variable.name = 'Candidate')
sd.plot.FWER.int <- melt(df.FWER.sd.int, id.vars = "intervention.strength", variable.name = 'Candidate')


# compute t-test confidence intervals
df.plot.jaccard.int$ci <- sd.plot.jaccard.int$value * qt(0.95, df = nreps - 1) / sqrt(nreps)
df.plot.FWER.int$ci <- sd.plot.FWER.int$value * qt(0.95, df = nreps - 1) / sqrt(nreps)





color.palette <- c(scales::hue_pal()(6), "black")
names(color.palette) <- icp.methods



alpha <- 0.3
size <- 10


p.FWER.int <- ggplot(df.plot.FWER.int, aes(intervention.strength, value)) +
  geom_line(aes(colour = Candidate)) +
  #geom_ribbon(aes(ymin = value - ci, ymax = value + ci, fill = Candidate), alpha = alpha, show.legend = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "black", linewidth = 0.75) +
  scale_color_manual(values = color.palette, labels=c('DeLong (RF)', 'DeLong (GLM)', 'TRAM-GCM (RF)', 'TRAM-GCM (GLM)', 'Correlation', 'Residual', 'Oracle')) +
  scale_fill_manual(values = color.palette) +
  xlab("Intervention strength") +
  ylab("Average FWER") +
  ggtitle("Family-Wise Error Rate") +
  theme_bw(base_size = size)

#p.FWER.int



p.jaccard.int <- ggplot(df.plot.jaccard.int, aes(intervention.strength, value)) +
  geom_line(aes(colour = Candidate)) +
  #geom_ribbon(aes(ymin = value - ci, ymax = value + ci, fill = Candidate), alpha = alpha, show.legend = F) +
  scale_color_manual(values = color.palette, labels=c('DeLong (RF)', 'DeLong (GLM)', 'TRAM-GCM (RF)', 'TRAM-GCM (GLM)', 'Correlation', 'Residual', 'Oracle')) +
  scale_fill_manual(values = color.palette) +
  xlab("Intervention strength") +
  ylab("Average Jaccard index") +
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Jaccard Index") +
  theme_bw(base_size = size)

#p.jaccard.int








combined.int <- p.FWER.int + p.jaccard.int & theme(legend.position = "bottom", legend.title=element_blank(), legend.text = element_text(size=size))

combined.int + plot_layout(guides = "collect")


ggsave(filename = file.path(script_dir, "saved_plots/ICP_level_power_random_int.pdf"), width = 6, height = 4)

#-------------------------------------------------------------------------------
















################################################################################
# Experiment 2: Investigate effect of sample size on ICP
################################################################################





#-------------------------------------------------------------------------------
# parameters for experiment 2
#-------------------------------------------------------------------------------

# how many different sample sizes we want to consider
nt <- 10

# this is the vector of different samples per environment => in total, we have 5 times that number
sample.size <- floor(seq(from = 20, to = 250, length = nt))


# fixed intervention strength:
t <- 1

#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# run experiment 2
#-------------------------------------------------------------------------------


set.seed(1)


# store the jaccard scores and family-wise errors for each intervention strength
df.init <- as.data.frame(matrix(0, nrow = nreps, ncol = length(icp.methods)))
names(df.init) <- icp.methods

list.jaccard.samp <- list.FWE.samp <- replicate(nt, df.init, simplify=FALSE)



for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  for(i in 1:nt){
    
    # set sample size
    n <- sample.size[i]
    
    # generate data from the random SCM
    s <- generate.samples.random(n = n, n.test = 1000, d = d, max.pa = max.pa, num.int = num.int, t = t)
    
    # extract training data
    sample <- s$sample_train
    
    # compute the pvalues for all sets using all invariance tests
    pvals.sample <- pvalues.all(sample)
    
    # find variable names of parents
    pars <- parents(x = s$dag.cov, v = "Y")
    
    # convert to a vector of variable indeces
    parents <- numeric(0)
    if(length(pars)>0){
      parents <- unname(sapply(X = pars, function(x) { as.numeric(substr(x, start = 2, stop = 2))}))
    }
    
    # compute the ICP oracle output in terms of covariate names for comparison
    s.ICP <- ICP.oracle(s$dag.cov, s$dag.full, num.int)
    
    # convert to a vector of variable indeces
    S.ICP <- numeric(0)
    if(length(s.ICP) > 0){
      S.ICP <- unname(sapply(X = s.ICP, function(x) { as.numeric(substr(x, start = 2, stop = 2))}))
    }
    
    # compute the ICP output from the p-values (i.e., intersection of all invariant sets)
    icp.output <- apply(X = pvals.sample, MARGIN = 2, FUN = calc.ICP.output, simplify = F)
    
    # append the oracle output
    icp.output$oracle <- S.ICP
    
    # compute family-wise errors
    FWE.i <- sapply(X = icp.output, FUN = FWE, parents = parents)
    list.FWE.samp[[i]][b, ] <- FWE.i
    
    # compute jaccard index
    jaccard.i <- sapply(X = icp.output, FUN = jaccard.ICP, parents = parents)
    list.jaccard.samp[[i]][b, ] <- jaccard.i
  }
}


save(list.FWE.samp, list.jaccard.samp, file = file.path(script_dir, "saved_data/ICP_level_power_random_samp.rdata"))


#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# plotting
#-------------------------------------------------------------------------------



# compute the family-wise error RATE and standard deviations
FWER.list.samp <- lapply(X = list.FWE.samp, FUN = colMeans)
FWER.sd.list.samp <- lapply(X = list.FWE.samp, function(x) { apply(X = x, MARGIN = 2, FUN = sd) })

# convert to dataframe
df.FWER.samp <- data.frame(matrix(unlist(FWER.list.samp), nrow = nt, byrow = T))
names(df.FWER.samp) <- icp.methods
df.FWER.samp$sample.size <- sample.size

# convert to dataframe
df.FWER.sd.samp <- data.frame(matrix(unlist(FWER.sd.list.samp), nrow = nt, byrow = T))
names(df.FWER.sd.samp) <- icp.methods
df.FWER.sd.samp$sample.size <- sample.size




# compute the mean and standard devitation of Jaccard index
jaccard.means.list.samp <- lapply(X = list.jaccard.samp, FUN = colMeans)
jaccard.sd.list.samp <- lapply(X = list.jaccard.samp, function(x) { apply(X = x, MARGIN = 2, FUN = sd) })

# convert to dataframe
df.jaccard.samp <- data.frame(matrix(unlist(jaccard.means.list.samp), nrow = nt, byrow = T))
names(df.jaccard.samp) <- icp.methods
df.jaccard.samp$sample.size <- sample.size

# convert to dataframe
df.jaccard.sd.samp <- data.frame(matrix(unlist(jaccard.sd.list.samp), nrow = nt, byrow = T))
names(df.jaccard.sd.samp) <- icp.methods
df.jaccard.sd.samp$sample.size <- sample.size



# prepare data for plotting
df.plot.jaccard.samp <- melt(df.jaccard.samp,  id.vars = "sample.size", variable.name = 'Candidate')
df.plot.FWER.samp <- melt(df.FWER.samp,  id.vars = "sample.size", variable.name = 'Candidate')

sd.plot.jaccard.samp <- melt(df.jaccard.sd.samp,  id.vars = "sample.size", variable.name = 'Candidate')
sd.plot.FWER.samp <- melt(df.FWER.sd.samp,  id.vars = "sample.size", variable.name = 'Candidate')


# compute t-test confidence intervals
df.plot.jaccard.samp$ci <- sd.plot.jaccard.samp$value * qt(0.95, df = nreps - 1) / sqrt(nreps)
df.plot.FWER.samp$ci <- sd.plot.FWER.samp$value * qt(0.95, df = nreps - 1) / sqrt(nreps)

df.plot.jaccard.samp$sample.size <- 5*df.plot.jaccard.samp$sample.size
df.plot.FWER.samp$sample.size <- 5*df.plot.FWER.samp$sample.size





color.palette <- c(scales::hue_pal()(6), "black")
names(color.palette) <- icp.methods



alpha <- 0.3
size <- 10


p.FWER.samp <- ggplot(df.plot.FWER.samp, aes(sample.size, value)) +
  geom_line(aes(colour = Candidate)) +
  #geom_ribbon(aes(ymin = value - ci, ymax = value + ci, fill = Candidate), alpha = alpha, show.legend = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "black", linewidth = 0.75) +
  scale_color_manual(values = color.palette, labels=c('DeLong (RF)', 'DeLong (GLM)', 'TRAM-GCM (RF)', 'TRAM-GCM (GLM)', 'Correlation', 'Residual', 'Oracle')) +
  scale_fill_manual(values = color.palette) +
  xlab("Sample size") +
  ylab("Average FWER") +
  ggtitle("Family-Wise Error Rate") +
  theme_bw(base_size = size)

#p.FWER.samp



p.jaccard.samp <- ggplot(df.plot.jaccard.samp, aes(sample.size, value)) +
  geom_line(aes(colour = Candidate)) +
  #geom_ribbon(aes(ymin = value - ci, ymax = value + ci, fill = Candidate), alpha = alpha, show.legend = F) +
  scale_color_manual(values = color.palette, labels=c('DeLong (RF)', 'DeLong (GLM)', 'TRAM-GCM (RF)', 'TRAM-GCM (GLM)', 'Correlation', 'Residual', 'Oracle')) +
  scale_fill_manual(values = color.palette) +
  xlab("Sample size") +
  ylab("Average Jaccard index") +
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Jaccard Index") +
  theme_bw(base_size = size)

#p.jaccard.samp




combined.samp <- p.FWER.samp + p.jaccard.samp & theme(legend.position = "bottom", legend.title=element_blank(), legend.text = element_text(size=size))

combined.samp + plot_layout(guides = "collect")

ggsave(filename = file.path(script_dir, "saved_plots/ICP_level_power_random_samp.pdf"), width = 6, height = 4)

#-------------------------------------------------------------------------------


# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/ICP_level_power_random.txt"))

