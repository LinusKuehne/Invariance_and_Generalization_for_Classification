# In this script, we generate data from the standard DAG and compute the p-values
# using all of our invariance tests for the invariant sets S = {1} and S = {1,3}. 
# We then plot the empirical CDF (ECDF) of these p-values. If the tests are level,
# the ECDFs should lie on or below the diagonal.


library(ggplot2)
library(patchwork)


# get the path of this script
script_dir <- getwd()

# load in the functions needed
source(file.path(script_dir, "../code/code_simulations/invariance_tests.R"))
source(file.path(script_dir, "../code/code_simulations/data_generating_process.R"))
source(file.path(script_dir, "../code/code_simulations/utils.R"))




#-------------------------------------------------------------------------------
# parameters
#-------------------------------------------------------------------------------

set.seed(1)

# these are the possible subsets of predictors 1:3
sets <- list(c(1,2,3), c(1,3), c(1,2), c(2,3), c(1), c(2), c(3))


# number of simulation repetitions
# 500 nreps should take around one hour on my MB
nreps <- 3


# n observations per environment => 5*n observations in total
n <- 200

#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# run the simulation
#-------------------------------------------------------------------------------



# here we store the p-values for the set S = {1} computed by all invariance tests
df.S1 <- data.frame(delong.rf = numeric(nreps),
                    delong.glm = numeric(nreps),
                    tram.rf = numeric(nreps),
                    tram.glm = numeric(nreps),
                    corr = numeric(nreps),
                    residual = numeric(nreps))

# here we store the p-values for the set S = {1,3} computed by all invariance tests
df.S13 <- data.frame(delong.rf = numeric(nreps),
                     delong.glm = numeric(nreps),
                     tram.rf = numeric(nreps),
                     tram.glm = numeric(nreps),
                     corr = numeric(nreps),
                     residual = numeric(nreps))



set.seed(1)


for(b in 1:nreps){
  
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # obtain a sample of 5*n observations with n observations from each of five environments
  sample <- gen.sample.standard(n)
  
  # we compute TRAM-GCM (RF) and (GLM) separately
  icp.ranger <- rangerICP(Y ~ X1 + X2 + X3, data = sample, env = ~ Env, test = "gcm.test", verbose = F)
  res.ranger <- pvalues(icp.ranger, "set")
  
  icp.glm <- glmICP(Y ~ X1 + X2 + X3, data = sample, env = ~ Env, family = "binomial", verbose = F)
  res.glm <- pvalues(icp.glm, "set")
  
  # compute invariance tests for S = {1}
  df.S1$delong.rf[b] <- pval.delong.rf.set(sample, 1)
  df.S1$delong.glm[b] <- pval.delong.glm.set(sample, 1)
  df.S1$tram.rf[b] <- res.ranger[2]
  df.S1$tram.glm[b] <- res.glm[2]
  df.S1$corr[b] <- pval.correlation.set(sample, 1)
  df.S1$residual[b] <- pval.residual.set(sample, 1)
  
  # compute invariance tests for S = {1,3}
  df.S13$delong.rf[b] <- pval.delong.rf.set(sample, c(1,3))
  df.S13$delong.glm[b] <- pval.delong.glm.set(sample, c(1,3))
  df.S13$tram.rf[b] <- res.ranger[6]
  df.S13$tram.glm[b] <- res.glm[6]
  df.S13$corr[b] <- pval.correlation.set(sample, c(1,3))
  df.S13$residual[b] <- pval.residual.set(sample, c(1,3))
  
}

save(df.S1, df.S13, file = file.path(script_dir, "/results/ecdf_standard.rdata"))

#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# plotting
#-------------------------------------------------------------------------------



# prepare data for plotting
df.S1.ecdf <- data.frame(
  x = c(df.S1$delong.rf, df.S1$delong.glm, df.S1$tram.rf, df.S1$tram.glm, df.S1$corr, df.S1$residual),
  g = gl(n = ncol(df.S1), k = nreps)
)

df.S13.ecdf <- data.frame(
  x = c(df.S13$delong.rf, df.S13$delong.glm, df.S13$tram.rf, df.S13$tram.glm, df.S13$corr, df.S13$residual),
  g = gl(n = ncol(df.S13), k = nreps)
)


size <- 10

p1 <- ggplot(df.S1.ecdf, aes(x, colour = g)) +
  stat_ecdf() +
  geom_abline() +
  scale_color_hue(labels=c('DeLong (RF)', 'DeLong (GLM)', 'TRAM-GCM (RF)', 'TRAM-GCM (GLM)', 'Correlation', 'Residual')) +
  labs(color='') +
  xlab("p-value") +
  ylab("Empirical CDF") +
  ggtitle("Subset S = {1}") +
  theme_bw(base_size = size)
#p1




p13 <- ggplot(df.S13.ecdf, aes(x, colour = g)) +
  stat_ecdf() +
  geom_abline() +
  scale_color_hue(labels=c('DeLong (RF)', 'DeLong (GLM)', 'TRAM-GCM (RF)', 'TRAM-GCM (GLM)', 'Correlation', 'Residual')) +
  labs(color='') +
  xlab("p-value") +
  ylab("Empirical CDF") +
  ggtitle("Subset S = {1,3}") +
  theme_bw(base_size = size) 

#p13

combined <- p13 + p1 & theme(legend.position = "bottom", legend.title=element_blank(), legend.text = element_text(size=size)) 

combined + plot_layout(guides = "collect")

ggsave(filename = file.path(script_dir, "results/ecdf_standard.pdf"), width = 6, height = 4)
  


#-------------------------------------------------------------------------------




writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/ecdf_standard.txt"))
















