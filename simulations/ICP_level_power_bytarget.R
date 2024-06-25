# In this script, we compute level and power of ICP for the semirandom SCM using 
# the different invariance tests for the logistic regression model. Here, we now
# stratify according to which covariate is chosen as the response, to investigate 
# whether invariance is more difficult to test for some responses.


library(ggplot2)
library(rje)
library(patchwork)


source("utils.R")
source("DGP.R")
source("invariance_tests.R")



#-------------------------------------------------------------------------------
# parameters
#-------------------------------------------------------------------------------

# number of variables (including Y, excluding I)
d <- 6

# number of interventions
num.int <- 2 

# 0 encodes the empty set
sets <- powerSet(1:(d-1))
sets[[1]] <- 0

# number of simulation repetitions (nreps = 300 would take roughly 7h for both experiments combined)
nreps <- 20

# number of observations per environment => 5*n total observations
n <- 150

# intervention strength scaling
t <- 1


# names of the predictors
Xnames <- c("X1", "X2", "X3", "X4", "X5")

ps <- powerSet(Xnames)

# for plotting
size <- 10

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# run the simulation
#-------------------------------------------------------------------------------

set.seed(1)



# store the Jaccard index scores
scores.j.mod1 <- data.frame("delongRF" = numeric(nreps),
                            "delongGLM" = numeric(nreps),
                            "tramRF" = numeric(nreps),
                            "tramGLM" = numeric(nreps),
                            "correlation" = numeric(nreps),
                            "residual" = numeric(nreps),
                            "oracle" = numeric(nreps))

# store the family-wise errors
FWE.mod1 <- scores.j.mod1


# store which variable was chosen as the response/target
target <- numeric(nreps)




for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # generate data from the semirandom SCM
  s <- generate.samples.semirandom(n = n, n.test = 1000, t = t)
  
  # store the response
  target[b] <- s$response
  
  # extract the training sample
  sample <- s$sample_train
  
  # compute the p-values for all sets using all invariance tests
  pvals.sample <- pvalues.all(sample)
  
  # extract the parents
  pars <- parents(x = s$dag.cov, v = "Y")
  
  # convert pars to the corresponding vector of covariate indeces
  parents <- numeric(0)
  if(length(pars)>0){
    parents <- unname(sapply(X = pars, function(x) { as.numeric(substr(x, start = 2, stop = 2))}))
  }
  
  # compute oracle ICP output
  s.ICP <- ICP.oracle(s$dag.cov, s$dag.full, num.int)
  
  # convert oracle ICP output to the corresponding vector of covariate indeces
  S.ICP <- numeric(0)
  if(length(s.ICP) > 0){
    S.ICP <- unname(sapply(X = s.ICP, function(x) { as.numeric(substr(x, start = 2, stop = 2))}))
  }
  
  # compute ICP output (intersection over all invariant sets) based on the pvalues
  icp.output <- apply(X = pvals.sample, MARGIN = 2, FUN = calc.ICP.output, simplify = F)
  
  # append oracle ICP output
  icp.output$oracle <- S.ICP
  
  # compute family-wise errors
  FWE.i <- sapply(X = icp.output, FUN = FWE, parents = parents)
  FWE.mod1[b, ] <- FWE.i
  
  # compute jaccard index
  jaccard.i <- sapply(X = icp.output, FUN = jaccard.ICP, parents = parents)
  scores.j.mod1[b, ] <- jaccard.i
}




# initialize dataframe to store the family-wise error rate for all tests and choice of response variable
df.FWER <- data.frame("delongRF" = numeric(6),
                      "delongGLM" = numeric(6),
                      "tramRF" = numeric(6),
                      "tramGLM" = numeric(6),
                      "correlation" = numeric(6),
                      "residual" = numeric(6),
                      "oracle" = numeric(6))

# the same for the jaccard index 
df.JAC <- df.FWER


# stratify by response/target
for(tar in 1:6){
  ind.tar <- which(target == tar)
  
  df.FWER[tar, ] <- colMeans(FWE.mod1[ind.tar, ])
  df.JAC[tar, ] <- colMeans(scores.j.mod1[ind.tar, ])
  
}




# colors for each test
group_colors <- scales::hue_pal()(6)

response_labels <- c("Response: 1", "Response: 2", "Response: 3", "Response: 4", "Response: 5", "Response: 6")


#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# DeLong (RF)
#-------------------------------------------------------------------------------


points.mod1.delongRF <- data.frame(
  x = df.FWER$delongRF,
  y = df.JAC$delongRF,
  group = factor(1:6, labels = response_labels)
)



p.delongRF <- ggplot(points.mod1.delongRF, aes(x = x, y = y, color = group)) +
  geom_point(shape = 1, size = 2, stroke = 1.5) + 
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  theme_bw(base_size = size) +
  guides(
    color = guide_legend(override.aes = list(shape = 1, size = 2))  
  ) +
  labs(
    color = "Group",
    x = "FWER",
    y = "Jaccard index",
    title = "DeLong (RF)"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size)
  )+
  coord_cartesian(ylim = c(0, 1)) 

#p.delongRF



#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# DeLong (GLM)
#-------------------------------------------------------------------------------


points.mod1.delongGLM <- data.frame(
  x = df.FWER$delongGLM,
  y = df.JAC$delongGLM,
  group = factor(1:6, labels = response_labels)
  )



p.delongGLM <- ggplot(points.mod1.delongGLM, aes(x = x, y = y, color = group)) +
  geom_point(shape = 1, size = 2, stroke = 1.5) + 
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  theme_bw(base_size = size) +
  guides(
    color = guide_legend(override.aes = list(shape = 1, size = 2))  
  ) +
  labs(
    color = "Group",
    x = "FWER",
    y = "Jaccard index",
    title = "DeLong (GLM)"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size)
  )+
  coord_cartesian(ylim = c(0, 1)) 

#p.delongGLM



#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# tram-gcm (RF)
#-------------------------------------------------------------------------------


points.mod1.tramRF <- data.frame(
  x = df.FWER$tramRF,
  y = df.JAC$tramRF,
  group = factor(1:6, labels = response_labels)
)



p.tramRF <- ggplot(points.mod1.tramRF, aes(x = x, y = y, color = group)) +
  geom_point(shape = 1, size = 2, stroke = 1.5) + 
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  theme_bw(base_size = size) +
  guides(
    color = guide_legend(override.aes = list(shape = 1, size = 2))  
  ) +
  labs(
    color = "Group",
    x = "FWER",
    y = "Jaccard index",
    title = "TRAM-GCM (RF)"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size)
  )+
  coord_cartesian(ylim = c(0, 1)) 

#p.tramRF



#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# tram-gcm (GLM)
#-------------------------------------------------------------------------------


points.mod1.tramGLM <- data.frame(
  x = df.FWER$tramGLM,
  y = df.JAC$tramGLM,
  group = factor(1:6, labels = response_labels)
)



p.tramGLM <- ggplot(points.mod1.tramGLM, aes(x = x, y = y, color = group)) +
  geom_point(shape = 1, size = 2, stroke = 1.5) + 
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  theme_bw(base_size = size) +
  guides(
    color = guide_legend(override.aes = list(shape = 1, size = 2))  
  ) +
  labs(
    color = "Group",
    x = "FWER",
    y = "Jaccard index",
    title = "TRAM-GCM (GLM)"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size)
  )+
  coord_cartesian(ylim = c(0, 1)) 

#p.tramGLM



#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# correlation
#-------------------------------------------------------------------------------


points.mod1.corr <- data.frame(
  x = df.FWER$correlation,
  y = df.JAC$correlation,
  group = factor(1:6, labels = response_labels)
)



p.corr <- ggplot(points.mod1.corr, aes(x = x, y = y, color = group)) +
  geom_point(shape = 1, size = 2, stroke = 1.5) + 
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  theme_bw(base_size = size) +
  guides(
    color = guide_legend(override.aes = list(shape = 1, size = 2))  
  ) +
  labs(
    color = "Group",
    x = "FWER",
    y = "Jaccard index",
    title = "Correlation"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size)
  )+
  coord_cartesian(ylim = c(0, 1)) 

#p.corr



#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# residual
#-------------------------------------------------------------------------------


points.mod1.res <- data.frame(
  x = df.FWER$residual,
  y = df.JAC$residual,
  group = factor(1:6, labels = response_labels)
)



p.residual <- ggplot(points.mod1.res, aes(x = x, y = y, color = group)) +
  geom_point(shape = 1, size = 2, stroke = 1.5) + 
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  theme_bw(base_size = size) +
  guides(
    color = guide_legend(override.aes = list(shape = 1, size = 2))  
  ) +
  labs(
    color = "Group",
    x = "FWER",
    y = "Jaccard index",
    title = "Residual"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size)
  )+
  coord_cartesian(ylim = c(0, 1)) 

#p.residual



#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# oracle
#-------------------------------------------------------------------------------


points.mod1.oracle <- data.frame(
  x = df.FWER$oracle,
  y = df.JAC$oracle,
  group = factor(1:6, labels = response_labels)
)



p.oracle <- ggplot(points.mod1.oracle, aes(x = x, y = y, color = group)) +
  geom_point(shape = 1, size = 2, stroke = 1.5) + 
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  theme_bw(base_size = size) +
  guides(
    color = guide_legend(override.aes = list(shape = 1, size = 2))  
  ) +
  labs(
    color = "Group",
    x = "FWER",
    y = "Jaccard index",
    title = "Oracle"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size)
  )+
  coord_cartesian(ylim = c(0, 1)) 

#p.oracle



#-------------------------------------------------------------------------------

save(points.mod1.delongRF, points.mod1.delongGLM, points.mod1.tramRF, points.mod1.tramGLM, points.mod1.corr, points.mod1.res, points.mod1.oracle, file = "ICP_level_power_bytarget.rdata")


# Combine the plots with a common legend
combined <- p.delongRF + p.delongGLM + p.tramRF + p.tramGLM + p.corr + p.residual + p.oracle &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size),
    legend.box = "vertical"  
  )

combined_plot <- combined + plot_layout(guides = "collect") &
  theme(
    legend.box.just = "center",
    legend.margin = margin(t = 10)
  )


combined_plot


ggsave(filename = "ICP_level_power_bytarget.pdf", width = 7, height = 7.5)






