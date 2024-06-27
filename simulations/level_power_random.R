# In this script, we see how level and power changes for the different invariance tests
# when the structural equation (i.e. the "model") for the response Y changes for data
# from the semirandom SCM.


library(ggplot2)
library(patchwork)
library(rje)


# get the path of this script
script_dir <- getwd()

# load in the functions needed
source(file.path(script_dir, "../code/code_simulations/invariance_tests.R"))
source(file.path(script_dir, "../code/code_simulations/data_generating_process.R"))
source(file.path(script_dir, "../code/code_simulations/utils.R"))




#-------------------------------------------------------------------------------
# parameters
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

# number of simulation repetitions (nreps = 300 would take roughly 7h for both experiments combined)
nreps <- 300

# number of observations per environment => 5*n total observations
n <- 200

# intervention strength scaling
t <- 1


# names of the predictors
Xnames <- rep("A", 5)
for(w in 1:5){
  number <- as.character(w)
  name <- paste("X", number, sep="") 
  Xnames[w] <- name  
}

ps <- powerSet(Xnames)

# for plotting
size <- 10

#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# model 1 (logistic regression)
#-------------------------------------------------------------------------------

set.seed(1)

# here we store the Jaccard scores for the different invariance tests
scores.j.mod1 <- data.frame("delongRF" = numeric(nreps),
                                "delongGLM" = numeric(nreps),
                                "tramRF" = numeric(nreps),
                                "tramGLM" = numeric(nreps),
                                "correlation" = numeric(nreps),
                                "residual" = numeric(nreps))



# these vectors store the p-values associated with the ground truth invariant sets
inv.pvals.delongRF.mod1 <- numeric(0)
inv.pvals.delongGLM.mod1 <- numeric(0)
inv.pvals.tramRF.mod1 <- numeric(0)
inv.pvals.tramGLM.mod1 <- numeric(0)
inv.pvals.correlation.mod1 <- numeric(0)
inv.pvals.residual.mod1 <- numeric(0)



for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # generate data from the semirandom SCM
  s <- generate.samples.random(n = n, n.test = 1000, d = d, max.pa = max.pa, num.int = num.int, t = t, mod = "logreg")
  
  # compute p-values for all tests and subsets
  pvals.all <- pvalues.all(s$sample_train)
  
  # find the sets from the DAG satisfying the d-separation associated with invariance
  stable.sets <- int.stable.sets(dag.full = s$dag.full, d = d, num.int = num.int)
  
  # find indeces of ground truth invariant sets with reference to the list of sets
  ind.invariant <- (1:length(ps))[ps %in% stable.sets]
  
  # extract rows of pvalues matrix corresponding to invariant sets
  pvals.inv <- pvals.all[ind.invariant, , drop = F]
  
  # append p-vals of invariant sets
  inv.pvals.delongRF.mod1 <- c(inv.pvals.delongRF.mod1, pvals.inv$delong.rf)
  inv.pvals.delongGLM.mod1 <- c(inv.pvals.delongGLM.mod1, pvals.inv$delong.glm)
  inv.pvals.tramRF.mod1 <- c(inv.pvals.tramRF.mod1, pvals.inv$tram.rf)
  inv.pvals.tramGLM.mod1 <- c(inv.pvals.tramGLM.mod1, pvals.inv$tram.glm)
  inv.pvals.correlation.mod1 <- c(inv.pvals.correlation.mod1, pvals.inv$correlation)
  inv.pvals.residual.mod1 <- c(inv.pvals.residual.mod1, pvals.inv$residual)
  
  # compute jaccard index
  scores.j.mod1[b, ] <- apply(X = pvals.all, MARGIN = 2, FUN = jaccard, thresh = 0.05, stable.sets = stable.sets, num.covariates = (d-1))
  
}

# compute rejection rate for invariant subsets
r.rate.delongRF.mod1 <- mean(inv.pvals.delongRF.mod1 < 0.05)
r.rate.delongGLM.mod1 <- mean(inv.pvals.delongGLM.mod1 < 0.05)
r.rate.tramRF.mod1 <- mean(inv.pvals.tramRF.mod1 < 0.05)
r.rate.tramGLM.mod1 <- mean(inv.pvals.tramGLM.mod1 < 0.05)
r.rate.corr.mod1 <- mean(inv.pvals.correlation.mod1 < 0.05)
r.rate.residual.mod1 <- mean(inv.pvals.residual.mod1 < 0.05)

# average the Jaccard index across all simulation repetitions
avg.jaccard.mod1 <- colMeans(scores.j.mod1)






# create a data frame with the points and group information
points.mod1 <- data.frame(
  x = c(r.rate.delongRF.mod1, r.rate.delongGLM.mod1, r.rate.tramRF.mod1, r.rate.tramGLM.mod1, r.rate.corr.mod1, r.rate.residual.mod1),
  y = avg.jaccard.mod1,
  group = factor(1:6, labels = c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual"))
  #type = rep("circle", 6)
)



# colors for each test
group_colors <- scales::hue_pal()(6)

# labels for the tests and subsets
test_labels <- c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual")






p.mod1 <- ggplot(points.mod1, aes(x = x, y = y, color = group)) +
  geom_point(shape = 1, size = 2, stroke = 1.5) +  # Use shape 1 for hollow circles
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  
  theme_bw(base_size = size) +
  guides(
    color = guide_legend(override.aes = list(shape = 1, size = 2))  # Legend with hollow circles
  ) +
  labs(
    color = "Group",
    x = "Average rejection rate across\ninvariant subsets",
    y = "Average Jaccard index",
    title = "Logistic Regression Model"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size)
  )+
  coord_cartesian(ylim = c(0, 1)) 

#p.mod1



#-------------------------------------------------------------------------------













#-------------------------------------------------------------------------------
# model 2 (probit regression)
#-------------------------------------------------------------------------------

set.seed(1)

# here we store the Jaccard scores for the different invariance tests
scores.j.mod2 <- data.frame("delongRF" = numeric(nreps),
                            "delongGLM" = numeric(nreps),
                            "tramRF" = numeric(nreps),
                            "tramGLM" = numeric(nreps),
                            "correlation" = numeric(nreps),
                            "residual" = numeric(nreps))



# these vectors store the p-values associated with the ground truth invariant sets
inv.pvals.delongRF.mod2 <- numeric(0)
inv.pvals.delongGLM.mod2 <- numeric(0)
inv.pvals.tramRF.mod2 <- numeric(0)
inv.pvals.tramGLM.mod2 <- numeric(0)
inv.pvals.correlation.mod2 <- numeric(0)
inv.pvals.residual.mod2 <- numeric(0)



for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # generate data from the semirandom SCM
  s <- generate.samples.random(n = n, n.test = 1000, d = d, max.pa = max.pa, num.int = num.int, t = t, mod = "probit")
  
  # compute p-values for all tests and subsets
  pvals.all <- pvalues.all(s$sample_train)
  
  # find the sets from the DAG satisfying the d-separation associated with invariance
  stable.sets <- int.stable.sets(dag.full = s$dag.full, d = d, num.int = num.int)
  
  # find indeces of ground truth invariant sets with reference to the list of sets
  ind.invariant <- (1:length(ps))[ps %in% stable.sets]
  
  # extract rows of pvalues matrix corresponding to invariant sets
  pvals.inv <- pvals.all[ind.invariant, , drop = F]
  
  # append p-vals of invariant sets
  inv.pvals.delongRF.mod2 <- c(inv.pvals.delongRF.mod2, pvals.inv$delong.rf)
  inv.pvals.delongGLM.mod2 <- c(inv.pvals.delongGLM.mod2, pvals.inv$delong.glm)
  inv.pvals.tramRF.mod2 <- c(inv.pvals.tramRF.mod2, pvals.inv$tram.rf)
  inv.pvals.tramGLM.mod2 <- c(inv.pvals.tramGLM.mod2, pvals.inv$tram.glm)
  inv.pvals.correlation.mod2 <- c(inv.pvals.correlation.mod2, pvals.inv$correlation)
  inv.pvals.residual.mod2 <- c(inv.pvals.residual.mod2, pvals.inv$residual)
  
  # compute jaccard index
  scores.j.mod2[b, ] <- apply(X = pvals.all, MARGIN = 2, FUN = jaccard, thresh = 0.05, stable.sets = stable.sets, num.covariates = (d-1))
  
}

# compute rejection rate for invariant subsets
r.rate.delongRF.mod2 <- mean(inv.pvals.delongRF.mod2 < 0.05)
r.rate.delongGLM.mod2 <- mean(inv.pvals.delongGLM.mod2 < 0.05)
r.rate.tramRF.mod2 <- mean(inv.pvals.tramRF.mod2 < 0.05)
r.rate.tramGLM.mod2 <- mean(inv.pvals.tramGLM.mod2 < 0.05)
r.rate.corr.mod2 <- mean(inv.pvals.correlation.mod2 < 0.05)
r.rate.residual.mod2 <- mean(inv.pvals.residual.mod2 < 0.05)

# average the Jaccard index across all simulation repetitions
avg.jaccard.mod2 <- colMeans(scores.j.mod2)






# create a data frame with the points and group information
points.mod2 <- data.frame(
  x = c(r.rate.delongRF.mod2, r.rate.delongGLM.mod2, r.rate.tramRF.mod2, r.rate.tramGLM.mod2, r.rate.corr.mod2, r.rate.residual.mod2),
  y = avg.jaccard.mod2,
  group = factor(1:6, labels = c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual"))
  #type = rep("circle", 6)
)



# colors for each test
group_colors <- scales::hue_pal()(6)

# labels for the tests and subsets
test_labels <- c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual")






p.mod2 <- ggplot(points.mod2, aes(x = x, y = y, color = group)) +
  geom_point(shape = 1, size = 2, stroke = 1.5) +  # Use shape 1 for hollow circles
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  
  theme_bw(base_size = size) +
  guides(
    color = guide_legend(override.aes = list(shape = 1, size = 2))  # Legend with hollow circles
  ) +
  labs(
    color = "Group",
    x = "Average rejection rate across\ninvariant subsets",
    y = "Average Jaccard index",
    title = "Probit Regression Model"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size)
  )+
  coord_cartesian(ylim = c(0, 1)) 

#p.mod2



#-------------------------------------------------------------------------------













#-------------------------------------------------------------------------------
# model 3 ("non-linear logistic regression")
#-------------------------------------------------------------------------------

set.seed(1)

# here we store the Jaccard scores for the different invariance tests
scores.j.mod3 <- data.frame("delongRF" = numeric(nreps),
                            "delongGLM" = numeric(nreps),
                            "tramRF" = numeric(nreps),
                            "tramGLM" = numeric(nreps),
                            "correlation" = numeric(nreps),
                            "residual" = numeric(nreps))



# these vectors store the p-values associated with the ground truth invariant sets
inv.pvals.delongRF.mod3 <- numeric(0)
inv.pvals.delongGLM.mod3 <- numeric(0)
inv.pvals.tramRF.mod3 <- numeric(0)
inv.pvals.tramGLM.mod3 <- numeric(0)
inv.pvals.correlation.mod3 <- numeric(0)
inv.pvals.residual.mod3 <- numeric(0)



for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # generate data from the semirandom SCM
  s <- generate.samples.random(n = n, n.test = 1000, d = d, max.pa = max.pa, num.int = num.int, t = t, mod = "nonlin")
  
  # compute p-values for all tests and subsets
  pvals.all <- pvalues.all(s$sample_train)
  
  # find the sets from the DAG satisfying the d-separation associated with invariance
  stable.sets <- int.stable.sets(dag.full = s$dag.full, d = d, num.int = num.int)
  
  # find indeces of ground truth invariant sets with reference to the list of sets
  ind.invariant <- (1:length(ps))[ps %in% stable.sets]
  
  # extract rows of pvalues matrix corresponding to invariant sets
  pvals.inv <- pvals.all[ind.invariant, , drop = F]
  
  # append p-vals of invariant sets
  inv.pvals.delongRF.mod3 <- c(inv.pvals.delongRF.mod3, pvals.inv$delong.rf)
  inv.pvals.delongGLM.mod3 <- c(inv.pvals.delongGLM.mod3, pvals.inv$delong.glm)
  inv.pvals.tramRF.mod3 <- c(inv.pvals.tramRF.mod3, pvals.inv$tram.rf)
  inv.pvals.tramGLM.mod3 <- c(inv.pvals.tramGLM.mod3, pvals.inv$tram.glm)
  inv.pvals.correlation.mod3 <- c(inv.pvals.correlation.mod3, pvals.inv$correlation)
  inv.pvals.residual.mod3 <- c(inv.pvals.residual.mod3, pvals.inv$residual)
  
  # compute jaccard index
  scores.j.mod3[b, ] <- apply(X = pvals.all, MARGIN = 2, FUN = jaccard, thresh = 0.05, stable.sets = stable.sets, num.covariates = (d-1))
  
}

# compute rejection rate for invariant subsets
r.rate.delongRF.mod3 <- mean(inv.pvals.delongRF.mod3 < 0.05)
r.rate.delongGLM.mod3 <- mean(inv.pvals.delongGLM.mod3 < 0.05)
r.rate.tramRF.mod3 <- mean(inv.pvals.tramRF.mod3 < 0.05)
r.rate.tramGLM.mod3 <- mean(inv.pvals.tramGLM.mod3 < 0.05)
r.rate.corr.mod3 <- mean(inv.pvals.correlation.mod3 < 0.05)
r.rate.residual.mod3 <- mean(inv.pvals.residual.mod3 < 0.05)

# average the Jaccard index across all simulation repetitions
avg.jaccard.mod3 <- colMeans(scores.j.mod3)






# create a data frame with the points and group information
points.mod3 <- data.frame(
  x = c(r.rate.delongRF.mod3, r.rate.delongGLM.mod3, r.rate.tramRF.mod3, r.rate.tramGLM.mod3, r.rate.corr.mod3, r.rate.residual.mod3),
  y = avg.jaccard.mod3,
  group = factor(1:6, labels = c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual"))
  #type = rep("circle", 6)
)



# colors for each test
group_colors <- scales::hue_pal()(6)

# labels for the tests and subsets
test_labels <- c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual")






p.mod3 <- ggplot(points.mod3, aes(x = x, y = y, color = group)) +
  geom_point(shape = 1, size = 2, stroke = 1.5) +  # Use shape 1 for hollow circles
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  
  theme_bw(base_size = size) +
  guides(
    color = guide_legend(override.aes = list(shape = 1, size = 2))  # Legend with hollow circles
  ) +
  labs(
    color = "Group",
    x = "Average rejection rate across\ninvariant subsets",
    y = "Average Jaccard index",
    title = "Non-Linear Logistic Regression Model"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size)
  )+
  coord_cartesian(ylim = c(0, 1)) 

#p.mod3






#-------------------------------------------------------------------------------












#-------------------------------------------------------------------------------
# model 4 ("bump model")
#-------------------------------------------------------------------------------

set.seed(1)

# here we store the Jaccard scores for the different invariance tests
scores.j.mod4 <- data.frame("delongRF" = numeric(nreps),
                            "delongGLM" = numeric(nreps),
                            "tramRF" = numeric(nreps),
                            "tramGLM" = numeric(nreps),
                            "correlation" = numeric(nreps),
                            "residual" = numeric(nreps))



# these vectors store the p-values associated with the ground truth invariant sets
inv.pvals.delongRF.mod4 <- numeric(0)
inv.pvals.delongGLM.mod4 <- numeric(0)
inv.pvals.tramRF.mod4 <- numeric(0)
inv.pvals.tramGLM.mod4 <- numeric(0)
inv.pvals.correlation.mod4 <- numeric(0)
inv.pvals.residual.mod4 <- numeric(0)



for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # generate data from the semirandom SCM
  s <- generate.samples.random(n = n, n.test = 1000, d = d, max.pa = max.pa, num.int = num.int, t = t, mod = "bump")
  
  # compute p-values for all tests and subsets
  pvals.all <- pvalues.all(s$sample_train)
  
  # find the sets from the DAG satisfying the d-separation associated with invariance
  stable.sets <- int.stable.sets(dag.full = s$dag.full, d = d, num.int = num.int)
  
  # find indeces of ground truth invariant sets with reference to the list of sets
  ind.invariant <- (1:length(ps))[ps %in% stable.sets]
  
  # extract rows of pvalues matrix corresponding to invariant sets
  pvals.inv <- pvals.all[ind.invariant, , drop = F]
  
  # append p-vals of invariant sets
  inv.pvals.delongRF.mod4 <- c(inv.pvals.delongRF.mod4, pvals.inv$delong.rf)
  inv.pvals.delongGLM.mod4 <- c(inv.pvals.delongGLM.mod4, pvals.inv$delong.glm)
  inv.pvals.tramRF.mod4 <- c(inv.pvals.tramRF.mod4, pvals.inv$tram.rf)
  inv.pvals.tramGLM.mod4 <- c(inv.pvals.tramGLM.mod4, pvals.inv$tram.glm)
  inv.pvals.correlation.mod4 <- c(inv.pvals.correlation.mod4, pvals.inv$correlation)
  inv.pvals.residual.mod4 <- c(inv.pvals.residual.mod4, pvals.inv$residual)
  
  # compute jaccard index
  scores.j.mod4[b, ] <- apply(X = pvals.all, MARGIN = 2, FUN = jaccard, thresh = 0.05, stable.sets = stable.sets, num.covariates = (d-1))
  
}

# compute rejection rate for invariant subsets
r.rate.delongRF.mod4 <- mean(inv.pvals.delongRF.mod4 < 0.05)
r.rate.delongGLM.mod4 <- mean(inv.pvals.delongGLM.mod4 < 0.05)
r.rate.tramRF.mod4 <- mean(inv.pvals.tramRF.mod4 < 0.05)
r.rate.tramGLM.mod4 <- mean(inv.pvals.tramGLM.mod4 < 0.05)
r.rate.corr.mod4 <- mean(inv.pvals.correlation.mod4 < 0.05)
r.rate.residual.mod4 <- mean(inv.pvals.residual.mod4 < 0.05)

# average the Jaccard index across all simulation repetitions
avg.jaccard.mod4 <- colMeans(scores.j.mod4)






# create a data frame with the points and group information
points.mod4 <- data.frame(
  x = c(r.rate.delongRF.mod4, r.rate.delongGLM.mod4, r.rate.tramRF.mod4, r.rate.tramGLM.mod4, r.rate.corr.mod4, r.rate.residual.mod4),
  y = avg.jaccard.mod4,
  group = factor(1:6, labels = c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual"))
  #type = rep("circle", 6)
)



# colors for each test
group_colors <- scales::hue_pal()(6)

# labels for the tests and subsets
test_labels <- c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual")






p.mod4 <- ggplot(points.mod4, aes(x = x, y = y, color = group)) +
  geom_point(shape = 1, size = 2, stroke = 1.5) +  # Use shape 1 for hollow circles
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  
  theme_bw(base_size = size) +
  guides(
    color = guide_legend(override.aes = list(shape = 1, size = 2))  # Legend with hollow circles
  ) +
  labs(
    color = "Group",
    x = "Average rejection rate across\ninvariant subsets",
    y = "Average Jaccard index",
    title = "Bump Model"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size)
  )+
  coord_cartesian(ylim = c(0, 1)) 

#p.mod4






#-------------------------------------------------------------------------------



save(points.mod1, points.mod2, points.mod3, points.mod4, file = file.path(script_dir, "saved_data/level_power_random.rdata"))






# Combine the plots with a common legend
combined <- p.mod1 + p.mod2 + p.mod3 + p.mod4 &
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


ggsave(filename = file.path(script_dir, "saved_plots/level_power_random.pdf"), width = 7, height = 7.5)






# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/level_power_random.txt"))

















