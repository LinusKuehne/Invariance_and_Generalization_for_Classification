# In this script, we see how level and power changes for the different invariance tests
# when the structural equation (i.e. the "model") for the response Y changes for data
# from the standard SCM.



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



# All subsets of predictors
# 0 encodes the empty set
sets <- list(c(1,2,3), c(1,3), c(1,2), c(2,3), c(1), c(2), c(3), c(0))




# number of simulation repetitions (nreps = 500 would take roughly 3 h for both experiments combined)
nreps <- 500

# number of observations per environment => 5*n total observations
n <- 200

# intervention strength scaling
t <- 1

# for plotting
size <- 10

#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# model 1 (logistic regression)
#-------------------------------------------------------------------------------

set.seed(1)

# here we store the p-values corresponding to each set in sets for the different invariance tests
pvals.delongRF.mod1 <- data.frame("p123" = numeric(nreps),
                                  "p13" = numeric(nreps),
                                  "p12" = numeric(nreps),
                                  "p23" = numeric(nreps),
                                  "p1" = numeric(nreps),
                                  "p2" = numeric(nreps),
                                  "p3" = numeric(nreps),
                                  "p0" = numeric(nreps))


pvals.delongGLM.mod1 <- pvals.tramRF.mod1 <- pvals.tramGLM.mod1 <- pvals.corr.mod1 <- pvals.residual.mod1 <- pvals.delongRF.mod1


for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # generate sample from the standard SCM
  sample <- gen.sample.standard(n = n, t = t, mod = "logreg")
  
  # compute p-values for all tests and subsets
  pvals.all <- pvalues.all(sample)
  
  # store the p-values
  pvals.delongRF.mod1[b,] <- pvals.all$delong.rf
  pvals.delongGLM.mod1[b,] <- pvals.all$delong.glm
  pvals.tramRF.mod1[b,] <- pvals.all$tram.rf
  pvals.tramGLM.mod1[b,] <- pvals.all$tram.glm
  pvals.corr.mod1[b,] <- pvals.all$correlation
  pvals.residual.mod1[b,] <- pvals.all$residual

}

# compute the rejection rates with respect to all sets
r.rate.delongRF.mod1 <- colMeans(pvals.delongRF.mod1 < 0.05)
r.rate.delongGLM.mod1 <- colMeans(pvals.delongGLM.mod1 < 0.05)
r.rate.tramRF.mod1 <- colMeans(pvals.tramRF.mod1 < 0.05)
r.rate.tramGLM.mod1 <- colMeans(pvals.tramGLM.mod1 < 0.05)
r.rate.corr.mod1 <- colMeans(pvals.corr.mod1 < 0.05)
r.rate.residual.mod1 <- colMeans(pvals.residual.mod1 < 0.05)

# compute the Jaccard index for each test
jaccard.delongRF.mod1 <- mean(apply(X = pvals.delongRF.mod1, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.delongGLM.mod1 <- mean(apply(X = pvals.delongGLM.mod1, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.tramRF.mod1 <- mean(apply(X = pvals.tramRF.mod1, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.tramGLM.mod1 <- mean(apply(X = pvals.tramGLM.mod1, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.corr.mod1 <- mean(apply(X = pvals.corr.mod1, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.residual.mod1 <- mean(apply(X = pvals.residual.mod1, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))


# Create a data frame with the points and group information
points.mod1 <- data.frame(
  x = c(r.rate.delongRF.mod1["p13"], r.rate.delongRF.mod1["p1"], r.rate.delongGLM.mod1["p13"], r.rate.delongGLM.mod1["p1"], r.rate.tramRF.mod1["p13"], r.rate.tramRF.mod1["p1"], r.rate.tramGLM.mod1["p13"], r.rate.tramGLM.mod1["p1"], r.rate.corr.mod1["p13"], r.rate.corr.mod1["p1"], r.rate.residual.mod1["p13"], r.rate.residual.mod1["p1"]),
  y = c(jaccard.delongRF.mod1, jaccard.delongRF.mod1, jaccard.delongGLM.mod1, jaccard.delongGLM.mod1, jaccard.tramRF.mod1, jaccard.tramRF.mod1, jaccard.tramGLM.mod1, jaccard.tramGLM.mod1, jaccard.corr.mod1, jaccard.corr.mod1, jaccard.residual.mod1, jaccard.residual.mod1),
  group = rep(1:6, each = 2),
  type = rep(c("circle", "cross"), 6)
)



# colors for each test
group_colors <- scales::hue_pal()(6)

# custom labels for the tests and subsets
test_labels <- c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual")
set_labels <- c("S = {1,3}", "S = {1}")



p.mod1 <- ggplot(points.mod1, aes(x = x, y = y, group = group)) +
  geom_point(aes(shape = type, color = as.factor(group)), size = 2, stroke = 1.5) +
  geom_line(aes(color = as.factor(group))) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  scale_shape_manual(values = c("circle" = 1, "cross" = 4), labels = set_labels) +
  scale_color_manual(values = group_colors, labels = test_labels) +
  theme_bw(base_size = size) +
  guides(
    shape = guide_legend(override.aes = list(color = "black", size = 2)),
    color = guide_legend(override.aes = list(shape = NA))
  ) +
  labs(
    shape = "Invariant subset",
    color = "Invariance test",
    x = "Rejection rate of invariant S",    
    y = "Average Jaccard index",   
    title = "Logistic Regression Model"  
  ) +
  coord_cartesian(ylim = c(0, 1)) 

#p.mod1



#-------------------------------------------------------------------------------













#-------------------------------------------------------------------------------
# model 2 (probit regression)
#-------------------------------------------------------------------------------

set.seed(1)

# here we store the p-values corresponding to each set in sets for the different invariance tests
pvals.delongRF.mod2 <- data.frame("p123" = numeric(nreps),
                                  "p13" = numeric(nreps),
                                  "p12" = numeric(nreps),
                                  "p23" = numeric(nreps),
                                  "p1" = numeric(nreps),
                                  "p2" = numeric(nreps),
                                  "p3" = numeric(nreps),
                                  "p0" = numeric(nreps))

pvals.delongGLM.mod2 <- pvals.tramRF.mod2 <- pvals.tramGLM.mod2 <- pvals.corr.mod2 <- pvals.residual.mod2 <- pvals.delongRF.mod2





for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # generate sample from the standard SCM
  sample <- gen.sample.standard(n = n, t = t, mod = "probit")
  
  # compute p-values for all tests and subsets
  pvals.all <- pvalues.all(sample)
  
  # store the p-values
  pvals.delongRF.mod2[b,] <- pvals.all$delong.rf
  pvals.delongGLM.mod2[b,] <- pvals.all$delong.glm
  pvals.tramRF.mod2[b,] <- pvals.all$tram.rf
  pvals.tramGLM.mod2[b,] <- pvals.all$tram.glm
  pvals.corr.mod2[b,] <- pvals.all$correlation
  pvals.residual.mod2[b,] <- pvals.all$residual
  
}

# compute the rejection rates with respect to all sets
r.rate.delongRF.mod2 <- colMeans(pvals.delongRF.mod2 < 0.05)
r.rate.delongGLM.mod2 <- colMeans(pvals.delongGLM.mod2 < 0.05)
r.rate.tramRF.mod2 <- colMeans(pvals.tramRF.mod2 < 0.05)
r.rate.tramGLM.mod2 <- colMeans(pvals.tramGLM.mod2 < 0.05)
r.rate.corr.mod2 <- colMeans(pvals.corr.mod2 < 0.05)
r.rate.residual.mod2 <- colMeans(pvals.residual.mod2 < 0.05)

# compute the Jaccard index for each test
jaccard.delongRF.mod2 <- mean(apply(X = pvals.delongRF.mod2, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.delongGLM.mod2 <- mean(apply(X = pvals.delongGLM.mod2, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.tramRF.mod2 <- mean(apply(X = pvals.tramRF.mod2, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.tramGLM.mod2 <- mean(apply(X = pvals.tramGLM.mod2, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.corr.mod2 <- mean(apply(X = pvals.corr.mod2, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.residual.mod2 <- mean(apply(X = pvals.residual.mod2, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))


# create a data frame with the points and group information
points.mod2 <- data.frame(
  x = c(r.rate.delongRF.mod2["p13"], r.rate.delongRF.mod2["p1"], r.rate.delongGLM.mod2["p13"], r.rate.delongGLM.mod2["p1"], r.rate.tramRF.mod2["p13"], r.rate.tramRF.mod2["p1"], r.rate.tramGLM.mod2["p13"], r.rate.tramGLM.mod2["p1"], r.rate.corr.mod2["p13"], r.rate.corr.mod2["p1"], r.rate.residual.mod2["p13"], r.rate.residual.mod2["p1"]),
  y = c(jaccard.delongRF.mod2, jaccard.delongRF.mod2, jaccard.delongGLM.mod2, jaccard.delongGLM.mod2, jaccard.tramRF.mod2, jaccard.tramRF.mod2, jaccard.tramGLM.mod2, jaccard.tramGLM.mod2, jaccard.corr.mod2, jaccard.corr.mod2, jaccard.residual.mod2, jaccard.residual.mod2),
  group = rep(1:6, each = 2),
  type = rep(c("circle", "cross"), 6)
)



# colors for each test
group_colors <- scales::hue_pal()(6)

# labels for the tests and subsets
test_labels <- c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual")
set_labels <- c("S = {1,3}", "S = {1}")



p.mod2 <- ggplot(points.mod2, aes(x = x, y = y, group = group)) +
  geom_point(aes(shape = type, color = as.factor(group)), size = 2, stroke = 1.5) +
  geom_line(aes(color = as.factor(group))) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  scale_shape_manual(values = c("circle" = 1, "cross" = 4), labels = set_labels) +
  scale_color_manual(values = group_colors, labels = test_labels) +
  theme_bw(base_size = size) +
  guides(
    shape = guide_legend(override.aes = list(color = "black", size = 2)),
    color = guide_legend(override.aes = list(shape = NA))
  ) +
  labs(
    shape = "Invariant subset",
    color = "Invariance test",
    x = "Rejection rate of invariant S",    
    y = "Average Jaccard index",   
    title = "Probit Regression Model"  
  ) +
  coord_cartesian(ylim = c(0, 1)) 

#p.mod2







#-------------------------------------------------------------------------------













#-------------------------------------------------------------------------------
# model 3 ("non-linear logistic regression")
#-------------------------------------------------------------------------------

set.seed(1)

# here we store the p-values corresponding to each set in sets for the different invariance tests
pvals.delongRF.mod3 <- data.frame("p123" = numeric(nreps),
                                  "p13" = numeric(nreps),
                                  "p12" = numeric(nreps),
                                  "p23" = numeric(nreps),
                                  "p1" = numeric(nreps),
                                  "p2" = numeric(nreps),
                                  "p3" = numeric(nreps),
                                  "p0" = numeric(nreps))

pvals.delongGLM.mod3 <- pvals.tramRF.mod3 <- pvals.tramGLM.mod3 <- pvals.corr.mod3 <- pvals.residual.mod3 <- pvals.delongRF.mod3





for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # generate a sample from the standard SCM
  sample <- gen.sample.standard(n = n, t = t, mod = "nonlin")
  
  # compute p-values for all tests and subsets
  pvals.all <- pvalues.all(sample)
  
  # store the p-values
  pvals.delongRF.mod3[b,] <- pvals.all$delong.rf
  pvals.delongGLM.mod3[b,] <- pvals.all$delong.glm
  pvals.tramRF.mod3[b,] <- pvals.all$tram.rf
  pvals.tramGLM.mod3[b,] <- pvals.all$tram.glm
  pvals.corr.mod3[b,] <- pvals.all$correlation
  pvals.residual.mod3[b,] <- pvals.all$residual
  
}

# compute the rejection rates with respect to all sets
r.rate.delongRF.mod3 <- colMeans(pvals.delongRF.mod3 < 0.05)
r.rate.delongGLM.mod3 <- colMeans(pvals.delongGLM.mod3 < 0.05)
r.rate.tramRF.mod3 <- colMeans(pvals.tramRF.mod3 < 0.05)
r.rate.tramGLM.mod3 <- colMeans(pvals.tramGLM.mod3 < 0.05)
r.rate.corr.mod3 <- colMeans(pvals.corr.mod3 < 0.05)
r.rate.residual.mod3 <- colMeans(pvals.residual.mod3 < 0.05)

# compute the Jaccard index for each test
jaccard.delongRF.mod3 <- mean(apply(X = pvals.delongRF.mod3, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.delongGLM.mod3 <- mean(apply(X = pvals.delongGLM.mod3, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.tramRF.mod3 <- mean(apply(X = pvals.tramRF.mod3, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.tramGLM.mod3 <- mean(apply(X = pvals.tramGLM.mod3, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.corr.mod3 <- mean(apply(X = pvals.corr.mod3, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.residual.mod3 <- mean(apply(X = pvals.residual.mod3, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))


# Create a data frame with the points and group information
points.mod3 <- data.frame(
  x = c(r.rate.delongRF.mod3["p13"], r.rate.delongRF.mod3["p1"], r.rate.delongGLM.mod3["p13"], r.rate.delongGLM.mod3["p1"], r.rate.tramRF.mod3["p13"], r.rate.tramRF.mod3["p1"], r.rate.tramGLM.mod3["p13"], r.rate.tramGLM.mod3["p1"], r.rate.corr.mod3["p13"], r.rate.corr.mod3["p1"], r.rate.residual.mod3["p13"], r.rate.residual.mod3["p1"]),
  y = c(jaccard.delongRF.mod3, jaccard.delongRF.mod3, jaccard.delongGLM.mod3, jaccard.delongGLM.mod3, jaccard.tramRF.mod3, jaccard.tramRF.mod3, jaccard.tramGLM.mod3, jaccard.tramGLM.mod3, jaccard.corr.mod3, jaccard.corr.mod3, jaccard.residual.mod3, jaccard.residual.mod3),
  group = rep(1:6, each = 2),
  type = rep(c("circle", "cross"), 6)
)



# colors for each test
group_colors <- scales::hue_pal()(6)

# custom labels for the tests and subsets
test_labels <- c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual")
set_labels <- c("S = {1,3}", "S = {1}")



p.mod3 <- ggplot(points.mod3, aes(x = x, y = y, group = group)) +
  geom_point(aes(shape = type, color = as.factor(group)), size = 2, stroke = 1.5) +
  geom_line(aes(color = as.factor(group))) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  scale_shape_manual(values = c("circle" = 1, "cross" = 4), labels = set_labels) +
  scale_color_manual(values = group_colors, labels = test_labels) +
  theme_bw(base_size = size) +
  guides(
    shape = guide_legend(override.aes = list(color = "black", size = 2)),
    color = guide_legend(override.aes = list(shape = NA))
  ) +
  labs(
    shape = "Invariant subset",
    color = "Invariance test",
    x = "Rejection rate of invariant S",    
    y = "Average Jaccard index",   
    title = "Non-Linear Logistic Regression Model"  
  ) +
  coord_cartesian(ylim = c(0, 1)) 

#p.mod3







#-------------------------------------------------------------------------------












#-------------------------------------------------------------------------------
# model 4 ("bump model")
#-------------------------------------------------------------------------------

set.seed(1)

# here we store the p-values corresponding to each set in sets for the different invariance tests
pvals.delongRF.mod4 <- data.frame("p123" = numeric(nreps),
                                  "p13" = numeric(nreps),
                                  "p12" = numeric(nreps),
                                  "p23" = numeric(nreps),
                                  "p1" = numeric(nreps),
                                  "p2" = numeric(nreps),
                                  "p3" = numeric(nreps),
                                  "p0" = numeric(nreps))

pvals.delongGLM.mod4 <- pvals.tramRF.mod4 <- pvals.tramGLM.mod4 <- pvals.corr.mod4 <- pvals.residual.mod4 <- pvals.delongRF.mod4





for(b in 1:nreps){
  print(paste0("Simulation iteration ", b, " out of ", nreps))
  
  # compute p-values for all tests and subsets
  sample <- gen.sample.standard(n = n, t = t, mod = "bump")
  
  # compute p-values for all tests and subsets
  pvals.all <- pvalues.all(sample)
  
  # store the p-values
  pvals.delongRF.mod4[b,] <- pvals.all$delong.rf
  pvals.delongGLM.mod4[b,] <- pvals.all$delong.glm
  pvals.tramRF.mod4[b,] <- pvals.all$tram.rf
  pvals.tramGLM.mod4[b,] <- pvals.all$tram.glm
  pvals.corr.mod4[b,] <- pvals.all$correlation
  pvals.residual.mod4[b,] <- pvals.all$residual
  
}

# compute the rejection rates with respect to all sets
r.rate.delongRF.mod4 <- colMeans(pvals.delongRF.mod4 < 0.05)
r.rate.delongGLM.mod4 <- colMeans(pvals.delongGLM.mod4 < 0.05)
r.rate.tramRF.mod4 <- colMeans(pvals.tramRF.mod4 < 0.05)
r.rate.tramGLM.mod4 <- colMeans(pvals.tramGLM.mod4 < 0.05)
r.rate.corr.mod4 <- colMeans(pvals.corr.mod4 < 0.05)
r.rate.residual.mod4 <- colMeans(pvals.residual.mod4 < 0.05)

# compute the Jaccard index for each test
jaccard.delongRF.mod4 <- mean(apply(X = pvals.delongRF.mod4, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.delongGLM.mod4 <- mean(apply(X = pvals.delongGLM.mod4, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.tramRF.mod4 <- mean(apply(X = pvals.tramRF.mod4, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.tramGLM.mod4 <- mean(apply(X = pvals.tramGLM.mod4, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.corr.mod4 <- mean(apply(X = pvals.corr.mod4, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))
jaccard.residual.mod4 <- mean(apply(X = pvals.residual.mod4, FUN = jaccard.standard, MARGIN = 1, thresh = 0.05))


# Create a data frame with the points and group information
points.mod4 <- data.frame(
  x = c(r.rate.delongRF.mod4["p13"], r.rate.delongRF.mod4["p1"], r.rate.delongGLM.mod4["p13"], r.rate.delongGLM.mod4["p1"], r.rate.tramRF.mod4["p13"], r.rate.tramRF.mod4["p1"], r.rate.tramGLM.mod4["p13"], r.rate.tramGLM.mod4["p1"], r.rate.corr.mod4["p13"], r.rate.corr.mod4["p1"], r.rate.residual.mod4["p13"], r.rate.residual.mod4["p1"]),
  y = c(jaccard.delongRF.mod4, jaccard.delongRF.mod4, jaccard.delongGLM.mod4, jaccard.delongGLM.mod4, jaccard.tramRF.mod4, jaccard.tramRF.mod4, jaccard.tramGLM.mod4, jaccard.tramGLM.mod4, jaccard.corr.mod4, jaccard.corr.mod4, jaccard.residual.mod4, jaccard.residual.mod4),
  group = rep(1:6, each = 2),
  type = rep(c("circle", "cross"), 6)
)



# colors for each test
group_colors <- scales::hue_pal()(6)

# custom labels for the tests and subsets
test_labels <- c("DeLong (RF)", "DeLong (GLM)", "TRAM-GCM (RF)", "TRAM-GCM (GLM)", "Correlation", "Residual")
set_labels <- c("S = {1,3}", "S = {1}")



p.mod4 <- ggplot(points.mod4, aes(x = x, y = y, group = group)) +
  geom_point(aes(shape = type, color = as.factor(group)), size = 2, stroke = 1.5) +
  geom_line(aes(color = as.factor(group))) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "black") +
  scale_shape_manual(values = c("circle" = 1, "cross" = 4), labels = set_labels) +
  scale_color_manual(values = group_colors, labels = test_labels) +
  theme_bw(base_size = size) +
  guides(
    shape = guide_legend(override.aes = list(color = "black", size = 2)),
    color = guide_legend(override.aes = list(shape = NA))
  ) +
  labs(
    shape = "Invariant subset",
    color = "Invariance test",
    x = "Rejection rate of invariant S",    
    y = "Average Jaccard index",   
    title = "Bump Model"  
  ) +
  coord_cartesian(ylim = c(0, 1)) 

#p.mod4







#-------------------------------------------------------------------------------



save(points.mod1, points.mod2, points.mod3, points.mod4, file = file.path(script_dir, "saved_data/level_power_standard.rdata"))



# Combine the plots with a common legend
combined <- p.mod1 + p.mod2 + p.mod3 + p.mod4 &
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = size),
    legend.box = "vertical")

combined_plot <- combined + plot_layout(guides = "collect") &
  theme(
    legend.box.just = "center",
    legend.margin = margin(t = 0)
  )


#combined_plot


ggsave(filename = file.path(script_dir, "saved_plots/level_power_standard.pdf"), width = 7, height = 7.5)



# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/level_power_standard.txt"))

