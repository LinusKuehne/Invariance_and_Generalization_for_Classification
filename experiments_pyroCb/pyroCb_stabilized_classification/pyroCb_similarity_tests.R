# in this script, we evaluate to which degree the rankings of subsets according to 
# the p-values by the different tests are similar. We also consider the LOEO CV losses.


library(rje)
library(ggplot2)
library(reshape2)
library(gridExtra)



# get the path of this script
script_dir <- getwd()



# load all computed p-values and LOEO CV losses
load(file.path(script_dir, "../saved_data/pyroCb_ICP_corr.rdata"))
load(file.path(script_dir, "../saved_data/pyroCb_ICP_delong_5.rdata"))
load(file.path(script_dir, "../saved_data/pyroCb_ICP_delong_9.rdata"))
load(file.path(script_dir, "../saved_data/pyroCb_ICP_delong_cont.rdata"))
load(file.path(script_dir, "../saved_data/pyroCb_ICP_res.rdata"))
load(file.path(script_dir, "../saved_data/pyroCb_ICP_tram_5.rdata"))
load(file.path(script_dir, "../saved_data/pyroCb_ICP_tram_9.rdata"))
load(file.path(script_dir, "../saved_data/pyroCb_subsets_generalization.rdata"))
 


# from the variable screening script
# using glm group lasso to get 13 variables
varincl <- c(3, 5, 8, 9, 10, 11, 12, 13, 14, 23, 28, 29, 30)
varincl <- varincl[order(varincl)]

# sets to check invariance
sets <- powerSet(varincl)
sets[[1]] <- c(0)





#-------------------------------------------------------------------------------
# generate Kendall's tau plot for the p-values of all subsets
#-------------------------------------------------------------------------------



# put all vectors used for a ranking into a dataframe
df.p <- data.frame("DeLong_cont" = pvals.delong.cont.1tail,
                   "DeLong_5" = pvals.delong.5.1tail,
                   "DeLong_9" = pvals.delong.9.1tail,
                   "GCM_5" = pvals.tram.5, 
                   "GCM_9" = pvals.tram.9,
                   "Corr_5" = pvals.corr.5,
                   "Corr_9" = pvals.corr.9,
                   "Residual_5" = pvals.res.5,
                   "Residual_9" = pvals.res.5,
                   "Neg. LOEO CV loss" = -wbce.worst)


# initialize a matrix for the Kendall's tau value for all pairwise comparisons
tau.mat <- matrix(NA, nrow = ncol(df.p), ncol = ncol(df.p))



# iterate over all comparisons
for(j in 1:(ncol(df.p)-1)){
  for(i in (j+1):ncol(df.p)){
    # compute kendall's tau
    tau.mat[i,j] <- cor(x = df.p[,i], y = df.p[,j], method = "kendall")
  }
}



# store result
mat <- tau.mat

# fill the NA's with zeros such that we can add the transpose
mat[is.na(mat)] <- 0


# convert to a symmetric matrix by adding the transpose
mat <- mat + t(mat)

# fill in NA's on diagonal (important for the plotting)
diag(mat) <- NA


# turn into a better shape for the plot
df <- melt(mat)


# generate grid with comparisons
x <- c("DeLong (cont.)", "DeLong (5)", "DeLong (9)", "TRAM-GCM (5)", "TRAM-GCM (9)", "Correlation (5)", "Correlation (9)", "Residual (5)", "Residual (9)", "Neg. LOEO CV loss")
y <- c("DeLong (cont.)", "DeLong (5)", "DeLong (9)", "TRAM-GCM (5)", "TRAM-GCM (9)", "Correlation (5)", "Correlation (9)", "Residual (5)", "Residual (9)", "Neg. LOEO CV loss")

data <- expand.grid(X=x, Y=y)
data$Z <- df$value

# font size
size <- 10

# assign colors to kendall tau values and plot
plt <- ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", na.value = "black")  +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Consistency Between Different Rankings", 
       subtitle = "All subsets",
       fill = "Kendall's Tau") +
  theme_bw(base_size = size) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

plt





#-------------------------------------------------------------------------------









#-------------------------------------------------------------------------------
# focus just on subsets of size >= 7
#-------------------------------------------------------------------------------


# extract the sets with size >= 7 (TRUE <=> subset contains at least 7 predictors)
big.sets <- rep(T, length(sets))

big.sets[1] <- F

for(s in 2:length(sets)){
  set <- sets[[s]]
  if(length(set)<7){
    big.sets[s] <- F
  }
}




# put the corresponding entries into a dataframe
df.p.big <- data.frame("DeLong_cont" = pvals.delong.cont.1tail[big.sets],
                       "DeLong_5" = pvals.delong.5.1tail[big.sets],
                       "DeLong_9" = pvals.delong.9.1tail[big.sets],
                       "GCM_5" = pvals.tram.5[big.sets], 
                       "GCM_9" = pvals.tram.9[big.sets],
                       "Corr_5" = pvals.corr.5[big.sets],
                       "Corr_9" = pvals.corr.9[big.sets],
                       "Residual_5" = pvals.res.5[big.sets],
                       "Residual_9" = pvals.res.5[big.sets],
                       "Neg. LOEO CV loss" = -wbce.worst[big.sets])


# initialize matrix for Kendall's tau values
tau.mat.big <- matrix(NA, nrow = ncol(df.p.big), ncol = ncol(df.p.big))



# iterate over all comparisons for the lower triangular matrix
for(j in 1:(ncol(df.p.big)-1)){
  for(i in (j+1):ncol(df.p.big)){
    # compute kendall's tau
    tau.mat.big[i,j] <- cor(x = df.p.big[,i], y = df.p.big[,j], method = "kendall")
  }
}


# store in matrix
mat <- tau.mat.big

# fill in NA's with 0 to be able to add the transpose
mat[is.na(mat)] <- 0


# convert to a symmetric matrix
mat <- mat + t(mat)

# put NA's on diagonal for plotting
diag(mat) <- NA


# convert to a different shape
df <- melt(mat)


# generate grid for comparisons
x <- c("DeLong (cont.)", "DeLong (5)", "DeLong (9)", "TRAM-GCM (5)", "TRAM-GCM (9)", "Correlation (5)", "Correlation (9)", "Residual (5)", "Residual (9)", "Neg. LOEO CV loss")
y <- c("DeLong (cont.)", "DeLong (5)", "DeLong (9)", "TRAM-GCM (5)", "TRAM-GCM (9)", "Correlation (5)", "Correlation (9)", "Residual (5)", "Residual (9)", "Neg. LOEO CV loss")

data <- expand.grid(X=x, Y=y)
data$Z <- df$value

# font size
size <- 10

# assign colors to kendall tau values and plot
plt.big <- ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", na.value = "black")  +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Consistency Between Different Rankings", 
       subtitle = "Subsets of size > 6",
       fill = "Kendall's Tau") +
  theme_bw(base_size = size) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

plt.big


#-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------
# combine both plots
#-------------------------------------------------------------------------------


# combine the plot into one
combined_plot <- grid.arrange(plt.big, plt.big, ncol = 2)

plot_big_with_legend <- arrangeGrob(plt.big)
plot_with_legend <- arrangeGrob(plt)

# Combine both plots side by side
combined_plot <- grid.arrange(plot_with_legend, plot_big_with_legend, ncol = 2)



 
ggsave(file.path(script_dir, "../saved_plots/pyroCb_similarity_tests.pdf"), plot = combined_plot, width = 8.5, height = 5.55)



#-------------------------------------------------------------------------------













#-------------------------------------------------------------------------------

# store sessionInfo()
writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/pyroCb_similarity_tests.txt"))





