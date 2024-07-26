# in this script, we evaluate to which degree the ranking of subsets according to 
# the p-values by the different tests is similar


library(rje)
library(ggplot2)
library(reshape2)


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

# sets to check stability
sets <- powerSet(varincl)
sets[[1]] <- c(0)



df.p <- data.frame("Neg. LOEO CV loss" = -wbce.worst,
                   "DeLong_cont" = pvals.delong.cont.1tail,
                   "DeLong_5" = pvals.delong.5.1tail,
                   "DeLong_9" = pvals.delong.9.1tail,
                   "GCM_5" = pvals.tram.5, 
                   "GCM_9" = pvals.tram.9,
                   "Corr_5" = pvals.corr.5,
                   "Corr_9" = pvals.corr.9,
                   "Residual_5" = pvals.res.5,
                   "Residual_9" = pvals.res.5)


tau.mat <- matrix(NA, nrow = ncol(df.p), ncol = ncol(df.p))




for(j in 1:(ncol(df.p)-1)){
  for(i in (j+1):ncol(df.p)){
    print(i)
    tau.mat[i,j] <- cor(x = df.p[,i], y = df.p[,j], method = "kendall")
  }
}




mat <- tau.mat

mat[is.na(mat)] <- 0


# convert to a symmetric matrix
mat <- mat + t(mat)

diag(mat) <- NA

df <- melt(mat)


# generate grid
x <- c("Neg. LOEO CV loss", "DeLong_cont", "DeLong_5", "DeLong_9", "GCM_5", "GCM_9", "Corr_5", "Corr_9", "Residual_5", "Residual_9")
y <- c("Neg. LOEO CV loss", "DeLong_cont", "DeLong_5", "DeLong_9", "GCM_5", "GCM_9", "Corr_5", "Corr_9", "Residual_5", "Residual_9")

data <- expand.grid(X=x, Y=y)
data$Z <- df$value

size <- 10

# assign colors to kendall tau values and plot
plt <- ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue", na.value = "black")  +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Kendall's Tau between different rankings", fill = "Kendall's Tau") +
  theme_bw(base_size = size) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_blank())

plt



ggsave(file.path(script_dir, "saved_plots/pyroCb_similarity_tests.pdf"), width = 6, height = 5.07)






#-------------------------------------------------------------------------------

# store sessionInfo()
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/pyroCb_similarity_tests.txt"))





