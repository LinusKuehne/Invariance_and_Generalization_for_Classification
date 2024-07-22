# In this script, we visualize a sample from the standard SCM


library(ggplot2)
library(GGally)



# get the path of this script
script_dir <- getwd()

# load in the functions needed
source(file.path(script_dir, "../code/code_simulations/data_generating_process.R"))



set.seed(1)


sample <- gen.sample.standard(n=100)
sample$Y <- factor(sample$Y, levels = c(0,1))
levels(sample$Env) <- c("0", "1", "2", "-1", "-2")
sample$E <- sample$Env
sample$Env <- NULL


plt <- ggpairs(sample, 
        diag = list(mapping = aes(color = E, alpha = 0.4)),
        lower = list(
          continuous = "density",
          combo = "box",
          mapping = aes(color = E)
        ),
        upper = list(
          continuous = "autopoint",
          combo = "autopoint",
          discrete = "autopoint",
          mapping = aes(color = E)
        ))
plt + theme_bw(base_size = 10)






ggsave(filename = file.path(script_dir, "saved_plots/visualization_standard_SCM.pdf"), width = 6, height = 6)

# store the sessionInfo:
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo/visualization_standard_SCM.txt"))



