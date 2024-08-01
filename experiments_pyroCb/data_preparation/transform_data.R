# this script transforms the preprocessed pyroCb data from the numpy format .npy 
# and the pickle format .pkl in python to dataframes in R, stored as a .rdata file


# we use this library to be able to read the files
library(reticulate)
pd <- import("pandas")
np <- import("numpy")



# get the path of this script
script_dir <- getwd()

# load in the data and convert to dataframes and factors
cube <- as.data.frame(np$load(file.path(script_dir, "../../data/pyroCb_data_python/cube_dat.npy")))

labels <- as.factor(np$load(file.path(script_dir, "../../data/pyroCb_data_python/labels_dat.npy")))

envVars <- as.data.frame(np$load(file.path(script_dir, "../../data/pyroCb_data_python/envVars_dat.npy")))

event_df <- as.data.frame(pd$read_pickle(file.path(script_dir, "../../data/pyroCb_data_python/event_df_dat.pkl")))

# in python, the posts reference the start index and end index of the summary statistics
# for each variable. In R, we need to shift this by 1 
posts <- 1 + np$load(file.path(script_dir, "../../data/pyroCb_data_python/posts_dat.npy"))

# variable list
varss <- c('ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'u10', 'v10', 'fg10',
           'blh', 'cape', 'cin', 'z', 'slhf', 'sshf', 'w', 'u', 'v', 'cvh',
           'cvl', 'tvh', 'tvl', 'r650', 'r750', 'r850', 'uv10', 'uv250', 'alt',
           'typeH', 'typeL')

# these are the included variables (excluding tvl and tvh from the paper)
# again, needs to be shifted by 1
indxIncl <- 1+c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27, 28, 29)




save(cube, envVars, event_df, indxIncl, labels, posts, varss, file = file.path(script_dir, "../../data/exported_pyrocb.rdata"))









#-------------------------------------------------------------------------------


writeLines(capture.output(sessionInfo()), file.path(script_dir, "../sessionInfo/transform_data.txt"))

