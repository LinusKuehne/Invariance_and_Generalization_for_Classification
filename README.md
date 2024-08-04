# Invariance and Generalization for Classification
 This repo contains the code associated with the master's thesis "Invariance and Generalization for Classification" written by Linus Kühne at the Seminar for Statistics, Department of Mathematics, ETH Zurich, in 2024. The thesis was supervised by Prof. Jonas Peters and co-supervised by Sorawit (James) Saengkyongam. 

This repo contains all R code used to create parts of the thesis. It also contains pre-processed data from the pyrocast database presented in (Tazi et al., [2022](https://arxiv.org/abs/2211.13052v1)). The preprocessing is explained in the thesis and in (Salas-Porras et al., [2022](https://arxiv.org/abs/2211.08883v3)). I have explicit permission by Kenza Tazi to publish these files here. 

When running a script, make sure the working directory is set to the script's location within the folder structure of this repository. For every script we run, the output of the `sessionInfo()` command is saved with the name of the script (and a .txt ending). These files contain the loaded packages, their versions, the operating system, etc. We explain the location of these files below. 

First, we provide tables explicitly stating which experiment in the thesis corresponds to which script in this repository. The references correspond to the thesis uploaded in this repository. 

Then, we explain the folder structure of the repository. This contains essentially the same information as the previous tables, but organized in a different format.

Finally, we explain the contents of the pyroCb data stored in the file data/exported_pyrocb.rdata. 

## Explanation of the content of this repository with respect to the thesis

### Code used across different experiments

| Method | Script location/name |
| --- | --- |
| Utilities (loss functions, DAG functions, etc.) | code/code_simulations/utils.R |
| Invariance tests from Chapter 3 | code/code_simulations/invariance_tests.R |
| Stabilized classification from Subsection 4.4.3 | code/code_simulations/stabilized_classification.R |
| Data generating process (Section 5.1 and Subsection B.3.1) | code/code_simulations/data_generating_process.R |
| Invariance tests for pyroCb data | code/code_pyroCb/pyroCb_invariance_tests.R |
| Utilities for stabilized classification on pyroCb data | code/code_pyroCb/pyroCb_stabilized_classification_utils.R |


### Specific experiments

| Object in thesis | Script location/name |
| --- | --- |
| Table 4.1 | simulations/example_predicting_on_inv_sets.R |
| Fig. 5.3 | simulations/ecdf_standard.R |
| Fig. 5.4 | simulations/level_power_standard.R |
| Fig. 5.5 | simulations/power_standard.R |
| Fig. 5.6 and 5.7 | simulations/ICP_level_power_standard.R |
| Fig. 5.8 and 5.9 | simulations/stable_stabclass_standard.R |
| Fig. 6.3 and 6.4 | experiments_pyroCb/data_preparation/create_environments.R |
| Table 6.2 | experiments_pyroCb/data_preparation/variable_screening.R |
| Table 6.3 | scripts in folder experiments_pyroCb/pyroCb_ICP |
| Section 6.3 | scripts in folder experiments_pyroCb/pyroCb_ICP |
| Table 6.4 | scripts in folder experiments_pyroCb/pyroCb_stabilized_classification |
| Subsection 6.4.2 | experiments_pyroCb/pyroCb_stabilized_classification/pyroCb_stabilized_classification_oracle.R |
| Fig. 6.5 | experiments_pyroCb/pyroCb_stabilized_classification/pyroCb_similarity_tests.R |
| Fig. B.2 | simulations/ecdf_randomSCM.R |
| Fig. B.3 | simulations/level_power_random.R |
| Fig. B.4 | simulations/power_random.R |
| Fig. B.5 | simulations/ICP_level_power_bytarget.R |
| Fig. B.6 | simulations/stable_stabclass_random.R |
| Fig. B.7 | simulations/delong_not_level.R |
| Fig. B.8 | simulations/visualization_standard_SCM.R |
| Table B.1 | experiments_pyroCb/pyroCb_stabilized_classification/pyroCb_small_invariant_subsets.R |
| Fig. B.9 | simulations/similarity_tests.R |


### Data

| Kind of data | Data location/name |
| --- | --- |
| Stored data for simulation experiments for the experiments in Chapter 5, and Sections B.3 and B.4 | simulations/saved_data |
| Stored sessionInfo files for simulation experiments for the experiments in Chapter 5, and Sections B.3 and B.4 | simulations/sessionInfo |
| Generated plots for simulation experiments for the experiments in Chapter 5, and Sections B.3 and B.4 | simulations/saved_plots |
| PyroCb data provided by authors of (Salas-Porras et al., [2022](https://arxiv.org/abs/2211.08883v3)) | data/pyroCb_data_python |
| Exported pyroCb data into .rdata format with script experiments_pyroCb/data_preparation/transform_data.R | data/exported_pyrocb.rdata |
| Stored data for the experiments in Chapter 6 and Subsection B.4.3 | experiments_pyroCb/saved_data |
| Stored sessionInfo files for the experiments in Chapter 6 and Subsection B.4.3 | experiments_pyroCb/sessionInfo |
| Generated plots from Chapter 6 | experiments_pyroCb/saved_plots |




 ## Folder structure of this repository
 

* "code" contains the functions used by the different experiments on simulated data and the pyroCb data.
    * "code_pyroCb" contains scripts implementing the invariance tests (pyroCb_invariance_tests.R) and utility functions for Stabilized Classification (pyroCb_stabilized_classification_utils.R) *for the pyroCb dataset*
    * "code_simulation" contains a script to generate data from the standard, semi-random, and random SCM (data_generating_process.R), a script implementing the different invariance tests (invariance_tests.R), a script containing the code for Stabilized Classification (stabilized_classification.R), and utility functions (utils.R). These scripts are written for the data generated from the script data_generating_process.R. For example, they assume that the first d columns of a dataframe containing a generated dataset correspond to the observations of the predictors 1, ..., d. The next column (d+1) is the response Y.
* "data" contains the pyroCb data. The folder "pyroCb_data_python" contains the data I have received from Emiliano Dìaz Salas-Porras, co-author of (Salas-Porras et al., [2022](https://arxiv.org/abs/2211.08883v3)). It is in python formats .npy and .pkl. The file exported_pyrocb.rdata contains the data transformed into a .rdata file with the script experiments_pyroCb/data_preparation/transform_data.R in this repo.
* "experiments_pyroCb" contains the scripts used for the analyses in Chapter 6 of the thesis (and a small part of Appendix B).
    * "data_preparation" contains the scripts to divide the observations into 5 or 9 environments (create_environments.R), the aforementioned script to transform the python files into .rdata files (transform_data.R), and the file with the variable screening procedure using the group lasso for logistic regression (variable_screening.R).
    * "pyroCb_ICP" contains scripts to run the ICP algorithm on the pyroCb data with different configurations and divisions into different environments. The files are called pyroCb_ICP_<test>.R, where <test> stands for the name of the invariance test used.
    * "pyroCb_stabilized_classification" contains scripts for Stabilized Classification and related methods for the pyroCb dataset. pyroCb_similarity_tests.R compares Kendall's Tau of the vectors of p-values computed by the scripts from the folder "pyroCb_ICP". pyroCb_subsets_generalization.R computes the LOEO CV loss for every subset of the 13 screened variables. The other scripts compute the LOEO CV loss for the models presented in Table 6.4.
    * "saved_data" contains the data saved by all experiments from the folder "experiments_pyroCb" except those from "pyroCb_stabilized_classification_appendix". The files have the same name as the corresponding script (with the .rdata ending instead of .R).
    * "saved_plots" contains all generated plots by the experiments from the folder "experiments_pyroCb". The plots have the same name as the corresponding script (with the .pdf ending instead of .R).
    * "sessionInfo" contains the output of the R command sessionInfo() run after each experiment from the folder "experiments_pyroCb" except those from "pyroCb_stabilized_classification_appendix". The text files have the same name as the corresponding script (with the .txt ending instead of .R).


* "simulations" contains all scripts and outputs for the experiments conducted on synthetically generated data.
    * "saved_data" contains the data generated by each script. The data files have the same name as the corresponding script (with the .rdata ending instead of .R).
    * "saved_plots" contains the plots generated by each script. The plots have the same name as the corresponding script (with the .pdf ending instead of .R).
    * "sessionInfo" contains the output of the R command sessionInfo() run at the end of each script. The text files have the same name as the corresponding script (with the .txt ending instead of .R).
    * The rest of the files in this folder are R scripts for the experiments on synthetic data.
 




## Content of the pyroCb data in data/exported_pyrocb.rdata


The file data/exported_pyrocb.rdata contains the following items:

* `cube`: dataframe with 6919 rows and 318 columns. For each of 30 variables (the 28 variables from Table 6.1 plus two further variables which (Salas-Porras et al., [2022](https://arxiv.org/abs/2211.08883v3)) do not consider), there are 11 summary statistics computed over a 200 km x 200 km grid. Two variables are only associated with 6 and 4 numbers, respectively. Each row corresponds to one observation of all summary statistics of the variables.
*  `envVars`: Contains the latitude, longitude, and date for all 6919 observations.
*  `event_df`: dataframe with 6919 rows with further information corresponding to the observations. We use the column `event_df$wildfire_id`, which specifies from which wildfire a certain observation is collected. This is useful to run cross-validation on the pyroCb dataset. Since there are multiple observations collected from each wildfire at different points in time, the information in `event_df$wildfire_id` allows us to ensure that every observation from a certain wildfire is contained in the same cross-validation fold.
*  `indxIncl`: Indices of the 28 variables among the 30 available variables (Salas-Porras et al., [2022](https://arxiv.org/abs/2211.08883v3)) include in their analysis. The indices are ordered according to the occurrence of the variables in the columns of `cube`.
*  `labels`: factor of length 6919 with a value of 1 if the corresponding observation corresponds to pyroCb occurrence, and a value of 0 is no pyroCb is present.
*  `posts`: Vector of the start indices of each group of summary statistics belonging to one variable. Note that we use the R language, where indexing starts at 1. For example, `posts[1] = 1`, `posts[2] = 12`, `posts[3] = 23`. This means that the 11 summary statistics belonging to the first variable ("ch1") correspond to the rows 1 up to and including 11 in "cube". The columns for the second variable are 12 up to and including 22. Except for two variables ("typeH" and "typeL"), there are 11 summary statistics for each variable. If one wants to extract the columns of "cube" corresponding to the variables `set <- c(2,3,4,5,9,30)`, for example, the following useful snippet returns the desired column indices:

         ind.set <- as.vector(unlist(sapply(X = set, function(i) posts[i]:(posts[i+1]-1))))

* `varss`: a vector of length 30 with the names of the corresponding variables. For example, the ith variable would have the name `varss[i]`, and corresponds to the columns `cube[, posts[i]:(posts[i+1]-1)]`.







