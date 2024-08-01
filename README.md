# Invariance and Generalization for Classification
 This repo contains the code associated with the master's thesis "Invariance and Generalization for Classification" written by Linus Kühne at the Seminar for Statistics, Department of Mathematics, ETH Zurich, in 2024. The thesis was supervised by Prof. Jonas Peters and co-supervised by Sorawit (James) Saengkyongam. 

 ## Content of the repository
 This repo contains all code used to create parts of the thesis. It also contains pre-processed data from the pyrocast database presented in (Tazi et al., [2022](https://arxiv.org/abs/2211.13052v1)). The preprocessing is explained in the thesis and in (Salas-Porras et al., [2022](https://arxiv.org/abs/2211.08883v3)). I have explicit permission by Kenza Tazi to publish these files here. 

We now explain the folder structure.

* "code" contains the functions used by the different experiments on simulated data and the pyroCb data.
    * "code_pyroCb" contains scripts implementing the invariance tests (pyroCb_invariance_tests.R) and utility functions for Stabilized Classification (pyroCb_stabilized_classification_utils.R) *for the pyroCb dataset*
    * "code_simulation" contains a script to generate data from the standard, semi-random, and random SCM (data_generating_process.R), a script implementing the different invariance tests (invariance_tests.R), a script containing the code for Stabilized Classification (stabilized_classification.R), and utility functions (utils.R). These scripts are written for the data generated from the script data_generating_process.R. For example, they assume that the first d columns of a dataframe containing a generated dataset correspond to the observations of the predictors 1, ..., d. The next column (d+1) is the response Y.
* "data" contains the pyroCb data. The folder "pyroCb_data_python" contains the data I have received from Emiliano Dìaz Salas-Porras, co-author of (Salas-Porras et al., [2022](https://arxiv.org/abs/2211.08883v3)). It is in python formats .npy and .pkl. The file exported_pyrocb.rdata contains the data transformed into a .rdata file with the script experiments_pyroCb/data_preparation/transform_data.R in this repo.
* "experiments_pyroCb" contains the scripts used for the analyses in Chapter 6 of the thesis (and a small part of Appendix B).
    * "data_preparation" contains the scripts to divide the observations into 5 or 9 environments (create_environments.R), the aforementioned script to transform the python files into .rdata files (transform_data.R), and the file with the variable screening procedure using the group lasso for logistic regression (variable_screening.R).
    * "pyroCb_ICP" contains scripts to run the ICP algorithm on the pyroCb data with different configurations and divisions into different environments. The files are called pyroCb_ICP_<test>.R, where <test> stands for the name of the invariance test used.
    * "pyroCb_stabilized_classification" contains scripts for Stabilized Classification and related methods for the pyroCb dataset. pyroCb_similarity_tests.R compares Kendall's Tau of the vectors of p-values computed by the scripts from the folder "pyroCb_ICP". pyroCb_subsets_generalization.R 

 
 Run every script with the working directory set to it's location directory. Create table with figure/table <-> script
