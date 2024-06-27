This folder contains the code for the first of two steps to reproduce the results of Section 5  of the article: 
"Large-scale Bayesian Structure Learning for Gaussian Graphical Models using Marginal Pseudo-likelihood"

This folder contains the following files:
- run_and_save.R, this is the main file. For a selected p, n, graph type, density, algorithm and replication number, it 
  i) creates Gaussian data,
  ii) runs the algorithm, and
  iii) saves the solution in an .Rdata file
- metric_functions.R, this supporting file is used by run_and_save.R to calculate the metric_functions
- MPLRJ_functions.R, this supporting file is used by run_and_save.R and contains the code to run the MPLRJ algorithm
- MPLBD_functions.R, this supporting file is used by run_and_save.R and contains the code to run the MPLBD algorithm
- SSO_functions.R, this supporting file is used by run_and_save.R and contains the code to run the SSO algorithm
- BDA_functions.R, this supporting file is used by run_and_save.R and contains the code to run the BDA algorithm
- CONCORD_functions.R, this supporting file is used by run_and_save.R and contains the code to run the CONCORD algorithm
- BConcord.cpp, this file contains C++ code that is used by CONCORD_functions.R


