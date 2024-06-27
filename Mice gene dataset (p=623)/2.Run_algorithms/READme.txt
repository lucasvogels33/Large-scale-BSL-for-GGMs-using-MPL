
This folder contains the code to run the algorithms. This is the second of three steps needed to reproduce the results 
of Subsection 6.2 of the article "Large-scale Bayesian Structure Learning for Gaussian Graphical Models using Marginal Pseudo-likelihood". 

This folder contains the following files:
- Cleaned_data.Rdata (cleaned data containing the n=653 observations of the 2.5% (p=623) variables with the highest variance. This file is created in the first step)
- Solve_data.R (the main Rscript to run all the algorithms. It loads the data, runs the algorithms, and saves the results in .Rdata files)
- Supporting_functions.R (Rscript containing functions used by solve_data.R)
- BConcord.Cpp (C++ file containing the script for the CONCORD method)

  








