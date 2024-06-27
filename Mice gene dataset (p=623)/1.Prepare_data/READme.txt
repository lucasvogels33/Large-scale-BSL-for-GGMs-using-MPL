
This folder contains the code to load and prepare the data. This is the first of three steps needed to reproduce the results 
of Subsection 6.2 ("Large-scale data application") of the article "Large-scale Bayesian Structure Learning for Gaussian Graphical 
Models using Marginal Pseudo-likelihood". 

This folder contains the following files:
- Prepare_data.R (R script). This file loads the raw data and creates three output files.
- Output files:
	- Cleaned_data.Rdata (cleaned data containing the n=653 observations of the 2.5% (p=623) variables with the highest variance)
	- Cleaned_data_small_p.Rdata (cleaned data containing the n=653 observations of the 0.5% (p=125) variables with the highest variance)
	- Gene_ID_names.Rdata (List of 24922 gene IDs with corresponding gene names)

This folder does not contain the raw data, because it is too big. To access it:
1. Go to http://rstats.immgen.org/DataPage/
2. Search for GSE15907
3. Click "GSE15907 Normalized Data"
4. Save the file as Data.xlsx in this folder
  








