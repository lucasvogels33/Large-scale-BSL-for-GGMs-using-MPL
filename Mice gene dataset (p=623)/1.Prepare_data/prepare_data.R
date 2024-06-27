#############
#Step 1: download the excel file containign the data
#       - Go to http://rstats.immgen.org/DataPage/
#       - Search for GSE15907
#       - Click "GSE15907 Normalized Data"
#       - Save excel file in the same folder as this R script

library(readxl) #for loading excel file to R
library(bioRad) #for checking NaN values in dataframe
library(tidyverse) #for ordering dataframe
library(huge) #for Gaussianization of the data

##################################################################
#####load data, logtranform and rank by variance #################
##################################################################

#load data from excel file to R

# setwd("C:/Users/lvogels/OneDrive - UvA/Documents/PhD-Project-Lucas/1_Projects/2_MPL_paper/Code_MPL_paper/Data application/Large-scale/mice/1.Prepare_data")
setwd("/Users/a.mohammadiuva.nl/Library/CloudStorage/OneDrive-UvA/PhD-Project-Lucas/1_Projects/2_MPL_paper/Code_MPL_paper/Data application/Large-scale/mice_github/1.Prepare_data")
data_orig <- read_excel("data.xlsx")

#save data_frame
data = data_orig

#test for NA or NaN
sum(is.na(data)) # test for NA (if 0, then there are no NA's)
sum(is.na.data.frame(data)) #test for NaN (if 0, then there are no NaN's)

#save names of genes (variables) and cells (obserations)
gene_ID = data$ProbeSetID #save the IDs of the 24922 genes
gene_names = data$GeneSymbol #save the names of the 24922 genes

#remove non-number columns (we want our matrix to contain only data, for the log transform)
data$ProbeSetID = NULL 
data$GeneSymbol = NULL

#save the names of the cells
cell_names = colnames(data) #save the names of the 653 cells

#check whether dimensions, max, min make sense
dim(data)

n = dim(data)[2] #n should be 653
p = dim(data)[1] #p should be 24922
max(data)
min(data)

#observe heteroskedacity
var_vec = apply(data,1,function(x) var(x))
mean_vec = apply(data,1,function(x) mean(x))
plot(NA,xlim=c(0,10000),ylim=c(0,20000000),xlab="mean of the variable",ylab="variance of the variable",main="heteroskedacity plot")
points(x=mean_vec,y=var_vec)

#log 2 transformation to avoid heteroskedacity
data = log(data,2)

#add a column with the row variance
var_vec = apply(data,1,function(x) var(x))
data = cbind(data,variance = var_vec)

#make histogram of variable variances
hist(data$variance,breaks=seq(0,13,0.25),xlab="variance")

#add removed columns back
data = cbind(data,gene_ID = gene_ID) #add the gene IDs 
data = cbind(data,gene_name = gene_names) #add the gene names

#order the rows by descending variance (top row has the highest variance)
data_ordered = data %>%
  arrange(desc(variance))  # arrange in descending order

##################################################################
#####save gene_IDs with gene_names in Rdata file  #################
##################################################################

#save ranked gene IDs and gene names in a new dataframe
gene_IDs_ordered = data_ordered$gene_ID
gene_names_ordered = data_ordered$gene_name
gene_ID_names = data.frame(gene_ID = gene_IDs_ordered,gene_name = gene_names_ordered, rank = 1:p)
filename = paste0("gene_ID_names.Rdata")
save(gene_ID_names,file=filename)

#########################################################################
############ Select 2.5% most variable genes, Gaussianize and save ######
#########################################################################

#remove non-number columns from the data 
data_ordered$gene_name = NULL 
data_ordered$gene_ID = NULL 
data_ordered$variance = NULL 

#select the top 2.5% variables
select = round(0.025*p)
data_select = data_ordered[1:select,] 

#transpose data (both normalisation (next step) and algorithms require n (rows) x p (columns) data )
data_select = t(data_select)

#verify that dimension are n x p = 653 x 623
nrow(data_select) #n should be 653
ncol(data_select) #p should be 623

#observe data is not normal
hist(data_select[,1]) 
hist(data_select[,2]) 

#Gaussianize data
data_normal = huge.npn(x=data_select,npn.func="shrinkage") 

#observe data is normal
hist(data_normal[,1]) #observe data is normal
hist(data_normal[,2]) #observe data is normal

#save cleaned data in Rdata file
output = list(data = data_normal)
filename = paste0("cleaned_data.Rdata")
save(output,file=filename)

#########################################################################
############ Select 0.5% most variable genes, Gaussianize and save ######
#########################################################################

#select the top 0.5% variables
p = dim(data_ordered)[1]
select = round(0.005*p)
data_select = data_ordered[1:select,] 

#transpose data (both normalisation (next step) and algorithms require n (rows) x p (columns) data )
data_select = t(data_select)

#verify that dimension are n x p = 653 x 125
nrow(data_select) #n should be 653
ncol(data_select) #p should be 125

#observe data is not normal
hist(data_select[,1]) 
hist(data_select[,2]) 

#Gaussianize data
data_normal = huge.npn(x=data_select,npn.func="shrinkage") 

#observe data is normal
hist(data_normal[,1]) #observe data is normal
hist(data_normal[,2]) #observe data is normal

#save cleaned data in Rdata file
output = list(data = data_normal)
filename = paste0("cleaned_data_small_p.Rdata")
save(output,file=filename)





