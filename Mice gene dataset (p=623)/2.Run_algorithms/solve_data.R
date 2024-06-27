#--set main parameters
algorithm = "MPLBD" #other options are "MPLRJ", "SS", "BDA" or "CONCORD"
g.prior = 0.01 

#--download the necessary libraries ---#
library(BDgraph) #for MPL-RJ, MPL-BD and BD
library(ssgraph) #for ss
library(Rcpp) #for CONCORD
library(RcppArmadillo) #for CONCORD

#--load supporting functions---#
source("supporting_functions.R")
Rcpp::sourceCpp('BConcord.cpp') #code for CONCORD method

#--load data ----------#
filename=paste0("cleaned_data.Rdata")
load(file=filename)
data = output$data

#--run algorithms ----#

#MPLBD
if (algorithm=="MPLBD"){
  iter_vec_thin = c(10,50,75,100,200,300,400,500,600,750,1000,2000,3000,5000,7500,10000,15000,20000,25000,30000,40000,50000,60000,75000,100000,150000,200000,250000,300000,400000,500000,750000,1000000,1250000,1500000,1750000,2000000)
  burnin_iter_vec_thin = iter_vec_thin
  MCMCstart = "empty"
  output = MPLBD_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart,g.prior=g.prior)
  filename=paste0("final_",g.prior,"_",algorithm,".Rdata")
  save(output,file = filename)
}

#MPLRJ
if (algorithm=="MPLRJ"){
  iter_vec_thin = c(10,50,75,100,150,200,300,400,500,600,750,1000,2000,3000,5000,7500,10000,25000,50000,75000,100000,250000,500000,1000000,2000000,4000000,7000000,10000000,15000000,20000000,25000000,30000000)
  burnin_iter_vec_thin = iter_vec_thin
  MCMCstart = "empty"
  output = MPLRJ_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart,g.prior=g.prior)
  filename=paste0("final_",g.prior,"_",algorithm,".Rdata")
  save(output,file = filename)
}

#BDA
if (algorithm=="BDA"){
  iter_vec_thin = c(10,20,30,40,50,60,70,80,90,100,1000,10000)
  burnin_iter_vec_thin = iter_vec_thin
  MCMCstart = "empty"
  output = BDA_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart,g.prior=g.prior)
  filename=paste0("final_",g.prior,"_",algorithm,".Rdata")
  save(output,file = filename)
}

#SS
if (algorithm=="SS"){
  iter_vec_thin = c(2,4,6,8,10,15,20,30,40,50,75,100,150,250,500,750,1000,1250,1500)
  burnin_iter_vec_thin = iter_vec_thin
  MCMCstart = "empty"
  output = SS_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart,g.prior=g.prior)
  filename=paste0("final_",g.prior,"_",algorithm,".Rdata")
  save(output,file = filename)
}

#CONCORD
if (algorithm=="CONCORD"){
  iter_vec_thin = c(10,20,30,50,100,150,250,500,750,1000,2000,3000,5000,7500,10000,15000,20000,25000,30000,40000,50000,60000,75000,100000,150000,200000,250000)
  burnin_iter_vec_thin = iter_vec_thin
  MCMCstart = "empty"
  output = CONCORD_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart,g.prior=g.prior)
  filename=paste0("final_",g.prior,"_",algorithm,".Rdata")
  save(output,file = filename)
}










