#set working directory (this is where the output files will be pasted)
#setwd()

#select what algorithm you want to run on what instance
p = 10 #set number of variables (either 10,100,1000)
n = 20 #set number of observations 
graph = "random" #set graph type (either "random", "cluster" or "scale-free")
density = "sparse" #set density (either "sparse" or "dense")
algorithm = "CONCORD" #set algorithm (either "MPLRJ", "MPLBD", "SSO", "BDA", or "CONCORD")
rep = 1 #set the replication number 


#--download the necessary libraries and load the run_experiments file --
library(PRROC) #to calculate the AUC PR
library(BDgraph) #for MPL-RJ, MPL-BD and BD
library(ssgraph) #for ss
library(Rcpp) #for CONCORD
library(RcppArmadillo) #for CONCORD

# download the supporting functions
source("1.MPLRJ_functions.R")
source("2.MPLBD_functions.R")
source("3.BDA_functions.R")
source("4.SSO_functions.R")
source("5.CONCORD_functions.R")
source("metric_functions.R")
Rcpp::sourceCpp('BConcord.cpp')


#--set value of "size" ---#
size = NULL
if ((p==10)*(graph=="random")*(density=="sparse")){size = 5}
if ((p==10)*(graph=="random")*(density=="dense")){size = 20}
if ((p==10)*(graph=="cluster")*(density=="sparse")){size = c(2,3)}
if ((p==10)*(graph=="cluster")*(density=="dense")){size = c(10,10)}

if ((p==100)*(graph=="random")*(density=="sparse")){size = 50}
if ((p==100)*(graph=="random")*(density=="dense")){size = 248}
if ((p==100)*(graph=="cluster")*(density=="sparse")){size = c(25,25)}
if ((p==100)*(graph=="cluster")*(density=="dense")){size = c(124,124)}

if ((p==1000)*(graph=="random")*(density=="sparse")){size = 2498}
if ((p==1000)*(graph=="random")*(density=="dense")){size = 24975}
if ((p==1000)*(graph=="cluster")*(density=="sparse")){size = rep(312,8)}
if ((p==1000)*(graph=="cluster")*(density=="dense")){size = rep(3122,8)}

#set remaining parameters
type = "Gaussian"
prob = 0.1 #this value does not matter but needs to be specified
vis = FALSE

#=--set the iterations at which the metrics need to be outputted
if ((p==10)*(algorithm=="MPLRJ")){iter_vec_thin =c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,35000,40000,45000,50000,60000,70000,80000,90000,100000)}
if ((p==10)*(algorithm=="MPLBD")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,7500,10000,15000,20000,25000,30000) }
if ((p==10)*(algorithm=="BDA")){iter_vec_thin = c(10,25,50,75,100,250,500,750,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000) }
if ((p==10)*(algorithm=="SSO")){iter_vec_thin =c(10,20,30,40,50,75,100,250,500,750,1000,1500,2000,2500,3000)}
if ((p==10)*(algorithm=="CONCORD")){iter_vec_thin =c(10,20,30,40,50,75,100,250,500,750,1000,1500,2000,2500,3000,3500,4000,4500,5000,6000,7000,8000,9000,10000)}

if ((p==100)*(algorithm=="MPLRJ")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,150000,200000,300000,400000,500000,750000,1000000,1500000,2000000,2500000,3000000,4000000,5000000,7500000,10000000,20000000,30000000,40000000,50000000,60000000,70000000,80000000,90000000,100000000,125000000)}
if ((p==100)*(algorithm=="MPLBD")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,250000,500000,1000000,1500000,2000000,2500000)}
if ((p==100)*(algorithm=="BDA")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,7500,10000,15000,20000,25000,30000)}
if ((p==100)*(algorithm=="SSO")){iter_vec_thin = c(2,4,6,8,10,20,30,40,50,75,100,250,500,750,1000,1500,2000,2500,3000,4000,5000,7500,10000,15000,20000,25000,30000,35000,40000,45000)}
if ((p==100)*(algorithm=="CONCORD")){iter_vec_thin =c(10,20,30,40,50,75,100,250,500,750,1000,1500,2000,2500,3000,3500,4000,4500,5000,6000,7000,8000,9000,10000,15000,20000,30000,40000,50000,75000,100000,150000,200000,250000,300000,350000,400000)}

#p=1000, sparse, n = 400
if ((p==1000)*(n==400)*(algorithm=="MPLRJ")*(density=="sparse")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,150000,200000,250000,500000,750000,1000000,1250000,1500000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000,12500000,15000000,16000000)}
if ((p==1000)*(n==400)*(algorithm=="MPLBD")*(density=="sparse")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,125000,150000,175000,200000,250000,300000)}
if ((p==1000)*(n==400)*(algorithm=="SSO")*(density=="sparse")){iter_vec_thin = c(2,4,6,8,10,15,20,30,40,50,75,100,150,200,250,300,350,400,450,500,550,600)}
if ((p==1000)*(n==400)*(algorithm=="CONCORD")*(density=="sparse")){iter_vec_thin = c(5,10,50,100,150,200,250,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500,3000,4000,5000,6000,7000,8000,9000,10000,12500,15000,17500,20000,25000,30000,35000,40000)}
if ((p==1000)*(n==400)*(algorithm=="BDA")*(density=="sparse")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,125000,150000,175000,200000,250000,300000)}

#p=1000, sparse, n = 1050
if ((p==1000)*(n==1050)*(algorithm=="MPLRJ")*(density=="sparse")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,150000,200000,250000,500000,750000,1000000,1250000,1500000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000,12500000,15000000,16000000)}
if ((p==1000)*(n==1050)*(algorithm=="MPLBD")*(density=="sparse")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,125000,150000,175000,200000,250000,300000)}
if ((p==1000)*(n==1050)*(algorithm=="SSO")*(density=="sparse")){iter_vec_thin = c(2,4,6,8,10,15,20,30,40,50,75,100,150,200,250,300,350,400)}
if ((p==1000)*(n==1050)*(algorithm=="CONCORD")*(density=="sparse")){iter_vec_thin = c(5,10,50,100,150,200,250,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500,3000,4000,5000,6000,7000,8000,9000,10000,12500,15000,17500,20000,25000,30000,35000,40000)}
if ((p==1000)*(n==1050)*(algorithm=="BDA")*(density=="sparse")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,125000,150000,175000,200000,250000,300000)}


#p=1000, dense and n=400
if ((p==1000)*(n==400)*(algorithm=="MPLRJ")*(density=="dense")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,150000,200000,250000,500000,750000,1000000,1250000,1500000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000,12500000,15000000,20000000,25000000,30000000)}
if ((p==1000)*(n==400)*(algorithm=="MPLBD")*(density=="dense")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,125000,150000,175000,200000,250000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000)}
if ((p==1000)*(n==400)*(algorithm=="SSO")*(density=="dense")){iter_vec_thin = c(2,4,6,8,10,15,20,30,40,50,75,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,1100,1200,1300,1400,1500)}
if ((p==1000)*(n==400)*(algorithm=="CONCORD")*(density=="dense")){iter_vec_thin = c(5,10,50,100,150,200,250,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500,3000,4000,5000,6000,7000,8000,9000,10000,12500,15000,17500,20000,25000,30000,40000,50000,60000,70000,80000,90000,100000,125000,150000,175000,200000,225000,250000)}
if ((p==1000)*(n==400)*(algorithm=="BDA")*(density=="dense")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,125000,150000,175000,200000,250000,300000)}

#p=1000, dense and n=1050
if ((p==1000)*(n==1050)*(algorithm=="MPLRJ")*(density=="dense")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,150000,200000,250000,500000,750000,1000000,1250000,1500000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000,12500000,15000000,17500000,20000000,25000000,30000000)}
if ((p==1000)*(n==1050)*(algorithm=="MPLBD")*(density=="dense")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,125000,150000,175000,200000,250000,300000,400000,500000)}
if ((p==1000)*(n==1050)*(algorithm=="SSO")*(density=="dense")){iter_vec_thin = c(2,4,6,8,10,15,20,30,40,50,75,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,1100,1200,1300,1400,1500)}
if ((p==1000)*(n==1050)*(algorithm=="CONCORD")*(density=="dense")){iter_vec_thin = c(5,10,50,100,150,200,250,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500,3000,4000,5000,6000,7000,8000,9000,10000,12500,15000,17500,20000,25000,30000,35000,40000,45000,50000,60000,70000,80000)}
if ((p==1000)*(n==1050)*(algorithm=="BDA")*(density=="dense")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,125000,150000,175000,200000,250000,300000)}


#p=1000, scale-free
if ((graph=="scale-free")*(p==1000)*(algorithm=="MPLRJ")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,150000,200000,250000,500000,750000,1000000,1250000,1500000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000,12500000,15000000,17500000,20000000,25000000,30000000)}
if ((graph=="scale-free")*(p==1000)*(algorithm=="MPLBD")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,125000,150000,175000,200000)}
if ((graph=="scale-free")*(p==1000)*(algorithm=="SSO")){iter_vec_thin = c(2,4,6,8,10,15,20,30,40,50,75,100,150,200,250,300,350,400,450,500,550,600)}
if ((graph=="scale-free")*(p==1000)*(algorithm=="CONCORD")){iter_vec_thin = c(5,10,50,100,150,200,250,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000,17500,20000,25000,30000,35000,40000,45000,50000)}
if ((graph=="scale-free")*(p==1000)*(algorithm=="BDA")){iter_vec_thin = c(10,50,100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,50000,75000,100000,125000,150000,175000,200000,250000,300000)}

#create data
set.seed(rep)
data.sim = bdgraph.sim( p = p, n = n, graph = graph, prob = prob, size = size, type = type, vis = vis )
data.sim$density = density
data.sim$rep = rep
data_output = list(G = data.sim$G,p=p,n=n,graph = graph,density=density,rep=rep)

#save true graph (all algorithms use the same data, hence we only save it for glasso to avoid duplicates)
#if (algorithm=="MPLRJ"){
#    filename = paste0("truth_p",p,"_n",n,"_",graph,"_",density,"_rep",rep,".Rdata")
#    save(data_output,file = filename)
#}

#solve data

if (algorithm=="MPLRJ"){output = MPLRJ_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}
if (algorithm=="MPLBD"){output = MPLBD_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}
if (algorithm=="BDA"){output = BDA_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}
if (algorithm=="SSO"){output = SSO_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}
if (algorithm=="CONCORD"){output = CONCORD_solve(data=data.sim,iter_vec_thin=iter_vec_thin)}

#save solution
filename = paste0("results_p",p,"_n",n,"_",graph,"_",density,"_",algorithm,"_rep",rep,".Rdata")
save(output,file = filename)











