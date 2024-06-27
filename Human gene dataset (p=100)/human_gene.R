
#################################################
#####download libraries and supporting files#####
#################################################

#set working directory
setwd("C:/Users/lvogels/OneDrive - UvA/Documents/PhD-Project-Lucas/1_Projects/2_MPL_paper/Code_MPL_paper/Data application/Medium-scale")

#download libraries
library(BDgraph) #for BD-MPL, RJ-MPL, BD
library(ssgraph) #for SS
library(huge) #for Gaussianization of the data
library(Rcpp) #for CONCORD
library(RcppArmadillo) #for CONCORD
Rcpp::sourceCpp('BConcord.cpp') #load c++ file containing the CONCORD code

#######################################
#####load data  ########################
######################################

data(geneExpression) #load data from BDgraph
gene_data = geneExpression 

#############################################################
#####transform non-normal data into normal data using npn ###
#############################################################

#observe that the data is not normally distributed
hist(geneExpression[,1])
hist(geneExpression[,2])

#transform data
new_data = huge.npn(x=geneExpression,npn.func="shrinkage") #new data is n x p matrix containing normalized data

#observe that the transformed data is normally distribured
hist(new_data[,1])
hist(new_data[,2])

#############################################################
################define all the functions################
#############################################################

#- we create 5 functions (one for each algorithm)

#Input:
#- the normalized data 
#- the amount of burnin iterations, specified as a vector at which to output the running time and #edges in the Markov state
#- the amount of MCMC iterations after burnin, specified as a vector at which to output the running time and #edges in the Markov state
#- the initital graph

#Output:
#- burnin_iter_vec_thin: the amount of burnin iterations, specified as a vector at which to output the running time and #edges in the Markov state
#- iter_vec_thin: the amount of MCMC iterations after burnin, specified as a vector at which to output the running time and #edges in the Markov state
#- edge_vec: contains #edges at every iteration mentioned in burnin_iter_vec_thin and iter_vec_thin
#- time_vec: contains runtime at every iteration mentioned in burnin_iter_vec_thin and iter_vec_thin
#- plinks: matrix containing the edge inclusion probabilities


MPLBD_solve = function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart){
  
  #parameters
  burnin = 0
  jump = 1
  save = FALSE
  cores = 1
  g.prior = 0.2
  verbose = TRUE
  
  #intialize outputs
  len_iter = length(iter_vec_thin)
  len_burnin = length(burnin_iter_vec_thin)
  edge_vec = rep(0,len_iter+len_burnin)
  time_vec = rep(0,len_iter+len_burnin)
  
  #run the burnin iterations
  for (j in 1:len_burnin){
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = burnin_iter_vec_thin[j-1]
    }
    newiter = burnin_iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_mpl_bd  = bdgraph.mpl( data = data, algorithm = "bdmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #calculate the runtime
    time_init = 0
    time_end = 0
    if (j==1){time_init = sample_mpl_bd$time_init}
    if (j==len_burnin){time_end = sample_mpl_bd$time_end}
    runtime = time_init + sample_mpl_bd$time_all_iterations + time_end
    
    #save metrics
    edge_vec[j] = sum(sample_mpl_bd$last_graph)/2
    time_vec[j] = runtime
    
    #save starting point next run
    MCMCstart = sample_mpl_bd$last_graph
    
  }
  
  #run the MCMC iterations (after burnin)
  weights_old = 0
  plinks_old = 0
  for (j in 1:len_iter){
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_mpl_bd  = bdgraph.mpl( data = data, algorithm = "bdmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #calculate new p matrix
    weights_new = sample_mpl_bd$sum_weights
    plinks_new = (weights_old*plinks_old + weights_new*sample_mpl_bd$p_links)/(weights_old+weights_new)
    
    #calculate the runtime
    time_init = 0
    time_end = 0
    if (j==1){time_init = sample_mpl_bd$time_init}
    if (j==len_iter){time_end = sample_mpl_bd$time_end}
    runtime = time_init + sample_mpl_bd$time_all_iterations + time_end
    
    #save metrics
    edge_vec[len_burnin+j] = sum(sample_mpl_bd$last_graph)/2
    time_vec[len_burnin+j] = runtime
    
    #save data for next run
    plinks_old = plinks_new
    weights_old = weights_old+weights_new
    MCMCstart = sample_mpl_bd$last_graph
    
  }
  
  
  plinks = plinks_new
  time_vec = cumsum(time_vec)
  time_vec = time_vec/(10^6)
  return(list(edge_vec=edge_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,burnin_iter_vec_thin=burnin_iter_vec_thin,plinks=plinks_new))
} 
MPLRJ_solve = function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart){
  
  #parameters
  burnin = 0
  jump = 1
  save = FALSE
  cores = 1
  g.prior = 0.2
  verbose = TRUE
  
  #intialize outputs
  len_iter = length(iter_vec_thin)
  len_burnin = length(burnin_iter_vec_thin)
  edge_vec = rep(0,len_iter+len_burnin)
  time_vec = rep(0,len_iter+len_burnin)
  
  for (j in 1:len_burnin){
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = burnin_iter_vec_thin[j-1]
    }
    newiter = burnin_iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_mpl_rj = bdgraph.mpl( data = data, algorithm = "rjmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #calculate the runtime
    time_init = 0
    time_end = 0
    if (j==1){time_init = sample_mpl_rj$time_init}
    if (j==len_burnin){time_end = sample_mpl_rj$time_end}
    runtime = time_init + sample_mpl_rj$time_all_iterations + time_end
    
    #save metrics
    edge_vec[j] = sum(sample_mpl_rj$last_graph)/2
    time_vec[j] = runtime
    
    #save data for next run
    MCMCstart = sample_mpl_rj
    
  }
  
  #run MCMC iterations (after burnin)
  plinks_old = 0
  for (j in 1:len_iter){
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_mpl_rj = bdgraph.mpl( data = data, algorithm = "rjmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #calculate new p matrix
    plinks_new = (olditer*plinks_old + iter*sample_mpl_rj$p_links)/(newiter)
    
    #calculate the runtime
    time_init = 0
    time_end = 0
    if (j==1){time_init = sample_mpl_rj$time_init}
    if (j==len_iter){time_end = sample_mpl_rj$time_end}
    runtime = time_init + sample_mpl_rj$time_all_iterations + time_end
    
    #save metrics
    edge_vec[len_burnin+j] = sum(sample_mpl_rj$last_graph)/2
    time_vec[len_burnin+j] = runtime
    
    #save data for next run
    plinks_old = plinks_new
    MCMCstart = sample_mpl_rj
    
  }
  
  time_vec = cumsum(time_vec)
  time_vec = time_vec/(10^6)
  return(list(edge_vec=edge_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,burnin_iter_vec_thin=burnin_iter_vec_thin,plinks=plinks_new))
} 
SS_solve= function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart){
  
  #obtain p and n
  p = ncol(data)
  n = nrow(data)
  
  #parameters
  burnin = 0
  var1 = 0.02
  var2 = 2
  lambda = 1
  save = FALSE
  cores = 1
  sig.start = NULL
  if (p > n){sig.start = diag(1,p)}
  g.prior = 0.2
  verbose = FALSE
  
  #intialize outputs
  len_iter = length(iter_vec_thin)
  len_burnin = length(burnin_iter_vec_thin)
  edge_vec = rep(0,len_iter+len_burnin)
  time_vec = rep(0,len_iter+len_burnin)
  
  #run burnin_iterations
  for (j in 1:len_burnin){
    print(burnin_iter_vec_thin[j]) #output progress
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = burnin_iter_vec_thin[j-1]
    }
    newiter = burnin_iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_ss  = ssgraph( data = data,iter = iter, burnin = burnin, var1 = var1, var2=var2,lambda=lambda,g.start=MCMCstart,sig.start=sig.start,save=save,cores=cores,g.prior=g.prior,verbose=verbose)
    
    #calculate the runtime
    time_init = 0
    if (j==1){time_init = sample_ss$time_init}
    runtime = time_init + sample_ss$time_all_iterations
    
    #save metrics
    edge_vec[j] = sum(sample_ss$last_graph)/2
    time_vec[j] = runtime
    
    #save data for next run
    MCMCstart = sample_ss
    
  }
  
  #run MCMC iterations after burnin
  plinks_old = 0
  for (j in 1:len_iter){
    print(iter_vec_thin[j]) #output progress
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_ss  = ssgraph( data = data,iter = iter, burnin = burnin, var1 = var1, var2=var2,lambda=lambda,g.start=MCMCstart,sig.start=sig.start,save=save,cores=cores,g.prior=g.prior,verbose=verbose)
    
    #calculate new p matrix
    plinks_new = (olditer*plinks_old + iter*sample_ss$p_links)/(newiter)
    
    #calculate the runtime
    time_init = 0
    if (j==1){time_init = sample_ss$time_init}
    runtime = time_init + sample_ss$time_all_iterations
    
    #save metrics
    edge_vec[len_burnin + j] = sum(sample_ss$last_graph)/2
    time_vec[len_burnin + j] = runtime
    
    #save data for next run
    plinks_old = plinks_new
    MCMCstart = sample_ss
    
  }
  
  
  time_vec = cumsum(time_vec)
  time_vec = time_vec/(10^6)
  return(list(edge_vec=edge_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,burnin_iter_vec_thin=burnin_iter_vec_thin,plinks=plinks_new))
}
CONCORD_solve = function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart){
  
  #obtain amount of variables
  n = nrow(data)
  p = ncol(data)
  S = crossprod(data)/n
  
  #parameters
  burnin = 0
  omega_start = matrix(rnorm(p*p,mean=0,sd=1), p, p)
  
  #intialize outputs
  len_iter = length(iter_vec_thin)
  len_burnin = length(burnin_iter_vec_thin)
  edge_vec = rep(0,len_iter+len_burnin)
  time_vec = rep(0,len_iter+len_burnin)
  
  for (j in 1:len_burnin){
    print(burnin_iter_vec_thin[j])
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = burnin_iter_vec_thin[j-1]
    }
    newiter = burnin_iter_vec_thin[j]
    iter = newiter-olditer
    
    #run CONCORD
    start = proc.time()[3]
    sample_concord  = BCONCORDSS(S=S, nmc=iter, burnin=burnin, n=n, omega_start=omega_start)
    end = proc.time()[3]
    runtime = as.numeric(end-start)
    
    #save metrics
    last_omega = sample_concord$last_omega
    last_omega_vec = last_omega[upper.tri(last_omega)]
    edge_vec[j] = length(which(last_omega_vec!=0))
    time_vec[j] = runtime
    
    #save starting point for next run
    omega_start = sample_concord$last_omega
    
  
  }
  
  #run MCMC iterations after burnin
  plinks_old = 0
  for (j in 1:len_iter){
    print(iter_vec_thin[j])
    
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run CONCORD
    start = proc.time()[3]
    sample_concord  = BCONCORDSS(S=S, nmc=iter, burnin=burnin, n=n, omega_start=omega_start)
    end = proc.time()[3]
    runtime = as.numeric(end-start)
    
    #calculate new p matrix
    plinks_new = (olditer*plinks_old + iter*sample_concord$p_links)/(newiter)
    
    #save metrics
    last_omega = sample_concord$last_omega
    last_omega_vec = last_omega[upper.tri(last_omega)]
    edge_vec[len_burnin+j] = length(which(last_omega_vec!=0))
    time_vec[len_burnin+j] = runtime
    
    #store data for next run
    plinks_old = plinks_new
    omega_start = sample_concord$last_omega
    
    
  }
  
  time_vec = cumsum(time_vec)
  return(list(edge_vec=edge_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,burnin_iter_vec_thin=burnin_iter_vec_thin,plinks=plinks_new))
} 
BDA_solve = function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart){
  
  #parameters
  burnin = 0
  jump = 1
  save = FALSE
  cores = 1
  g.prior = 0.2
  verbose = TRUE
  
  #run the algorithm, produce output at every iteration in the iter_vec_thin vector
  len_iter = length(iter_vec_thin)
  len_burnin = length(burnin_iter_vec_thin)
  edge_vec = rep(0,len_iter+len_burnin)
  time_vec = rep(0,len_iter+len_burnin)
  
  for (j in 1:len_burnin){
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = burnin_iter_vec_thin[j-1]
    }
    newiter = burnin_iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_bd  = bdgraph( data = data,algorithm = "bdmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #calculate the runtime
    time_init = 0
    time_end = 0
    if (j==1){time_init = sample_bd$time_init}
    if (j==len_burnin){time_end = sample_bd$time_end}
    runtime = time_init + sample_bd$time_all_iterations + time_end
    
    #save metrics
    edge_vec[j] = sum(sample_bd$last_graph)/2
    time_vec[j] = runtime
    
    #save data for next run
    MCMCstart = sample_bd
    
  }

  #run MCMC iterations (after burnin)
  weights_old = 0
  plinks_old = 0
  for (j in 1:len_iter){
    
    #determine amount of iterations
    olditer = 0
    if (j > 1){
      olditer = iter_vec_thin[j-1]
    }
    newiter = iter_vec_thin[j]
    iter = newiter-olditer
    
    #run bdgraph for maxiter iterations
    sample_bd  = bdgraph( data = data,algorithm = "bdmcmc", iter = iter, burnin = burnin, jump = jump, save = save,cores=cores,g.start=MCMCstart,g.prior=g.prior,verbose=verbose)  
    
    #calculate new p matrix
    weights_new = sample_bd$sum_weights
    plinks_new = (weights_old*plinks_old + weights_new*sample_bd$p_links)/(weights_old+weights_new)
    
    #calculate the runtime
    time_init = 0
    time_end = 0
    if (j==1){time_init = sample_bd$time_init}
    if (j==len_iter){time_end = sample_bd$time_end}
    runtime = time_init + sample_bd$time_all_iterations + time_end
    
    #save metrics
    edge_vec[len_burnin+j] = sum(sample_bd$last_graph)/2
    time_vec[len_burnin+j] = runtime
    
    #save data for next run
    plinks_old = plinks_new
    weights_old = weights_old+weights_new
    MCMCstart = sample_bd
    
  }

  plinks = plinks_new
  time_vec = cumsum(time_vec)
  time_vec = time_vec/(10^6)
  return(list(edge_vec=edge_vec,time_vec=time_vec,iter_vec_thin=iter_vec_thin,burnin_iter_vec_thin=burnin_iter_vec_thin,plinks=plinks_new))
} 

#############################################################
################run algorithms ##############################
#############################################################

#settings
data = new_data
MCMCstart = "empty" #CONCORD will always start at the full graph

#run algorithms
algorithm = "MPLRJ"
iter_vec_thin=c(20,50,100,200,500,750,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000,30000,40000,50000,75000,100000,125000,150000,175000,200000,250000,300000,350000,400000,500000,600000,700000,800000,900000,1000000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000,12500000,15000000)
burnin_iter_vec_thin = iter_vec_thin
output = MPLRJ_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart)
filename=paste0("gene_solution_",algorithm,".Rdata")
save(output,file = filename)

algorithm = "MPLBD"
iter_vec_thin=c(20,50,100,200,500,750,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,12500,15000,17500,20000,22500,25000,27500,30000,35000,40000,45000,50000,60000,70000,80000,90000,100000,150000,200000,250000,300000,350000,400000,450000,500000)
burnin_iter_vec_thin = iter_vec_thin
output = MPLBD_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart)
filename=paste0("gene_solution_",algorithm,".Rdata")
save(output,file = filename)

algorithm = "BD" #this algorithm is slow. Better to run it on a cluster computer.
iter_vec_thin= c(10,20,50,100,200,500,750,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,12500,15000,17500,20000,22500,25000,27500,30000,35000,40000,45000,50000)
burnin_iter_vec_thin = iter_vec_thin
output = BDA_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart)
filename=paste0("gene_solution_",algorithm,".Rdata")
save(output,file = filename)

algorithm = "SS"
iter_vec_thin=c(1,2,3,4,6,8,10,15,20,30,40,50,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,2500,5000,7500,10000)
burnin_iter_vec_thin=iter_vec_thin
output = SS_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart)
filename=paste0("gene_solution_",algorithm,".Rdata")
save(output,file = filename)

algorithm = "CONCORD"
iter_vec_thin=c(5,10,20,30,40,50,100,200,500,750,1000,1500,2000,2500,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,30000,40000,50000,60000,70000,80000,90000,100000)
burnin_iter_vec_thin=iter_vec_thin
output = CONCORD_solve(data=data,burnin_iter_vec_thin=burnin_iter_vec_thin,iter_vec_thin=iter_vec_thin,MCMCstart=MCMCstart)
filename=paste0("gene_solution_",algorithm,".Rdata")
save(output,file = filename)

#############################################################
################make edge convergence plot ##################
#############################################################

par(mfrow=c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1) )
algorithm_vec = c("MPLRJ","MPLBD","BD","SS","CONCORD")
col_vec = c(1,2,3,4,5)
lty_vec = c(1,2,3,4,5)

#standard x-axis
log_on = FALSE
if (log_on){
  log = 'x'
  xlim = c(0.001,2000)
  xlab = " log time (seconds)"
} else {
  log = ""
  xlim = c(0,1000)
  xlab = "time (seconds)"
}
plot(NA,xlim=xlim,ylim=c(0,600),xlab=xlab,log=log,ylab="edges",main="")
count = 0
for (algorithm in algorithm_vec){
  
  
  #load data
  filename = paste0("gene_solution_",algorithm,".Rdata")
  load(file=filename)
  edge_vec=output$edge_vec
  time_vec = output$time_vec
  
  #add the point at time t=0 to the vectors
  if (algorithm=="CONCORD"){
    edge_vec=c(4950,edge_vec)
    time_vec=c(0,time_vec)
  }
  else {
    edge_vec=c(0,edge_vec)
    time_vec=c(0,time_vec)
  }
  
  #set layout
  count = count + 1
  col = col_vec[count]
  lty = lty_vec[count]
  
  #add line
  points(x=time_vec,y=edge_vec,col=col,lty=lty,lw = 2,type="l")
}
legend( "bottomright", inset=.03, algorithm_vec, col = col_vec,lw = c(2,2,2,2,2),lty = lty_vec, cex = 1, box.lty = 0 )

#############################
#############networks########
#############################

library(igraph)
algorithm_vec = c("MPLBD","MPLRJ","BD","SS","CONCORD")
count = 0

top_links = 100
cutoff = 0.9 
edge_select = "cutoff" #either cutoff or top_links
G_sum = matrix(0,100,100)

#load data and create undirected graphs for all algorithms 
for (algorithm in algorithm_vec){
  filename = paste0("gene_solution_",algorithm,".Rdata")
  load(file=filename)
  plinks = output$plinks
  
  #fill the lower triangle of the matrix
  p = 100
  count_2 = 0
  for (i in 2:p){ #fill lower matrix
    for (j in 1:(i-1)){
      count_2 = count_2 + 1
      plinks[i,j] = plinks[j,i]
    } 
  }
  
  #determine x links with highest edge incl. prob.
  if (edge_select == "top_links"){
    nr_links_show = top_links
    plinks_sort = sort(plinks,decreasing=TRUE)
    cutoff = 0.5*plinks_sort[nr_links_show*2]+0.5*plinks_sort[nr_links_show*2+1]
  }
  index = which(plinks>cutoff)
  
  #create undirected graph with those x links
  G = matrix(0,100,100)
  G[index] = 1
  G_sum = G_sum + G
  
  #create graph
  if (algorithm == "MPLBD"){G_ig_mplbd = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")}
  if (algorithm == "MPLRJ"){G_ig_mplrj = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")}
  if (algorithm == "BD"){G_ig_bd = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")}
  if (algorithm == "SS"){G_ig_ss = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")}
  if (algorithm == "CONCORD"){G_ig_cc = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")}

}

#position the nodes
index_incl = which(G_sum > 0)
G_all = matrix(0,100,100)
G_all[index_incl] = 1
G_ig_all = graph_from_adjacency_matrix(adjmatrix=G_all,mode="undirected")
V(G_ig_all)$label = 1:100
layout = layout_nicely(G_ig_all,dim=2)

#plot all the graphs
dev.off()
par(mfcol = c(2,3))
par(mar=c(2,1,1,1))
plot.igraph(G_ig_mplbd,mode="undirected",layout=layout,main="MPL-BD",diag=FALSE,edge.width=3,vertex.size=0.1,edge.color = "blue",vertex.label.color ="black",vertex.label.dist=0.5)
plot.igraph(G_ig_mplrj,mode="undirected",layout=layout,main="MPL-RJ",diag=FALSE,edge.width=3,vertex.size=0.1,edge.color = "blue",vertex.label.color ="black",vertex.label.dist=0.5)
plot.igraph(G_ig_bd,mode="undirected",layout=layout,main="BD",diag=FALSE,edge.width=3,vertex.size=0.1,edge.color = "blue",vertex.label.color ="black",vertex.label.dist=0.5)
plot.igraph(G_ig_ss,mode="undirected",layout=layout,main="SS",diag=FALSE,edge.width=3,vertex.size=0.1,edge.color = "blue",vertex.label.color ="black",vertex.label.dist=0.5)
plot.igraph(G_ig_cc,mode="undirected",layout=layout,main="CONCORD",diag=FALSE,edge.width=3,vertex.size=0.1,edge.color = "blue",vertex.label.color ="black",vertex.label.dist=0.5)



########################################
#############difference matrix #########
########################################

algorithm = "MPLBD"
filename = paste0("gene_solution_",algorithm,".Rdata")
load(file=filename)
plinks = output$plinks
p_links_mplbd = plinks[upper.tri(plinks)]

algorithm = "MPLRJ"
filename = paste0("gene_solution_",algorithm,".Rdata")
load(file=filename)
plinks = output$plinks
p_links_mplrj = plinks[upper.tri(plinks)]

algorithm = "BD"
filename = paste0("gene_solution_",algorithm,".Rdata")
load(file=filename)
plinks = output$plinks
p_links_bd = plinks[upper.tri(plinks)]

algorithm = "SS"
filename = paste0("gene_solution_",algorithm,".Rdata")
load(file=filename)
plinks = output$plinks
p_links_ss = plinks[upper.tri(plinks)]

algorithm = "CONCORD"
filename = paste0("gene_solution_",algorithm,".Rdata")
load(file=filename)
plinks = output$plinks
p_links_cc = plinks[upper.tri(plinks)]

#difference via average absolute difference
sum(abs(p_links_mplbd-p_links_mplrj))/4950
sum(abs(p_links_mplbd-p_links_bd))/4950
sum(abs(p_links_mplbd-p_links_ss))/4950
sum(abs(p_links_mplbd-p_links_cc))/4950

sum(abs(p_links_mplrj-p_links_bd))/4950
sum(abs(p_links_mplrj-p_links_ss))/4950
sum(abs(p_links_mplrj-p_links_cc))/4950

sum(abs(p_links_bd-p_links_ss))/4950
sum(abs(p_links_bd-p_links_cc))/4950
sum(abs(p_links_ss-p_links_cc))/4950

########################################
############# overlap tables ###########
########################################

#overlap in inclusion
overlap_func = function(vector1=c(1),vector2=c(1),threshold=0.5){
  index_1 = which(vector1 > threshold)
  index_2 = which(vector2 > threshold)
  overlap = length(intersect(index_1,index_2))
  total = length(index_1)
  return(overlap/total)
}

threshold = 0.9 

length(which(p_links_mplbd > threshold))
length(which(p_links_mplrj > threshold))
length(which(p_links_bd > threshold))
length(which(p_links_ss > threshold))
length(which(p_links_cc > threshold))

vector1 = p_links_mplbd
vector2 = p_links_mplbd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_mplbd
vector2 = p_links_mplrj
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_mplbd
vector2 = p_links_bd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_mplbd
vector2 = p_links_ss
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_mplbd
vector2 = p_links_cc
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)

vector1 = p_links_mplrj
vector2 = p_links_mplbd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_mplrj
vector2 = p_links_mplrj
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_mplrj
vector2 = p_links_bd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_mplrj
vector2 = p_links_ss
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_mplrj
vector2 = p_links_cc
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)

vector1 = p_links_bd
vector2 = p_links_mplbd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_bd
vector2 = p_links_mplrj
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_bd
vector2 = p_links_bd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_bd
vector2 = p_links_ss
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_bd
vector2 = p_links_cc
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)

vector1 = p_links_ss
vector2 = p_links_mplbd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_ss
vector2 = p_links_mplrj
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_ss
vector2 = p_links_bd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_ss
vector2 = p_links_ss
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_ss
vector2 = p_links_cc
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)


vector1 = p_links_cc
vector2 = p_links_mplbd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_cc
vector2 = p_links_mplrj
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_cc
vector2 = p_links_bd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_cc
vector2 = p_links_ss
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_cc
vector2 = p_links_cc
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)


#################################################################################################################################
###########find similarities between our network and networks in literature ####################################
#################################################################################################################################

algorithm = "MPLBD"
filename = paste0("gene_solution_",algorithm,".Rdata")
load(file=filename)
plinks = output$plinks

#fill lower part of matrix
p= 100
count = 0
for (i in 2:p){ #fill lower matrix
  for (j in 1:(i-1)){
    count = count + 1
    plinks[i,j] = plinks[j,i]
  } 
}

#####find similarities with Li et al. 2019 ################
gene_names = colnames(plinks)
degree_cont = rowSums(plinks)
degree_cont_rank = rank(degree_cont)
nr = 20
top_vec = (p-nr+1):100
for (top_nr in top_vec){
  index = which(degree_cont_rank==top_nr)
  print(gene_names[index])
}
  
#####find similarities with the Bhadra and Mallick 2013 ####
#the below tree structure also appear in Bhadra and Mallick 2013
gene_names = colnames(plinks)

gene_names[c(86,3,5,23,57,27)]
gene_names[c(38,13,10,16,18,8,4)]

gene_names[c(90,63,53,47,97,58,26,95,51)]
gene_names[c(73,40)]
gene_names[c(38,13,18,4,8,6,9,10,16)]
gene_names[c(25,93,1,2,12,17)]

####create table containing the conversion of labels 1:100 to the gene names
install.packages("lazyWeave")
library(lazyWeave)

table_r = matrix(0,nrow=100,ncol=2)
table_r[,1] = 1:100
table_r[,2] = gene_names
table_latex = lazy.matrix(table_r,align="left")







  
