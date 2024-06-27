#################################
#####supporting functions ###############
#################################

#- we create 5 functions (one for each algorithm)

#Input:
#- the normalized data 
#- the amount of burnin iterations, specified as a vector at which to output the running time and #edges in the Markov state
#- the amount of MCMC iterations after burnin, specified as a vector at which to output the running time and #edges in the Markov state
#- the initital graph
#- the prior on the graph

#Output:
#- burnin_iter_vec_thin: the amount of burnin iterations, specified as a vector at which to output the running time and #edges in the Markov state
#- iter_vec_thin: the amount of MCMC iterations after burnin, specified as a vector at which to output the running time and #edges in the Markov state
#- edge_vec: contains #edges at every iteration mentioned in burnin_iter_vec_thin and iter_vec_thin
#- time_vec: contains runtime at every iteration mentioned in burnin_iter_vec_thin and iter_vec_thin
#- plinks: matrix containing the edge inclusion probabilities
 
MPLBD_solve = function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart,g.prior){
  
  #parameters
  burnin = 0
  jump = 1
  save = FALSE
  cores = 1
  g.prior = g.prior
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
MPLRJ_solve = function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart,g.prior){
  
  #parameters
  burnin = 0
  jump = 1
  save = FALSE
  cores = 1
  g.prior = g.prior
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
SS_solve= function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart,g.prior){
  
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
  g.prior = g.prior
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
CONCORD_solve = function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart,g.prior){
  
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
BDA_solve = function(data,burnin_iter_vec_thin,iter_vec_thin,MCMCstart,g.prior){
  
  #parameters
  burnin = 0
  jump = 1
  save = FALSE
  cores = 1
  verbose = TRUE
  
  #initialize output vectors
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
    
    #save starting point next run
    MCMCstart = sample_bd
        
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

