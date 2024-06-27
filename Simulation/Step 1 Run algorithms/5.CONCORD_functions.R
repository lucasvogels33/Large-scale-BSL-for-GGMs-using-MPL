
############################
### CONCORD function#######
############################

#the function CONCORD_solve takes as input:
#- data from BDgraph. This "data" type contains data, a true graph 
#- iter_vec_thin. The iteration at which it needs to produce output

#It outputs:
#- the AUC ROC at every iteration number indicated in iter_vec_thin
#- the p+ at every iteration number indicated in iter_vec_thin
#- the p- at every iteration number indicated in iter_vec_thin
#- the AUC PR at every iteration number indicated in iter_vec_thin
#- the F1 score at every iteration number indicated in iter_vec_thin (at a threshold of 0.5)
#- the running time at every iteration number indicated in iter_vec_thin

CONCORD_solve = function(data,iter_vec_thin){
  
  #obtain true graph
  G_true = data$G
  response = G_true[upper.tri(G_true)]
  
  #obtain amount of observations and variables
  n = nrow(data$data)
  p = ncol(data$data)
  density = data$density
  graph = data$graph
  rep = data$rep
  algorithm_name = "CONCORD"
  S = crossprod(data$data)/n
  
  #parameters
  burnin = 0
  plinks_old = 0
  omega_start = matrix(rnorm(p*p,mean=0,sd=1), p, p)
 
  #run the algorithm, produce output at every iteration in the iter_vec_thin vector
  len = length(iter_vec_thin)
  AUC_vec = rep(0,len)
  pplus_vec = rep(0,len)
  pmin_vec = rep(0,len)
  AUC_PR_vec = rep(0,len)
  F1_vec = rep(0,len)
  time_vec = rep(0,len)
  
  for (j in 1:len){
    
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
    predictor = plinks_new[upper.tri(plinks_new)]
    
    #calculate the AUC, pplus, pmin, AUC_PR, F1
    AUC = calc_AUC_ROC(predictor=predictor,response=response)
    obj = calc_pplus_pmin(predictor=predictor,response=response)
    p_plus = obj$p_plus
    p_min = obj$p_min
    AUC_PR = calc_AUC_PR(predictor=predictor,response=response)
    F1 = calc_F1(predictor=predictor,response=response)
    
    #save metrics
    AUC_vec[j] = AUC
    pplus_vec[j] = p_plus
    pmin_vec[j] = p_min
    time_vec[j] = runtime
    AUC_PR_vec[j] = AUC_PR
    F1_vec[j] = F1
        
    #store data for next run
    plinks_old = plinks_new
    omega_start = sample_concord$last_omega
    
    #save output per iteration in Rdata file (only for large scale problems)
    if (p > 250){
      output_per_iter = list()
      output_per_iter$AUC = AUC
      output_per_iter$pplus = p_plus
      output_per_iter$pmin = p_min
      output_per_iter$AUC_PR = AUC_PR
      output_per_iter$F1 = F1
      output_per_iter$time = sum(time_vec)
      output_per_iter$iter = newiter
      filename = paste0("results_p",p,"_n",n,"_",graph,"_",density,"_",algorithm_name,"_rep",rep,"_iter",j,".Rdata")
      save(output_per_iter,file = filename)
    }
  }
  time_vec = cumsum(time_vec)
  time_vec = time_vec
  
  return(list(AUC_vec=AUC_vec,
              pplus_vec=pplus_vec,
              pmin_vec=pmin_vec,
              AUC_PR_vec=AUC_PR_vec,
              F1_vec=F1_vec,
              time_vec=time_vec,
              iter_vec_thin=iter_vec_thin,
              predictor=predictor))
} 




