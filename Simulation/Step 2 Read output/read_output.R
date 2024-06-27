#set working directory. Make sure that the directory contains all the output files created in step 1

#setwd()

#########################
#####true graphs #####
#########################

p = 1000
n = 400
density = "dense"
graph = "random"
mean_degree_vec=c()
max_degree_vec=c()
sparsity_vec = c()
total_edge_vec = c()

for (rep in 1:8){
  filename = paste0("truth_p",p,"_n",n,"_",graph,"_",density,"_rep",rep,".Rdata")
  load(file=filename)
  G = data_output$G
  
  mean_degree_vec = c(mean_degree_vec,mean(colSums(G)))
  max_degree_vec = c(max_degree_vec,max(colSums(G)))
  sparsity_vec = c(sparsity_vec,sum(G)/2/(p*(p-1)/2))
  total_edge_vec = c(total_edge_vec,sum(G)/2)
  
}

mean(max_degree_vec)
mean(sparsity_vec)

##################################################
#####convergence of AUC, Pr+, Pr-, AUC_PR and F1#####
##################################################
p = 1000
n = 1050
graph = "scale-free"
density = "sparse"
algorithm_vec = c("MPLRJ","MPLBD","BDA","SSO","CONCORD")
yaxis = "AUC_PR" #possible values are AUC, p_plus, or p_min, AUC_PR, F1
col_vec = c(1,2,3,4,5)
lty_vec = c(1,2,3,4,5)
unit = "min" #select either "sec", "min" or "hour"

#calculate sparsity for AUC_PR initial setting
rep = 1
filename = paste0("truth_p",p,"_n",n,"_",graph,"_",density,"_rep",rep,".Rdata")
load(file=filename)
G = data_output$G
sparsity = (sum(G)/2)/(p*(p-1)/2)

title = paste0("p =",p,", n = ",n,", ",graph," graph with ",density,"density")
title = ""
if (yaxis=="AUC"){ylim = c(0.5,1)}
if (yaxis=="p_plus"){ylim = c(0,1)}
if (yaxis=="p_min"){ylim = c(0,0.2)}
if (yaxis=="AUC_PR"){ylim = c(0,1)}
if (yaxis=="F1"){ylim = c(0,1)}

if (yaxis=="AUC"){ylab = "AUC-ROC"}
if (yaxis=="p_plus"){ylab = "Pr+"}
if (yaxis=="p_min"){ylab = "Pr-"}
if (yaxis=="AUC_PR"){ylab = "AUC-PR"}
if (yaxis=="F1"){ylab = "F1"}

if (unit=="sec"){xlab="Time (in seconds)"}
if (unit=="min"){xlab="Time (in minutes)"}
if (unit=="hour"){xlab="Time (in hours)"}

plot(NA,xlim=c(0,15000),ylim=ylim,xlab=xlab,ylab=ylab,main = title)

count = 0
count_lty = 0
for (method in algorithm_vec ){
  
  ##################
  ####load data#####
  ##################
  #initialize metrics
  AUC = 0
  p_plus = 0
  p_min = 0
  time = 0
  AUC_matrix = c()
  p_plus_matrix = c()
  p_min_matrix = c()
  AUC_PR_matrix = c()
  F1_matrix = c()
  time_matrix = c()
  
  rep_list = 1:16
  if (p==1000){rep_list=1:8}
  if ((p==1000)*(method=="BDA")*(graph=="cluster")*(density=="dense")){rep_list=c(1,2,4,5,6,7)}
  replications = length(rep_list)
  
  
  for (rep in rep_list){
    
    #load file
    filename = paste0("results_p",p,"_n",n,"_",graph,"_",density,"_",method,"_rep",rep,".Rdata")
    load(file=filename)
    
    #load metrics
    time_matrix = rbind(time_matrix,output$time_vec)
    AUC_matrix = rbind(AUC_matrix,output$AUC_vec)
    p_plus_matrix = rbind(p_plus_matrix,output$pplus_vec)
    p_min_matrix = rbind(p_min_matrix,output$pmin_vec)
    AUC_PR_matrix = rbind(AUC_PR_matrix,output$AUC_PR_vec)
    F1_matrix = rbind(F1_matrix,output$F1_vec)
    
  }
  
  
  #######################
  ####plot averages #####
  #######################
  
  #select colour
  count = count +1
  col = col_vec[count]
  
  #select line type
  count_lty = count_lty + 1
  lty = lty_vec[count_lty]
  
  #average each measurement
  time_vec = c(0,colMeans(time_matrix)) #we add the point time t=0)
  AUC_vec = c(0.5,colMeans(AUC_matrix)) #we add the point AUC=0.5 at time t=0)
  AUC_PR_vec = c(sparsity,colMeans(AUC_PR_matrix)) #we add the point AUC_PR = sparsity% at time t=0)
  F1_vec = c(0,colMeans(F1_matrix)) #we add the point F1 = 0 at time t=0)
  p_plus_vec = c(0,colMeans(p_plus_matrix)) #we add the point p+ = 0 at time t=0)
  p_min_vec = c(0,colMeans(p_min_matrix)) #we add the point p- = 0 at time t=0)
  
  if (method=="CONCORD"){ #concord starts at the full graph
    p_plus_vec = c(1,colMeans(p_plus_matrix)) #we add the point p+ = 1 at time t=0)
    p_min_vec = c(1,colMeans(p_min_matrix)) #we add the point p- = 1 at time t=0)
  }
  
  #converge time to correct measurement
  if (unit=="min"){time_vec=time_vec/60}
  if (unit=="hour"){time_vec=time_vec/3600}
  
  if (yaxis == "AUC"){points(x=time_vec,y=AUC_vec,type="l",col=col,lw=2,lty=lty) }
  if (yaxis == "p_plus"){points(x=time_vec,y=p_plus_vec,type="l",col=col,lw=2,lty=lty) }
  if (yaxis == "p_min"){points(x=time_vec,y=p_min_vec,type="l",col=col,lw=2,lty=lty) }
  if (yaxis == "AUC_PR"){points(x=time_vec,y=AUC_PR_vec,type="l",col=col,lw=2,lty=lty) }
  if (yaxis == "F1"){points(x=time_vec,y=F1_vec,type="l",col=col,lw=2,lty=lty) }

}

if (p < 250){
  algorithm_vec = algorithm_vec = c("MPLRJ","MPLBD","BDA","SSO")
  legend( "topright", inset=.03, algorithm_vec, col = 1:4,lw = c(2,2,2,2),lty = c(1,2,3,4), cex = 1, box.lty = 0 )
}
if (p > 250){
  algorithm_vec = algorithm_vec = c("MPL-RJ","MPL-BD","BD","SS","B-CONCORD")
  legend(x=75000,y=0.45, inset=.03, algorithm_vec, col = c(1,2,3,4,5),lw = c(2,2,2,2,2),lty = c(1,2,3,4,5), cex = 1, box.lty = 0 )
}


##################################################
#####convergence of metric per method ###############
##################################################
par(mfrow=c(1,1))
p = 1000
n = 400
graph = "cluster"
density = "dense"
method_list = c("SSO") #c("MPLRJ","MPLBD","SSO","CONCORD")
ylab = "AUC_PR"  #select either AUC or AUC_PR or p_plus


title = paste0("p =",p,", n = ",n,", ",graph," graph,",density)
if (ylab=="AUC"){ylim = c(0.5,1)}
if (ylab=="AUC_PR"){ylim = c(0,1)}
if (ylab=="p_plus"){ylim = c(0,1)}
plot(NA,xlim=c(0,250000),ylim=ylim,xlab="Time (in seconds)",ylab=ylab,main = title)

rep_list = 1:8
col_vec = c(1,2,4,5)
count = 0
for (method in method_list){
  count = count + 1
  for (rep in rep_list){
  
    ##################
    ####load data#####
    ##################
    filename = paste0("results_p",p,"_n",n,"_",graph,"_",density,"_",method,"_rep",rep,".Rdata")
    load(file=filename)
    time_vec = output$time_vec
    if (ylab=="AUC"){output_vec = output$AUC_vec}
    if (ylab=="AUC_PR"){output_vec = output$AUC_PR_vec}
    if (ylab=="p_plus"){output_vec = output$pplus_vec}
    
    #######################
    ####plot lines #####
    #######################
    lty = 1
    col = col_vec[count]
    points(x=time_vec,y=output_vec,type="l",col=col,lw=2,lty=lty)
  }
}

if (p < 250){
  algorithm_vec = algorithm_vec = c("MPLRJ","MPLBD","BDA","SSO")
  legend( "topright", inset=.03, algorithm_vec, col = 1:4,lw = c(2,2,2,2),lty = c(1,2,3,4), cex = 1, box.lty = 0 )
}
if (p > 250){
  algorithm_vec = algorithm_vec = c("MPLRJ","MPLBD","SSO","CONCORD")
  legend( "bottomright", inset=.03, algorithm_vec, col = c(1,2,4,5),lw = c(2,2,2,2),lty = c(1,2,4,5), cex = 1, box.lty = 0 )
}




#########################
#####latex tables #####
#########################

p_list = c(10,100,1000)
graph_list = c("random","cluster","scale-free")
method_list = c("MPLRJ","MPLBD","BDA","SSO","CONCORD")
sd_on = FALSE #set to true if standard deviations need to be printed

final_values = function(p=10,n=20,graph="cluster",density="sparse",rep_list=1:16,method="RJWWA"){
  
  #initiate vectors
  AUC_vec = c()
  pplus_vec = c()
  pmin_vec = c()
  AUC_PR_vec = c()
  F1_vec = c()
  
  #load data into the vectors
  for (rep in rep_list){
    
    #load file
    filename = paste0("results_p",p,"_n",n,"_",graph,"_",density,"_",method,"_rep",rep,".Rdata")
    load(file = filename)
    
    #save metrics
    AUC_vec = c(AUC_vec,tail(output$AUC_vec, n=1))
    pplus_vec = c(pplus_vec,tail(output$pplus_vec, n=1))
    pmin_vec = c(pmin_vec,tail(output$pmin_vec, n=1))
    AUC_PR_vec = c(AUC_PR_vec,tail(output$AUC_PR_vec, n=1))
    F1_vec = c(F1_vec,tail(output$F1_vec, n=1))  
    
    
  }
  
  output_list = list()
  output_list$AUC_vec = AUC_vec
  output_list$pplus_vec = pplus_vec
  output_list$pmin_vec = pmin_vec
  output_list$AUC_PR_vec = AUC_PR_vec
  output_list$F1_vec = F1_vec
  
  return(output_list)
}

make_latex_table = function(graph_list=c("cluster"),p_list = c(10),method_list=c("mpl_bd"),
                            file_p_plus="pplus.txt",file_p_min="pmin.txt",
                            file_AUC = "AUC.txt",file_AUC_PR = "AUC_PR.txt",file_F1 = "F1.txt", round_nr=2,sd_on=TRUE){
  
  method_string = paste(method_list,collapse=" & ")
  
  for (paste_filename in c(file_p_plus,file_p_min,file_AUC,file_AUC_PR,file_F1)){
    cat("p & graph & density & n & ",method_string," \\\\",file=paste_filename,"\n", append = TRUE )
    cat("\\hline",file=paste_filename,"\n", append = TRUE )
  }
  
  for (p in p_list){
    print(p)
    for (graph in graph_list){
      print(graph)
      if (graph=="random"){density_list=c("sparse","dense")}
      if (graph=="cluster"){density_list=c("sparse","dense")}
      if (graph=="scale-free"){density_list=c("sparse")}
      for (density in density_list){
      
        print(density)
      
        if (p==10){n_list=c(20,350)}
        if (p==100){n_list=c(40,700)}
        if (p==1000){n_list=c(400,1050)}
        
        for (n in n_list){
          for (paste_filename in c(file_p_plus,file_p_min,file_AUC,file_AUC_PR,file_F1)){
            cat(p,"&",graph," & ",density," & ",n," &",  file=paste_filename, append = TRUE )
          }
          print(n)
          for (method in method_list){
              print(method)
        
              rep_list = 1:16
              if (p==1000){rep_list=1:8}
              if ((p==1000)*(method=="BDA")*(graph=="cluster")*(density=="dense")){rep_list=c(1,2,4,5,6,7)}
  
              output_list = final_values(p=p,n=n,graph=graph,density=density,rep_list=rep_list,method=method)
              sd_text = ""
              
              #p_plus
              mean = (round(mean(output_list$pplus_vec),round_nr))
              if (sd_on){
                sd_val = (round(sd(output_list$pplus_vec),round_nr))
                sd_text = paste0("(",sd_val,")")
              }
              text = paste0(mean," ",sd_text," & ")
              cat(text,file = file_p_plus, append = TRUE )
              
              #p_min
              mean = (round(mean(output_list$pmin_vec),round_nr))
              if (sd_on){
                sd_val = (round(sd(output_list$pmin_vec),round_nr))
                sd_text = paste0("(",sd_val,")")
              }
              text = paste0(mean," ",sd_text," & ")
              cat(text,file = file_p_min, append = TRUE )
              
              #AUC
              mean = (round(mean(output_list$AUC_vec),round_nr))
              if (sd_on){
                sd_val = (round(sd(output_list$AUC_vec),round_nr))
                sd_text = paste0("(",sd_val,")")
              }
              text = paste0(mean," ",sd_text," & ")
              cat(text,file = file_AUC, append = TRUE )
              
              #AUC_PR
              mean = (round(mean(output_list$AUC_PR_vec),round_nr))
              if (sd_on){
                sd_val = (round(sd(output_list$AUC_PR_vec),round_nr))
                sd_text = paste0("(",sd_val,")")
              }
              text = paste0(mean," ",sd_text," & ")
              cat(text,file = file_AUC_PR, append = TRUE )
              
              #F1
              mean = (round(mean(output_list$F1_vec),round_nr))
              if (sd_on){
                sd_val = (round(sd(output_list$F1_vec),round_nr))
                sd_text = paste0("(",sd_val,")")
              }
              text = paste0(mean," ",sd_text," & ")
              cat(text,file = file_F1, append = TRUE )
          
          } #closes method
          
          for (paste_filename in c(file_p_plus,file_p_min,file_AUC,file_AUC_PR,file_F1)){
            cat("\\\\",file = paste_filename,"\n", append = TRUE )
          }
          
        }#closes n
      }  #closes density
    }#closes graph
}#closes p
  
}#closes function

make_latex_table(graph_list=graph_list,p_list = p_list,method_list=method_list,sd_on=sd_on)

########################################
############## computational cost ############
########################################

#supporting functions
convergence_values = function(graph="random",p = 10,method="GMPLBD",n=100,metric="AUC",density=0.1,rep_list=1:16,unit="sec"){
  
  cost_vec = c()
  for (rep in rep_list){
    
    #load AUC_vec and time_vec
    filename = paste0("results_p",p,"_n",n,"_",graph,"_",density,"_",method,"_rep",rep,".Rdata")
    load(filename)
    
    if (metric=="AUC"){output_vec = output$AUC_vec}
    if (metric=="AUC_PR"){output_vec = output$AUC_PR_vec}
    
    time_vec = output$time_vec
    
    #calculate computational cost
    output_final = tail(output_vec, n=1)
    diff_output_vec = output_final - output_vec
    epsilon = 0.01
    conv_index = min(which(diff_output_vec < epsilon))
    conv_time = time_vec[conv_index]
    
    #save computaional cost
    cost_vec = c(cost_vec,conv_time)
    
  }
  
  if (unit=="min"){
    cost_vec = cost_vec/60
  }
  if (unit=="hour"){
    cost_vec = cost_vec/3600
  }
  return(cost_vec)
}

make_latex_table_comp_cost = function(graph_list=c("cluster"),metric="AUC",p_list = c(10),method_list=c("mpl_bd"),
                            paste_filename="cost.txt",round_nr=2,sd_on=TRUE,unit="sec"){
  
  method_string = paste(method_list,collapse=" & ")
  cat("p & graph & density & n & ",method_string," \\\\",file=paste_filename,"\n", append = TRUE )
  cat("\\hline",file=paste_filename,"\n", append = TRUE )
  
  for (p in p_list){
    print(p)
    for (graph in graph_list){
      print(graph)
      if (graph=="random"){density_list=c("sparse","dense")}
      if (graph=="cluster"){density_list=c("sparse","dense")}
      if (graph=="scale-free"){density_list=c("sparse")}
      for (density in density_list){
        
        print(density)
        
        if (p==10){n_list=c(20,350)}
        if (p==100){n_list=c(40,700)}
        if (p==1000){n_list=c(400,1050)}
        
        for (n in n_list){
          print(n)
          cat(p,"&",graph," & ",density," & ",n," &",  file=paste_filename, append = TRUE )
          
          for (method in method_list){
            print(method)
            if ((p==1000)*(method=="BDA")){
              cat("-&",file = paste_filename, append = TRUE )
            }
            else {
              rep_list = 1:16
              if(p==1000){rep_list=1:8}
              cost_vec = convergence_values(graph=graph,p = p,method=method,n=n,metric=metric,density=density,rep_list=rep_list,unit=unit)
              
              mean = round(mean(cost_vec),round_nr)
              sd_text = ""
              if (sd_on){
                sd_val = round(sd(cost_vec),round_nr)
                sd_text = paste0("(",sd_val,")")
              }
              text = paste0(mean," ",sd_text," & ")
              cat(text,file = paste_filename, append = TRUE )
            }
            
          } #closes method
          
          cat("\\\\",file = paste_filename,"\n", append = TRUE )
          
        }#closes n
      }  #closes density
    }#closes graph
  }#closes p
  
}#closes function

#creating the table
graph_list = c("random","cluster","scale-free")
p_list = c(10,100,1000)
method_list = c("MPLRJ","MPLBD","BDA","SSO","CONCORD")
metric = "AUC_PR" #options are AUC or AUC_PR
paste_filename = "cost_AUC_PR.txt"
round_nr = 0
sd_on = FALSE
unit = "min" #select either "sec" or "min" or "hour"
make_latex_table_comp_cost(graph_list=graph_list,metric=metric,p_list = p_list,method_list=method_list,paste_filename=paste_filename,round_nr =round_nr,sd_on=sd_on,unit=unit)

##############################################
############## MCMC iter and time ############
##############################################
make_iter_time_table= function(graph_list=c("cluster"),p_list = c(10),method_list=c("mpl_bd"),
                                      iter_filename="iter.txt",time_filename="time.txt",round_nr=2){
  
  method_string = paste(method_list,collapse=" & ")
  
  for (paste_filename in c(iter_filename,time_filename)){
    cat("p & graph & density & n & ",method_string," \\\\",file=paste_filename,"\n", append = TRUE )
    cat("\\hline",file=paste_filename,"\n", append = TRUE )
  }
  
  for (p in p_list){
    print(p)
    for (graph in graph_list){
      print(graph)
      if (graph=="random"){density_list=c("sparse","dense")}
      if (graph=="cluster"){density_list=c("sparse","dense")}
      if (graph=="scale-free"){density_list=c("sparse")}
      for (density in density_list){
        
        print(density)
        
        if (p==10){n_list=c(20,350)}
        if (p==100){n_list=c(40,700)}
        if (p==1000){n_list=c(400,1050)}
        
        for (n in n_list){
          print(n)
          for (paste_filename in c(iter_filename,time_filename)){
            cat(p,"&",graph," & ",density," & ",n," &",  file=paste_filename, append = TRUE )
          }
          
          for (method in method_list){
            print(method)
                          
            #amount of MCMC iterations
            rep = 1
            filename = paste0("results_p",p,"_n",n,"_",graph,"_",density,"_",method,"_rep",rep,".Rdata")
            load(filename)
            MCMCiter = tail(output$iter_vec_thin,1)
            text = paste0(MCMCiter," & ")
            cat(text,file = iter_filename, append = TRUE )
            
            #total running time
            rep_list = 1:16
            if (p==1000){rep_list=1:8}
            if ((p==1000)*(method=="BDA")*(graph=="cluster")*(density=="dense")){rep_list=c(1,2,4,5,6,7)}
            time_vec = rep(0,length(rep_list))
            for (rep in rep_list){
              filename = paste0("results_p",p,"_n",n,"_",graph,"_",density,"_",method,"_rep",rep,".Rdata")
              load(filename)
              time_vec[rep] = tail(output$time_vec,1)
            }
            mean = round(mean(time_vec),round_nr)
            sd = round(sd(time_vec),round_nr)
            text = paste0(mean," (",sd,") & ")
            cat(text,file = time_filename, append = TRUE )
            
            
          
            
          } #closes method
          for (paste_filename in c(iter_filename,time_filename)){
            cat("\\\\",file = paste_filename,"\n", append = TRUE )
          }
          
        }#closes n
      }  #closes density
    }#closes graph
  }#closes p
  
}#closes function

graph_list = c("random","cluster","scale-free")
p_list = c(10,100,1000)
method_list = c("MPLRJ","MPLBD","BDA","SSO","CONCORD")
options(scipen=999) #ensures correct notation of large numbers
make_iter_time_table(graph_list=graph_list,p_list = p_list,method_list=method_list,
                                iter_filename="iter.txt",time_filename="time.txt",round_nr=2)






