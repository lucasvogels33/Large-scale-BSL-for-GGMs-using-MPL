
# set working directory. Make sure that the working directory contains the files:
# final_0.01_BD.Rdata
# final_0.01_CONCORD.Rdata
# final_0.01_MPLBD.Rdata
# final_0.01_MPLRJ.Rdata
# final_0.01_SS.Rdata
# gene_ID_names.Rdata

#setwd(....)

#############################################################
################make edge convergence plots ##############################
#############################################################

#set layout
par(mfrow=c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1) )
algorithm_vec = c("MPLRJ","MPLBD","BD","SS","CONCORD")
col_vec = c(1,2,3,4,5)
lty_vec = c(1,2,3,4,5)
xlim = c(0,10000)
xlab = "time (seconds)"

#make plot
plot(NA,xlim=xlim,ylim=c(0,20000),xlab=xlab,ylab="edges",main="")

#fill plot with lines
count = 0
g.prior = 0.01
for (algorithm in algorithm_vec){

  #load data
  filename = paste0("final_",g.prior,"_",algorithm,".Rdata")
  load(file=filename)
  edge_vec=output$edge_vec
  time_vec = output$time_vec
  
  #add the point at time t=0 to the vectors
  if (algorithm=="CONCORD"){ #the full graph for the CONCORD algorithm
    p = dim(output$plinks)[1]
    total_links = p*(p-1)/2
    edge_vec=c(total_links,edge_vec)
    time_vec=c(0,time_vec)
  }
  else { #the empty graph for the other algorithms
    edge_vec=c(0,edge_vec)
    time_vec=c(0,time_vec)
  }
  
  #set colour and line tye
  count = count + 1
  col = col_vec[count]
  lty = lty_vec[count]
  
  #add line
  points(x=time_vec,y=edge_vec,col=col,lty=lty,lw = 2,type="l")
}
#add legend
legend(x=8000,y=14000, algorithm_vec, col = col_vec,lw = c(2,2,2,2,2),lty = lty_vec, cex = 1, box.lty = 0 )

##################################
#############gene networks########
##################################

#download necessary library
library(igraph)

#set parameters 
algorithm_vec = c("MPLBD","MPLRJ","SS","CONCORD")
cutoff = 0.9 
count = 0
G_sum = matrix(0,623,623)
g.prior = 0.01

#create graphs for all algorithms 
for (algorithm in algorithm_vec){
  filename = paste0("final_",g.prior,"_",algorithm,".Rdata")
  load(file=filename)
  plinks = output$plinks
  
  #fill the lower triangle of the matrix
  p = 623
  for (i in 2:p){ #fill lower matrix
    for (j in 1:(i-1)){
      plinks[i,j] = plinks[j,i]
    } 
  }
  
  #determine indices of included links
  index = which(plinks >= cutoff)
  
  #create undirected graph with included links
  G = matrix(0,623,623)
  G[index] = 1
  G_sum = G_sum + G
  
  #create graph
  if (algorithm == "MPLBD"){G_ig_mplbd = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")}
  if (algorithm == "MPLRJ"){G_ig_mplrj = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")}
  if (algorithm == "SS"){G_ig_ss = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")}
  if (algorithm == "CONCORD"){G_ig_cc = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")}
  
}

#make a layout for the positioning of the nodes
index_incl = which(G_sum > 0)
G_all = matrix(0,623,623)
G_all[index_incl] = 1
G_ig_all = graph_from_adjacency_matrix(adjmatrix=G_all,mode="undirected")
V(G_ig_all)$label = 1:623
layout = layout_nicely(G_ig_all)

#plot all the graphs
dev.off()
par(mfcol = c(2 ,2))
par(mar=c(2,1,1,1))
plot.igraph(G_ig_mplbd,mode="undirected",layout=layout,main="MPL-BD",diag=FALSE,edge.width=0.1,vertex.size=0.1,edge.color = "blue",vertex.label = NA,vertex.label.dist=0.5)
plot.igraph(G_ig_mplrj,mode="undirected",layout=layout,main="MPL-RJ",diag=FALSE,edge.width=0.1,vertex.size=0.1,edge.color = "blue",vertex.label = NA,vertex.label.dist=0.5)
plot.igraph(G_ig_ss,mode="undirected",layout=layout,main="SS",diag=FALSE,edge.width=0.1,vertex.size=0.1,edge.color = "blue",vertex.label = NA,vertex.label.dist=0.5)
plot.igraph(G_ig_cc,mode="undirected",layout=layout,main="CONCORD",diag=FALSE,edge.width=0.1,vertex.size=0.1,edge.color = "blue",vertex.label = NA,vertex.label.dist=0.5)


####################################################
#############difference and overlap matrices #########
####################################################
g.prior=0.01

algorithm = "MPLBD"
filename = paste0("final_",g.prior,"_",algorithm,".Rdata")
load(file=filename)
plinks = output$plinks
plinks_matrix_mplbd = plinks
p_links_mplbd = plinks[upper.tri(plinks)]

algorithm = "MPLRJ"
filename = paste0("final_",g.prior,"_",algorithm,".Rdata")
load(file=filename)
plinks = output$plinks
plinks_matrix_mplrj = plinks
p_links_mplrj = plinks[upper.tri(plinks)]

algorithm = "SS"
filename = paste0("final_",g.prior,"_",algorithm,".Rdata")
load(file=filename)
plinks = output$plinks
plinks_matrix_ss = plinks
p_links_ss = plinks[upper.tri(plinks)]

algorithm = "CONCORD"
filename = paste0("final_",g.prior,"_",algorithm,".Rdata")
load(file=filename)
plinks = output$plinks
plinks_matrix_cc = plinks
p_links_cc = plinks[upper.tri(plinks)]

#average absolute difference between algorithms
total_links = length(p_links_mplbd)
sum(abs(p_links_mplbd-p_links_mplrj))/total_links
sum(abs(p_links_mplbd-p_links_ss))/total_links
sum(abs(p_links_mplbd-p_links_cc))/total_links
sum(abs(p_links_mplrj-p_links_ss))/total_links
sum(abs(p_links_mplrj-p_links_cc))/total_links
sum(abs(p_links_ss-p_links_cc))/total_links

#overlap between algorithms
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
length(which(p_links_ss > threshold))
length(which(p_links_cc > threshold))

vector1 = p_links_mplbd
vector2 = p_links_mplbd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_mplbd
vector2 = p_links_mplrj
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
vector2 = p_links_ss
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_mplrj
vector2 = p_links_cc
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)

vector1 = p_links_ss
vector2 = p_links_mplbd
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_ss
vector2 = p_links_mplrj
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
vector2 = p_links_ss
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)
vector1 = p_links_cc
vector2 = p_links_cc
overlap_func(vector1=vector1,vector2=vector2,threshold = threshold)


###################################################
#########literature comparison ####################
###################################################

#load plinks matrix of an algorithm 
g.prior = 0.01 
algorithm = "MPLRJ" #choose either MPLRJ, MPLBD, SS or CONCORD
filename = paste0("final_",g.prior,"_",algorithm,".Rdata")
load(file=filename)
plinks = output$plinks
p_links_vec = plinks[upper.tri(plinks)]

#load file containing gene IDs and names
filename = paste0("gene_ID_names.Rdata")
load(file=filename)
p = dim(plinks)[1]
gene_ID_names = gene_ID_names[1:p,]

#set margins and layout settings
library(igraph)
dev.off()
par(mfrow=c(2,2))
par(mar=c(2,1,2,1)) #bottom, left, top and right margins
edge_width = 2

#create network for histone genes
first_four = substr(gene_ID_names$gene_name,start=1,stop=4)
index_histone = which(first_four=="Hist")
plinks_histone = plinks[index_histone,index_histone]
plinks_histone = plinks_histone + t(plinks_histone)
index_include = which(plinks_histone > 0.5)
nr_nodes = length(index_histone)
G = matrix(0,nr_nodes,nr_nodes)
G[index_include] = 1
G_histone = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")
V(G_histone)$label = gene_ID_names$gene_name[index_histone]
layout = layout_nicely(G_histone)
plot.igraph(G_histone,mode="undirected",layout=layout,main="Histone genes",diag=FALSE,edge.width=edge_width,vertex.size=10,edge.color = "black",vertex.label.dist=2)

#create network for ly6c genes
index_ly6c1 = which(gene_ID_names$gene_name=="Ly6c1")
index_ly6c2 = which(gene_ID_names$gene_name=="Ly6c2")
index_ly = c(index_ly6c1,index_ly6c2)
plinks_ly = plinks[index_ly,index_ly]
plinks_ly = plinks_ly + t(plinks_ly)
index_include = which(plinks_ly > 0.5)
nr_nodes = length(index_ly)
G = matrix(0,nr_nodes,nr_nodes)
G[index_include] = 1
G_ly = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")
V(G_ly)$label = c("Ly6c1","Ly6c2")
layout = layout_nicely(G_ly)
plot.igraph(G_ly,mode="undirected",layout=layout,main="Leukocyte antigen genes",diag=FALSE,edge.width=edge_width,vertex.size=10,edge.color = "black",vertex.label.dist=4)

#create network for Bc12a genes
index_Bcl2a1a = which(gene_ID_names$gene_name=="Bcl2a1a")
index_Bcl2a1b = which(gene_ID_names$gene_name=="Bcl2a1b")
index_Bcl2a1c = which(gene_ID_names$gene_name=="Bcl2a1c")
index_Bcl2a1d = which(gene_ID_names$gene_name=="Bcl2a1d")
index_Bc12a = c(index_Bcl2a1a,index_Bcl2a1b,index_Bcl2a1c,index_Bcl2a1d)
plinks_Bc12a = plinks[index_Bc12a,index_Bc12a]
plinks_Bc12a = plinks_Bc12a + t(plinks_Bc12a)
index_include = which(plinks_Bc12a > 0.5)
nr_nodes = length(index_Bc12a)
G = matrix(0,nr_nodes,nr_nodes)
G[index_include] = 1
G_Bc12a = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")
V(G_Bc12a)$label = c("Bcl2a1a","Bcl2a1b","Bcl2a1c","Bcl2a1d")
layout = layout_nicely(G_Bc12a)
plot.igraph(G_Bc12a,mode="undirected",layout=layout,main="B-cell leukemia genes",diag=FALSE,edge.width=edge_width,vertex.size=10,edge.color = "black",vertex.label.dist=2)

#create network for Ms4a genes
first_four = substr(gene_ID_names$gene_name,start=1,stop=4)
index_membrane4A = which(first_four=="Ms4a")
plinks_membrane4A = plinks[index_membrane4A,index_membrane4A]
plinks_membrane4A = plinks_membrane4A + t(plinks_membrane4A)
index_include = which(plinks_membrane4A > 0.5)
nr_nodes = length(index_membrane4A)
G = matrix(0,nr_nodes,nr_nodes)
G[index_include] = 1
G_membrane4A = graph_from_adjacency_matrix(adjmatrix=G,mode="undirected")
V(G_membrane4A)$label = gene_ID_names$gene_name[index_membrane4A]
layout = layout_nicely(G_membrane4A)
plot.igraph(G_membrane4A,mode="undirected",layout=layout,main="membrane spanning 4A genes",diag=FALSE,edge.width=edge_width,vertex.size=10,edge.color = "black",vertex.label.dist=2)

