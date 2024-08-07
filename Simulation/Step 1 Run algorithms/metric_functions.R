####################################################
############install packages and functions##########
###################################################

calc_AUC_PR = function(predictor,response){
  if (length(response) != length(predictor)) {
    stop("response and predictor vector must be of same length")
  }
  else {
    index_ones = which(response==1)
    index_zeroes = which(response==0)
    predictor_ones = predictor[index_ones]
    predictor_zeroes = predictor[index_zeroes]
    AUC_PR_obj = pr.curve(scores.class0 = predictor_ones, scores.class1 = predictor_zeroes)
    AUC_PR = AUC_PR_obj$auc.davis.goadrich
  }
  return(AUC_PR)
}

calc_F1=function(predictor,response,threshold=0.5){
  
  if (length(response) != length(predictor)) {
    stop("response and predictor vector must be of same length")
  }
  
  #calculate number of true instances
  ones = sum(response)
  
  #make binary predictor
  index_ones = which(predictor>=threshold)
  predictor_binary =  rep(0,length(predictor))
  predictor_binary[index_ones] = 1
  
  #calculate TP, FP and FN
  TP = sum(predictor_binary[which(response==1)])
  FP = sum(predictor_binary[which(response==0)])
  FN = ones - TP
  
  #calculate F1
  F1 = (2*TP)/(2*TP + FP + FN)
  return(F1) 
}

calc_pplus_pmin = function(predictor,response){
  if (length(response) != length(predictor)) {
    stop("response and predictor vector must be of same length")
  }
  
  #calculate the weighed CE and weighed MSE
  ones_index = which(response==1)
  zeroes_index = which(response==0)
  
  predictor_ones =  predictor[ones_index]
  predictor_zeroes = predictor[zeroes_index]
  
  p_plus = mean(predictor_ones)
  p_min = mean(predictor_zeroes)
  
  return(list(p_plus=p_plus,p_min=p_min))
}

calc_AUC_ROC = function(predictor,response){
  
  if (length(response) != length(predictor)) {
    stop("response and predictor vector must be of same length")
  }
  
  #order the vectors so that the predictor is increasing
  predictor.order = order(predictor,decreasing=FALSE)
  predictor.sorted = predictor[predictor.order]
  response.sorted = response[predictor.order]
  
  #determine amount of zeroes and ones
  ones = sum(response)
  zeroes = length(response)-ones
  
  #if there are duplicates
  if (sum(duplicated(predictor.sorted))>0){
    #create a vector with one index for every group of duplicates
    dup_index = cumsum(duplicated(predictor.sorted)==0)
    
    #create a vector sum_vec that sums the true positives in each group of duplicates   
    df <- data.frame(duplicates=dup_index,response.sorted=response.sorted)
    sum_vec = aggregate(response.sorted ~ duplicates, data=df, sum)[,2]
    
    #create a vector that averages the maximum amount of false positives of the current group with the previous group
    fp = cumsum(response.sorted==0)
    df <- data.frame(duplicates=dup_index,fp=fp)
    max_vec = aggregate(fp ~ duplicates, data=df, max)[,2]
    top = c(0,max_vec)
    bottom = c(max_vec,0)
    average_vec = head((top+bottom)/2,-1)
    
    #AUC is the dot product of the two vectors divided by the normalizing constant
    AUC = (sum_vec%*%average_vec)/(ones*zeroes)
    
  }
  
  #if there are no duplicates
  if (sum(duplicated(predictor.sorted))==0){
    fp = cumsum(response.sorted==0)
    AUC = sum(fp * response.sorted)
    AUC = AUC/(zeroes*ones)  
    
  }
  
  return(AUC)
}

