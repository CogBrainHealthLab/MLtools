## CONNECTOME-BASED PREDICTION MODEL WITH CV-TUNED P-VALUES
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

##################################################################################################################
##################################################################################################################
##to train CPM models using a fixed p value threshold

cpm.train=function(data,outcome,p=0.05)
{
  ##checks
  #pvalues
  data=data.matrix(data)
  if(length(p)==1)  {p.posneg=c(p,p)} 
  else if (length(p)==2)  {p.posneg=p} 
  else {stop("Only one or two pvalues should be entered")}
  
  #number of rows match
  if(NROW(data)!=NROW(outcome))  {stop(paste("\nThe number of rows in the data (",NROW(data),") and outcome variable (",NROW(outcome),") do not match!\n",sep=""))}
  
  #missing data
  idx.missing=which(is.na(outcome)==T)
  if(NROW(idx.missing)>0)
  {
    cat(paste("\n",NROW(idx.missing)," missing values are detected in the outcome variable. Subjects with missing values will be excluded in the training procedure\n",sep=""))
    data=data[-idx.missing,]
    outcome=outcome[-idx.missing]
  }
  
  ##feature selection
  #critical values 
  pos.tcrit=qt((p.posneg[1]/2), NROW(outcome)-2, lower=FALSE)
  neg.tcrit=qt((p.posneg[2]/2), NROW(outcome)-2, lower=FALSE)
  pos.rcrit=sqrt(pos.tcrit^2/(pos.tcrit^2+NROW(outcome)-2))
  neg.rcrit=sqrt(neg.tcrit^2/(neg.tcrit^2+NROW(outcome)-2))
  
  #binarizing pearsons r values
  r.mat=cor(data,outcome)
  weights=rep(0,NCOL(data))
  weights[r.mat> pos.rcrit]=1
  weights[r.mat< -neg.rcrit]=-1
  
  ##network models
  ##positive model
  if(NROW(which(weights==1))>1) ## proceed only if at least 2 edges are selected
  {
    pos.netstr=rowSums(data[,which(weights==1)])  
    pos.netstr.coef=lm(outcome~pos.netstr)$coefficients
  } else 
  {
    pos.netstr.coef=NA
    cat("\nNone of edges are significantly and positively correlated with the outcome. The positive network model cannot be constructed\n")
  }
  
  ##negative model  
  if(NROW(which(weights==-1))>1) ## proceed only if at least 2 edges are selected
  {
    neg.netstr=rowSums(data[,which(weights==-1)])  
    neg.netstr.coef=lm(outcome~neg.netstr)$coefficients
  } else
  {
    neg.netstr.coef=NA
    cat("\nNone of edges are significantly and negatively correlated with the outcome. The negative network model cannot be constructed\n")
  } 
  
  ## positive + negative model  
  if(NROW(which(weights==-1))>1 & NROW(which(weights==1))>1) ## proceed only if at least 2 edges are selected in each of the earlier models
  {both.netstr.coef=lm(outcome~pos.netstr+neg.netstr)$coefficients} 
  else {both.netstr.coef=NA}
    
  ##listing objects to return 
  model=list(weights,pos.netstr.coef,neg.netstr.coef,both.netstr.coef)
  names(model)=c("weights","pos.network.coef","neg.network.coef","both.network.coef")
  return(model)
}
##################################################################################################################
##################################################################################################################
##to predict scores from previously generated cpm.train() models
  
  cpm.predict=function(model,test.data, network="both")
  {
    test.data=data.matrix(test.data)
    ##checks
    #number of rows match
    if(NROW(model$weights)!=NCOL(test.data))
    {stop(paste("\nThe number of predictors in the training data (",NROW(model$weights),") and testing data(",NCOL(test.data),")do not match!\n",sep=""))}
    
    ##select model {compute predscore}
    if(network=="positive")  {predscore=rowSums(test.data[,which(model$weights==1)])*model$pos.network.coef[2]+model$pos.network.coef[1]} 
    else if(network=="negative")  {predscore=rowSums(test.data[,which(model$weights==-1)])*model$neg.network.coef[2]+model$neg.network.coef[1]} 
    else if(network=="both")  
    {  
      #check if positive and negative models are valid, and proceed accordingly
      if(is.na(model$pos.network.coef)[1])  {positive=rep(NA,NROW(test.data))} 
      else {positive=rowSums(test.data[,which(model$weights==1)])*model$pos.network.coef[2]+model$pos.network.coef[1]}
      
      if(is.na(model$neg.network.coef)[1])  {negative=rep(NA,NROW(test.data))} 
      else {negative=rowSums(test.data[,which(model$weights==-1)])*model$neg.network.coef[2]+model$neg.network.coef[1]}
      
      if(is.na(model$pos.network.coef)[1] | is.na(model$neg.network.coef)[1])  {both=rep(NA,NROW(test.data))} 
      else  {both=rowSums(test.data[,which(model$weights==1)])*model$both.network.coef[2]+rowSums(test.data[,which(model$weights==-1)])*model$both.network.coef[3]+model$both.network.coef[1]}
      
      predscore=data.frame(positive,negative,both)
    }
    return(predscore)
  }
##################################################################################################################
##################################################################################################################
##cross-validation procedure to identify optimal p values for subsequent use of cpm.train()
  
  cpm.train.cv=function(data,outcome,p,nfolds=5)
  { 
    data=data.matrix(data)
    
    ##checks
    #require packages
    list.of.packages = c("caret")
    new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    
    if(length(new.packages)) 
    {
      cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
      install.packages(new.packages)
    }  
    
    #number of rows match
    if(NROW(data)!=NROW(outcome))  {stop(paste("\nThe number of rows in the data (",NROW(data),") and outcome variable (",NROW(outcome),") do not match!\n",sep=""))}
    
    #missing data
    idx.missing=which(is.na(outcome)==T)
    if(NROW(idx.missing)>0)
    {
      data=data[-idx.missing,]
      outcome=outcome[-idx.missing]
    }
    #check pvalues
    if(length(p)==1)  {stop("At least 2 p-values should be entered")} 
    else if (NCOL(p)==1)  {p=cbind(p,p)}
    
    ##generate geometric sequence of p-values from distribution of p values in the data
    if(missing(p))
    {
      #converting r to p
      r_to_p=function(r)
      {
        t=(r*sqrt(NROW(outcome)-2))/(sqrt(1-(r^2)))
        return(2 * (1 - pt(abs(t), NROW(outcome)-2)))
      }
      
      r.thresh=0.15 ##lower cutoff
      r.mat=cor(data,outcome)
      pos.r.mat=r.mat[r.mat>r.thresh]
      neg.r.mat=r.mat[r.mat< (-r.thresh)]
      pos.r.mat=pos.r.mat[order(pos.r.mat)]
      neg.r.mat=neg.r.mat[order(-neg.r.mat)]
      
      n.int=47
      
      #generating p values    
      p=matrix(NA,nrow =n.int-1, ncol=2)
      
      #positive model: iterating p values across fixed intervals
      pos.interval=NROW(pos.r.mat)/n.int
      for(iter in 1:(n.int-2))
      {
        p[iter+1,1]=r_to_p(pos.r.mat[round(pos.interval*iter)])  
      }
      p[1,1]=r_to_p(r.thresh) #lower cutoff r=0.15
      
      #negative model: iterating p values across fixed intervals
      neg.interval=NROW(neg.r.mat)/n.int
      for(iter in 1:(n.int-2))
      {
        p[iter+1,2]=r_to_p(neg.r.mat[round(neg.interval*iter)])  
      }
      p[1,2]=r_to_p(-r.thresh) #lower cutoff r=-0.15  
      
      #selecting p values using a geometric step progression from the fixed-interval p values; steps become progressively smaller towards the end
      p=p[c(1,9,17,24,30,35,39,42,44,45),]
    }
    
    ##setup CV folds
    folds=caret::createFolds(outcome,k=nfolds)
    
    ##training
    
    rvals=matrix(NA,nrow=NROW(p),ncol=2)
    for(iter in 1:NROW(p))
    {
      p.iter=p[iter,]    
      for (fold in 1:length(folds))
      {
        #leaving a fold out
        train_outcome=outcome[-folds[[fold]]]
        train_dat=data[-folds[[fold]],]
        
        #training the CPM and applying model to predict scores
        train.mod=cpm.train(data = train_dat,outcome = train_outcome,p=p.iter)
        predscores.fold=cpm.predict(train.mod, test.dat=data[folds[[fold]],])
        
        #saving the scores from each fold to a single vector
        if(fold==1)  {predscores.all=predscores.fold} 
        else  {predscores.all=rbind(predscores.all,predscores.fold)}
      }      
      #check if positive and negative models are valid, and proceed accordingly
      if(anyNA(predscores.all[,1]))  {pos.net.r=NA} 
      else {pos.net.r=cor(outcome[unlist(folds)],predscores.all[,1])}
      
      if(anyNA(predscores.all[,2]))  {neg.net.r=NA} 
      else {neg.net.r=cor(outcome[unlist(folds)],predscores.all[,2])}
      
      #saving the r values for the iteration
      rvals[iter,]=c(pos.net.r,neg.net.r)
    }
    
    #inform user if p values are too small such that no edges are selected; only if user-defined p values are entered  
    idx.NA.pos=which(is.na(rvals[,1]))
    if(length(idx.NA.pos)>0)
    {cat(paste("\nNote: No positive edges were selected when the p value of ",p[min(idx.NA.pos)]," or smaller was used", sep=""))}
    
    idx.NA.neg=which(is.na(rvals[,2]))
    if(length(idx.NA.neg)>0)
    {cat(paste("\nNote: No negative edges were selected when the p value of ",p[min(idx.NA.neg)]," or smaller was used", sep=""))}
    
    #identify the indices for p-values that produce the most accurate predictions
    r.pos.min.idx=which(rvals[,1]==max(rvals[,1],na.rm = T))
    r.neg.min.idx=which(rvals[,2]==max(rvals[,2],na.rm = T))
    
    #listing out objects to return
    results=list(c(p[r.pos.min.idx,1],p[r.neg.min.idx,2]),rvals,p)
    names(results)=c("opt.pvals","results","pvals")
    
    return(results)
  }
  ##################################################################################################################
  ##################################################################################################################
  ##Using the lesion approach (leave-one-network out) for CPM
  cpm.lesion=function(train.data,test.data,train.outcome, test.outcome,p)
  {
    data=data.matrix(dat_FC)
    
    ##atlas selection
    labels.url=c("https://github.com/CogBrainHealthLab/VizConnectome/blob/main/labels/labelsSC_AAL90.csv?raw=TRUE",
                 "https://github.com/CogBrainHealthLab/VizConnectome/blob/main/labels/labelsFC_schaefer119.csv?raw=TRUE",
                 "https://github.com/CogBrainHealthLab/VizConnectome/blob/main/labels/labelsFC_schaefer219.csv?raw=TRUE",
                 "https://github.com/CogBrainHealthLab/VizConnectome/blob/main/labels/labelsFC_brainnetome_yeo.csv?raw=TRUE")
    
    edge_lengths=c(4005,7021,23871,30135)
    
    if(is.na(match(NCOL(data),edge_lengths)))
    {stop("The number of columns in the input matrix is not consistent with any of the recognized parcellation schemes. The input matrix should contain 4005, 7021, 23871 or 30135 columns")} 
    else  {atlas=match(NCOL(data),edge_lengths)}
    
    ##preparing atlas labels for removing networks of edges
    label=read.csv(labels.url[atlas])
    nnode=NROW(label)
    networks.list=data.frame(unique(cbind(as.numeric(label$region),label$regionlabel)))
    names(networks.list)=c("netno","network.name")
    networks.list=networks.list[order(networks.list$netno),]
    
    ##CPM 
    results=matrix(NA,nrow=NROW(networks.list)+1,ncol=4)
    
    #training the CPM (without any network exclusions) and applying model to predict scores
    model.allnetworks=cpm.train(data=train.data, outcome=train.outcome, p=p)
    pred.allnetwork=cpm.predict(model = model.allnetworks, test.data=test.data)
    results[1,2:4]=cor(test.outcome,pred.allnetwork)
    
    #CPM with one network removed each time
    for (netno in 1:NROW(networks.list))
    {
      #identifying indice of edges to remove
      FC_matrix=array(rep(NA,nnode^2),dim=c(nnode,nnode))
      FC_matrix[upper.tri(FC_matrix, diag=FALSE)] = 1:edge_lengths[atlas]
      FC_matrix.1net=FC_matrix
      FC_matrix.1net[which(label$region==netno),which(label$region==netno)]=NA
      edge.column=FC_matrix.1net[upper.tri(FC_matrix.1net, diag=FALSE)]
      remove.idx=which(is.na(edge.column)==T)
      
      #training the CPM and applying model to predict scores
      model.1net=cpm.train(data=train.data[,-remove.idx], outcome=train.outcome, p=p)
      pred.1net=cpm.predict(model = model.1net, test.data=test.data[,-remove.idx])
      results[netno+1,2:4]=cor(test.outcome,pred.1net)
    }
    
    ##saving results in a data.frame object for returning
    results[1,1]="none removed"
    results[c(2:(NROW(networks.list)+1)),1]=paste("removed",networks.list$network.name, sep=" ")
    results=data.frame(results)
    results[, c(2:4)]=sapply(results[, c(2:4)], as.numeric)
    names(results)=c("lesion.model","positive","negative","both")
    return(results)
  }
  ##################################################################################################################
  ##################################################################################################################
  ## EXAMPLE:
  
  #source("https://github.com/CogBrainHealthLab/MLtools/blob/main/cpm.r?raw=TRUE")
  #cv.model=cpm.train.cv(data=dat, outcome=outcome)
  #model=cpm.train(data=dat_FC, outcome=dat_beh$age, p=model$opt.pvals)
  #predicted.score=cpm.predict(model = model, test.data=test.dat_FC)
