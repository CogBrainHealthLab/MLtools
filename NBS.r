## A QUICKER AND MORE EFFICIENT NETWORK-BASED-STATISTICS
## ADAPTED FROM NBR::nbr_lm()
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY
############################################################################################################################
############################################################################################################################
### extract t-stats efficiently
extract.t=function(mod,row)
{
  p = mod$rank
  df.residual=NROW(mod$residuals)-NROW(mod$coefficients)
  rdf = df.residual
  Qr = mod$qr
  p1 = 1L:p
  r = mod$residuals
  rss = colSums(r^2)
  resvar = rss/rdf
  R = chol2inv(Qr[p1, p1, drop = FALSE])  
  se = (sqrt(diag(R) %*% t(resvar)))[row,]
  est = mod$coefficients[row,]
  tval = est/se 
  return(tval)
}
############################################################################################################################
############################################################################################################################
### extract cluster-stats from graphs
cluster.stat=function(data,nnodes,tcrit)
{
  ##thresholding
  tstat.thresholded=data
  tstat.thresholded[abs(data)<tcrit]=0
  tstat.thresholded.bin=tstat.thresholded
  tstat.thresholded.bin[abs(tstat.thresholded.bin)>0]=1
  
  ##setting up FCmatrices
  nnodes=(0.5 + sqrt(0.5^2 - 4 * 0.5 * -length(data))) / (2 * 0.5)
  FC_mat.unweighted=matrix(0,nrow=nnodes,ncol=nnodes)
  FC_mat.weighted=matrix(0,nrow=nnodes,ncol=nnodes)
  
  ##thresholding
  FC_mat.weighted[upper.tri(FC_mat.weighted,diag = F)]=tstat.thresholded
  FC_mat.unweighted[upper.tri(FC_mat.unweighted,diag = F)]=tstat.thresholded.bin
  FC_mat.weighted=abs(FC_mat.weighted)-(FC_mat.unweighted*tcrit)
  #clustering
  com=igraph::components(igraph::graph_from_adjacency_matrix(FC_mat.unweighted, mode='undirected', weighted=NULL))
  
  #count edges in clusters
  if(length(which(com$csize>2)>0)) #proceed only if there is at least one cluster with 3 nodes (i.e 2 edges). Isolated/unconnected edges are removed
  {
    cluster.idx=which(com$csize>2)
    clust.results=matrix(NA,nrow=length(cluster.idx), ncol=2)
    
    for (cluster.no in 1:length(cluster.idx))
    {
      idx=which(com$membership==cluster.idx[cluster.no])
      clust.results[cluster.no,1]=strength.unweighted=sum(FC_mat.unweighted[idx,idx])
      clust.results[cluster.no,2]=strength.weighted=sum(FC_mat.weighted[idx,idx])
    }
  } else  { clust.results=c(0,0)} #if no clusters (with 3 nodes) are detected
  return(clust.results)
}

############################################################################################################################
############################################################################################################################
NBS=function(all_predictors,IV_of_interest, FC_data, nperm=50, nthread=1, p=0.001)
{
  all_predictors=data.matrix(all_predictors)
  FC_data=data.matrix(FC_data)
  ##checks
  #check required packages
  list.of.packages = c("parallel", "doParallel","igraph","doSNOW","foreach")
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) 
  {
    cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
    install.packages(new.packages)
  }  
  
  #check if nrow is consistent for all_predictors and FC_data
  if(NROW(FC_data)!=NROW(all_predictors))  {stop(paste("The number of rows for FC_data (",NROW(FC_data),") and all_predictors (",NROW(all_predictors),") are not the same",sep=""))}
  
  #incomplete data check
  idxF=which(complete.cases(all_predictors)==F)
  if(length(idxF)>0)
  {
    cat(paste("all_predictors contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis\n"))
    all_predictors=all_predictors[-idxF,]
    IV_of_interest=IV_of_interest[-idxF]
    FC_data=FC_data[-idxF,]
  }
  
  #check IV_of_interest
  if(NCOL(all_predictors)>1)
  {
    for(colno in 1:(NCOL(all_predictors)+1))
    {
      if(colno==(NCOL(all_predictors)+1))  {stop("IV_of_interest is not contained within all_predictors")}
      
      if(class(IV_of_interest) != "integer" & class(IV_of_interest) != "numeric") 
      {
        if(identical(IV_of_interest,all_predictors[,colno]))  {break} 
      } else 
      {
        if(identical(as.numeric(IV_of_interest),as.numeric(all_predictors[,colno])))  {break}
      }
    }
  }  else
  {
    if(class(IV_of_interest) != "integer" & class(IV_of_interest) != "numeric") 
    {
      if(identical(IV_of_interest,all_predictors))  {colno=1} 
      else  {stop("IV_of_interest is not contained within all_predictors")}
    } else
    {
      if(identical(as.numeric(IV_of_interest),as.numeric(all_predictors)))  {colno=1}
      else  {stop("IV_of_interest is not contained within all_predictors")}
    }
  }
  
  #check categorical variable
  if(NCOL(all_predictors)>1)
  {
    for (column in 1:NCOL(all_predictors))
    {
      if(class(all_predictors[,column]) != "integer" & class(all_predictors[,column]) != "numeric")
      {
        if(length(unique(all_predictors[,column]))==2)
        {
          cat(paste("The binary variable '",colnames(all_predictors)[column],"' will be recoded with ",unique(all_predictors[,column])[1],"=0 and ",unique(all_predictors[,column])[2],"=1 for the analysis\n",sep=""))
          
          recode=rep(0,NROW(all_predictors))
          recode[all_predictors[,column]==unique(all_predictors[,column])[2]]=1
          all_predictors[,column]=recode
          IV_of_interest=all_predictors[,colno]
        } else if(length(unique(all_predictors[,column]))>2)    {stop(paste("The categorical variable '",colnames(all_predictors)[column],"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
      }      
    }
  } else
  {
    if (!suppressWarnings(all(!is.na(as.numeric(as.character(all_predictors)))))) 
    {
      if(length(unique(all_predictors))==2)
      {
        cat(paste("The binary variable '",colnames(all_predictors),"' will be recoded such that ",unique(all_predictors)[1],"=0 and ",unique(all_predictors)[2],"=1 for the analysis\n",sep=""))
        
        recode=rep(0,NROW(all_predictors))
        recode[all_predictors==unique(all_predictors)[2]]=1
        all_predictors=recode
        IV_of_interest=all_predictors
      } else if(length(unique(all_predictors))>2)    {stop(paste("The categorical variable '",colnames(all_predictors),"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
    }      
  }
  
  #collinearity check
  if(NCOL(all_predictors)>1)
  {
    cormat=cor(all_predictors,use = "pairwise.complete.obs")
    cormat.0=cormat
    cormat.0[cormat.0==1]=NA
    if(max(abs(cormat.0),na.rm = T) >0.5)
    {
      warning(paste("correlations among variables in all_predictors are observed to be as high as ",round(max(abs(cormat.0),na.rm = T),2),", suggesting potential collinearity among predictors.\nAnalysis will continue...",sep=""))
    }  
  }
  
  ##unpermuted model
  mod=.lm.fit(y = FC_data,x=data.matrix(cbind(1,all_predictors)))
  
  #define/init variables
  t.orig=extract.t(mod,colno+1)
  nnodes=(0.5 + sqrt(0.5^2 - 4 * 0.5 * -NCOL(FC_data))) / (2 * 0.5)
  tcrit=qt(p/2, NROW(all_predictors)-NCOL(all_predictors)-1, lower=FALSE)
  orig.clust=cluster.stat(t.orig,nnodes,tcrit)
  
  remove(mod)
  
  ##permuted models
  #generating permutation sequences  
  permseq=matrix(NA, nrow=NROW(all_predictors), ncol=nperm)
  for (perm in 1:nperm)  {permseq[,perm]=sample(1:NROW(all_predictors))}

  #activate parallel processing
  unregister_dopar = function() {
    env = foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }
  unregister_dopar()
  
  cl=parallel::makeCluster(nthread)
  doParallel::registerDoParallel(nthread)
  `%dopar%` = foreach::`%dopar%`
  
  #progress bar
  doSNOW::registerDoSNOW(cl)
  pb=txtProgressBar(max = nperm, style = 3)
  progress=function(n) setTxtProgressBar(pb, n)
  opts=list(progress = progress)
  
  
  start=Sys.time()
  cat("\nEstimating permuted network strengths...\n")
  
  max.netstr=foreach::foreach(perm=1:nperm, .combine="rbind",.export=c("extract.t","cluster.stat"), .options.snow = opts)  %dopar%
    {
      #fitting permuted regression model and extracting max netstr in parallel streams
      mod.permuted=.lm.fit(y = FC_data,x=data.matrix(cbind(1,all_predictors))[permseq[,perm],])
      t.perm=extract.t(mod.permuted,colno+1)
      netstr=cluster.stat(t.perm,nnodes,tcrit)
      
      remove(t.perm,mod.permuted)
      
      if(length(netstr)>2)  {max.netstr=c(max(netstr[,1]),max(netstr[,2]))} 
      else {max.netstr=netstr} #if only one row of results is obtained, there is no need to use the max() function
      
      return(max.netstr)
    }
  end=Sys.time()
  cat(paste("\nCompleted in :",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))
  
  ##processing results
  #saving cluster-related results into a data.frame object
  orig.clust=data.frame(orig.clust)
  orig.clust$p.unweighted=NA
  orig.clust$p.weighted=NA
  
  #thresholding clusters using permuted null distribution
  for(row in 1:nrow(orig.clust))
  {
    orig.clust[row,3]=sum(max.netstr[,1] > orig.clust[row,1])/nperm
    orig.clust[row,4]=sum(max.netstr[,2] > orig.clust[row,2])/nperm
  }
  
  #formatting results table
  orig.clust[,c(3,4)][orig.clust[,c(3,4)]==0]=paste("<",1/nperm,sep="") #if p=0
  orig.clust=cbind(c(1:nrow(orig.clust)),orig.clust)
  colnames(orig.clust)=c("network no.","strength.unweighted","strength.weighted","pFWE.unweighted","pFWE.weighted")
  
  #objects to return
  returnobj=list(orig.clust,t.orig, tcrit,max.netstr)
  names(returnobj)=c("results","t.orig","tcrit","max.netstr")
  return(returnobj)
}
############################################################################################################################
############################################################################################################################
extract.edges=function(NBS.obj,clust.no=1)
{
  nnodes=(0.5 + sqrt(0.5^2 - 4 * 0.5 * -length(NBS.obj$t.orig))) / (2 * 0.5)
  ##recode all p="<0.**" into 0 for subsequent thresholding
  if(is.character(NBS.obj$results[,4]))
  {
    NBS.obj$results[,4]=suppressWarnings(as.numeric(NBS.obj$results[,4]))
    NBS.obj$results[,4][is.na(NBS.obj$results[,4])]=0
  }
  if(is.character(NBS.obj$results[,5]))
  {
    NBS.obj$results[,5]=suppressWarnings(as.numeric(NBS.obj$results[,5]))
    NBS.obj$results[,5][is.na(NBS.obj$results[,5])]=0
  }
  
  ##thresholding tstats
  tstat.thresholded=NBS.obj$t.orig
  tstat.thresholded[abs(NBS.obj$t.orig)<NBS.obj$tcrit]=0
  tstat.thresholded.bin=tstat.thresholded
  tstat.thresholded.bin[abs(tstat.thresholded.bin)>0]=1
  
  ##reshaping 1D tstat vector to 2D matrices
  FC_mat.unweighted=matrix(0,nrow=nnodes,ncol=nnodes)
  FC_mat.weighted=matrix(0,nrow=nnodes,ncol=nnodes)
  
  FC_mat.weighted[upper.tri(FC_mat.weighted,diag = F)]=tstat.thresholded-(tstat.thresholded.bin*NBS.obj$) ## subtracting tcrit values to be consist with NBR::nbr_lm()
  FC_mat.unweighted[upper.tri(FC_mat.unweighted,diag = F)]=tstat.thresholded.bin
  
  ##clustering
  com=igraph::components(igraph::graph_from_adjacency_matrix(FC_mat.unweighted, mode='undirected', weighted=NULL))
  idx=which(com$membership==clust.no)
  
  ##masking out edges from other networks
  FC_mat.mask=matrix(0,nrow=nnodes,ncol=nnodes)
  FC_mat.mask[idx,idx]=1
  mask=FC_mat.mask[upper.tri(FC_mat.mask,diag = F)]
  clust.tstat=tstat.thresholded*mask
  
  #positive mask
  clust.pos.mask=clust.tstat
  clust.pos.mask[clust.pos.mask>0]=1
  clust.pos.mask[clust.pos.mask<0]=0
  
  #negative mask
  clust.neg.mask=clust.tstat
  clust.neg.mask[clust.neg.mask>0]=0
  clust.neg.mask[clust.neg.mask<0]=1
  
  ##objects to return
  returnobj=list(as.numeric(clust.tstat),as.numeric(clust.pos.mask),as.numeric(clust.neg.mask))
  names(returnobj)=c("clust.tstat","clust.pos.mask","clust.neg.mask")
  return(returnobj)
}
############################################################################################################################
############################################################################################################################
##source("https://github.com/CogBrainHealthLab/MLtools/blob/main/NBS.r?raw=TRUE")
