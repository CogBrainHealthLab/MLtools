
##perform regression
extract.t=function(mod,row)
{
  p = mod$rank
  rdf = mod$df.residual
  Qr = mod$qr
  p1 = 1L:p
  r = mod$residuals
  rss = colSums(r^2)
  resvar = rss/rdf
  R = chol2inv(Qr$qr[p1, p1, drop = FALSE])  
  se = (sqrt(diag(R) %*% t(resvar)))[row,]
  est = mod$coefficients[row,]
  tval = est/se                          
}
############################################################################################################################
############################################################################################################################
cluster.stat=function(data,nnodes,tcrit)
{
  #thresholding
  tstat.thresholded=data
  tstat.thresholded[abs(data)<tcrit]=0
  tstat.thresholded.bin=tstat.thresholded
  tstat.thresholded.bin[abs(tstat.thresholded.bin)>0]=1
  
  #setting up FCmatrices
  nnodes=(0.5 + sqrt(0.5^2 - 4 * 0.5 * -NCOL(FC_data))) / (2 * 0.5)
  FC_mat.unweighted=matrix(0,nrow=nnodes,ncol=nnodes)
  FC_mat.weighted=matrix(0,nrow=nnodes,ncol=nnodes)
  
  #thresholding
  FC_mat.weighted[upper.tri(FC_mat.weighted,diag = F)]=tstat.thresholded
  FC_mat.unweighted[upper.tri(FC_mat.unweighted,diag = F)]=tstat.thresholded.bin
  FC_mat.weighted=abs(FC_mat.weighted)-(FC_mat.unweighted*tcrit)
  #clustering
  com=igraph::components(igraph::graph_from_adjacency_matrix(FC_mat.unweighted, mode='undirected', weighted=NULL))
  
  #count edges in clusters
  if(length(which(com$csize>2)>0))
  {
    cluster.idx=which(com$csize>2)
    clust.results=matrix(NA,nrow=length(cluster.idx), ncol=2)
    
    for (cluster.no in 1:length(cluster.idx))
    {
      idx=which(com$membership==cluster.idx[cluster.no])
      clust.results[cluster.no,1]=strength.unweighted=sum(FC_mat.unweighted[idx,idx])
      clust.results[cluster.no,2]=strength.weighted=sum(FC_mat.weighted[idx,idx])
    }
  } else
  {
    clust.results=c(0,0)
  }
  return(clust.results)
}

############################################################################################################################
############################################################################################################################
NBS=function(all_predictors,IV_of_interest, FC_data, nperm=50, nthread, p=0.001)
{
  ##unpermuted model
  mod=lm(FC_data~data.matrix(all_predictors))
  
  #identify contrast
  for(colno in 1:(NCOL(all_predictors)+1))
  {
    if(colno==(NCOL(all_predictors)+1))
    {stop("IV_of_interest is not contained within all_predictors")}
    if(identical(IV_of_interest,all_predictors[,colno]))
    {break}
  }
  
  #define variables
  t.orig=extract.t(mod,colno+1)
  nnodes=(0.5 + sqrt(0.5^2 - 4 * 0.5 * -NCOL(FC_data))) / (2 * 0.5)
  tcrit=qt(p/2, NROW(all_predictors)-2, lower=FALSE)
  orig.clust=cluster.stat(t.orig,nnodes,tcrit)
  
  ##permuted models
  ##generating permutation sequences  
  permseq=matrix(NA, nrow=NROW(all_predictors), ncol=nperm)
  for (perm in 1:nperm)  {permseq[,perm]=sample.int(NROW(all_predictors))}
  
  #activate parallel processing
  cl=parallel::makeCluster(nthread)
  doParallel::registerDoParallel(cl)
  `%dopar%` = foreach::`%dopar%`
  
  #progress bar
  doSNOW::registerDoSNOW(cl)
  pb=txtProgressBar(max = nperm, style = 3)
  progress=function(n) setTxtProgressBar(pb, n)
  opts=list(progress = progress)
  
  ##fitting permuted regression model and extracting max netstr in parallel streams
  start=Sys.time()
  cat("\nEstimating permuted network strengths...\n")
  
  max.netstr=foreach::foreach(perm=1:nperm, .combine="rbind",.export=c("extract.t","cluster.stat"), .options.snow = opts)  %dopar%
    {
      all_predictors.permuted=all_predictors
      mod.permuted=lm(FC_data~data.matrix(all_predictors.permuted)[permseq[,perm],])
      t.perm=extract.t(mod.permuted,colno+1)
      netstr=cluster.stat(t.perm,nnodes,tcrit)
      if(length(netstr)>2)
      {
        max.netstr=c(max(netstr[,1]),max(netstr[,1]))  
      } else 
      {
        max.netstr=netstr
      }
      return(max.netstr)
    }
  end=Sys.time()
  cat(paste("\nCompleted in :",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))
  
  orig.clust=data.frame(orig.clust)
  orig.clust$p.unweighted=NA
  orig.clust$p.weighted=NA
  
  for(row in 1:nrow(orig.clust))
  {
    orig.clust[row,3]=length(which(max.netstr[,1]>orig.clust[row,1]))/nperm
    orig.clust[row,4]=length(which(max.netstr[,2]>orig.clust[row,2]))/nperm
  }
  
  orig.clust[,c(3,4)][orig.clust[,c(3,4)]==0]=paste("<",1/nperm,sep="")
  
  orig.clust=cbind(c(1:1:nrow(orig.clust)),orig.clust)
  colnames(orig.clust)=c("network no.","strength.unweighted","strength.weighted","pFWE.unweighted","pFWE.weighted")
  returnobj=list(orig.clust,t.orig, tcrit)
  names(returnobj)=c("results","t.orig","tcrit")
  return(returnobj)
}
############################################################################################################################
############################################################################################################################
extract.edges=function(NBS.obj,clust.no=1)
{
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
  
  tstat.thresholded=NBS.obj$t.orig
  tstat.thresholded[abs(NBS.obj$t.orig)<NBS.obj$tcrit]=0
  tstat.thresholded.bin=tstat.thresholded
  tstat.thresholded.bin[abs(tstat.thresholded.bin)>0]=1
  

  #setting up FCmatrices
  nnodes=(0.5 + sqrt(0.5^2 - 4 * 0.5 * -length(NBS.obj$t.orig))) / (2 * 0.5)
  FC_mat.unweighted=matrix(0,nrow=nnodes,ncol=nnodes)
  FC_mat.weighted=matrix(0,nrow=nnodes,ncol=nnodes)
  
  #thresholding
  FC_mat.weighted[upper.tri(FC_mat.weighted,diag = F)]=tstat.thresholded
  FC_mat.unweighted[upper.tri(FC_mat.unweighted,diag = F)]=tstat.thresholded.bin
  FC_mat.weighted=abs(FC_mat.weighted)-(FC_mat.unweighted*NBS.obj$tcrit)
  #clustering
  com=igraph::components(igraph::graph_from_adjacency_matrix(FC_mat.unweighted, mode='undirected', weighted=NULL))
  
  idx=which(com$membership==clust.no)
  
  #masking out edges from other networks
  FC_mat.mask=matrix(0,nrow=nnodes,ncol=nnodes)
  FC_mat.mask[idx,idx]=1
  mask=FC_mat.mask[upper.tri(FC_mat.mask,diag = F)]
  clust.tstat=tstat.thresholded*mask
  clust.pos.mask=clust.tstat
  clust.pos.mask[clust.pos.mask>0]=1
  clust.pos.mask[clust.pos.mask<0]=0
  
  clust.neg.mask=clust.tstat
  clust.neg.mask[clust.neg.mask>0]=0
  clust.neg.mask[clust.neg.mask<0]=1
  
  returnobj=list(as.numeric(clust.tstat),as.numeric(clust.pos.mask),as.numeric(clust.neg.mask))
  names(returnobj)=c("clust.tstat","clust.pos.mask","clust.neg.mask")
  return(returnobj)
}
############################################################################################################################
############################################################################################################################
##source("https://github.com/CogBrainHealthLab/MLtools/blob/main/NBS.r?raw=TRUE")
