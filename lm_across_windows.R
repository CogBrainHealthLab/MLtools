## Linear model analyses across FC time windows
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

##################################################################################################################
##################################################################################################################

lm_across_window=function(window_duration,nframes,nthread, timeseries,model)
{
  ##check and remove subjects with incomplete data
  idx.na=which(!complete.cases(model))
  if(length(idx.na)>0)
  {
    model=model[-idx.na,]
    timeseries=timeseries[-idx.na]
    cat(paste0(length(idx.na)," subjects with missing data detected"))
  }
  ## correlation function to be used within lapply()
  window.corr=function(ts,start,end)
  {
    cor.mat=cor(ts[start:end,])
    cor.mat[upper.tri(cor.mat,diag = F)]
  }

  ##activate parallel clusters
  cl=parallel::makeCluster(nthread)
  doParallel::registerDoParallel(nthread)
  `%dopar%` = foreach::`%dopar%`

  ##looping across different window lengths
  av.coef.vector=list()
  for (window_length in 1:length(window_duration))
  {
    nwindows=nframes-window_duration[window_length]+1

    av.coef.vector[[window_length]]=foreach::foreach(window=1:nwindows, .combine="c")  %dopar%
      {
        #perform correlation analyses within each list element
        all.cor.vector=lapply(timeseries, window.corr,start=1+window-1,end=window_duration[window_length]+window-1) 
        #flatten all list elements into a 2D matrix
        all.cor.vector=scale(do.call(rbind,all.cor.vector))
        
        mod=.lm.fit(y=all.cor.vector,x=data.matrix(cbind(1,model)))
        coef=mod$coefficients[NCOL(model)+1,]
        return(mean(abs(coef)))
      }
  }
  return(av.coef.vector)  
  closeAllConnections()
}
