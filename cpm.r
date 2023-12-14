## CONNECTOME-BASED PREDICTION MODEL WITH CV-TUNED P-VALUES
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

##################################################################################################################
##################################################################################################################

cpm.train=function(data,outcome,p=0.05)
{
  ##checks
  #pvalues
  data=data.matrix(data)
  if(length(p)==1)
  {
    p.posneg=c(p,p)
  } else if (length(p)==2)
  {
    p.posneg=p
  } else 
  {
    stop("Only one or two pvalues should be entered")
  }
  #number of rows match
  if(NROW(data)!=NROW(outcome))
  {
    stop(paste("\nThe number of rows in the data (",NROW(data),") and outcome variable (",NROW(outcome),") do not match!\n",sep=""))
  }
  #missing data
  idx.missing=which(is.na(outcome)==T)
  if(NROW(idx.missing)>0)
  {
    cat(paste("\n",NROW(idx.missing)," missing values are detected in the outcome variable. Subjects with missing values will be excluded in the training procedure\n",sep=""))
    data=data[-idx.missing,]
    outcome=outcome[-idx.missing]
  }
  
  ## feature selection
  weights=rep(NA,NCOL(data))

  pos.tcrit=qt((p.posneg[1]/2), NROW(outcome)-2, lower=FALSE)
  neg.tcrit=qt((p.posneg[2]/2), NROW(outcome)-2, lower=FALSE)
  pos.rcrit=sqrt(pos.tcrit^2/(pos.tcrit^2+NROW(outcome)-2))
  neg.rcrit=sqrt(neg.tcrit^2/(neg.tcrit^2+NROW(outcome)-2))
  
  r.mat=cor(data,outcome)
  
  weights[r.mat> pos.rcrit]=1
  weights[r.mat< -neg.rcrit]=-1

  ## network models
  if(NROW(which(weights==1))>1)
  {
    pos.netstr=rowSums(data[,which(weights==1)])  
    pos.netstr.coef=lm(outcome~pos.netstr)$coefficients
  } else 
  {
    pos.netstr.coef=NA
    cat("\nNone of edges are significantly and positively correlated with the outcome. The positive network model cannot be constructed\n")
  }
  if(NROW(which(weights==-1))>1)
  {
    neg.netstr=rowSums(data[,which(weights==-1)])  
    neg.netstr.coef=lm(outcome~neg.netstr)$coefficients
  } else
  {
    neg.netstr.coef=NA
    cat("\nNone of edges are significantly and negatively correlated with the outcome. The negative network model cannot be constructed\n")
  }
  if(NROW(which(weights==-1))>1 & NROW(which(weights==1))>1)
  {
    both.netstr.coef=lm(outcome~pos.netstr+neg.netstr)$coefficients
  } else {both.netstr.coef=NA}
  model=list(weights,pos.netstr.coef,neg.netstr.coef,both.netstr.coef)
  names(model)=c("weights","pos.network.coef","neg.network.coef","both.network.coef")
  return(model)
}
##################################################################################################################
##################################################################################################################

cpm.predict=function(model,test.data, network="both")
{
  ##checks
  if(NROW(model$weights)!=NCOL(test.data))
  {
    stop(paste("\nThe number of predictors in the training data (",NROW(model$weights),") and testing data(",NCOL(test.data),")do not match!\n",sep=""))
  }
  if(network=="positive")
  {
    predscore=rowSums(test.data[,which(model$weights==1)])*model$pos.network.coef[2]+model$pos.network.coef[1]
  } else if(network=="negative")
  {
    predscore=rowSums(test.data[,which(model$weights==-1)])*model$neg.network.coef[2]+model$neg.network.coef[1]      
  } else if(network=="both")
  {
    if(is.na(model$pos.network.coef)[1])
    {
      positive=rep(NA,NROW(test.data))
    } else {positive=rowSums(test.data[,which(model$weights==1)])*model$pos.network.coef[2]+model$pos.network.coef[1]}
    
    if(is.na(model$neg.network.coef)[1])
    {
      negative=rep(NA,NROW(test.data))
    } else {negative=rowSums(test.data[,which(model$weights==-1)])*model$neg.network.coef[2]+model$neg.network.coef[1]}
    
    if(is.na(model$pos.network.coef)[1] | is.na(model$neg.network.coef)[1])
    {
      both=rep(NA,NROW(test.data))
    } else
    {
      both=rowSums(test.data[,which(model$weights==1)])*model$both.network.coef[2]+rowSums(test.data[,which(model$weights==-1)])*model$both.network.coef[3]+model$both.network.coef[1]
    }
    predscore=data.frame(positive,negative,both)
  }
  return(predscore)
}
##################################################################################################################
##################################################################################################################

cpm.train.cv=function(data,outcome,p,nfolds=5)
{ 
  data=data.matrix(data)
  ##checks
  #require packages
  list.of.packages = c("caret","foreach", "doParallel", "parallel")
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  
  if(length(new.packages)) 
  {
    cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
    install.packages(new.packages)
  }  
  #number of rows match
  if(NROW(data)!=NROW(outcome))
  {
    stop(paste("\nThe number of rows in the data (",NROW(data),") and outcome variable (",NROW(outcome),") do not match!\n",sep=""))
  }
  
  #missing data
  idx.missing=which(is.na(outcome)==T)
  if(NROW(idx.missing)>0)
  {
    data=data[-idx.missing,]
    outcome=outcome[-idx.missing]
  }
  
  #generate sequence of p-values to tune
  if(missing(p))
  {
    r_to_p=function(r)
    {
      t=(r*sqrt(NROW(outcome)-2))/(sqrt(1-(r^2)))
      return(2 * (1 - pt(abs(t), NROW(outcome)-2)))
    }
    r.thresh=0.15
    r.mat=cor(data,outcome)
    pos.r.mat=r.mat[r.mat>r.thresh]
    neg.r.mat=r.mat[r.mat< (-r.thresh)]
    pos.r.mat=pos.r.mat[order(pos.r.mat)]
    neg.r.mat=neg.r.mat[order(-neg.r.mat)]
   
    n.int=11
    
    p=matrix(NA,nrow =n.int-1, ncol=2)
    pos.interval=NROW(pos.r.mat)/n.int
    for(iter in 1:(n.int-2))
    {
      p[iter+1,1]=r_to_p(pos.r.mat[round(pos.interval*iter)])  
    }
    p[1,1]=r_to_p(r.thresh)
    
    neg.interval=NROW(neg.r.mat)/n.int
    for(iter in 1:(n.int-2))
    {
      p[iter+1,2]=r_to_p(neg.r.mat[round(neg.interval*iter)])  
    }
    p[1,2]=r_to_p(-r.thresh)  
    
  }
  #check pvalues
  
  if(length(p)==1)
   {
     stop("At least 2 p-values should be entered")
   } else if (NCOL(p)==1)
   {
     p=cbind(p,p)
   }
  
  ##setup environment
  folds=caret::createFolds(outcome,k=nfolds)

  
  ##training
  rvals=matrix(NA,nrow=NROW(p),ncol=2)
  for(iter in 1:NROW(p))
    {
      p.iter=p[iter,]    
      for (fold in 1:length(folds))
      {
        train_outcome=outcome[-folds[[fold]]]
        train_dat=data[-folds[[fold]],]
        train.mod=cpm.train(data = train_dat,outcome = train_outcome,p=p.iter)
        predscores.fold=cpm.predict(train.mod, test.dat=data[folds[[fold]],])
        if(fold==1)
        {
          predscores.all=predscores.fold
        } else
        {
          predscores.all=rbind(predscores.all,predscores.fold)
        }
      }      
      if(anyNA(predscores.all[,1]))
      {
        pos.net.r=NA  
      } else {pos.net.r=cor(outcome[unlist(folds)],predscores.all[,1])}
      
      if(anyNA(predscores.all[,2]))
      {
        neg.net.r=NA  
      } else {neg.net.r=cor(outcome[unlist(folds)],predscores.all[,2])}

      rvals[iter,]=c(pos.net.r,neg.net.r)
  }
  r.pos.min.idx=which(rvals[,1]==max(rvals[,1],na.rm = T))
  r.neg.min.idx=which(rvals[,2]==max(rvals[,2],na.rm = T))
  results=list(c(p[r.pos.min.idx,1],p[r.neg.min.idx,2]),rvals,p)
  names(results)=c("opt.pvals","results","pvals")
  
  idx.NA.pos=which(is.na(rvals[,1]))
  if(length(idx.NA.pos)>0)
  {cat(paste("\nNote: No positive edges were selected when the p value of ",p[min(idx.NA.pos)]," or smaller was used", sep=""))}
  
  idx.NA.neg=which(is.na(rvals[,2]))
  if(length(idx.NA.neg)>0)
  {cat(paste("\nNote: No negative edges were selected when the p value of ",p[min(idx.NA.neg)]," or smaller was used", sep=""))}

  return(results)
}
##################################################################################################################
##################################################################################################################
## EXAMPLE:
#p=0.05-(1:9)*0.005
#cv.model=cpm.train.cv(data=dat, outcome=outcome, p=p,nthread=10)
#model=cpm.train(data=dat_FC, outcome=dat_beh$age, p=model$opt.pvals)
#predicted.score=cpm.predict(model = model, test.data=test.dat_FC,network="positive")
