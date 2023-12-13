## CONNECTOME-BASED PREDICTION MODEL WITH CV-TUNED P-VALUES
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

##################################################################################################################
##################################################################################################################

cpm.train=function(data,outcome,p=0.05)
{
  ##checks
  #pvalues
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
  
  weights=rep(0,NCOL(data))
  ## feature selection
  for (k in 1:NCOL(data))
  {
    corr=cor.test(data[,k],outcome)
    if(corr$p.value<p.posneg[1] & corr$estimate>0)
    {
      weights[k]=1
    } else if(corr$p.value<p.posneg[2] & corr$estimate<0)
    {
      weights[k]=-1
    }
  }
  ## network models
  
  if(NROW(which(weights==1))>0)
  {
    pos.netstr=rowSums(data[,which(weights==1)])  
    pos.netstr.coef=lm(outcome~pos.netstr)$coefficients
  } else 
  {
    pos.netstr.coef=NA
    cat("\nNone of edges are significantly and positively correlated with the outcome. The positive network model cannot be constructed\n")
  }
  
  if(NROW(which(weights==-1))>0)
  {
    neg.netstr=rowSums(data[,which(weights==-1)])  
    neg.netstr.coef=lm(outcome~neg.netstr)$coefficients
  } else
  {
    neg.netstr.coef=NA
    cat("\nNone of edges are significantly and negatively correlated with the outcome. The negative network model cannot be constructed\n")
  }
  if(NROW(which(weights==-1))>0 & NROW(which(weights==1))>0)
  {
    both.netstr.coef=lm(outcome~pos.netstr+neg.netstr)$coefficients
  }
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
    positive=rowSums(test.data[,which(model$weights==1)])*model$pos.network.coef[2]+model$pos.network.coef[1]
    negative=rowSums(test.data[,which(model$weights==-1)])*model$neg.network.coef[2]+model$neg.network.coef[1]      
    both=rowSums(test.data[,which(model$weights==1)])*model$both.network.coef[2]+rowSums(test.data[,which(model$weights==-1)])*model$both.network.coef[3]+model$both.network.coef[1]
    
    predscore=data.frame(positive,negative,both)
  }
  return(predscore)
}
##################################################################################################################
##################################################################################################################

cpm.train.cv=function(data,outcome,p,nfolds=5,nthread)
{
  ## check require packages and load them
  list.of.packages = c("foreach", "doParallel")
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  
  if(length(new.packages)) 
  {
    cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
    install.packages(new.packages)
  }  
  
  `%dopar%` = foreach::`%dopar%`
  
  folds=caret::createFolds(outcome,k=nfolds)
  doParallel::registerDoParallel(nthread)
  
  rvals=foreach::foreach (iter=1:length(p), .combine=rbind,.export=c("cpm.train","cpm.predict")) %dopar% 
    {
      p.iter=p[iter]    
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
      return(c(cor(outcome[unlist(folds)],predscores.all[,1]),
               cor(outcome[unlist(folds)],predscores.all[,2])))
    }
  
  r.pos.min.idx=which(rvals[,1]==max(rvals[,1]))
  r.neg.min.idx=which(rvals[,2]==max(rvals[,2]))
  results=list(c(p[r.pos.min.idx],p[r.neg.min.idx]),rvals)
  names(results)=c("opt.pvals","results")
  return(results)
}
##################################################################################################################
##################################################################################################################
## EXAMPLE:
#p=0.05-(1:9)*0.005
#cv.model=cpm.train.cv(data=dat, outcome=outcome, p=p,nthread=10)
#model=cpm.train(data=dat_FC, outcome=dat_beh$age, p=model$opt.pvals)
#predicted.score=cpm.predict(model = model, test.data=test.dat_FC,network="positive")


