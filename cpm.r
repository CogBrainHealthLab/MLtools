## CONNECTOME-BASED PREDICTION MODEL
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

##################################################################################################################
##################################################################################################################

cpm.train=function(data,outcome,p=0.05)
{
  ##checks
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
    if(corr$p.value<p & corr$estimate>0)
      {
      weights[k]=1
      } else if(corr$p.value<p & corr$estimate<0)
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
  predscore=rowSums(test.data[,which(model$weights==1)])*model$both.network.coef[2]+rowSums(dat[,which(model$weights==-1)])*model$both.network.coef[3]+model$both.network.coef[1]
  } else if(network=="all")
  {
  positive=rowSums(test.data[,which(model$weights==1)])*model$pos.network.coef[2]+model$pos.network.coef[1]
  negative=rowSums(test.data[,which(model$weights==-1)])*model$neg.network.coef[2]+model$neg.network.coef[1]      
  both=rowSums(test.dat[,which(model$weights==1)])*model$both.network.coef[2]+rowSums(test.dat[,which(model$weights==-1)])*model$both.network.coef[3]+model$both.network.coef[1]
  predscore=data.frame(positive,negative,both)
  }
  return(predscore)
}
##################################################################################################################
##################################################################################################################
## EXAMPLE:
#model=cpm.train(data=dat_FC, outcome=dat_beh$age, p=0.05)
#predicted.score=cpm.predict(model = model, test.data=test.dat_FC,network="positive")

  
  
