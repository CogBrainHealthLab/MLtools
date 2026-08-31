extractmetric.bysex=function(model,test_feat, test_outcome)
{
  #PLSR model has to be treated differently
  if(class(model)[1]=="mvr")
  {
    rmsep=RMSEP(model)
    pred_outcome=predict(model,test_feat, ncomp=which.min(data.frame(rmsep$val)[2,2:21]))
    return(list(c(cor(pred_outcome,test_outcome),mean(abs(pred_outcome-test_outcome)),cor((pred_outcome-test_outcome),test_outcome)),pred_outcome))
  } else
  {
    #all other regression models
    pred_outcome=kernlab::predict(model,test_feat)
    return(list(c(cor(pred_outcome,test_outcome),mean(abs(pred_outcome-test_outcome)),cor((pred_outcome-test_outcome),test_outcome)),pred_outcome))
  }
}

##Runs regression models
pred.allmodels.bysex=function(train_outcome, train_feat,train_sex,test_outcome, test_feat,test_sex, xgb=F, harm="combat", eb=F, train_site, test_site, model='lm', formula=y~outcome)
{
  #check if train_feat contains columns of 0s, if so, these columns are removed
  col0_idx=union(which(colSums(train_feat)==0),which(colSums(test_feat)==0))
  if(length(col0_idx)>0)
  {
    train_feat=train_feat[,-col0_idx]
    test_feat=test_feat[,-col0_idx]
  }
  
  #split datasets by sex
  train.M.idx=which(train_sex==0)
  train.F.idx=which(train_sex==1)
  
  test.M.idx=which(test_sex==0)
  test.F.idx=which(test_sex==1)
  
  train_outcome.bysex=list(train_outcome[train.M.idx],train_outcome[train.F.idx])
  train_feat.bysex=list(train_feat[train.M.idx,],train_feat[train.F.idx,])
  
  if(missing("train_site"))
  {
    train_site=rep("train",length(train_outcome))
  }
  train_site.bysex=list(train_site[train.M.idx],train_site[train.F.idx])
  
  if(missing("test_site"))
  {
    test_site=rep("test",length(test_outcome))
  }
  test_site.bysex=list(test_site[test.M.idx],test_site[test.F.idx])
  
  test_outcome.bysex=list(test_outcome[test.M.idx],test_outcome[test.F.idx])
  test_feat.bysex=list(test_feat[test.M.idx,],test_feat[test.F.idx,])

  #set default formula
  if(model=='lm' && missing("formula")) {formula=y~outcome}
  if(model=='gam' && missing("formula")) {formula=y~s(outcome)}
  
  remove(test_outcome,train_outcome,test_feat,train_feat,train.M.idx,train.F.idx)

  #activate parallel processing
  unregister_dopar = function() {
    env = foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }
  unregister_dopar()

# Detect if the operating system is macOS (Darwin)
if (Sys.info()[["sysname"]] == "Darwin") {
  cluster_type <- "PSOCK"
} else {
  cluster_type <- "FORK"
}
# Initialize the cluster safely based on the OS
  registerDoParallel(cl)
  cl=parallel::makeCluster(2, type = cluster_type)
  doParallel::registerDoParallel(cl)
  `%dopar%` = foreach::`%dopar%`
  
  results=foreach::foreach(sex=1:2, .combine="c",.packages = c("glmnet","pls","kernlab","ComBatFamily"), .export ="extractmetric.bysex")  %dopar%
    {
      ##harmonization by sex
      dat.all=rbind(data.matrix(train_feat.bysex[[sex]]),data.matrix(test_feat.bysex[[sex]]))
        ##combat
        if(harm=="combat") 
        {
          covar=data.frame(outcome=c(train_outcome.bysex[[sex]],test_outcome.bysex[[sex]]))
          dat.harmonized =comfam(data=dat.all,covar=covar,model=model,eb=eb,
                          bat=as.factor(c(train_site.bysex[[sex]],test_site.bysex[[sex]])),
                          formula=formula)
          
          train_feat.bysex[[sex]]=dat.harmonized$dat.combat[1:length(train_outcome.bysex[[sex]]),]
          test_feat.bysex[[sex]]=dat.harmonized$dat.combat[(length(train_outcome.bysex[[sex]])+1):(length(train_outcome.bysex[[sex]])+length(test_outcome.bysex[[sex]])),]  
        }
        
        if(harm=="covbat") 
        {
          covar=data.frame(outcome=c(train_outcome.bysex[[sex]],test_outcome.bysex[[sex]]))
          dat.harmonized =covfam(data=dat.all,covar=covar,model=model,eb=eb,
                                 bat=as.factor(c(train_site.bysex[[sex]],test_site.bysex[[sex]])),
                                 formula=formula)
          
          train_feat.bysex[[sex]]=dat.harmonized$dat.covbat[1:length(train_outcome.bysex[[sex]]),]
          test_feat.bysex[[sex]]=dat.harmonized$dat.covbat[(length(train_outcome.bysex[[sex]])+1):(length(train_outcome.bysex[[sex]])+length(test_outcome.bysex[[sex]])),]  
        }
        remove(dat.harmonized)
      

      #setting up results matrix
      predmetrics=matrix(NA,nrow=10, ncol=4)
      predmetrics[,1]=c("RidgeR", "LassoR","ElasticNet","PLSR","GPR (Linear)","SVM (Linear)","KQR (Linear)", "GPR (RBF)", "SVM (RBF)","KQR (RBF)")
      
      predscores=matrix(NA,nrow=length(test_outcome.bysex[[sex]]),ncol=10)
      #start of training/testing
      #1) Fitting regression models on training dataset
      #2) applying models to testing dataset
      #3) calculate prediction metrics
      #4) calculate predicted scores
      
      count=1
      
      CV.RR.CT = glmnet::cv.glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 0,nfolds = 5)
      RR.mod=glmnet::glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 0, lambda = CV.RR.CT$lambda.min)
      predmetrics[count,2:4]=extractmetric.bysex(RR.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,count]=extractmetric.bysex(RR.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(RR.mod,CV.RR.CT)
      
      count=count+1
      
      CV.lasso.CT = glmnet::cv.glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 1,nfolds = 5)
      lasso.mod=glmnet::glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 1, lambda = CV.lasso.CT$lambda.min)
      predmetrics[count,2:4]=extractmetric.bysex(lasso.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,count]=extractmetric.bysex(lasso.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(lasso.mod,CV.lasso.CT)
      
      count=count+1
      
      CV.Enet.CT = glmnet::cv.glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 1,nfolds = 5)
      Enet.mod=glmnet::glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 0.5, lambda = CV.Enet.CT$lambda.min)
      predmetrics[count,2:4]=extractmetric.bysex(Enet.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,count]=extractmetric.bysex(Enet.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(Enet.mod,CV.Enet.CT)
      
      count=count+1
      
      pls.mod = pls::plsr(train_outcome.bysex[[sex]]~train_feat.bysex[[sex]],ncomp=20,segments=5, validation="CV")
      predmetrics[count,2:4]=extractmetric.bysex(pls.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,count]=extractmetric.bysex(pls.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(pls.mod)
      
      count=count+1
      
      gprlinear.mod=kernlab::gausspr(train_feat.bysex[[sex]],train_outcome.bysex[[sex]], kernel="vanilladot")
      predmetrics[count,2:4]=extractmetric.bysex(gprlinear.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,count]=extractmetric.bysex(gprlinear.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(gprlinear.mod)
      
      count=count+1

      svmlinear.mod=kernlab::ksvm(train_feat.bysex[[sex]],train_outcome.bysex[[sex]], kernel="vanilladot")
      predmetrics[count,2:4]=extractmetric.bysex(svmlinear.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,count]=extractmetric.bysex(svmlinear.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(svmlinear.mod)
      
      count=count+1

      kqrlinear.mod=kernlab::kqr(train_feat.bysex[[sex]],train_outcome.bysex[[sex]], kernel="vanilladot")
      predmetrics[count,2:4]=extractmetric.bysex(kqrlinear.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,count]=extractmetric.bysex(kqrlinear.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(kqrlinear.mod)
      
      count=count+1
      
      gprrbf.mod=kernlab::gausspr(train_feat.bysex[[sex]],train_outcome.bysex[[sex]], kernel="rbfdot")
      predmetrics[count,2:4]=extractmetric.bysex(gprrbf.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,count]=extractmetric.bysex(gprrbf.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(gprrbf.mod)
      
      count=count+1
      
      svmrbf.mod=kernlab::ksvm(train_feat.bysex[[sex]],train_outcome.bysex[[sex]], kernel="rbfdot")
      predmetrics[count,2:4]=extractmetric.bysex(svmrbf.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,count]=extractmetric.bysex(svmrbf.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(svmrbf.mod)
      
      count=count+1
      
      kqrrbf.mod=kernlab::kqr(train_feat.bysex[[sex]],train_outcome.bysex[[sex]], kernel="rbfdot")
      predmetrics[count,2:4]=extractmetric.bysex(kqrrbf.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,count]=extractmetric.bysex(kqrrbf.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(kqrrbf.mod)
      
      #formatting results matrix
      predmetrics=data.frame(predmetrics)
      colnames(predmetrics)=c("model","r","MAE","bias")
      predmetrics$r=as.numeric(predmetrics$r)
      predmetrics$MAE=as.numeric(predmetrics$MAE)
      
      return(list(predmetrics,predscores))
      closeAllConnections()  
    }
  
  ## XGB needs to be executed outside the foreach loops
 if(xgb==T)
  {
    source("https://github.com/CogBrainHealthLab/MLtools/blob/main/xgb.R?raw=TRUE")
    xgbresults=list()
    for (sex in 1:2)
    {
      #results matrix
      xgbpredmetrics=matrix(NA,nrow=1, ncol=4)
      xgbpredscores=matrix(NA,nrow=length(test_outcome.bysex[[sex]]),ncol=1)
      
      #training models
      xgb.mod=XGBlinear(train_feat.bysex[[sex]], train_outcome.bysex[[sex]])
      xgbpredmetrics[1,2:4]=extractmetric.bysex(xgb.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      xgbpredscores[,1]=extractmetric.bysex(xgb.mod,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(xgb.mod)
      
      gc()
      
      #formatting results matrix
      xgbpredmetrics=data.frame(xgbpredmetrics)
      colnames(xgbpredmetrics)=c("model","r","MAE","bias")
      xgbpredmetrics$r=as.numeric(xgbpredmetrics$r)
      xgbpredmetrics$MAE=as.numeric(xgbpredmetrics$MAE)
      
      xgbresults[[1+((sex-1)*2)]]=xgbpredmetrics
      xgbresults[[2+((sex-1)*2)]]=xgbpredscores
    }
    results[[1]]=rbind(results[[1]],xgbresults[[1]])
    results[[3]]=rbind(results[[3]],xgbresults[[3]])
    results[[2]]=cbind(results[[2]],xgbresults[[2]])
    results[[4]]=cbind(results[[4]],xgbresults[[4]])
  } 
  
  #recombine sex partitions  
  pred_outcome.recomb=rbind(results[[2]],results[[4]])
  test_outcome.recomb=c(test_outcome.bysex[[1]],test_outcome.bysex[[2]])
  
  #prediction metrics in recombined data
  predmetrics.recomb=matrix(NA,nrow=NCOL(results[[4]]), ncol=4)
  
  if(xgb==F)
  {
    predmetrics.recomb[,1]=c("RidgeR", "LassoR","ElasticNet","PLSR","GPR (Linear)","SVM (Linear)","KQR (Linear)", "GPR (RBF)", "SVM (RBF)","KQR (RBF)")
  } else
  {
    predmetrics.recomb[,1]=c("RidgeR", "LassoR","ElasticNet","PLSR","GPR (Linear)","SVM (Linear)","KQR (Linear)", "GPR (RBF)", "SVM (RBF)","KQR (RBF)", "XGB (linear)")
  }
  
  predmetrics.recomb=data.frame(predmetrics.recomb)
  colnames(predmetrics.recomb)=c("model","r","MAE","bias")
  predmetrics.recomb[,2]=cor(pred_outcome.recomb,test_outcome.recomb)
  predmetrics.recomb[,3]=colMeans(abs(pred_outcome.recomb-test_outcome.recomb))
  predmetrics.recomb[,4]=cor((pred_outcome.recomb-test_outcome.recomb),test_outcome.recomb)
  
  max.idx=which(as.numeric(predmetrics.recomb$r)==max(as.numeric(predmetrics.recomb$r),na.rm = T))
  min.idx=which(as.numeric(predmetrics.recomb$MAE)==min(as.numeric(predmetrics.recomb$MAE),na.rm = T))
  
  cat(paste("\nModel with highest r: ",predmetrics.recomb$model[max.idx],"; r=",round(max(as.numeric(predmetrics.recomb$r),na.rm=T),3),"\n",sep=""))
  cat(paste("Model with lowest MAE: ",predmetrics.recomb$model[min.idx],"; MAE=",round(min(as.numeric(predmetrics.recomb$MAE),na.rm=T),3),sep=""))
  
  pred_outcome.recomb.ordered=pred_outcome.recomb[order(c(test.M.idx,test.F.idx)),]
  
  ##bias correction
  pred_outcome.recomb.ordered.biascorrected=pred_outcome.recomb.ordered
  for(model in 1:NCOL(pred_outcome.recomb.ordered))
  {
    agegap=pred_outcome.recomb.ordered[,model]-test_outcome.recomb
    biasmod=lm(agegap~test_outcome.recomb)
    offset=(coef(biasmod)[2]*test_outcome.recomb)+coef(biasmod)[1]
    pred_outcome.recomb.ordered.biascorrected[,model]=pred_outcome.recomb.ordered[,model]-offset
    remove(agegap, biasmod, offset)
  }
  
  returnobj=list(results[[1]],results[[3]],predmetrics.recomb,pred_outcome.recomb.ordered,pred_outcome.recomb.ordered.biascorrected)
  names(returnobj)=c("predmetrics.M","predmetrics.F","predmetrics.all","predscores","predscores.biascorrected")
  return(returnobj)
}

plot.metrics=function(results)
{
  results$modelno=1:NROW(results)
  a=ggplot2::ggplot(results,ggplot2::aes(x=modelno,y=as.numeric(r), group=1))+
    ggplot2::geom_point()+
    ggplot2::geom_line()+
    ggplot2::scale_x_continuous(breaks=1:NROW(results))+
    ggplot2::labs(x=NULL, y="r")+
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
  b=ggplot2::ggplot(results,ggplot2::aes(x=modelno,y=as.numeric(MAE), group=1))+
    ggplot2::geom_point()+
    ggplot2::geom_line()+
    ggplot2::scale_x_continuous(breaks=1:NROW(results))+
    ggplot2::labs(x=NULL, y="MAE")+
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
  c=ggplot2::ggplot(results,ggplot2::aes(x=modelno,y=as.numeric(bias), group=1))+
    ggplot2::geom_point()+
    ggplot2::geom_line()+
    ggplot2::scale_x_continuous(breaks=1:NROW(results),labels=results$model)+
    ggplot2::labs(x=NULL, y="bias")+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1),plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
  return(cowplot::plot_grid(a,b,c,nrow = 3,rel_heights = c(0.3,0.3,0.45),align="v",axis="lr"))
}
