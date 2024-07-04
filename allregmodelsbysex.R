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
pred.allmodels.bysex=function(train_outcome, train_feat,train_sex,test_outcome, test_feat,test_sex, xgb=F, harm=1, eb=F, train_site)
{
  
  #check if train_feat contains columns of 0s, if so, these columns are removed
  col0_idx=which(colSums(train_feat)==0)
  if(length(col0_idx)>1)
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
  
  test_outcome.bysex=list(test_outcome[test.M.idx],test_outcome[test.F.idx])
  test_feat.bysex=list(test_feat[test.M.idx,],test_feat[test.F.idx,])
  
  remove(test_outcome,train_outcome,test_feat,train_feat,train.M.idx,train.F.idx)
  
  ##harmonize different sexes separately

  for (sex in 1:2)
  {
    dat.all=rbind(data.matrix(train_feat.bysex[[sex]]),data.matrix(test_feat.bysex[[sex]]))
    if(harm==1) 
    {
      dat.harmonized =neuroCombat::neuroCombat(dat=t(dat.all), eb=eb,
                                               batch=c(train_site.bysex[[sex]],rep("test",length(test_outcome.bysex[[sex]]))),
                                               mod=c(train_outcome.bysex[[sex]],test_outcome.bysex[[sex]]))  
      
      train_feat.bysex[[sex]]=t(dat.harmonized$dat.combat)[1:length(train_outcome.bysex[[sex]]),]
      test_feat.bysex[[sex]]=t(dat.harmonized$dat.combat)[(length(train_outcome.bysex[[sex]])+1):(length(train_outcome.bysex[[sex]])+length(test_outcome.bysex[[sex]])),]  
    }
    if(harm==2) 
    {
      dat.harmonized =CovBat::covbat(dat=t(dat.all), eb=eb,
                                     bat=c(rep("train",length(train_outcome.bysex[[sex]])),rep("test",length(test_outcome.bysex[[sex]]))),
                                     mod=c(train_outcome.bysex[[sex]],test_outcome.bysex[[sex]]))  
      
      train_feat.bysex[[sex]]=t(dat.harmonized$dat.covbat)[1:length(train_outcome.bysex[[sex]]),]
      test_feat.bysex[[sex]]=t(dat.harmonized$dat.covbat)[(length(train_outcome.bysex[[sex]])+1):(length(train_outcome.bysex[[sex]])+length(test_outcome.bysex[[sex]])),]  
    }
    remove(dat.harmonized)
  }
 
  cat("completed harmonization\n")
  #activate parallel processing
  unregister_dopar = function() {
    env = foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }
  unregister_dopar()
  
  cl=parallel::makeCluster(2)
  doParallel::registerDoParallel(2)
  `%dopar%` = foreach::`%dopar%`
  
  results=foreach::foreach(sex=1:2, .combine="c",.packages = c("glmnet","pls","kernlab"), .export ="extractmetric.bysex")  %dopar%
    {
      #setting up results matrix
      predmetrics=matrix(NA,nrow=11, ncol=4)
      predmetrics[,1]=c("RidgeR", "LassoR","PLSR","GPR (Linear)","SVM (Linear)", "RVM (Linear)","KQR (Linear)", "GPR (RBF)", "SVM (RBF)", "RVM (RBF)", "KQR (RBF)")
      
      predscores=matrix(NA,nrow=length(test_outcome.bysex[[sex]]),ncol=11)
      #start of training/testing
      #1) Fitting regression models on training dataset
      #2) applying models to testing dataset
      #3) calculate prediction metrics
      #4) calculate predicted scores
      
      set.seed(123)
      CV.RR.CT = glmnet::cv.glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 0,nfolds = 5)
      model1=glmnet::glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 0, lambda = CV.RR.CT$lambda.1se)
      predmetrics[1,2:4]=extractmetric.bysex(model1,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,1]=extractmetric.bysex(model1,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model1,CV.RR.CT)
      
      CV.RR.CT = glmnet::cv.glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 1,nfolds = 5)
      model2=glmnet::glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 1, lambda = CV.RR.CT$lambda.1se)
      predmetrics[2,2:4]=extractmetric.bysex(model2,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,2]=extractmetric.bysex(model2,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model2,CV.RR.CT)
      
      model3 = pls::plsr(train_outcome.bysex[[sex]]~train_feat.bysex[[sex]],ncomp=20,segments=5, validation="CV")
      predmetrics[3,2:4]=extractmetric.bysex(model3,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,3]=extractmetric.bysex(model3,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model3)
      
      model4=kernlab::gausspr(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="vanilladot")
      predmetrics[4,2:4]=extractmetric.bysex(model4,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,4]=extractmetric.bysex(model4,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model4)
      
      model5=kernlab::ksvm(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="vanilladot")
      predmetrics[5,2:4]=extractmetric.bysex(model5,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,5]=extractmetric.bysex(model5,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model5)
      
      model6=kernlab::rvm(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="vanilladot")
      predmetrics[6,2:4]=extractmetric.bysex(model6,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,6]=extractmetric.bysex(model6,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model6)
      
      model7=kernlab::kqr(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="vanilladot")
      predmetrics[7,2:4]=extractmetric.bysex(model7,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,7]=extractmetric.bysex(model7,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model7)
      
      model8=kernlab::gausspr(x=train_feat.bysex[[sex]], y=as.numeric(train_outcome.bysex[[sex]]), kernel="rbfdot")
      predmetrics[8,2:4]=extractmetric.bysex(model8,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,8]=extractmetric.bysex(model8,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model8)
      
      model9=kernlab::ksvm(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="rbfdot")
      predmetrics[9,2:4]=extractmetric.bysex(model9,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,9]=extractmetric.bysex(model9,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model9)
      
      model10=kernlab::rvm(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="rbfdot")
      predmetrics[10,2:4]=extractmetric.bysex(model10,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,10]=extractmetric.bysex(model10,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model10)
      
      model11=kernlab::kqr(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="rbfdot")
      predmetrics[11,2:4]=extractmetric.bysex(model11,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      predscores[,11]=extractmetric.bysex(model11,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model11)
      
      #formatting results matrix
      predmetrics=data.frame(predmetrics)
      colnames(predmetrics)=c("model","r","MAE","bias")
      predmetrics$r=as.numeric(predmetrics$r)
      predmetrics$MAE=as.numeric(predmetrics$MAE)
      
      return(list(predmetrics,predscores))
    }
  
  ## XGB needs to be executed outside the foreach loops
  if(xgb==T)
  {
    source("https://github.com/CogBrainHealthLab/MLtools/blob/main/xgb.R?raw=TRUE")
    xgbresults=list()
    for (sex in 1:2)
    {
      #results matrix
      xgbpredmetrics=matrix(NA,nrow=2, ncol=4)
      xgbpredscores=matrix(NA,nrow=length(test_outcome.bysex[[sex]]),ncol=2)
      
      #training models
      model12=XGBlinear(train_feat.bysex[[sex]], train_outcome.bysex[[sex]])
      xgbpredmetrics[1,2:4]=extractmetric.bysex(model12,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      xgbpredscores[,1]=extractmetric.bysex(model12,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model12)
      
      model13=XGBtree(train_feat.bysex[[sex]], train_outcome.bysex[[sex]])
      xgbpredmetrics[2,2:4]=extractmetric.bysex(model13,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[1]]
      xgbpredscores[,2]=extractmetric.bysex(model13,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])[[2]]
      remove(model13)
      
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
    predmetrics.recomb[,1]=c("RidgeR", "LassoR","PLSR","GPR (Linear)","SVM (Linear)", "RVM (Linear)","KQR (Linear)", "GPR (RBF)", "SVM (RBF)", "RVM (RBF)", "KQR (RBF)")  
  } else
  {
    predmetrics.recomb[,1]=c("RidgeR", "LassoR","PLSR","GPR (Linear)","SVM (Linear)", "RVM (Linear)","KQR (Linear)", "GPR (RBF)", "SVM (RBF)", "RVM (RBF)", "KQR (RBF)", "XGB (linear)","XGB (Tree)")  
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
  
  returnobj=list(results[[1]],results[[3]],predmetrics.recomb,pred_outcome.recomb.ordered)
  names(returnobj)=c("predmetrics.M","predmetrics.F","predmetrics.all","predscores")
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
