extractmetric.bysex=function(model,test_feat, test_outcome)
{
  #PLSR model has to be treated differently
  if(class(model)[1]=="mvr")
  {
    pred_outcome=predict(model,test_feat, ncomp=pls::selectNcomp(model,  method = "randomization"))
    return(list(c(cor(pred_outcome,test_outcome),mean(abs(pred_outcome-test_outcome)),cor((pred_outcome-test_outcome),test_outcome)),pred_outcome))
  } else
  {
    #all other regression models
    pred_outcome=kernlab::predict(model,test_feat)
    return(list(c(cor(pred_outcome,test_outcome),mean(abs(pred_outcome-test_outcome)),cor((pred_outcome-test_outcome),test_outcome)),pred_outcome))
  }
}

##Runs 11 regression models
pred.allmodels.bysex=function(train_outcome, train_feat,train_sex,test_outcome, test_feat,test_sex, xgb=F)
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
  
  test_outcome.bysex=list(test_outcome[test.M.idx],test_outcome[test.F.idx])
  test_feat.bysex=list(test_feat[test.M.idx,],test_feat[test.F.idx,])
  
  remove(test_outcome,train_outcome,test_feat,train_feat,train.M.idx,train.F.idx)
  
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
      predmetrics=matrix(NA,nrow=13, ncol=4)
      predmetrics[,1]=c("RidgeR", "LassoR","PLSR","GPR (Linear)","SVM (Linear)", "RVM (Linear)","KQR (Linear)", "GPR (RBF)", "SVM (RBF)", "RVM (RBF)", "KQR (RBF)","XGB (linear)", "XGB (tree)")
      
      predscores=matrix(NA,nrow=length(test_outcome.bysex[[sex]]),ncol=11)
      #start of training/testing
      #1) Fitting regression models on training dataset
      #2) applying models to testing dataset
      #3) calculate prediction metrics
      
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
      
      model3 = pls::plsr(train_outcome.bysex[[sex]]~train_feat.bysex[[sex]],ncomp=20,segments=5, validation="CV",)
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
      
      #optional XGB models
      if(xgb==T)
      {
        source("https://github.com/CogBrainHealthLab/MLtools/blob/main/xgb.R?raw=TRUE")
        model12=XGBlinear(train_feat.bysex[[sex]], train_outcome.bysex[[sex]])
        predmetrics[12,2:4]=extractmetric.bysex(model12,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
        remove(model12)
        
        model13=XGBtree(train_feat.bysex[[sex]], train_outcome.bysex[[sex]])
        predmetrics[13,2:4]=extractmetric.bysex(model13,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
        remove(model13)
        
      } else
      {
        predmetrics=predmetrics[1:11,]
      }
      #formatting results matrix
      predmetrics=data.frame(predmetrics)
      colnames(predmetrics)=c("model","r","MAE","bias")
      predmetrics$r=as.numeric(predmetrics$r)
      predmetrics$MAE=as.numeric(predmetrics$MAE)
      
      return(list(predmetrics,predscores))
    }
  
    #recombine sex partitions  
    pred_outcome.recomb=rbind(results[[2]],results[[4]])
    test_outcome.recomb=c(test_outcome.bysex[[1]],test_outcome.bysex[[2]])

    #prediction metrics in recombined data
    predmetrics.recomb=matrix(NA,nrow=13, ncol=4)
    predmetrics.recomb[,1]=c("RidgeR", "LassoR","PLSR","GPR (Linear)","SVM (Linear)", "RVM (Linear)","KQR (Linear)", "GPR (RBF)", "SVM (RBF)", "RVM (RBF)", "KQR (RBF)","XGB (linear)", "XGB (tree)")
    predmetrics.recomb=data.frame(predmetrics.recomb)
    colnames(predmetrics.recomb)=c("model","r","MAE","bias")
    predmetrics.recomb[1:11,2]=cor(pred_outcome.recomb,test_outcome.recomb)
    predmetrics.recomb[1:11,3]=colMeans(abs(pred_outcome.recomb-test_outcome.recomb))
    predmetrics.recomb[1:11,4]=cor((pred_outcome.recomb-test_outcome.recomb),test_outcome.recomb)
    
    colnames(pred_outcome.recomb)=c("RidgeR", "LassoR","PLSR","GPR (Linear)","SVM (Linear)", "RVM (Linear)","KQR (Linear)", "GPR (RBF)", "SVM (RBF)", "RVM (RBF)","KQR (RBF)")   
    if(xgb==F)
    {
      predmetrics.recomb=predmetrics.recomb[1:11,]
      pred_outcome.recomb=pred_outcome.recomb[,1:11]
    }
    cat(paste("\nModel with highest r: ",predmetrics.recomb$model[which.max(as.numeric(predmetrics.recomb$r))],"; r=",round(max(as.numeric(predmetrics.recomb$r)),3),"\n",sep=""))
    cat(paste("Model with lowest MAE: ",predmetrics.recomb$model[which.min(as.numeric(predmetrics.recomb$MAE))],"; MAE=",round(min(as.numeric(predmetrics.recomb$MAE)),3),sep=""))
    
    
    return(list(results[[1]],results[[3]],predmetrics.recomb,pred_outcome.recomb,c(test.M.idx,test.F.idx)))
  }
  
