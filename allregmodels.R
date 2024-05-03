## 11 regression-based ML models in a single function
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

##################################################################################################################
##################################################################################################################
##output prediction metrics
extractmetric=function(model,test_feat, test_outcome)
{
  #PLSR model has to be treated differently
  if(class(model)[1]=="mvr")
  {
    pred_outcome=predict(model,test_feat, ncomp=pls::selectNcomp(model,  method = "randomization"))
    return(list(c(cor(pred_outcome,test_outcome),
                 mean(abs(pred_outcome-test_outcome)),
                 cor((pred_outcome-test_outcome),test_outcome)),pred_outcome))
  } else
  {
    #all other regression models
    pred_outcome=kernlab::predict(model,test_feat)
    return(list(c(cor(pred_outcome,test_outcome),
             mean(abs(pred_outcome-test_outcome)),
             cor((pred_outcome-test_outcome),test_outcome)),pred_outcome))
  }
}

##Runs 11 regression models
pred.allmodels=function(train_outcome, train_feat,test_outcome, test_feat, xgb=F)
{
  #check if train_feat contains columns of 0s, if so, these columns are removed
  col0_idx=which(colSums(train_feat)==0)
  if(length(col0_idx)>1)
  {
    train_feat=train_feat[,-col0_idx]
    test_feat=test_feat[,-col0_idx]
  }
  
  #setting up results matrix
  predmetrics=matrix(NA,nrow=13, ncol=4)
  predmetrics[,1]=c("RidgeR", "LassoR","PLSR","GPR (Linear)","SVM (Linear)", "RVM (Linear)","KQR (Linear)", "GPR (RBF)", "SVM (RBF)", "RVM (RBF)", "KQR (RBF)","XGB (linear)", "XGB (tree)")
  predscores=matrix(NA,nrow=length(test_outcome,ncol=13)
  #start of training/testing
  #1) Fitting regression models on training dataset
  #2) applying models to testing dataset
  #3) calculate prediction metrics
  
  CV.RR.CT = glmnet::cv.glmnet(train_feat, train_outcome, alpha = 0,nfolds = 5,parallel = T)
  model1=glmnet::glmnet(train_feat, train_outcome, alpha = 0, lambda = CV.RR.CT$lambda.1se)
  predmetrics[1,2:4]=extractmetric(model1,test_feat,test_outcome)[[1]]
  predscores[,1]=extractmetric(model1,test_feat,test_outcome)[[2]]
  remove(model1,CV.RR.CT)
  
  cat("1 Ridge regression completed\n")
  
  CV.RR.CT = glmnet::cv.glmnet(train_feat, train_outcome, alpha = 1,nfolds = 5,parallel = T)
  model2=glmnet::glmnet(train_feat, train_outcome, alpha = 1, lambda = CV.RR.CT$lambda.1se)
  predmetrics[2,2:4]=extractmetric(model2,test_feat,test_outcome)[[1]]
  predscores[,2]=extractmetric(model2,test_feat,test_outcome)[[2]]
  remove(model2,CV.RR.CT)
  
  cat("2 Lasso regression completed\n")
  
  model3 = pls::plsr(train_outcome~train_feat,ncomp=20,segments=5, validation="CV",)
  predmetrics[3,2:4]=extractmetric(model3,test_feat,test_outcome)[[1]]
  predscores[,3]=extractmetric(model3,test_feat,test_outcome)[[2]]
  remove(model3)
  
  cat("3 Partial least squares regression completed\n")
  
  model4=kernlab::gausspr(x=train_feat, y=train_outcome, kernel="vanilladot")
  predmetrics[4,2:4]=extractmetric(model4,test_feat,test_outcome)[[1]]
  predscores[,4]=extractmetric(model4,test_feat,test_outcome)[[2]]
  remove(model4)
  
  cat("4 Gaussian process regression completed\n")
  
  model5=kernlab::ksvm(x=train_feat, y=train_outcome, kernel="vanilladot")
  predmetrics[5,2:4]=extractmetric(model5,test_feat,test_outcome)[[1]]
  predscores[,5]=extractmetric(model5,test_feat,test_outcome)[[2]]
  remove(model5)
  
  cat("5 Support vector machine (Vanilla) completed\n")
  
  model6=kernlab::rvm(x=train_feat, y=train_outcome, kernel="vanilladot")
  predmetrics[6,2:4]=extractmetric(model6,test_feat,test_outcome)[[1]]
  predscores[,6]=extractmetric(model6,test_feat,test_outcome)[[2]]
  remove(model6)
  
  cat("6 Relevance vector machine (Vanilla) completed\n")
  
  model7=kernlab::kqr(x=train_feat, y=train_outcome, kernel="vanilladot")
  predmetrics[7,2:4]=extractmetric(model7,test_feat,test_outcome)[[1]]
  predscores[,7]=extractmetric(model7,test_feat,test_outcome)[[2]]
  remove(model7)
  
  cat("7 Kernel quantile regression (vanilla) completed\n")
  
  model8=kernlab::gausspr(x=train_feat, y=as.numeric(train_outcome), kernel="rbfdot")
  predmetrics[8,2:4]=extractmetric(model8,test_feat,test_outcome)[[1]]
  predscores[,8]=extractmetric(model8,test_feat,test_outcome)[[2]]
  remove(model8)
  
  cat("8 Gaussian process regression (RBF) completed\n")
  
  model9=kernlab::ksvm(x=train_feat, y=train_outcome, kernel="rbfdot")
  predmetrics[9,2:4]=extractmetric(model9,test_feat,test_outcome)[[1]]
  predscores[,9]=extractmetric(model9,test_feat,test_outcome)[[2]]
  remove(model9)
  
  cat("9 Support vector machine (RBF) completed\n")
  
  model10=kernlab::rvm(x=train_feat, y=train_outcome, kernel="rbfdot")
  predmetrics[10,2:4]=extractmetric(model10,test_feat,test_outcome)[[1]]
  predscores[,10]=extractmetric(model10,test_feat,test_outcome)[[2]]
  remove(model10)
  
  cat("10 Relevance vector machine (RBF) completed\n")
  
  model11=kernlab::kqr(x=train_feat, y=train_outcome, kernel="rbfdot")
  predmetrics[11,2:4]=extractmetric(model11,test_feat,test_outcome)[[1]]
  predscores[,11]=extractmetric(model11,test_feat,test_outcome)[[2]]
  remove(model11)
  
  cat("11 Kernel quantile regression (RBF) completed\n")
  
  #optional XGB models
  if(xgb==T)
  {
    source("https://github.com/CogBrainHealthLab/MLtools/blob/main/xgb.R?raw=TRUE")
    model12=XGBlinear(train_feat, train_outcome)
    predmetrics[12,2:4]=extractmetric(model12,test_feat,test_outcome)[[1]]
    predscores[,12]=extractmetric(model12,test_feat,test_outcome)[[2]]
    remove(model12)
    
    cat("12 XGBlinear completed\n")
    
    model13=XGBtree(train_feat, train_outcome)
    predmetrics[13,2:4]=extractmetric(model13,test_feat,test_outcome)[[1]]
    predscores[,13]=extractmetric(model13,test_feat,test_outcome)[[2]]
    remove(model13)
    
    cat("13 XGBtree completed\n")
  } else
  {
    predmetrics=predmetrics[1:11,]
  }
  #formatting results matrix
  predmetrics=data.frame(predmetrics)
  colnames(predmetrics)=c("model","r","MAE","bias")
  predmetrics$r=as.numeric(predmetrics$r)
  predmetrics$MAE=as.numeric(predmetrics$MAE)
  
  cat(paste("\nModel with highest r: ",predmetrics$model[which.max(predmetrics$r)],"; r=",round(max(predmetrics$r),3),"\n",sep=""))
  cat(paste("Model with lowest MAE: ",predmetrics$model[which.min(predmetrics$MAE)],"; MAE=",round(min(predmetrics$MAE),3),sep=""))

  predscores=data.frame(predscores)
  colnames(predscores)="RidgeR", "LassoR","PLSR","GPR (Linear)","SVM (Linear)", "RVM (Linear)","KQR (Linear)", "GPR (RBF)", "SVM (RBF)", "RVM (RBF)", "KQR (RBF)","XGB (linear)", "XGB (tree)")
                    
  return(list(predmetrics,predscores)
}

## plot out results using ggplot
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
##################################################################################################################
##################################################################################################################
## EXAMPLE:

#source("https://github.com/CogBrainHealthLab/MLtools/blob/main/allregmodels.R?raw=TRUE")
#pred.allmodels(train_outcome = HCP_age,train_feat =HCP_dat,test_outcome = CC_age,test_feat = CC_dat)
