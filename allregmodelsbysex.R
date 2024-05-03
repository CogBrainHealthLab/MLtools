CV.RR.CT = glmnet::cv.glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 0,nfolds = 5)
model1=glmnet::glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 0, lambda = CV.RR.CT$lambda.1se)
predmetrics[1,2:4]=extractmetric(model1,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
remove(model1,CV.RR.CT)

CV.RR.CT = glmnet::cv.glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 1,nfolds = 5)
model2=glmnet::glmnet(train_feat.bysex[[sex]], train_outcome.bysex[[sex]], alpha = 1, lambda = CV.RR.CT$lambda.1se)
predmetrics[2,2:4]=extractmetric(model2,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
remove(model2,CV.RR.CT)

model3 = pls::plsr(train_outcome.bysex[[sex]]~train_feat.bysex[[sex]],ncomp=20,segments=5, validation="CV",)
predmetrics[3,2:4]=extractmetric(model3,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
remove(model3)

model4=kernlab::gausspr(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="vanilladot")
predmetrics[4,2:4]=extractmetric(model4,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
remove(model4)

model5=kernlab::ksvm(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="vanilladot")
predmetrics[5,2:4]=extractmetric(model5,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
remove(model5)

model6=kernlab::rvm(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="vanilladot")
predmetrics[6,2:4]=extractmetric(model6,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
remove(model6)

model7=kernlab::kqr(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="vanilladot")
predmetrics[7,2:4]=extractmetric(model7,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
remove(model7)

model8=kernlab::gausspr(x=train_feat.bysex[[sex]], y=as.numeric(train_outcome.bysex[[sex]]), kernel="rbfdot")
predmetrics[8,2:4]=extractmetric(model8,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
remove(model8)

model9=kernlab::ksvm(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="rbfdot")
predmetrics[9,2:4]=extractmetric(model9,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
remove(model9)

model10=kernlab::rvm(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="rbfdot")
predmetrics[10,2:4]=extractmetric(model10,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
remove(model10)

model11=kernlab::kqr(x=train_feat.bysex[[sex]], y=train_outcome.bysex[[sex]], kernel="rbfdot")
predmetrics[11,2:4]=extractmetric(model11,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
remove(model11)

#optional XGB models
if(xgb==T)
{
  source("https://github.com/CogBrainHealthLab/MLtools/blob/main/xgb.R?raw=TRUE")
  model12=XGBlinear(train_feat.bysex[[sex]], train_outcome.bysex[[sex]])
  predmetrics[12,2:4]=extractmetric(model12,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
  remove(model12)
  
  model13=XGBtree(train_feat.bysex[[sex]], train_outcome.bysex[[sex]])
  predmetrics[13,2:4]=extractmetric(model13,test_feat.bysex[[sex]],test_outcome.bysex[[sex]])
  remove(model13)
  
} else
{
  predmetrics=predmetrics[1:11,]
}
