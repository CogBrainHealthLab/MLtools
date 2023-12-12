## WRAPPERS FOR BAYESIAN OPTIMIZED XGBOOSTING
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

##################################################################################################################
##################################################################################################################
#XGBlinear
XGBlinear=function(train_dat,outcome, alpharange=c(0,1), lambdarange=c(1,5),etarange=c(0.001, 0.2),nthread=4, nround=500)
{

  ## checked required packages
  list.of.packages = c("doParallel","parallel","ParBayesianOptimization","xgboost","caret")
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) 
  {
    cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
    install.packages(new.packages)
  }
  doParallel::registerDoParallel(nthread)
  folds=caret::createFolds(outcome,k=5)
  bounds = list(lambda = lambdarange,alpha = alpharange,eta = etarange)
  dtrain = xgboost::xgb.DMatrix(data = train_dat, label = outcome)
  cl = parallel::makeCluster(nthread)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,c('train_dat','outcome','folds'),envir = environment())

  obj_func = function(eta,lambda,alpha,nround) 
  {
    param = list(
      eta=eta,
      lambda = lambda,
      alpha = alpha,
      booster = "gblinear",
      objective = "reg:squarederror",
      eval_metric = "mae")
    
    xgbcv = xgboost::xgb.cv(params = param,
                   data = xgboost::xgb.DMatrix(data = train_dat, label = outcome),
                   nround = 500,
                   folds = folds,
                   prediction = TRUE,
                   early_stopping_rounds = 5,
                   verbose = 0,
                   maximize = F)
    
    lst = list(Score = -min(xgbcv$evaluation_log$test_mae_mean),nrounds = xgbcv$best_iteration)
    return(lst)
  }
  
  bayes_out = ParBayesianOptimization::bayesOpt(FUN = obj_func, bounds = bounds,initPoints = length(bounds) + 2,parallel = T,iters.n = nthread,iters.k=nthread)
  
  opt_params = append(list(booster = "gblinear",objective = "reg:squarederror",eval_metric = "mae"),ParBayesianOptimization::getBestPars(bayes_out))
  
  xgbcv = xgboost::xgb.cv(params = opt_params,data=xgboost::xgb.DMatrix(data = train_dat, label = outcome),nround = nround,folds = folds,prediction = TRUE,early_stopping_rounds = 5,verbose = 0,maximize = F)
  
  #final XGB model with optimal hyperparameters and nround
  xg_mod = xgboost::xgboost(data = xgboost::xgb.DMatrix(data = train_dat, label = outcome), params = opt_params, nround = xgbcv$best_iteration, verbose = F)
  return(xg_mod)
}
##################################################################################################################
##################################################################################################################
#XGBtree
XGBtree=function(train_dat,outcome, alpharange=c(0,1), lambdarange=c(1,5),etarange=c(0.001, 0.2),max_depth = c(1L, 10L),min_child_weight = c(1, 50),subsample = c(0.1, 1),nthread=4,nround=500)
{
  list.of.packages = c("doParallel", "ParBayesianOptimization","xgboost")
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) 
  {
    cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
    install.packages(new.packages)
  }
  doParallel::registerDoParallel(nthread)
  folds=caret::createFolds(outcome,k=5)
  bounds = list(lambda = lambdarange,alpha = alpharange,eta = etarange, max_depth=max_depth, min_child_weight=min_child_weight,subsample=subsample)
  
  dtrain = xgboost::xgb.DMatrix(data = train_dat, label = outcome)
  cl = parallel::makeCluster(nthread)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,c('train_dat','outcome','folds'),envir = environment())

  #optimizer
  obj_func = function(eta,lambda,alpha,nround,max_depth,min_child_weight,subsample) 
  {
    param = list(
      eta=eta,
      lambda = lambda,
      alpha = alpha,
      max_depth=max_depth,
      min_child_weight=min_child_weight,
      subsample=subsample,
      booster = "gbtree",
      objective = "reg:squarederror",
      eval_metric = "mae")
    
    xgbcv =xgboost::xgb.cv(params = param,
                   data = xgboost::xgb.DMatrix(data = train_dat, label = outcome),
                   nround = 500,
                   folds = folds,
                   prediction = TRUE,
                   early_stopping_rounds = 5,
                   verbose = 0,
                   maximize = F)
    
    lst = list(Score = -min(xgbcv$evaluation_log$test_mae_mean),nrounds = xgbcv$best_iteration)
    return(lst)
  }
  
  bayes_out = ParBayesianOptimization::bayesOpt(FUN = obj_func, bounds = bounds,initPoints = length(bounds) + 2,parallel = T,iters.n = nthread,iters.k=nthread)
  
  opt_params = append(list(booster = "gbtree",objective = "reg:squarederror",eval_metric = "mae"), ParBayesianOptimization::getBestPars(bayes_out))
  
  xgbcv = xgboost::xgb.cv(params = opt_params,data=xgboost::xgb.DMatrix(data = train_dat, label = outcome),nround = nround,folds = folds,prediction = TRUE,early_stopping_rounds = 5,verbose = 0,maximize = F)
  
  #final XGB model with optimal hyperparameters and nround
  xg_mod = xgboost::xgboost(data = xgboost::xgb.DMatrix(data = train_dat, label = outcome), params = opt_params, nround = xgbcv$best_iteration, verbose = F)
  return(xg_mod)
}

