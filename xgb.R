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

  folds=caret::createFolds(outcome,k=5)
  bounds = list(lambda = lambdarange,
                alpha = alpharange,
                eta = etarange)
  
  dtrain = xgboost::xgb.DMatrix(data = train_dat, label = outcome)
  xgboost::xgb.DMatrix.save(dtrain, "train.buffer")

  cl <- parallel::makeCluster(nthread)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(
    cl,
    c("dtrain", "folds"),
    envir = environment()
  )
  
  obj_func = function(eta,lambda,alpha,nround) 
  {
    param = list(
      learning_rate=eta,
      reg_lambda = lambda,
      reg_alpha = alpha,
      booster = "gblinear",
      objective = "reg:squarederror",
      eval_metric = "mae")
    
    xgbcv = xgboost::xgb.cv(params = param,
                            data =  xgboost::xgb.DMatrix("train.buffer"),
                            nround = 500,
                            folds = folds,
                            prediction = TRUE,
                            early_stopping_rounds = 5,
                            verbose = 0,
                            maximize = F)
    
    best_row <- which.min(xgbcv$evaluation_log$test_mae_mean)
    actual_best_iter <- xgbcv$evaluation_log$iter[best_row]
    
    lst = list(
      Score = -min(xgbcv$evaluation_log$test_mae_mean),
      nrounds = as.integer(actual_best_iter) # Guaranteed to be a length-1 integer
    )
    return(lst)
  }
  
  bayes_out = ParBayesianOptimization::bayesOpt(FUN = obj_func, 
                                                bounds = bounds,
                                                initPoints = length(bounds) + 2,
                                                parallel = T,
                                                iters.n = nthread,
                                                iters.k=nthread)
  
  opt_params = append(list(booster = "gblinear",
                           objective = "reg:squarederror",
                           eval_metric = "mae",
                           nthread=nthread),
                      ParBayesianOptimization::getBestPars(bayes_out))
  
  xgbcv = xgboost::xgb.cv(params = opt_params,
                          data= xgboost::xgb.DMatrix("train.buffer"),
                          nround = nround,
                          folds = folds,
                          prediction = TRUE,
                          early_stopping_rounds = 5,
                          verbose = 0,maximize = F)
  
  #final XGB model with optimal hyperparameters and nround
  
  xg_mod = xgboost::xgboost(y = outcome, x=train_dat, 
                            booster = opt_params$booster,
                            objective = opt_params$objective,
                            eval_metric=opt_params$eval_metric,
                            nthread=opt_params$nthread,
                            reg_lambda=opt_params$lambda, 
                            reg_alpha=opt_params$alpha, 
                            learning_rate=opt_params$eta,
                            nround = xgbcv$evaluation_log$iter[which.min(xgbcv$evaluation_log$test_mae_mean)])
  unlink("train.buffer")
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
  
  folds=caret::createFolds(outcome,k=5)
  bounds = list(lambda = lambdarange,
                alpha = alpharange,
                eta = etarange, 
                max_depth=max_depth,
                min_child_weight=min_child_weight,
                subsample=subsample)
  
  dtrain = xgboost::xgb.DMatrix(data = train_dat, label = outcome)
  xgboost::xgb.DMatrix.save(dtrain, "train.buffer")
  
  
  cl <- parallel::makeCluster(nthread)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(
    cl,
    c("dtrain", "folds"),
    envir = environment()
  )
  
  #optimizer
  obj_func = function(eta,lambda,alpha,nround,max_depth,min_child_weight,subsample) 
  {
    param = list(
      learning_rate=eta,
      reg_lambda = lambda,
      reg_alpha = alpha,
      max_depth=max_depth,
      min_child_weight=min_child_weight,
      subsample=subsample,
      booster = "gbtree",
      objective = "reg:squarederror",
      eval_metric = "mae")
    
    xgbcv =xgboost::xgb.cv(params = param,
                           data = xgboost::xgb.DMatrix("train.buffer"),
                           nround = 500,
                           folds = folds,
                           prediction = TRUE,
                           early_stopping_rounds = 5,
                           verbose = 0,
                           maximize = F)
    
    best_row <- which.min(xgbcv$evaluation_log$test_mae_mean)
    actual_best_iter <- xgbcv$evaluation_log$iter[best_row]
    
    lst = list(
      Score = -min(xgbcv$evaluation_log$test_mae_mean),
      nrounds = as.integer(actual_best_iter) # Guaranteed to be a length-1 integer
    )
    return(lst)
  }
  
  bayes_out = ParBayesianOptimization::bayesOpt(FUN = obj_func, 
                                                bounds = bounds,
                                                initPoints = length(bounds) + 2,
                                                parallel = T,
                                                iters.n = nthread,
                                                iters.k=nthread)
  
  opt_params = append(list(booster = "gbtree",
                           objective = "reg:squarederror",
                           eval_metric = "mae",
                           nthread=nthread), 
                      ParBayesianOptimization::getBestPars(bayes_out))
  
  xgbcv = xgboost::xgb.cv(params = opt_params,
                          data=xgboost::xgb.DMatrix("train.buffer"),
                          nround = nround,
                          folds = folds,
                          prediction = TRUE,
                          early_stopping_rounds = 5,
                          verbose = 0,
                          maximize = F)
  
  #final XGB model with optimal hyperparameters and nround
  dtrain
  xg_mod = xgboost::xgboost(y = outcome, x=train_dat, 
                            booster = opt_params$booster,
                            objective = opt_params$objective,
                            eval_metric=opt_params$eval_metric,
                            nthread=opt_params$nthread,
                            reg_lambda=opt_params$lambda, 
                            reg_alpha=opt_params$alpha, 
                            learning_rate=opt_params$eta,
                            
                            max_depth=opt_params$max_depth,
                            min_child_weight=opt_params$min_child_weight,
                            subsample=opt_params$subsample,
                            
                            nround = xgbcv$evaluation_log$iter[which.min(xgbcv$evaluation_log$test_mae_mean)])
  return(xg_mod)
}

