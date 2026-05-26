## WRAPPERS FOR BAYESIAN OPTIMIZED XGBOOSTING
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

##################################################################################################################
##################################################################################################################
#XGBlinear
XGBlinear <- function(train_dat, outcome, alpharange = c(0,10), lambdarange = c(0,10),
                      etarange = c(0.001, 0.3), nthread = 5, nround = 500) {
  
  # Check and install required packages
  list.of.packages <- c("doParallel", "parallel", "ParBayesianOptimization", "xgboost", "caret")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
  if (length(new.packages)) {
    cat(paste("Installing required packages:\n", paste(new.packages, collapse = ", "), "\n"))
    install.packages(new.packages)
  }
  
  folds <- caret::createFolds(outcome, k = 5)
  bounds <- list(
    eta        = c(0.001, 0.1),
    log_lambda = c(log(0.1), log(100)),
    log_alpha  = c(log(0.1), log(100))
  )
  
  # Save DMatrix to disk immediately, don't keep in memory
  dtrain <- xgboost::xgb.DMatrix(data = train_dat, label = outcome)
  xgboost::xgb.DMatrix.save(dtrain, "train.buffer")
  rm(dtrain, train_dat)   # <-- also drop train_dat reference early
  gc()
  
  cl <- parallel::makeCluster(nthread)
  doParallel::registerDoParallel(cl)
  on.exit({
    parallel::stopCluster(cl)
    unlink("train.buffer")  # <-- moved here so it always cleans up even on error
    gc()
  }, add = TRUE)
  
  obj_func <- function(eta, log_lambda, log_alpha, nround) {
    # Back-transform first
    lambda <- exp(log_lambda)
    alpha  <- exp(log_alpha)
    
    param <- list(
      nthread = 1,
      learning_rate  = eta,
      reg_lambda     = lambda,
      reg_alpha      = alpha,
      booster        = "gblinear",
      objective      = "reg:squarederror",
      eval_metric    = "mae"
    )
    
    xgbcv <- xgboost::xgb.cv(
      params               = param,
      data                 = xgboost::xgb.DMatrix("train.buffer"),
      nround               = 1000,
      folds                = folds,
      prediction           = FALSE,
      early_stopping_rounds = 30,
      verbose              = 0,
      maximize             = FALSE,
      save_models           = FALSE 
    )
    
    # Extract only the two scalars needed — drop the whole log object
    best_row  <- which.min(xgbcv$evaluation_log$test_mae_mean)
    best_score  <- xgbcv$evaluation_log$test_mae_mean[best_row]
    best_iter   <- xgbcv$evaluation_log$iter[best_row]
    rm(xgbcv); gc()
    
    list(Score = -best_score, nrounds = as.integer(best_iter))
  }
  
  bayes_out <- ParBayesianOptimization::bayesOpt(
    FUN        = obj_func,
    bounds     = bounds,
    initPoints = 20,
    parallel   = TRUE,
    iters.n    = 50,
    iters.k    = nthread,
    acqThresh  = 0.0,
    eps = 0.5
  )
  
  opt_params <- c(
    list(booster = "gblinear", objective = "reg:squarederror",
         eval_metric = "mae", nthread = nthread),
    ParBayesianOptimization::getBestPars(bayes_out)
  )
  rm(bayes_out); gc()   # <-- free Bayesian optimisation object before final CV
  
  xgbcv <- xgboost::xgb.cv(
    params                = opt_params,
    data                  = xgboost::xgb.DMatrix("train.buffer"),
    nround                = nround,
    folds                 = folds,
    prediction            = FALSE,
    early_stopping_rounds = 5,
    verbose               = 0,
    maximize              = FALSE
  )
  
  opt.nround <- xgbcv$evaluation_log$iter[which.min(xgbcv$evaluation_log$test_mae_mean)]
  rm(xgbcv, folds); gc()   # <-- folds no longer needed after this point
  
  # Reconstruct outcome vector from the buffer to avoid holding train_dat
  xg_mod <- xgboost::xgb.train(
    params  = opt_params,
    data    = xgboost::xgb.DMatrix("train.buffer"),
    nrounds = opt.nround,
    verbose = 0
  )
  
  return(xg_mod)
}
##################################################################################################################
##################################################################################################################
#XGBtree
XGBtree = function(train_dat, outcome, alpharange=c(0,1), lambdarange=c(1,5),
                   etarange=c(0.001, 0.2), max_depth=c(1L, 10L),
                   min_child_weight=c(1, 50), subsample=c(0.1, 1),
                   nthread=2, nround=500)
{
  buffer_path <- file.path(getwd(), "train.buffer")
  on.exit(unlink(buffer_path), add = TRUE)          # always clean up
  
  list.of.packages <- c("doParallel", "ParBayesianOptimization", "xgboost")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) {
    cat(paste("Installing required packages:\n", paste(new.packages, collapse=", "), "\n"))
    install.packages(new.packages)
  }
  
  # --- Save to disk, then FREE the large objects immediately ---
  folds   <- caret::createFolds(outcome, k=5)
  dtrain  <- xgboost::xgb.DMatrix(data=train_dat, label=outcome)
  xgboost::xgb.DMatrix.save(dtrain, buffer_path)
  rm(dtrain, train_dat, outcome)
  gc()
  
  bounds <- list(
    lambda           = lambdarange,
    alpha            = alpharange,
    eta              = etarange,
    max_depth        = max_depth,
    min_child_weight = min_child_weight,
    subsample        = subsample
  )
  
  cl <- parallel::makeCluster(nthread)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # Export only what workers actually need — not the raw data
  parallel::clusterExport(cl, c("folds", "buffer_path"), envir = environment())
  
  # --- Objective: load DMatrix from disk, never hold it across iterations ---
  obj_func <- function(eta, lambda, alpha, nround, max_depth, min_child_weight, subsample)
  {
    param <- list(
      learning_rate    = eta,
      reg_lambda       = lambda,
      reg_alpha        = alpha,
      max_depth        = max_depth,
      min_child_weight = min_child_weight,
      subsample        = subsample,
      booster          = "gbtree",
      objective        = "reg:squarederror",
      eval_metric      = "mae",
      nthread          = 1         # workers share CPU; avoid nested parallelism
    )
    
    xgbcv <- xgboost::xgb.cv(
      params               = param,
      data                 = xgboost::xgb.DMatrix(buffer_path),  # load fresh, local scope
      nround               = 500,
      folds                = folds,
      prediction           = FALSE,
      early_stopping_rounds = 5,
      verbose              = 0,
      maximize             = FALSE
    )
    
    score     <- -min(xgbcv$evaluation_log$test_mae_mean)
    best_iter <- as.integer(xgbcv$evaluation_log$iter[which.min(xgbcv$evaluation_log$test_mae_mean)])
    rm(xgbcv); gc()
    
    list(Score = score, nrounds = best_iter)
  }
  
  bayes_out <- ParBayesianOptimization::bayesOpt(
    FUN        = obj_func,
    bounds     = bounds,
    initPoints = length(bounds) + 2,
    parallel   = TRUE,
    iters.n    = nthread,
    iters.k    = nthread
  )
  gc()
  
  best_pars <- ParBayesianOptimization::getBestPars(bayes_out)
  rm(bayes_out); gc()
  
  opt_params <- c(
    list(booster     = "gbtree",
         objective   = "reg:squarederror",
         eval_metric = "mae",
         nthread     = nthread),
    best_pars
  )
  
  # --- Final CV to find optimal nround ---
  xgbcv <- xgboost::xgb.cv(
    params                = opt_params,
    data                  = xgboost::xgb.DMatrix(buffer_path),
    nround                = nround,
    folds                 = folds,
    prediction            = FALSE,
    early_stopping_rounds = 5,
    verbose               = 0,
    maximize              = FALSE
  )
  opt.nround <- xgbcv$evaluation_log$iter[which.min(xgbcv$evaluation_log$test_mae_mean)]
  rm(xgbcv); gc()
  
  # --- Train final model directly from the buffer ---
  xg_mod <- xgboost::xgb.train(
    params  = opt_params,
    data    = xgboost::xgb.DMatrix(buffer_path),   # no train_dat in memory
    nrounds = opt.nround,
    verbose = 0
  )
  
  return(xg_mod)
}

