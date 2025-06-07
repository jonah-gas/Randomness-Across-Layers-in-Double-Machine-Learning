# this list stores tuning, training, and prediction methods for each learner

# --- define learner-specific tuning, training, and prediction methods ---
learner_methods <- list(
  random_forest = list(
    tune = function(x, y, type) tune_rf(x, y, param_grids$random_forest, type), # random forest tuning
    train = function(x, y, params) regression_forest( # random forest training
      x, y,
      num.trees = params$num.trees,
      sample.fraction = params$sample.fraction,
      min.node.size = params$min.node.size,
      mtry = params$mtry,
      honesty = FALSE, ci.group.size = 1
    ),
    predict = function(model, x) predict(model, x)$predictions # random forest prediction
  ),
  
  xgboost = list(
    tune = function(x, y, type) tune_xgb(x, y, param_grids$xgboost, type), # xgboost tuning
    train = function(x, y, params, obj) xgboost( # xgboost training
      data = xgb.DMatrix(x, label = y),
      max_depth = params$max_depth, nrounds = params$nrounds, eta = params$eta,
      objective = obj, verbose = 0
    ),
    predict = function(model, x) predict(model, x) # xgboost prediction
  ),
  
  lasso = list(
    tune = function(x, y, type) tune_lasso(x, y, type), # lasso tuning
    train = function(x, y, params, family) glmnet( # lasso training
      x, y, alpha = 1, lambda = params$lambda, family = family, standardize = TRUE
    ),
    predict = function(model, x, params, type) as.vector(predict(model, x, s = params$lambda, type = type)) # lasso prediction
  )
)
