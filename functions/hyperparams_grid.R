# this list defines parameter grids for hyperparameter tuning of each learner

# --- parameter grids for tuning ---
param_grids <- list(
  random_forest = expand.grid(
    num.trees = c(100, 500),             # number of trees
    sample.fraction = c(0.5, 0.75),      # fraction of data used per tree
    min.node.size = c(5, 10, 20),        # minimum samples per leaf
    mtry = c(2, 5)                       # number of variables randomly sampled at each split
  ),
  xgboost = expand.grid(
    eta = c(0.01, 0.05, 0.1, 0.3),      # learning rate
    max_depth = c(2, 3, 4),             # tree depth
    nrounds = c(300)                    # number of boosting rounds
  ),
  lasso = NULL  # no grid needed for glmnet
)
