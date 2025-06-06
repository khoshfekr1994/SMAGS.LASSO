#' Cross-Validation for SMAGS-LASSO
#'
#' Performs k-fold cross-validation to select optimal lambda parameter for
#' SMAGS-LASSO by minimizing sensitivity mean squared error.
#'
#' @param X Feature matrix of dimensions n x p (samples x features)
#' @param y Binary response vector of length n (0 = negative, 1 = positive)
#' @param SP Target specificity level, a value between 0 and 1
#' @param lambda_seq Sequence of lambda values to test. If NULL, automatically generated.
#' @param nfolds Number of cross-validation folds (default: 10)
#' @param seed Random seed for reproducibility (default: 123)
#' @return A list containing:
#' \itemize{
#'   \item lambda_min: Optimal lambda value
#'   \item lambda_seq: Sequence of tested lambda values
#'   \item mean_sensitivity_mse: Mean squared error for sensitivity across folds
#'   \item mean_sensitivity: Mean sensitivity across folds for each lambda
#'   \item mean_norm_ratio: Mean coefficient norm ratio across folds
#'   \item final_model: SMAGS-LASSO model fitted with optimal lambda
#'   \item cv_results: Data frame with detailed CV results
#'   \item plot: Function to create CV plot
#' }
#' @export
#' @importFrom stats quantile
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 200
#' p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("feature_", 1:p)
#' y <- rbinom(n, 1, 0.5)
#'
#' # Run cross-validation
#' cv_result <- cv_SMAGS_LASSO(X, y, SP = 0.95)
#'
#' # Plot results
#' cv_result$plot()
#'
#' # Get optimal lambda
#' print(cv_result$lambda_min)
#' }
cv_SMAGS_LASSO <- function(X, y, SP, lambda_seq = NULL, nfolds = 10, seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)

  # Convert inputs to matrix/vector if needed
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)

  # Validate inputs
  if (SP <= 0 || SP >= 1) {
    stop("SP must be between 0 and 1")
  }
  if (nfolds < 2) {
    stop("nfolds must be at least 2")
  }
  if (nrow(X) != length(y)) {
    stop("Number of rows in X must equal length of y")
  }

  # Create lambda sequence if not provided
  if (is.null(lambda_seq)) {
    lambda_max <- max(abs(t(X) %*% y)) / n
    lambda_min <- lambda_max * 1e-10
    lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out = 100))
  }

  # Fit full model to calculate the norm ratio
  full_model <- glm(y ~ ., data = data.frame(y = y, X), family = binomial())
  full_coefs <- coef(full_model)[-1]  # Excluding intercept
  beta_norm_full <- sum(abs(full_coefs))  # ||β̂||₁

  # Create folds
  folds <- createFolds(y, k = nfolds, list = TRUE, returnTrain = FALSE)

  # Matrix to store MSE for each fold and lambda
  mse_matrix <- matrix(NA, nrow = nfolds, ncol = length(lambda_seq))

  # For sensitivity tracking
  sensitivity_matrix <- matrix(NA, nrow = nfolds, ncol = length(lambda_seq))

  # For tracking coefficient norm ratios
  norm_ratio_matrix <- matrix(NA, nrow = nfolds, ncol = length(lambda_seq))

  # Loop through folds
  for (i in 1:nfolds) {
    # Split data into training and validation
    train_indices <- setdiff(1:n, folds[[i]])
    val_indices <- folds[[i]]

    X_train <- X[train_indices, , drop = FALSE]
    y_train <- y[train_indices]
    X_val <- X[val_indices, , drop = FALSE]
    y_val <- y[val_indices]

    # Loop through lambda values
    for (j in 1:length(lambda_seq)) {
      tryCatch({
        # Fit model on training data
        model <- SMAGS_LASSO(X_train, y_train, SP, lambda_seq[j])

        # Extract coefficients (excluding metadata)
        coefs <- as.numeric(model[1:(ncol(X)+1)])

        # Calculate coefficient norm ratio ||β̂_λ||₁/||β̂||₁
        beta_lambda <- coefs[1:ncol(X)]  # Excluding intercept
        beta_norm_lambda <- sum(abs(beta_lambda))
        norm_ratio <- beta_norm_lambda / beta_norm_full
        norm_ratio_matrix[i, j] <- norm_ratio

        # Make predictions on validation set
        X_val1 <- cbind(X_val, I = 1)  # Add intercept column
        pred <- X_val1 %*% coefs

        # Apply sigmoid to get probabilities
        prob <- sigmoid(pred)

        # Determine threshold based on SP
        threshold <- as.numeric(model["threshold"])

        if(is.na(threshold)) threshold <- SP

        # Get predicted classes
        pred_class <- as.numeric(prob > threshold)

        # Calculate sensitivity
        if (sum(y_val) > 0) {
          sensitivity <- sum(pred_class * y_val) / sum(y_val)
        } else {
          sensitivity <- NA
        }
        sensitivity_matrix[i, j] <- sensitivity

        # Calculate MSE for sensitivity (error from optimal sensitivity of 1.0)
        mse_matrix[i, j] <- (1 - sensitivity)^2

      }, error = function(e) {
        # If model fitting fails, record NA
        mse_matrix[i, j] <<- NA
        sensitivity_matrix[i, j] <<- NA
        norm_ratio_matrix[i, j] <<- NA
      })
    }
  }

  # Calculate means across folds for each lambda
  mean_mse <- colMeans(mse_matrix, na.rm = TRUE)
  mean_sensitivity <- colMeans(sensitivity_matrix, na.rm = TRUE)
  mean_norm_ratio <- colMeans(norm_ratio_matrix, na.rm = TRUE)

  # Create a data frame with results for easier analysis
  cv_results_df <- data.frame(
    lambda = lambda_seq,
    sensitivity_mse = mean_mse,
    sensitivity = mean_sensitivity,
    norm_ratio = mean_norm_ratio
  )

  # Find lambda with minimum MSE
  valid_indices <- !is.na(mean_mse) & !is.infinite(mean_mse)
  if (sum(valid_indices) == 0) {
    stop("All cross-validation attempts failed")
  }

  min_index <- which.min(mean_mse[valid_indices])
  lambda_min <- lambda_seq[valid_indices][min_index]

  # Fit final model with optimal lambda
  final_model <- SMAGS_LASSO(X, y, SP, lambda_min)

  # Return results
  return(list(
    lambda_min = lambda_min,
    lambda_seq = lambda_seq,
    mean_sensitivity_mse = mean_mse,
    mean_sensitivity = mean_sensitivity,
    mean_norm_ratio = mean_norm_ratio,
    final_model = final_model,
    cv_results = cv_results_df,

    # Add a plot function
    plot = function() {
      plot(cv_results_df$norm_ratio, cv_results_df$sensitivity_mse,
           type = "l",
           xlab = expression("||"*hat(beta)[lambda]*"||"[1]*"/"*"||"*hat(beta)*"||"[1]),
           ylab = "Cross-Validation Error",
           main = "Sensitivity MSE vs. Coefficient Norm Ratio")
      # Add dotted vertical line at minimum MSE
      min_ratio <- cv_results_df$norm_ratio[min_index]
      abline(v = min_ratio, lty = 2)
    }
  ))
}

#' Create Cross-Validation Folds
#'
#' Creates stratified folds for cross-validation, ensuring balanced class
#' distribution across folds when possible.
#'
#' @param y Response vector for stratification
#' @param k Number of folds (default: 10)
#' @param list Return as list (default: TRUE)
#' @param returnTrain Return training indices instead of test indices (default: FALSE)
#' @return List of fold indices (test sets) or vector of fold assignments
#' @export
#' @examples
#' \dontrun{
#' # Create 5-fold CV indices
#' y <- c(rep(0, 50), rep(1, 50))
#' folds <- createFolds(y, k = 5)
#'
#' # Check fold sizes
#' sapply(folds, length)
#' }
createFolds <- function(y, k = 10, list = TRUE, returnTrain = FALSE) {
  # Check if caret is available for stratified sampling
  if (requireNamespace("caret", quietly = TRUE)) {
    return(caret::createFolds(y, k = k, list = list, returnTrain = returnTrain))
  } else {
    # Simple implementation if caret is not available
    n <- length(y)

    # For binary outcomes, try to stratify
    if (length(unique(y)) == 2) {
      y_factor <- as.factor(y)
      class_0_indices <- which(y_factor == levels(y_factor)[1])
      class_1_indices <- which(y_factor == levels(y_factor)[2])

      # Randomly assign each class to folds
      folds_0 <- rep(1:k, length.out = length(class_0_indices))
      folds_1 <- rep(1:k, length.out = length(class_1_indices))

      # Shuffle within each class
      folds_0 <- sample(folds_0)
      folds_1 <- sample(folds_1)

      # Combine fold assignments
      fold_indices <- numeric(n)
      fold_indices[class_0_indices] <- folds_0
      fold_indices[class_1_indices] <- folds_1
    } else {
      # Non-stratified for non-binary outcomes
      fold_indices <- sample(rep(1:k, length.out = n))
    }

    if (list) {
      folds <- vector("list", k)
      for (i in 1:k) {
        if (returnTrain) {
          folds[[i]] <- which(fold_indices != i)
        } else {
          folds[[i]] <- which(fold_indices == i)
        }
      }
      return(folds)
    } else {
      return(fold_indices)
    }
  }
}
