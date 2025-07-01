#' Sigmoid activation function
#'
#' Applies the sigmoid (logistic) function to transform values to (0,1) range.
#'
#' @param x Numeric vector or matrix of input values
#' @return Numeric vector or matrix of transformed values between 0 and 1
#' @export
#' @examples
#' sigmoid(0)    # Returns 0.5
#' sigmoid(c(-2, 0, 2))  # Returns values between 0 and 1
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

#' SMAGS: Sensitivity Maximization at a Given Specificity
#'
#' Performs sensitivity maximization at a specified specificity level using
#' logistic regression with custom optimization across multiple algorithms.
#'
#' @param X Feature matrix of dimensions n x p (samples x features)
#' @param y Binary response vector of length n (0 = negative, 1 = positive)
#' @param SP Target specificity level, a value between 0 and 1
#' @return A list containing:
#' \itemize{
#'   \item Optimized coefficients (including intercept)
#'   \item Optimization method used
#'   \item Tolerance level
#'   \item Threshold value
#'   \item Achieved sensitivity
#' }
#' @export
#' @importFrom stats glm coef quantile optim binomial
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %:% %dopar% %dopar%
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 100
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("feature_", 1:p)
#' y <- rbinom(n, 1, 0.5)
#'
#' # Run SMAGS with 95% specificity
#' result <- SMAGS(X, y, SP = 0.95)
#' print(result)
#' }
SMAGS <- function(X, y, SP) {
  # Convert inputs to matrix/vector if needed
  X <- as.matrix(X)
  y <- as.numeric(y)

  # Validate inputs
  if (SP <= 0 || SP >= 1) {
    stop("SP must be between 0 and 1")
  }
  if (length(unique(y)) != 2) {
    stop("y must be binary (0/1)")
  }
  if (nrow(X) != length(y)) {
    stop("Number of rows in X must equal length of y")
  }

  # Initial logistic regression
  init_model <- glm(y ~ ., data = data.frame(y = y, X), family = binomial())
  all_coefs <- c(coef(init_model)[-1], coef(init_model)[1])

  # Add intercept column
  X1 <- cbind(X, I = 1)

  # Custom loss function
  custom_loss <- function(coefs) {
    m <- as.matrix(X1) %*% coefs
    z <- sigmoid(m)
    threshold <- quantile(((1 - y) * z)[(1 - y) * z != 0], prob = SP)
    if(is.na(threshold)) threshold <- 0

    y_hat <- as.numeric(z > threshold)
    loss1 <- sum(y_hat * y) / sum(y)  # Sensitivity
    return(-loss1)  # Return negative sensitivity for minimization
  }

  # Set up optimization parameters
  opt_methods <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")
  tolerance <- c(1e-20, 1e-50, 1e-100)

  # Process parameters function
  process_params <- function(method, tol) {
    tryCatch({
      result <- optim(par = all_coefs,
                      fn = custom_loss,
                      method = method,
                      control = list(trace = FALSE, maxit = 1000, reltol = tol))

      # Calculate performance metrics
      pred <- X1 %*% result$par
      prob <- sigmoid(pred)
      threshold <- quantile(((1 - y) * z)[(1 - y) * z != 0], prob = SP)
      if(is.na(threshold)) threshold <- 0
      pred_class <- as.numeric(prob > threshold)

      sensitivity <- sum(pred_class * y) / sum(y)

      # Return results
      c(result$par,
        method = method,
        tolerance = tol,
        threshold = SP,
        Sensitivity = sensitivity)
    }, error = function(e) {
      rep(NA, length(all_coefs) + 4)
    })
  }

  # Set up parallel processing
  n_cores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  # Run optimization in parallel
  results <- foreach::foreach(method = opt_methods, .combine = rbind) %:%
    foreach::foreach(tol = tolerance, .combine = rbind) %dopar% {
      process_params(method, tol)
    }

  parallel::stopCluster(cl)

  # Clean results and select best model
  results <- as.data.frame(results)
  results <- results[complete.cases(results),]

  if (nrow(results) == 0) {
    stop("All optimization attempts failed")
  }

  best_model <- results[which.max(results$Sensitivity),]

  return(best_model)
}

#' SMAGS with LASSO Regularization
#'
#' Combines SMAGS optimization with LASSO regularization for feature selection
#' while maximizing sensitivity at a given specificity.
#'
#' @param X Feature matrix of dimensions n x p (samples x features)
#' @param y Binary response vector of length n (0 = negative, 1 = positive)
#' @param SP Target specificity level, a value between 0 and 1
#' @param lambda Regularization parameter (>= 0). Higher values increase sparsity.
#' @return A list containing:
#' \itemize{
#'   \item Optimized coefficients (including intercept)
#'   \item Optimization method used
#'   \item Tolerance level
#'   \item Threshold value
#'   \item Achieved specificity
#'   \item Achieved sensitivity
#'   \item AIC value
#' }
#' @export
#' @importFrom stats glm coef quantile optim binomial
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %:%
#' @examples
#' \dontrun{
#' # Generate example data with some noise features
#' set.seed(123)
#' n <- 200
#' p <- 50
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("feature_", 1:p)
#' # Make only first 5 features informative
#' y <- as.numeric(rowSums(X[, 1:5]) + rnorm(n) > 0)
#'
#' # Run SMAGS-LASSO
#' result <- SMAGS_LASSO(X, y, SP = 0.90, lambda = 0.1)
#'
#' # Check which features were selected (non-zero coefficients)
#' selected_features <- abs(as.numeric(result[1:p])) > 1e-6
#' print(which(selected_features))
#' }
SMAGS_LASSO <- function(X, y, SP, lambda) {
  # Convert inputs to matrix/vector if needed
  X <- as.matrix(X)
  y <- as.numeric(y)

  # Validate inputs
  if (SP <= 0 || SP >= 1) {
    stop("SP must be between 0 and 1")
  }
  if (lambda < 0) {
    stop("lambda must be non-negative")
  }
  if (length(unique(y)) != 2) {
    stop("y must be binary (0/1)")
  }
  if (nrow(X) != length(y)) {
    stop("Number of rows in X must equal length of y")
  }

  # Initial logistic regression
  init_model <- glm(y ~ ., data = data.frame(y = y, X), family = binomial())
  all_coefs <- c(coef(init_model)[-1], coef(init_model)[1])

  # Add intercept column
  X1 <- cbind(X, I = 1)

  # Custom loss function with LASSO penalty
  custom_loss <- function(coefs, lambda) {
    m <- X1 %*% coefs
    z <- sigmoid(m)
    threshold <- quantile(((1 - y) * z)[(1 - y) * z != 0], prob = SP)
    if(is.na(threshold)) threshold <- 0

    y_hat <- as.numeric(z > threshold)
    loss1 <- sum(y_hat * y) / sum(y * y)  # Sensitivity
    l1_penalty <- lambda * sum(abs(coefs[1:(length(coefs)-1)]))  # L1 penalty (exclude intercept)
    return(-loss1 + l1_penalty)  # Return negative sensitivity + penalty for minimization
  }

  # Set up optimization parameters
  opt_methods <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")
  tolerance <- c(1e-20, 1e-50, 1e-100)

  # Process parameters function
  process_params <- function(method, tol) {
    tryCatch({
      result <- optim(par = all_coefs,
                      fn = custom_loss,
                      method = method,
                      lambda = lambda,
                      control = list(trace = FALSE, maxit = 9000, reltol = tol))

      # Calculate performance metrics
      pred <- X1 %*% result$par
      prob <- sigmoid(pred)
      threshold <- quantile(((1 - y) * z)[(1 - y) * z != 0], prob = SP)
      if(is.na(threshold)) threshold <- 0
      pred_class <- as.numeric(prob > threshold)

      specificity <- sum((1 - pred_class) * (1 - y)) / sum(1 - y)
      sensitivity <- sum(pred_class * y) / sum(y)

      # Calculate AIC
      n_features <- ncol(X)
      loglik <- -sum(y * log(prob + 1e-15) + (1 - y) * log(1 - prob + 1e-15))
      aic <- 2 * ((n_features + 1) - loglik)

      # Return results
      c(result$par,
        method = method,
        tolerance = tol,
        threshold = threshold[[1]],
        Specificity = specificity,
        Sensitivity = sensitivity,
        AIC = aic)
    }, error = function(e) {
      rep(NA, length(all_coefs) + 7)
    })
  }

  # Set up parallel processing
  n_cores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  # Run optimization in parallel
  results <- foreach::foreach(method = opt_methods, .combine = rbind) %:%
    foreach::foreach(tol = tolerance, .combine = rbind) %dopar% {
      process_params(method, tol)
    }

  parallel::stopCluster(cl)

  # Clean results and select best model
  results <- as.data.frame(results)
  results <- results[complete.cases(results),]

  if (nrow(results) == 0) {
    stop("All optimization attempts failed")
  }

  best_model <- results[which.max(results$Sensitivity),]

  return(best_model)
}
