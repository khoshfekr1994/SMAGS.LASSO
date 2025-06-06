#' ROC Analysis for Synthetic Data
#'
#' Compares LASSO, SMAGS, and SMAGS-LASSO performance on synthetic datasets
#' by generating ROC curves and performance metrics for both training and test data.
#'
#' @param X_train Training feature matrix
#' @param y_train Training binary response vector
#' @param X_test Test feature matrix
#' @param y_test Test binary response vector
#' @param SP Target specificity level
#' @param lambda_lasso Lambda parameter for standard LASSO
#' @param lambda_smags Lambda parameter for SMAGS-LASSO
#' @param save_path_a Optional path to save training ROC plot
#' @param save_path_b Optional path to save test ROC plot
#' @return List containing performance metrics, plots, and model coefficients
#' @export
#' @importFrom glmnet glmnet
#' @importFrom ROCR prediction performance
#' @importFrom ggplot2 ggplot aes geom_line geom_vline scale_color_manual labs theme_classic theme element_text unit ggsave
#' @importFrom stats glm predict
#' @examples
#' \dontrun{
#' # Generate synthetic data
#' set.seed(123)
#' n <- 200
#' p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' y <- rbinom(n, 1, 0.5)
#'
#' # Split data
#' split <- train_test_split(X, y)
#'
#' # Run analysis
#' results <- synthetic_data_results(
#'   split$X_train, split$y_train, split$X_test, split$y_test,
#'   SP = 0.95, lambda_lasso = 0.01, lambda_smags = 0.1
#' )
#' }
synthetic_data_results <- function(X_train, y_train, X_test, y_test, SP,
                                   lambda_lasso, lambda_smags,
                                   save_path_a = NULL, save_path_b = NULL) {

  # Helper function to calculate sensitivity at a given SP using quantile method
  calculate_metrics_at_sp <- function(predictions, y_true, SP) {
    # Get quantile threshold based on predictions for negative cases
    threshold <- quantile((1 - y_true) * predictions[(1 - y_true) * predictions != 0], prob = SP)
    if(is.na(threshold)) threshold <- 0

    # Make binary predictions based on threshold
    pred_class <- as.numeric(predictions > threshold)

    # Calculate metrics
    sensitivity <- sum(pred_class * y_true) / sum(y_true)
    specificity <- sum((1 - pred_class) * (1 - y_true)) / sum(1 - y_true)

    return(list(
      threshold = threshold,
      sensitivity = sensitivity,
      specificity = specificity,
      pred_class = pred_class
    ))
  }

  # Function to create a ggplot ROC curve
  create_roc_plot <- function(
    perf_objects,
    metrics_list,
    nonzero_counts,
    title,
    SP,
    lambda_lasso,
    lambda_smags,
    include_specificity = FALSE
  ) {

    # Create data frames for each curve
    roc_data <- data.frame()

    # Add LASSO curve data
    lasso_data <- data.frame(
      x = unlist(perf_objects[[1]]@x.values),
      y = unlist(perf_objects[[1]]@y.values),
      Method = "LASSO"
    )
    roc_data <- rbind(roc_data, lasso_data)

    # Add SMAGS curve data
    smags_data <- data.frame(
      x = unlist(perf_objects[[2]]@x.values),
      y = unlist(perf_objects[[2]]@y.values),
      Method = "SMAGS"
    )
    roc_data <- rbind(roc_data, smags_data)

    # Add SMAGS-LASSO curve data
    smags_lasso_data <- data.frame(
      x = unlist(perf_objects[[3]]@x.values),
      y = unlist(perf_objects[[3]]@y.values),
      Method = "SMAGS-LASSO"
    )
    roc_data <- rbind(roc_data, smags_lasso_data)

    # Create legend labels
    if (include_specificity) {
      legend_labels <- c(
        sprintf("LASSO (λ=%.4f) - Sens = %.2f, Spec = %.2f (NonZero = %d)",
                lambda_lasso, metrics_list[[1]]$sensitivity, metrics_list[[1]]$specificity, nonzero_counts[1]),
        sprintf("SMAGS - Sens = %.2f, Spec = %.2f (NonZero = %d)",
                metrics_list[[2]]$sensitivity, metrics_list[[2]]$specificity, nonzero_counts[2]),
        sprintf("SMAGS-LASSO (λ=%.4f) - Sens = %.2f, Spec = %.2f (NonZero = %d)",
                lambda_smags, metrics_list[[3]]$sensitivity, metrics_list[[3]]$specificity, nonzero_counts[3])
      )
    } else {
      legend_labels <- c(
        sprintf("LASSO (λ=%.4f) - Sens at SP = %.2f (NonZero = %d)",
                lambda_lasso, metrics_list[[1]]$sensitivity, nonzero_counts[1]),
        sprintf("SMAGS - Sens at SP = %.2f (NonZero = %d)",
                metrics_list[[2]]$sensitivity, nonzero_counts[2]),
        sprintf("SMAGS-LASSO (λ=%.4f) - Sens at SP = %.2f (NonZero = %d)",
                lambda_smags, metrics_list[[3]]$sensitivity, nonzero_counts[3])
      )
    }

    # Create the plot
    p <- ggplot() +
      # Add ROC curves
      geom_line(data = roc_data, aes(x = x, y = y, color = Method), linewidth = 1.2) +
      # Add vertical line at 1-SP
      geom_vline(xintercept = 1-SP, linetype = "dashed", color = "gray") +
      # Set colors
      scale_color_manual(values = c("LASSO" = "green", "SMAGS" = "orange", "SMAGS-LASSO" = "red"),
                         labels = legend_labels) +
      # Set labels
      labs(
        title = title,
        x = "False positive rate",
        y = "True positive rate"
      ) +
      # Set theme with larger font
      theme_classic(base_size = 16) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = "gray90"),
        legend.margin = ggplot2::margin(6, 6, 6, 6),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        legend.spacing.y = unit(0.05, "cm"),
        legend.box = "vertical",
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14)
      )
    return(p)
  }

  # Fit logistic regression
  log_model <- glm(y_train ~ ., data = data.frame(y_train, X_train), family = binomial())
  pred_train_lr <- predict(log_model, newdata = data.frame(X_train), type = "response")

  # Get number of non-zero coefficients for LR
  lr_coef <- coef(log_model)[-1]  # Exclude intercept
  lr_nonzero <- sum(abs(lr_coef) > 1e-3)

  # Fit LASSO model
  lasso_model <- glmnet(as.matrix(X_train), y_train, family = "binomial", lambda = lambda_lasso)
  pred_train_lasso <- predict(lasso_model, newx = as.matrix(X_train), type = "response")

  # Get number of non-zero coefficients for LASSO
  lasso_coef <- coef(lasso_model, s = lambda_lasso)[-1]  # Exclude intercept
  lasso_nonzero <- sum(abs(lasso_coef) > 0)

  # Fit SMAGS models
  smags_result <- SMAGS_LASSO(X_train, y_train, SP, lambda = lambda_smags)
  M_smags_result <- SMAGS(X_train, y_train, SP)

  # Get predictions
  smags_pred_train <- sigmoid(as.matrix(cbind(X_train, 1)) %*% as.numeric(smags_result[1:(ncol(X_train) + 1)]))
  M_smags_pred_train <- sigmoid(as.matrix(cbind(X_train, 1)) %*% as.numeric(M_smags_result[1:(ncol(X_train) + 1)]))

  # Get number of non-zero coefficients for SMAGS models
  smags_coef <- smags_result[1:ncol(X_train)]  # Exclude intercept
  nonzero_param_coef_smags <- 0.05 * max(abs(as.numeric(smags_coef)))
  smags_nonzero <- sum(abs(as.numeric(smags_coef)) > nonzero_param_coef_smags)

  M_smags_coef <- M_smags_result[1:ncol(X_train)]  # Exclude intercept
  M_smags_nonzero <- sum(abs(as.numeric(M_smags_coef)) > 0)

  # Calculate metrics using quantile method
  lr_metrics <- calculate_metrics_at_sp(pred_train_lr, y_train, SP)
  lasso_metrics <- calculate_metrics_at_sp(pred_train_lasso, y_train, SP)
  smags_metrics <- calculate_metrics_at_sp(smags_pred_train, y_train, SP)
  M_smags_metrics <- calculate_metrics_at_sp(M_smags_pred_train, y_train, SP)

  # Calculate ROC curves for plotting
  pred_obj_lr <- prediction(pred_train_lr, y_train)
  perf_obj_lr <- performance(pred_obj_lr, "tpr", "fpr")

  pred_obj_lasso <- prediction(pred_train_lasso, y_train)
  perf_obj_lasso <- performance(pred_obj_lasso, "tpr", "fpr")

  smags_pred_obj <- prediction(smags_pred_train, y_train)
  smags_perf_obj <- performance(smags_pred_obj, "tpr", "fpr")

  M_smags_pred_obj <- prediction(M_smags_pred_train, y_train)
  M_smags_perf_obj <- performance(M_smags_pred_obj, "tpr", "fpr")

  # Create training plot
  train_perf_objects <- list(perf_obj_lasso, M_smags_perf_obj, smags_perf_obj)
  train_metrics_list <- list(lasso_metrics, M_smags_metrics, smags_metrics)
  train_nonzero_counts <- c(lasso_nonzero, M_smags_nonzero, smags_nonzero)

  train_plot <- create_roc_plot(
    perf_objects = train_perf_objects,
    metrics_list = train_metrics_list,
    nonzero_counts = train_nonzero_counts,
    title = "ROC Curve (Training Data)",
    SP = SP,
    lambda_lasso = lambda_lasso,
    lambda_smags = lambda_smags
  )

  # Test data predictions
  pred_test_lr <- predict(log_model, newdata = data.frame(as.matrix(X_test)), type = "response")
  pred_test_lasso <- predict(lasso_model, newx = as.matrix(X_test), type = "response")
  smags_pred_test <- sigmoid(as.matrix(cbind(X_test, 1)) %*% as.numeric(smags_result[1:(ncol(X_test) + 1)]))
  M_smags_pred_test <- sigmoid(as.matrix(cbind(X_test, 1)) %*% as.numeric(M_smags_result[1:(ncol(X_test) + 1)]))

  # Calculate test metrics using SAME thresholds from training data
  test_lr_metrics <- list(
    sensitivity = sum((pred_test_lr > lr_metrics$threshold) * y_test) / sum(y_test),
    specificity = sum((pred_test_lr <= lr_metrics$threshold) * (1 - y_test)) / sum(1 - y_test)
  )

  test_lasso_metrics <- list(
    sensitivity = sum((pred_test_lasso > lasso_metrics$threshold) * y_test) / sum(y_test),
    specificity = sum((pred_test_lasso <= lasso_metrics$threshold) * (1 - y_test)) / sum(1 - y_test)
  )

  test_smags_metrics <- list(
    sensitivity = sum((smags_pred_test > smags_metrics$threshold) * y_test) / sum(y_test),
    specificity = sum((smags_pred_test <= smags_metrics$threshold) * (1 - y_test)) / sum(1 - y_test)
  )

  test_M_smags_metrics <- list(
    sensitivity = sum((M_smags_pred_test > M_smags_metrics$threshold) * y_test) / sum(y_test),
    specificity = sum((M_smags_pred_test <= M_smags_metrics$threshold) * (1 - y_test)) / sum(1 - y_test)
  )

  # Calculate ROC curves for test data
  test_pred_obj_lasso <- prediction(pred_test_lasso, y_test)
  test_perf_obj_lasso <- performance(test_pred_obj_lasso, "tpr", "fpr")

  smags_test_pred_obj <- prediction(smags_pred_test, y_test)
  smags_test_perf_obj <- performance(smags_test_pred_obj, "tpr", "fpr")

  M_smags_test_pred_obj <- prediction(M_smags_pred_test, y_test)
  M_smags_test_perf_obj <- performance(M_smags_test_pred_obj, "tpr", "fpr")

  # Create test plot
  test_perf_objects <- list(test_perf_obj_lasso, M_smags_test_perf_obj, smags_test_perf_obj)
  test_metrics_list <- list(test_lasso_metrics, test_M_smags_metrics, test_smags_metrics)
  test_nonzero_counts <- c(lasso_nonzero, M_smags_nonzero, smags_nonzero)

  test_plot <- create_roc_plot(
    perf_objects = test_perf_objects,
    metrics_list = test_metrics_list,
    nonzero_counts = test_nonzero_counts,
    title = "ROC Curve (Test Data)",
    SP = SP,
    lambda_lasso = lambda_lasso,
    lambda_smags = lambda_smags,
    include_specificity = FALSE
  )

  # Save plots if paths provided
  if (!is.null(save_path_a)) {
    ggsave(save_path_a, train_plot, width = 8, height = 6.5, dpi = 300)
  }
  if (!is.null(save_path_b)) {
    ggsave(save_path_b, test_plot, width = 8, height = 6.5, dpi = 300)
  }

  # Return comprehensive results
  return(list(
    plots = list(train = train_plot, test = test_plot),
    train = list(
      lr = lr_metrics,
      lasso = lasso_metrics,
      smags = smags_metrics,
      M_smags = M_smags_metrics,
      nonzero = list(
        lr = lr_nonzero,
        lasso = lasso_nonzero,
        smags = smags_nonzero,
        M_smags = M_smags_nonzero
      )
    ),
    test = list(
      lr = test_lr_metrics,
      lasso = test_lasso_metrics,
      smags = test_smags_metrics,
      M_smags = test_M_smags_metrics
    ),
    thresholds = list(
      lr = lr_metrics$threshold,
      lasso = lasso_metrics$threshold,
      smags = smags_metrics$threshold,
      M_smags = M_smags_metrics$threshold
    ),
    coefficients = list(
      lr = lr_coef,
      lasso = lasso_coef,
      smags = smags_coef,
      M_smags = M_smags_coef
    )
  ))
}

#' ROC Analysis for CancerSEEK Data
#'
#' Comprehensive analysis comparing LASSO, SMAGS, and SMAGS-LASSO methods on
#' real cancer biomarker data with bootstrap confidence intervals.
#'
#' @param X_train Training feature matrix
#' @param y_train Training binary response vector
#' @param X_test Test feature matrix
#' @param y_test Test binary response vector
#' @param SP Target specificity level
#' @param lambda_lasso Lambda for LASSO (if NA, uses cross-validation)
#' @param lambda_smags Lambda for SMAGS-LASSO
#' @param save_path_a Optional path to save training ROC plot
#' @param save_path_b Optional path to save test ROC plot
#' @return List containing detailed performance metrics with confidence intervals
#' @export
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom pROC roc auc ggroc
#' @importFrom boot boot
#' @importFrom ggplot2 ggplot aes geom_line geom_vline scale_color_manual labs theme_classic theme element_text unit ggsave
#' @examples
#' \dontrun{
#' # Load CancerSEEK data
#' data_list <- prepare_cancerseek_data("path/to/data.xlsx")
#' colorectal <- data_list[["Colorectum"]]
#'
#' # Split data
#' split <- train_test_split(colorectal$X, colorectal$y)
#'
#' # Run analysis with cross-validation for lambda selection
#' cv_results <- cv_SMAGS_LASSO(split$X_train, split$y_train, SP = 0.985)
#'
#' # Generate comprehensive results
#' results <- cancerseek_results(
#'   split$X_train, split$y_train, split$X_test, split$y_test,
#'   SP = 0.985, lambda_smags = cv_results$lambda_min
#' )
#' }
cancerseek_results <- function(X_train, y_train, X_test, y_test, SP,
                               lambda_lasso = NA, lambda_smags,
                               save_path_a = NULL, save_path_b = NULL) {

  # Helper function to calculate sensitivity at a given SP using quantile method
  calculate_metrics_at_sp <- function(predictions, y_true, SP) {
    threshold <- quantile((1 - y_true) * predictions[(1 - y_true) * predictions != 0], prob = SP)
    if(is.na(threshold)) threshold <- 0

    pred_class <- as.numeric(predictions > threshold)

    sensitivity <- sum(pred_class * y_true) / sum(y_true)
    specificity <- sum((1 - pred_class) * (1 - y_true)) / sum(1 - y_true)

    return(list(
      threshold = threshold,
      sensitivity = sensitivity,
      specificity = specificity,
      pred_class = pred_class
    ))
  }

  # Function to calculate confidence intervals for sensitivity at a given SP
  bootstrap_sensitivity_ci <- function(predictions, y_true, SP, R = 1000) {
    original_sens <- calculate_metrics_at_sp(predictions, y_true, SP)$sensitivity

    bootstrap_sensitivities <- numeric(R)

    for (i in 1:R) {
      # Sample with replacement
      indices <- sample(length(y_true), replace = TRUE)
      boot_y <- y_true[indices]
      boot_pred <- predictions[indices]

      # Calculate sensitivity for bootstrap sample
      tryCatch({
        boot_sens <- calculate_metrics_at_sp(boot_pred, boot_y, SP)$sensitivity
        bootstrap_sensitivities[i] <- boot_sens
      }, error = function(e) {
        bootstrap_sensitivities[i] <<- NA
      })
    }

    # Remove NA values
    bootstrap_sensitivities <- bootstrap_sensitivities[!is.na(bootstrap_sensitivities)]

    if (length(bootstrap_sensitivities) == 0) {
      return(list(
        sensitivity = original_sens,
        ci_lower = original_sens,
        ci_upper = original_sens
      ))
    }

    # Calculate percentile confidence intervals
    ci_lower <- quantile(bootstrap_sensitivities, 0.025, na.rm = TRUE)
    ci_upper <- quantile(bootstrap_sensitivities, 0.975, na.rm = TRUE)

    return(list(
      sensitivity = original_sens,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    ))
  }

  # Calculate AUC with confidence interval using bootstrap
  calculate_auc_with_ci <- function(predictions, y_true, R = 1000) {
    original_auc <- auc(roc(y_true, predictions, quiet = TRUE))

    bootstrap_aucs <- numeric(R)

    for (i in 1:R) {
      indices <- sample(length(y_true), replace = TRUE)
      boot_y <- y_true[indices]
      boot_pred <- predictions[indices]

      tryCatch({
        boot_roc <- roc(boot_y, boot_pred, quiet = TRUE)
        bootstrap_aucs[i] <- auc(boot_roc)
      }, error = function(e) {
        bootstrap_aucs[i] <- NA
      })
    }

    bootstrap_aucs <- bootstrap_aucs[!is.na(bootstrap_aucs)]

    if (length(bootstrap_aucs) == 0) {
      return(list(
        auc = original_auc,
        ci_lower = original_auc,
        ci_upper = original_auc
      ))
    }

    ci_lower <- quantile(bootstrap_aucs, 0.025, na.rm = TRUE)
    ci_upper <- quantile(bootstrap_aucs, 0.975, na.rm = TRUE)

    return(list(
      auc = original_auc,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    ))
  }

  # Function to create a ggplot ROC curve with simplified legend
  create_roc_plot <- function(
    roc_objects,
    metrics_list,
    nonzero_counts,
    title,
    SP,
    lambda_lasso,
    lambda_smags
  ) {

    p <- ggroc(list(
      LASSO = roc_objects$lasso,
      SMAGS = roc_objects$smags,
      "SMAGS-LASSO" = roc_objects$smags_lasso
    ), legacy.axes = TRUE, linewidth = 1.8) +
      geom_vline(xintercept = 1-SP, linetype = "dashed", color = "gray") +
      scale_color_manual(values = c("LASSO" = "green", "SMAGS" = "orange", "SMAGS-LASSO" = "red"),
                         labels = c(
                           sprintf("LASSO (Selected Biomarkers = %d, λ = %.5f)",
                                   nonzero_counts[1], lambda_lasso),
                           sprintf("SMAGS (Selected Biomarkers = %d)",
                                   nonzero_counts[2]),
                           sprintf("SMAGS-LASSO (Selected Biomarkers = %d, λ = %.5f)",
                                   nonzero_counts[3], lambda_smags)
                         )) +
      labs(
        title = title,
        x = "False positive rate",
        y = "True positive rate"
      ) +
      theme_classic(base_size = 16) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.title = element_blank(),
        legend.position = c(0.55, 0.2),
        legend.background = element_rect(fill = "white", color = "gray90"),
        legend.margin = ggplot2::margin(6,6,6,6),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        legend.spacing.x = unit(0.5, "cm"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14)
      )

    return(p)
  }

  # Fit LASSO model
  if (is.na(lambda_lasso)) {
    lasso_cv <- cv.glmnet(as.matrix(X_train), y_train, family = "binomial")
    lambda_lasso <- lasso_cv$lambda.min
    pred_train_lasso <- predict(lasso_cv, newx = as.matrix(X_train), s = "lambda.min", type = "response")
  } else {
    lasso_model <- glmnet(as.matrix(X_train), y_train, family = "binomial", lambda = lambda_lasso)
    pred_train_lasso <- predict(lasso_model, newx = as.matrix(X_train), type = "response")
  }

  # Get number of non-zero coefficients for LASSO
  if (exists("lasso_cv")) {
    lasso_coef <- coef(lasso_cv, s = "lambda.min")[-1]
  } else {
    lasso_coef <- coef(lasso_model)[-1]
  }
  lasso_nonzero <- sum(abs(lasso_coef) > 0)

  # Fit SMAGS models
  smags_result <- SMAGS_LASSO(X_train, y_train, SP = SP, lambda = lambda_smags)
  M_smags_result <- SMAGS(X_train, y_train, SP)

  # Get predictions
  smags_pred_train <- sigmoid(as.matrix(cbind(X_train, 1)) %*% as.numeric(smags_result[1:(ncol(X_train) + 1)]))
  M_smags_pred_train <- sigmoid(as.matrix(cbind(X_train, 1)) %*% as.numeric(M_smags_result[1:(ncol(X_train) + 1)]))

  # Get number of non-zero coefficients for SMAGS models
  smags_coef <- smags_result[1:ncol(X_train)]
  nonzero_param_coef_smags <- 0.05 * max(abs(as.numeric(smags_coef)))
  smags_nonzero <- sum(abs(as.numeric(smags_coef)) > nonzero_param_coef_smags)

  M_smags_coef <- M_smags_result[1:ncol(X_train)]
  M_smags_nonzero <- length(M_smags_coef)

  # Print non-zero features
  lasso_features <- colnames(X_train)[abs(lasso_coef) > 0]
  cat("\n===== LASSO Non-Zero Features =====\n")
  cat(paste(lasso_features, collapse = "\n"))

  smags_features <- colnames(X_train)[abs(as.numeric(smags_coef)) > nonzero_param_coef_smags]
  cat("\n===== SMAGS-LASSO Non-Zero Features =====\n")
  cat(paste(smags_features, collapse = "\n"))

  # Calculate metrics and confidence intervals for training data
  lasso_metrics <- calculate_metrics_at_sp(pred_train_lasso, y_train, SP)
  smags_metrics <- calculate_metrics_at_sp(smags_pred_train, y_train, SP)
  M_smags_metrics <- calculate_metrics_at_sp(M_smags_pred_train, y_train, SP)

  train_lasso_ci <- bootstrap_sensitivity_ci(pred_train_lasso, y_train, SP)
  train_M_smags_ci <- bootstrap_sensitivity_ci(M_smags_pred_train, y_train, SP)
  train_smags_ci <- bootstrap_sensitivity_ci(smags_pred_train, y_train, SP)

  train_lasso_auc <- calculate_auc_with_ci(pred_train_lasso, y_train)
  train_M_smags_auc <- calculate_auc_with_ci(M_smags_pred_train, y_train)
  train_smags_auc <- calculate_auc_with_ci(smags_pred_train, y_train)

  # Print training results
  cat("\n=== TRAINING DATA METRICS ===\n")
  cat("LASSO:\n")
  cat(sprintf("  AUC = %.3f (95%% CI: %.3f-%.3f)\n",
              train_lasso_auc$auc, train_lasso_auc$ci_lower, train_lasso_auc$ci_upper))
  cat(sprintf("  Sensitivity at SP = %.2f: %.3f (95%% CI: %.3f-%.3f)\n",
              SP, train_lasso_ci$sensitivity, train_lasso_ci$ci_lower, train_lasso_ci$ci_upper))
  cat(sprintf("  Number of non-zero features: %d\n", lasso_nonzero))

  cat("SMAGS:\n")
  cat(sprintf("  AUC = %.3f (95%% CI: %.3f-%.3f)\n",
              train_M_smags_auc$auc, train_M_smags_auc$ci_lower, train_M_smags_auc$ci_upper))
  cat(sprintf("  Sensitivity at SP = %.2f: %.3f (95%% CI: %.3f-%.3f)\n",
              SP, train_M_smags_ci$sensitivity, train_M_smags_ci$ci_lower, train_M_smags_ci$ci_upper))
  cat(sprintf("  Number of non-zero features: %d\n", M_smags_nonzero))

  cat("SMAGS-LASSO:\n")
  cat(sprintf("  AUC = %.3f (95%% CI: %.3f-%.3f)\n",
              train_smags_auc$auc, train_smags_auc$ci_lower, train_smags_auc$ci_upper))
  cat(sprintf("  Sensitivity at SP = %.2f: %.3f (95%% CI: %.3f-%.3f)\n",
              SP, train_smags_ci$sensitivity, train_smags_ci$ci_lower, train_smags_ci$ci_upper))
  cat(sprintf("  Number of non-zero features: %d\n", smags_nonzero))

  # Create ROC objects for plotting
  train_roc_lasso <- roc(y_train, pred_train_lasso, quiet = TRUE)
  train_roc_smags <- roc(y_train, M_smags_pred_train, quiet = TRUE)
  train_roc_smags_lasso <- roc(y_train, smags_pred_train, quiet = TRUE)

  train_plot <- create_roc_plot(
    roc_objects = list(
      lasso = train_roc_lasso,
      smags = train_roc_smags,
      smags_lasso = train_roc_smags_lasso
    ),
    metrics_list = list(lasso_metrics, M_smags_metrics, smags_metrics),
    nonzero_counts = c(lasso_nonzero, M_smags_nonzero, smags_nonzero),
    title = "ROC Curve (Training Data)",
    SP = SP,
    lambda_lasso = lambda_lasso,
    lambda_smags = lambda_smags
  )

  # Test data predictions and analysis
  if (exists("lasso_cv")) {
    pred_test_lasso <- predict(lasso_cv, newx = as.matrix(X_test), s = "lambda.min", type = "response")
  } else {
    pred_test_lasso <- predict(lasso_model, newx = as.matrix(X_test), type = "response")
  }

  smags_pred_test <- sigmoid(as.matrix(cbind(X_test, 1)) %*% as.numeric(smags_result[1:(ncol(X_test) + 1)]))
  M_smags_pred_test <- sigmoid(as.matrix(cbind(X_test, 1)) %*% as.numeric(M_smags_result[1:(ncol(X_test) + 1)]))

  # Calculate test metrics
  test_lasso_metrics <- calculate_metrics_at_sp(pred_test_lasso, y_test, SP)
  test_smags_metrics <- calculate_metrics_at_sp(smags_pred_test, y_test, SP)
  test_M_smags_metrics <- calculate_metrics_at_sp(M_smags_pred_test, y_test, SP)

  test_lasso_ci <- bootstrap_sensitivity_ci(pred_test_lasso, y_test, SP)
  test_M_smags_ci <- bootstrap_sensitivity_ci(M_smags_pred_test, y_test, SP)
  test_smags_ci <- bootstrap_sensitivity_ci(smags_pred_test, y_test, SP)

  test_lasso_auc <- calculate_auc_with_ci(pred_test_lasso, y_test)
  test_M_smags_auc <- calculate_auc_with_ci(M_smags_pred_test, y_test)
  test_smags_auc <- calculate_auc_with_ci(smags_pred_test, y_test)

  # Print test results
  cat("\n=== TEST DATA METRICS ===\n")
  cat("LASSO:\n")
  cat(sprintf("  AUC = %.3f (95%% CI: %.3f-%.3f)\n",
              test_lasso_auc$auc, test_lasso_auc$ci_lower, test_lasso_auc$ci_upper))
  cat(sprintf("  Sensitivity at SP = %.2f: %.3f (95%% CI: %.3f-%.3f)\n",
              SP, test_lasso_ci$sensitivity, test_lasso_ci$ci_lower, test_lasso_ci$ci_upper))

  cat("SMAGS:\n")
  cat(sprintf("  AUC = %.3f (95%% CI: %.3f-%.3f)\n",
              test_M_smags_auc$auc, test_M_smags_auc$ci_lower, test_M_smags_auc$ci_upper))
  cat(sprintf("  Sensitivity at SP = %.2f: %.3f (95%% CI: %.3f-%.3f)\n",
              SP, test_M_smags_ci$sensitivity, test_M_smags_ci$ci_lower, test_M_smags_ci$ci_upper))

  cat("SMAGS-LASSO:\n")
  cat(sprintf("  AUC = %.3f (95%% CI: %.3f-%.3f)\n",
              test_smags_auc$auc, test_smags_auc$ci_lower, test_smags_auc$ci_upper))
  cat(sprintf("  Sensitivity at SP = %.2f: %.3f (95%% CI: %.3f-%.3f)\n",
              SP, test_smags_ci$sensitivity, test_smags_ci$ci_lower, test_smags_ci$ci_upper))

  # Create test plot
  test_roc_lasso <- roc(y_test, pred_test_lasso, quiet = TRUE)
  test_roc_smags <- roc(y_test, M_smags_pred_test, quiet = TRUE)
  test_roc_smags_lasso <- roc(y_test, smags_pred_test, quiet = TRUE)

  test_plot <- create_roc_plot(
    roc_objects = list(
      lasso = test_roc_lasso,
      smags = test_roc_smags,
      smags_lasso = test_roc_smags_lasso
    ),
    metrics_list = list(test_lasso_metrics, test_M_smags_metrics, test_smags_metrics),
    nonzero_counts = c(lasso_nonzero, M_smags_nonzero, smags_nonzero),
    title = "ROC Curve (Test Data)",
    SP = SP,
    lambda_lasso = lambda_lasso,
    lambda_smags = lambda_smags
  )

  # Save plots if paths provided
  if (!is.null(save_path_a)) {
    ggsave(save_path_a, train_plot, width = 10, height = 8, dpi = 300)
  }
  if (!is.null(save_path_b)) {
    ggsave(save_path_b, test_plot, width = 10, height = 8, dpi = 300)
  }

  # Return comprehensive results
  return(list(
    plots = list(train = train_plot, test = test_plot),
    train = list(
      lasso = lasso_metrics,
      smags = smags_metrics,
      M_smags = M_smags_metrics,
      nonzero = list(
        lasso = lasso_nonzero,
        smags = smags_nonzero,
        M_smags = M_smags_nonzero
      ),
      sensitivity_cis = list(
        lasso = train_lasso_ci,
        smags = train_M_smags_ci,
        smags_lasso = train_smags_ci
      ),
      auc_values = list(
        lasso = train_lasso_auc,
        smags = train_M_smags_auc,
        smags_lasso = train_smags_auc
      )
    ),
    test = list(
      lasso = test_lasso_metrics,
      smags = test_smags_metrics,
      M_smags = test_M_smags_metrics,
      sensitivity_cis = list(
        lasso = test_lasso_ci,
        smags = test_M_smags_ci,
        smags_lasso = test_smags_ci
      ),
      auc_values = list(
        lasso = test_lasso_auc,
        smags = test_M_smags_auc,
        smags_lasso = test_smags_auc
      )
    ),
    thresholds = list(
      lasso = lasso_metrics$threshold,
      smags = smags_metrics$threshold,
      M_smags = M_smags_metrics$threshold
    ),
    coefficients = list(
      lasso = lasso_coef,
      smags = smags_coef,
      M_smags = M_smags_coef
    )
  ))
}
