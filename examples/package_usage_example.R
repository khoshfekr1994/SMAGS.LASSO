# After installing the package with devtools::install_github("khoshfekr1994/SMAGS.LASSO")
library(SMAGS.LASSO)

# =====================================
# Example 1: Basic SMAGS Analysis
# =====================================

# Generate synthetic data
set.seed(123)
n <- 200
p <- 20
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("Feature_", 1:p)

# Create informative features (first 5 features are signal)
y <- as.numeric(rowSums(X[, 1:5]) + rnorm(n, 0, 0.5) > 0)

# Split into train/test
split_data <- train_test_split(X, y, train_prop = 0.8)

# Run basic SMAGS
smags_result <- SMAGS(split_data$X_train, split_data$y_train, SP = 0.95)
print("SMAGS Result:")
print(smags_result)

# Run SMAGS-LASSO with fixed lambda
smags_lasso_result <- SMAGS_LASSO(split_data$X_train, split_data$y_train,
                                  SP = 0.95, lambda = 0.1)
print("SMAGS-LASSO Result:")
print(smags_lasso_result)

# =====================================
# Example 2: Cross-Validation for Optimal Lambda
# =====================================

# Use cross-validation to find optimal lambda
cv_results <- cv_SMAGS_LASSO(split_data$X_train, split_data$y_train, SP = 0.95)

# Plot cross-validation results
cv_results$plot()

# Get optimal lambda
optimal_lambda <- cv_results$lambda_min
cat("Optimal lambda:", optimal_lambda, "\n")

# =====================================
# Example 3: Comprehensive Analysis
# =====================================

# Compare all methods with visualization
comparison_results <- synthetic_data_results(
  X_train = split_data$X_train,
  y_train = split_data$y_train,
  X_test = split_data$X_test,
  y_test = split_data$y_test,
  SP = 0.95,
  lambda_lasso = 0.01,  # Fixed lambda for LASSO
  lambda_smags = optimal_lambda  # Optimal lambda from CV
)

# Display the ROC plots
print(comparison_results$plots$train)
print(comparison_results$plots$test)

# Print performance summary
cat("\n=== PERFORMANCE SUMMARY ===\n")
cat("Training Sensitivity:\n")
cat("  LASSO:", comparison_results$train$lasso$sensitivity, "\n")
cat("  SMAGS:", comparison_results$train$M_smags$sensitivity, "\n")
cat("  SMAGS-LASSO:", comparison_results$train$smags$sensitivity, "\n")

cat("\nTest Sensitivity:\n")
cat("  LASSO:", comparison_results$test$lasso$sensitivity, "\n")
cat("  SMAGS:", comparison_results$test$M_smags$sensitivity, "\n")
cat("  SMAGS-LASSO:", comparison_results$test$smags$sensitivity, "\n")

# =====================================
# Example 4: CancerSEEK Data Analysis
# =====================================

# Note: This requires actual CancerSEEK data file
# cancerseek_data <- prepare_cancerseek_data("path/to/cancer_seek_protein.xlsx")

# Example with simulated cancer-like data
set.seed(42)
n_cancer <- 100
n_normal <- 100
p_biomarkers <- 39  # Similar to CancerSEEK

# Simulate cancer vs normal data
X_cancer <- matrix(rnorm(n_cancer * p_biomarkers, mean = 2, sd = 1), n_cancer, p_biomarkers)
X_normal <- matrix(rnorm(n_normal * p_biomarkers, mean = 0, sd = 1), n_normal, p_biomarkers)

X_combined <- rbind(X_cancer, X_normal)
y_combined <- c(rep(1, n_cancer), rep(0, n_normal))

colnames(X_combined) <- paste0("Biomarker_", 1:p_biomarkers)

# Split data
cancer_split <- train_test_split(X_combined, y_combined, train_prop = 0.8)

# Get recommended SP value for cancer analysis
sp_values <- get_cancer_sp_values()
target_sp <- 0.98  # High specificity for cancer screening

# Run comprehensive cancer analysis
cancer_results <- cancerseek_results(
  X_train = cancer_split$X_train,
  y_train = cancer_split$y_train,
  X_test = cancer_split$X_test,
  y_test = cancer_split$y_test,
  SP = target_sp,
  lambda_lasso = NA,  # Use cross-validation
  lambda_smags = 0.1
)

# Display results (will print detailed metrics with confidence intervals)

# =====================================
# Example 5: Feature Selection Analysis
# =====================================

# Examine which features were selected by each method
cat("\n=== FEATURE SELECTION COMPARISON ===\n")

# Get coefficients from SMAGS-LASSO
smags_coefs <- as.numeric(smags_lasso_result[1:p])
selected_features <- which(abs(smags_coefs) > 1e-6)

cat("SMAGS-LASSO selected features:", selected_features, "\n")
cat("True signal features: 1, 2, 3, 4, 5\n")

# Calculate overlap with true signal
overlap <- intersect(selected_features, 1:5)
cat("Correctly identified signal features:", overlap, "\n")
cat("Sensitivity for feature selection:", length(overlap) / 5, "\n")

# =====================================
# Example 6: Sensitivity to SP Values
# =====================================

# Test different SP values
sp_values_test <- c(0.90, 0.95, 0.98, 0.99)
sensitivity_results <- data.frame(
  SP = sp_values_test,
  SMAGS_Sensitivity = numeric(length(sp_values_test)),
  SMAGS_LASSO_Sensitivity = numeric(length(sp_values_test))
)

for (i in seq_along(sp_values_test)) {
  sp <- sp_values_test[i]

  # Run SMAGS
  smags_temp <- SMAGS(split_data$X_train, split_data$y_train, SP = sp)
  sensitivity_results$SMAGS_Sensitivity[i] <- as.numeric(smags_temp["Sensitivity"])

  # Run SMAGS-LASSO
  smags_lasso_temp <- SMAGS_LASSO(split_data$X_train, split_data$y_train,
                                  SP = sp, lambda = 0.1)
  sensitivity_results$SMAGS_LASSO_Sensitivity[i] <- as.numeric(smags_lasso_temp["Sensitivity"])
}

print("Sensitivity vs. Specificity Trade-off:")
print(sensitivity_results)

# Plot the trade-off
if (require(ggplot2)) {
  library(reshape2)
  plot_data <- melt(sensitivity_results, id.vars = "SP")

  p <- ggplot(plot_data, aes(x = SP, y = value, color = variable)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    labs(
      title = "Sensitivity vs. Specificity Trade-off",
      x = "Target Specificity (SP)",
      y = "Achieved Sensitivity",
      color = "Method"
    ) +
    theme_classic() +
    scale_color_manual(values = c("SMAGS_Sensitivity" = "orange",
                                  "SMAGS_LASSO_Sensitivity" = "red"))

  print(p)
}

cat("\nExample analysis complete!\n")
cat("The SMAGS.LASSO package provides comprehensive tools for:\n")
cat("1. Sensitivity maximization at given specificity\n")
cat("2. Feature selection with LASSO regularization\n")
cat("3. Cross-validation for parameter selection\n")
cat("4. Comprehensive performance evaluation\n")
cat("5. Visualization and comparison tools\n")
