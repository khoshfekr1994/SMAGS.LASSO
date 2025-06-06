# SMAGS.LASSO

> Sensitivity Maximization at a Given Specificity with LASSO Regularization

[![R](https://img.shields.io/badge/R-276DC3?style=flat&logo=r&logoColor=white)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**SMAGS.LASSO** implements algorithms for biomarker discovery and classification that maximize sensitivity at a user-specified specificity level. The package is particularly designed for cancer biomarker analysis and includes both standard SMAGS and LASSO-regularized versions for feature selection.

### Key Features

- 🎯 **SMAGS**: Core sensitivity maximization algorithm at given specificity
- 🔀 **SMAGS-LASSO**: SMAGS with LASSO regularization for automatic feature selection  
- 📊 **Cross-validation**: Automated lambda parameter selection with MSE optimization
- 📈 **Visualization**: ROC curves and comprehensive performance comparisons
- 🧬 **Cancer data support**: Built-in CancerSEEK dataset preprocessing tools
- 📊 **Bootstrap confidence intervals**: Statistical inference for real-world applications

## Installation

Install the package directly from GitHub:

```r
# Install devtools if you haven't already
if (!require(devtools)) install.packages("devtools")

# Install SMAGS.LASSO
devtools::install_github("yourusername/SMAGS.LASSO")
```

### Dependencies

The package automatically installs required dependencies:
- `glmnet`, `ROCR`, `pROC` for machine learning and ROC analysis
- `ggplot2` for visualization
- `dplyr`, `tidyr` for data manipulation
- `readxl` for Excel file processing
- `parallel`, `doParallel`, `foreach` for computational efficiency

## Quick Start

```r
library(SMAGS.LASSO)

# Generate example biomarker data
set.seed(123)
n <- 200
p <- 20
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("Biomarker_", 1:p)

# Create binary outcome (1 = disease, 0 = normal)
y <- rbinom(n, 1, 0.4)

# Basic SMAGS analysis (95% specificity)
result <- SMAGS(X, y, SP = 0.95)
print(paste("Achieved Sensitivity:", round(result$Sensitivity, 3)))

# SMAGS with LASSO regularization for feature selection
result_lasso <- SMAGS_LASSO(X, y, SP = 0.95, lambda = 0.1)

# Cross-validation for optimal lambda selection
cv_result <- cv_SMAGS_LASSO(X, y, SP = 0.95)
print(paste("Optimal lambda:", round(cv_result$lambda_min, 4)))
```

## Core Functions

| Function | Description | Key Parameters |
|----------|-------------|----------------|
| `SMAGS()` | Core sensitivity maximization | `X`, `y`, `SP` |
| `SMAGS_LASSO()` | SMAGS with LASSO regularization | `X`, `y`, `SP`, `lambda` |
| `cv_SMAGS_LASSO()` | Cross-validation for lambda selection | `X`, `y`, `SP`, `nfolds` |
| `synthetic_data_results()` | Compare methods on synthetic data | Training/test sets, lambdas |
| `cancerseek_results()` | Comprehensive real-data analysis | With confidence intervals |
| `prepare_cancerseek_data()` | Preprocess CancerSEEK datasets | Excel file path |
| `train_test_split()` | Stratified data splitting | `X`, `y`, `train_prop` |

## Example: Cancer Biomarker Analysis

### Synthetic Cancer Data Analysis

```r
# Simulate realistic cancer biomarker scenario
set.seed(42)
n_cases <- 150
n_controls <- 150
p_biomarkers <- 39

# Cancer samples: elevated biomarker levels
X_cancer <- matrix(rnorm(n_cases * p_biomarkers, mean = 1.2, sd = 1), 
                   n_cases, p_biomarkers)
# Normal samples: baseline levels  
X_normal <- matrix(rnorm(n_controls * p_biomarkers, mean = 0, sd = 1), 
                   n_controls, p_biomarkers)

# Combine datasets
X <- rbind(X_cancer, X_normal)
y <- c(rep(1, n_cases), rep(0, n_controls))
colnames(X) <- paste0("Protein_", 1:p_biomarkers)

# Split into training and testing sets
split_data <- train_test_split(X, y, train_prop = 0.8)

# Find optimal regularization parameter
cv_results <- cv_SMAGS_LASSO(split_data$X_train, split_data$y_train, SP = 0.98)

# Comprehensive analysis with visualization
results <- cancerseek_results(
  X_train = split_data$X_train,
  y_train = split_data$y_train,
  X_test = split_data$X_test,
  y_test = split_data$y_test,
  SP = 0.98,  # High specificity for cancer screening
  lambda_smags = cv_results$lambda_min
)

# Display ROC curves
print(results$plots$train)  # Training performance
print(results$plots$test)   # Test performance
```

### Real CancerSEEK Data Analysis

```r
# Load and preprocess CancerSEEK data (requires actual dataset)
# cancerseek_data <- prepare_cancerseek_data("path/to/cancer_seek_protein.xlsx")

# Get recommended specificity values for different cancer types
sp_values <- get_cancer_sp_values()
print(sp_values)
#> $Colorectum: 0.985
#> $Breast: 0.86  
#> $Lung: 0.983
#> $Pancreas: 0.985
#> $Liver: 0.988
#> $Stomach: 0.96
```

## Advanced Usage

### Cross-Validation and Model Selection

```r
# Generate high-dimensional sparse data
set.seed(123)
n <- 300
p <- 100
X <- matrix(rnorm(n * p), n, p)
# Only first 10 features are informative
y <- as.numeric(rowSums(X[, 1:10]) + rnorm(n, 0, 2) > 0)

# Comprehensive cross-validation
cv_results <- cv_SMAGS_LASSO(X, y, SP = 0.90, nfolds = 10)

# Plot cross-validation curve
cv_results$plot()

# Examine selected features
final_model <- cv_results$final_model
selected_features <- which(abs(as.numeric(final_model[1:p])) > 1e-6)
print(paste("Selected features:", paste(selected_features, collapse = ", ")))
```

### Performance Comparison

```r
# Compare LASSO vs SMAGS vs SMAGS-LASSO
comparison <- synthetic_data_results(
  X_train = split_data$X_train,
  y_train = split_data$y_train,
  X_test = split_data$X_test,
  y_test = split_data$y_test,
  SP = 0.95,
  lambda_lasso = 0.01,
  lambda_smags = cv_results$lambda_min
)

# Training vs Test Performance
cat("Training Sensitivity:\n")
cat("  LASSO:", round(comparison$train$lasso$sensitivity, 3), "\n")
cat("  SMAGS:", round(comparison$train$M_smags$sensitivity, 3), "\n")
cat("  SMAGS-LASSO:", round(comparison$train$smags$sensitivity, 3), "\n")

cat("\nTest Sensitivity:\n")
cat("  LASSO:", round(comparison$test$lasso$sensitivity, 3), "\n")
cat("  SMAGS:", round(comparison$test$M_smags$sensitivity, 3), "\n")
cat("  SMAGS-LASSO:", round(comparison$test$smags$sensitivity, 3), "\n")
```

## Methodology

### SMAGS Algorithm

SMAGS (Sensitivity Maximization at a Given Specificity) optimizes the following objective:

**Maximize**: Sensitivity = TP / (TP + FN)  
**Subject to**: Specificity ≥ SP

Where the threshold is determined by the quantile method to achieve the target specificity.

### SMAGS-LASSO

SMAGS-LASSO adds L1 regularization for automatic feature selection:

**Minimize**: -Sensitivity + λ||β||₁

This combines sensitivity optimization with sparsity-inducing regularization for interpretable biomarker panels.

### Cross-Validation

The package uses k-fold cross-validation to select optimal λ by minimizing sensitivity mean squared error:

**MSE** = (1 - Sensitivity)²

## Applications

- **Cancer Screening**: Early detection biomarker panels
- **Diagnostic Testing**: High-specificity clinical tests  
- **Genomics**: Gene expression analysis with thousands of features
- **Proteomics**: Protein biomarker discovery
- **Clinical Research**: Biomarker validation studies
- **Personalized Medicine**: Patient stratification

## Performance Metrics

The package provides comprehensive evaluation:

- **Sensitivity/Specificity** at user-defined operating points
- **AUC** with bootstrap confidence intervals
- **Feature selection** reporting and visualization
- **ROC curves** with publication-ready plots
- **Cross-validation** error curves
- **Model comparison** across different methods

## Citation

If you use SMAGS.LASSO in your research, please cite:

```bibtex
@misc{smagslasso2024,
  title={SMAGS.LASSO: Sensitivity Maximization at a Given Specificity with LASSO Regularization},
  author={[Your Name]},
  year={2024},
  url={https://github.com/yourusername/SMAGS.LASSO}
}
```

## Requirements

- R ≥ 3.5.0
- Required packages are automatically installed

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to:
- Report bugs via [Issues](https://github.com/yourusername/SMAGS.LASSO/issues)
- Suggest enhancements
- Submit pull requests

### Development

```r
# Clone and install development version
devtools::install_github("yourusername/SMAGS.LASSO", ref = "develop")

# Run tests
devtools::test()

# Check package
devtools::check()
```

## Support

- **Documentation**: Function help with `?function_name`
- **Issues**: [GitHub Issues](https://github.com/khoshfekr1994/SMAGS.LASSO/issues)
- **Email**: [khoshfekr1994@gmail.com] & [hkhoshfekr@mdanderson.org]

## Developed by

**[Hamid Khoshfekr Rudsari, PhD]** 


