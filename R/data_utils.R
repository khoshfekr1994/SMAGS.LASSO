#' Prepare CancerSEEK Data for Analysis
#'
#' Preprocesses CancerSEEK dataset for SMAGS analysis by creating separate
#' datasets for each cancer type versus normal controls. Handles data cleaning,
#' missing value imputation, and proper formatting for downstream analysis.
#'
#' @param file_path Character string specifying the path to the CancerSEEK Excel file
#' @return A named list where each element corresponds to a cancer type and contains:
#' \itemize{
#'   \item X: Feature matrix (biomarker measurements)
#'   \item y: Binary response vector (1 = cancer, 0 = normal)
#'   \item n: Total number of samples
#'   \item n_cases: Number of cancer cases
#'   \item n_controls: Number of normal controls
#' }
#' @details
#' This function performs the following preprocessing steps:
#' \itemize{
#'   \item Removes samples with missing tumor type information
#'   \item Selects relevant columns (features 8-46 are biomarker measurements)
#'   \item Cleans data by removing asterisks and handling NA values
#'   \item Converts measurement columns to numeric type
#'   \item For breast and ovary cancers, filters to female samples only
#'   \item Replaces missing values with 0
#'   \item Creates binary outcome variable (cancer vs. normal)
#' }
#' @export
#' @importFrom readxl read_xlsx
#' @importFrom dplyr select
#' @importFrom tidyr drop_na
#' @examples
#' \dontrun{
#' # Load and preprocess CancerSEEK data
#' cancerseek_data <- prepare_cancerseek_data("path/to/cancer_seek_protein.xlsx")
#'
#' # Examine available cancer types
#' names(cancerseek_data)
#'
#' # Look at colorectal cancer data
#' colorectal <- cancerseek_data[["Colorectum"]]
#' dim(colorectal$X)  # Feature matrix dimensions
#' table(colorectal$y)  # Case/control distribution
#' }
prepare_cancerseek_data <- function(file_path) {
  # Validate file path
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Read and preprocess the CancerSEEK dataset
  tryCatch({
    cancerseek <- readxl::read_xlsx(file_path)
  }, error = function(e) {
    stop("Error reading Excel file: ", e$message)
  })

  # Initial data processing
  data <- cancerseek

  # Clean up the data and handle missing values
  data <- data %>%
    drop_na(`Tumor type`) %>%
    dplyr::select(1, 3, 6, 8:46)  # Select relevant columns

  # Remove asterisks and handle NA values
  data[] <- lapply(data, function(x) gsub("\\*{1,3}", "", x))
  data[] <- lapply(data, function(x) ifelse(x == "NA", NA, x))

  # Convert numeric columns (columns 4:42) to numeric type
  numeric_cols <- names(data)[4:42]
  for (col in numeric_cols) {
    data[[col]] <- as.numeric(data[[col]])
  }

  # Get cancer types excluding Normal
  cancer_types <- unique(data$`Tumor type`)
  cancer_types <- cancer_types[cancer_types != "Normal"]

  if (length(cancer_types) == 0) {
    stop("No cancer types found in the data")
  }

  # Create list of datasets for each cancer type
  cancerseek_list <- list()

  for (cancer in cancer_types) {
    tryCatch({
      # Special handling for breast and ovary cancers (female only)
      if (!cancer %in% c("Breast", "Ovary")) {
        cancer_data <- data[data$`Tumor type` %in% c(cancer, "Normal"),]
        cancer_data$y <- ifelse(cancer_data$`Tumor type` != "Normal", 1, 0)
      } else {
        # Filter to female samples for breast and ovary cancers
        cancer_data <- data[data$`Tumor type` %in% c(cancer, "Normal") & data$Sex == "Female",]
        cancer_data$y <- ifelse(cancer_data$`Tumor type` != "Normal", 1, 0)
      }

      # Check if we have both cases and controls
      if (sum(cancer_data$y == 1) == 0) {
        warning("No cancer cases found for ", cancer, ". Skipping.")
        next
      }
      if (sum(cancer_data$y == 0) == 0) {
        warning("No normal controls found for ", cancer, ". Skipping.")
        next
      }

      # Extract features and handle missing values
      X <- cancer_data[,4:42]
      X[is.na(X)] <- 0  # Replace missing values with 0
      y <- cancer_data$y

      # Ensure X is a matrix with proper column names
      X <- as.matrix(X)
      if (is.null(colnames(X))) {
        colnames(X) <- paste0("biomarker_", 1:ncol(X))
      }

      # Store in list with metadata
      cancerseek_list[[cancer]] <- list(
        X = X,
        y = y,
        n = nrow(X),
        n_cases = sum(y),
        n_controls = sum(y == 0)
      )

    }, error = function(e) {
      warning("Error processing cancer type ", cancer, ": ", e$message)
    })
  }

  if (length(cancerseek_list) == 0) {
    stop("No valid cancer datasets were created")
  }

  # Print summary information
  cat("Successfully processed", length(cancerseek_list), "cancer types:\n")
  for (cancer in names(cancerseek_list)) {
    data_info <- cancerseek_list[[cancer]]
    cat(sprintf("  %s: %d cases, %d controls (%d total)\n",
                cancer, data_info$n_cases, data_info$n_controls, data_info$n))
  }

  return(cancerseek_list)
}

#' Generate Cancer-Specificity Parameter Values
#'
#' Returns the recommended specificity (SP) values for different cancer types
#' based on CancerSEEK publication standards.
#'
#' @return Named list of specificity values for each cancer type
#' @export
#' @examples
#' # Get recommended SP values
#' sp_values <- get_cancer_sp_values()
#' print(sp_values)
get_cancer_sp_values <- function() {
  list(
    "Colorectum" = 0.985,
    "Breast" = 0.86,
    "Lung" = 0.983,
    "Pancreas" = 0.985,
    "Liver" = 0.988,
    "Stomach" = 0.96
  )
}

#' Split Data into Training and Testing Sets
#'
#' Creates stratified train/test splits for binary classification problems,
#' ensuring balanced class representation in both sets.
#'
#' @param X Feature matrix
#' @param y Binary response vector
#' @param train_prop Proportion of data for training (default: 0.8)
#' @param seed Random seed for reproducibility (default: 42)
#' @return List containing X_train, X_test, y_train, y_test
#' @export
#' @examples
#' \dontrun{
#' # Generate example data
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' y <- rbinom(100, 1, 0.5)
#'
#' # Create train/test split
#' split_data <- train_test_split(X, y, train_prop = 0.7)
#'
#' # Access the splits
#' X_train <- split_data$X_train
#' y_train <- split_data$y_train
#' X_test <- split_data$X_test
#' y_test <- split_data$y_test
#' }
train_test_split <- function(X, y, train_prop = 0.8, seed = 42) {
  set.seed(seed)

  # Validate inputs
  if (nrow(X) != length(y)) {
    stop("Number of rows in X must equal length of y")
  }
  if (train_prop <= 0 || train_prop >= 1) {
    stop("train_prop must be between 0 and 1")
  }

  n <- length(y)

  # For binary outcomes, try stratified sampling
  if (length(unique(y)) == 2) {
    # Get indices for each class
    class_0_indices <- which(y == 0)
    class_1_indices <- which(y == 1)

    # Sample from each class
    n_train_0 <- round(length(class_0_indices) * train_prop)
    n_train_1 <- round(length(class_1_indices) * train_prop)

    train_0_indices <- sample(class_0_indices, n_train_0)
    train_1_indices <- sample(class_1_indices, n_train_1)

    train_indices <- c(train_0_indices, train_1_indices)
  } else {
    # Simple random sampling for non-binary outcomes
    n_train <- round(n * train_prop)
    train_indices <- sample(1:n, n_train)
  }

  test_indices <- setdiff(1:n, train_indices)

  return(list(
    X_train = X[train_indices, , drop = FALSE],
    X_test = X[test_indices, , drop = FALSE],
    y_train = y[train_indices],
    y_test = y[test_indices],
    train_indices = train_indices,
    test_indices = test_indices
  ))
}
