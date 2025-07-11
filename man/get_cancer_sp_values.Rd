% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_utils.R
\name{get_cancer_sp_values}
\alias{get_cancer_sp_values}
\title{Generate Cancer-Specificity Parameter Values}
\usage{
get_cancer_sp_values()
}
\value{
Named list of specificity values for each cancer type
}
\description{
Returns the recommended specificity (SP) values for different cancer types
based on CancerSEEK publication standards.
}
\examples{
# Get recommended SP values
sp_values <- get_cancer_sp_values()
print(sp_values)
}
