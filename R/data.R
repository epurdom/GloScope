#' @title example of data frame as the input for distMat
#'
#'
#' @description A data frame contains a subset of data
#' from Allen mouse data. The data includes the mouse
#' ID (donor_label), region of the brain where the data was
#' sampled, and PCA dimenstion reduction (PC_1 to PC_50)
#'
#'
#' @format A data frame with 5000 rows and 52 variables:
#' \describe{
#'   \item{donor_label}{mouse's ID}
#'   \item{joint_region_label}{Brain region where the cells were sampled}
#'   \item{PC_1 to PC_50}{PCA dimension reduction embedding}
#' }
#'
"example_data"

