

#' Berkeley Child Guidance Study Data for males
#'
#' @description A subset of the Berkeley Child Guidance Study data
#'   \strong{berkeley} that contains longitudinal growth data for 66 males (6 to
#'   20 years of age).
#'
#' @details A detailed description of the full data including the frequency of
#'   measurements per year is provided in the [bsitar::berkeley] dataset.
#'   
#' @name berkeley_mdata
#' @docType data
#' @format A data frame with 902 observations on the following 3 variables:
#' \describe{
#' \item{id}{factor with levels 201-278}
#' \item{age}{years, numeric vector}
#' \item{height}{cm, numeric vector}
#' }
#' @references
#'  \insertAllCited{}
#'  
#' @keywords datasets
#' 
#' @return A data frame with 3 columns.
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
"berkeley_mdata"