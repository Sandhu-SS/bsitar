
#' Model fit to serial measurements of height on 12 girls
#'
#' @description Model is fit to data comprise of height measurements of 12   
#' girls from the Chard Growth Study measured twice a year between 8 and 
#' 16 years of age.
#'
#' @details Data details are provided in the the *sitar* package
#'   \insertCite{R-sitar}{bsitar}. The data included in the **bsitar** package is 
#'   the same as \code{heights} data used in the *sitar* package.
#'
#'
#' @name heights_fit
#' @docType data
#' @format Data frame analysed have the following 3 variables:
#' \describe{
#' \item{id}{factor of subject ids (levels 1:12)}
#' \item{age}{vector of ages (years)}
#' \item{height}{vector of heights (cm)}
#' }
#' 
#' @references
#'  \insertAllCited{}
#' 
#' @keywords datasets
#' @return A object of class \code{bsitar} with posterior draws.
"heights_fit"