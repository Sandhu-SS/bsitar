
#' Model fit berkeley dataset
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
#' @name berkeley_fit
#' @docType data
#' @format A data frame with 124 observations on the following 4 variables:
#' \describe{
#' \item{id}{factor of subject ids (levels 1:12).}
#' \item{age}{vector of ages (years).}
#' \item{height}{vector of heights (cm).}
#' \item{men}{vector of ages at menarche (years), where negative values
#' are right censored.} }
#' 
#' @references
#'  \insertAllCited{}
#' 
#' @keywords datasets
#' @return A object of class \code{bsitar} with posterior draws.
"berkeley_fit"