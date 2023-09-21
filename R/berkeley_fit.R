
#' Model fit to the Berkeley datas et
#'
#' @description Model is fit to the Berkeley Child Guidance Study dataset
#'   containing longitudinal anthropometry data for 136 children from birth to
#'   21 years.
#'
#' @details Data details are provided in the the *sitar* package
#'   \insertCite{R-sitar}{bsitar}. 
#'
#'
#' @name berkeley_fit
#' @docType data
#' @format The Berkeley data analysed have the following 5 variables:
#' \describe{
#' \item{id}{factor with levels 201-278 male and 301-385 female}
#' \item{age}{years, numeric vector}
#' \item{height}{cm, numeric vector}
#' \item{sex}{factor with level 1 male and level 2 female}
#' }
#' 
#' @references
#'  \insertAllCited{}
#' 
#' @keywords datasets
#' @return A object of class \code{bgmfit} with posterior draws.
"berkeley_fit"