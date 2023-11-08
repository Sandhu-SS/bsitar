

#' Berkeley Child Guidance Study Data for males
#'
#' @description The Berkeley Child Guidance Study data for males is subset of
#'   the Berkeley Child Guidance Study data \code{berkeley} that contains
#'   longitudinal growth data for 136 children (66 males and 70 females) from
#'   birth to 21 years. The subset data \code{berkeley_mdata} includes growth
#'   measurements on males between 6 and 20 years of age.
#'
#' @details The data were originaly provided as an appendix to the book by
#'   \insertCite{Tuddenham1954;textual}{bsitar}. The same data has been used in
#'   the *sitar* \insertCite{R-sitar}{bsitar} and the *fda*
#'   \insertCite{R-fda}{bsitar} packages after correcting the transcription
#'   errors. The data included in the **bsitar** package is the same as
#'   \code{berkeley} data used in the *sitar* package.
#'
#'   A detailed description of the data including the frequency of measurements
#'   per year is provided in the the *sitar* package.
#'   \insertCite{R-sitar}{bsitar} and the [bsitar::berkeley_mdata] included in
#'   this \code{bsitar} package.
#'   
#' @name berkeley_mdata
#' @docType data
#' @format A data frame with 902 observations on the following 10 variables:
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
#' @return A data frame with 10 columns.
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
"berkeley_mdata"