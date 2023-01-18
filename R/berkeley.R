

#' The Berkeley Child Guidance Study
#'
#' @description The Berkeley Child Guidance Study dataset contains longitudinal
#'   anthropometry data for 136 children from birth to 21 years.
#'
#' @details The data were originaly provided as an appendix to the book by
#'   \insertCite{Tuddenham1954;textual}{bsitar}. The same data has been used in
#'   the *sitar* \insertCite{R-sitar}{bsitar} and the *fda*
#'   \insertCite{R-fda}{bsitar} packages after correcting the transcription
#'   errors. The data included in the **bsitar** package is 
#'   the same as \code{berkeley} data used in the *sitar* package.
#'
#'   A detailed description of the data including the frequency of measurements
#'   per year is provided in the the *sitar* package.
#'   \insertCite{R-sitar}{bsitar}. Briefly, the data comprise of repeated growth
#'   measurements made on 66 boys and 70 girls (birth to 21 years). Children
#'   were born in 1928-29 (Berkeley, California) and were of north European
#'   ancestry. Measurements were made at ages 0, 0.085, 0.25 to 2 (3-monthly), 2
#'   to 8 (annually), and 8 to 21 (6-monthly) years. The children were measured
#'   for height, weight (undressed), stem length, biacromial diameter, bi-iliac
#'   diameter, leg circumference, and dynamometric strength.
#' 
#' @name berkeley
#' @docType data
#' @format A data frame with 4884 observations on the following 10 variables:
#' \describe{
#' \item{id}{factor with levels 201-278 male and 301-385 female}
#' \item{age}{years, numeric vector}
#' \item{height}{cm, numeric vector}
#' \item{weight}{kg, numeric vector}
#' \item{stem.length}{cm, numeric vector}
#' \item{bi.acromial}{cm, numeric vector}
#' \item{bi.iliac}{cm, numeric vector}
#' \item{leg.circ}{cm, numeric vector}
#' \item{strength}{lb, numeric vector}
#' \item{sex}{factor with level 1 male and level 2 female}
#' }
#' @references
#'  \insertAllCited{}
#'  
#' @keywords datasets
#' 
#' @return A data frame with 10 columns.
"berkeley"