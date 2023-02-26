

#' Set up data for \pkg{bsitar} 
#'
#' @param data An input data frame
#' @param response_variable The name(s) of the outcome(s) variables. Must be a
#'   single name except when fitting a multivariate model.
#' @param uvarby An optional (default \code{NA}) to specify the indicator
#'   variable for fitting univariate-by-subgroup model. See \code{univariate_by}
#'   argument in the [bsitar::bsitar()] function. If not \code{NA}, then it
#'   should be a valid factor variable present in the \code{data}.
#' @param mvar A logical (default \code{FALSE}) to specify the the multivariate
#'   model. See \code{multivariate} argument in the [bsitar::bsitar()] function.
#' @param xfuns Optional name(s) of the transformation function(s) applied to
#'   the predictor variable (typically age). Default \code{NULL}.
#' @param yfuns Optional name(s) of the transformation function(s) applied to
#'   the outcome variable. Default \code{NULL}.
#'
#' @return A data frame with necessary information added a attributes.
#' 
#' @export
#'
#' @examples
#' \dontrun {
#' data(heights)
#' make_bsitar_data(data = heights, response_variable = heights)
#' }
#' 
make_bsitar_data <- function(data, response_variable, 
                             uvarby = NA, mvar = FALSE, 
                             xfuns = NULL, yfuns = NULL) {
  
  org.data <- data
  if (!(is.na(uvarby) | uvarby == "NA")) {
    if (!uvarby %in% colnames(data)) {
      stop(
        paste(
          "\nvariable",
          uvarby,
          "used for setting univariate submodels is missing"
        )
      )
    }
    if (!is.factor(data[[uvarby]])) {
      stop("subset by variable '",
           uvarby,
           "' should be a factor variable")
    }
    for (l in levels(data[[uvarby]])) {
      data[[l]] <- data[[response_variable[1]]]
    }
    unibyimat <-
      model.matrix( ~ 0 + eval(parse(text = uvarby)), data)
    subindicators <- paste0(uvarby, levels(data[[uvarby]]))
    colnames(unibyimat) <- subindicators
    response_variable <- levels(data[[uvarby]])
    data <- as.data.frame(cbind(data, unibyimat))
    attr(data, "ys") <- response_variable
    attr(data, "multivariate") <- FALSE
    attr(data, "uvarby") <- uvarby
    attr(data, "subindicators") <- subindicators
    # data_out <- data
  } else if(mvar) {
    for(myfunsi in 1:length(response_variable)) {
      mysi <- response_variable[[myfunsi]]
      myfunsi <- yfuns[[myfunsi]]
      if(grepl('.Primitive', myfunsi, fixed = T) & 
         grepl('log', myfunsi, fixed = T)) {
        myfunsi <- 'log'
      }
      if(grepl('.Primitive', myfunsi, fixed = T) & 
         grepl('sqrt', myfunsi, fixed = T)) {
        myfunsi <- 'sqrt'
      }
      if(myfunsi == 'log') data[[mysi]] <- log(data[[mysi]])
      if(myfunsi == 'sqrt') data[[mysi]] <- sqrt(data[[mysi]])
    }
    #  data_out <- org.data
    attr(data, "ys") <- response_variable
    attr(data, "multivariate") <- TRUE
    attr(data, "uvarby") <- NULL
    attr(data, "subindicators") <- NULL
  } else {
    # data_out <- org.data
    data <- org.data
    attr(data, "ys") <- response_variable
    attr(data, "multivariate") <- FALSE
    attr(data, "uvarby") <- NULL
    attr(data, "subindicators") <- NULL
  }
  return(data)
}

