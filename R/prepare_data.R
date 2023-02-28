

#' Set up data for bsitar model
#'
#' @param data An input data frame
#' @param x The predictor (typically age) variables.
#' @param y The outcome variables. Must be a single name except when fitting a
#'   multivariate model.
#' @param id The group identifier.
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
#' @param outliers An optional (default \code{NULL}) to remove velocity
#'   outliers. The argument should be a named list to pass options to the
#'   [bsitar::outliers] function. See [bsitar::outliers] for details.
#'
#' @return A data frame with necessary information added a attributes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(heights)
#' prepare_data(data = heights, y = heights)
#' }
#' 
prepare_data <- function(data, 
                         x,
                         y,
                         id,
                         uvarby = NA, 
                         mvar = FALSE, 
                         xfuns = NULL, 
                         yfuns = NULL,
                         outliers = NULL) {
  
  
  if(!is.null(outliers)) {
    if(is.null(outliers$remove))    outliers$remove <- TRUE
    if(is.null(outliers$icode))     outliers$icode <- c(4,5,6)
    if(is.null(outliers$limit))     outliers$limit <- 5
    if(is.null(outliers$velpower))  outliers$velpower <- 0.5
    if(is.null(outliers$lag))       outliers$lag <- 1
    if(is.null(outliers$linearise)) outliers$linearise <- FALSE
    if(is.null(outliers$verbose))   outliers$verbose <- FALSE
    
    remove_ <- outliers$remove
    icode_ <- outliers$icode
    icode_ <- deparse(substitute(icode_))
    limit_ <- outliers$limit
    velpower_ <- outliers$velpower
    lag_ <- outliers$lag
    linearise_ <- outliers$linearise
    verbose_ <- outliers$verbose
    data_ <- deparse(substitute(data))
    for(yi in 1:length(y)) {
      exe_outlier <- paste0("outliers(x = ", x, "," , 
                            "y = ", y[yi], "," , 
                            "id = ", id, "," , 
                            "data = ", data_, ", " ,
                            "icode = ", icode_, "," , 
                            "lag = ", lag_, "," , 
                            "velpower = ", velpower_,  "," , 
                            "limit = ", limit_, ", " , 
                            "linearise = ", linearise_,  ", " ,
                            "remove = ", remove_,  ", " ,
                            "verbose = ", verbose_, ")")
      
      data <- eval(parse(text = exe_outlier))
    }
  } # if(!is.null(outliers)) {
  
  
  
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
      data[[l]] <- data[[y[1]]]
    }
    unibyimat <-
      model.matrix( ~ 0 + eval(parse(text = uvarby)), data)
    subindicators <- paste0(uvarby, levels(data[[uvarby]]))
    colnames(unibyimat) <- subindicators
    y <- levels(data[[uvarby]])
    data <- as.data.frame(cbind(data, unibyimat))
    attr(data, "ys") <- y
    attr(data, "multivariate") <- FALSE
    attr(data, "uvarby") <- uvarby
    attr(data, "subindicators") <- subindicators
    # data_out <- data
  } else if(mvar) {
    for(myfunsi in 1:length(y)) {
      mysi <- y[[myfunsi]]
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
    attr(data, "ys") <- y
    attr(data, "multivariate") <- TRUE
    attr(data, "uvarby") <- NULL
    attr(data, "subindicators") <- NULL
  } else {
    # data_out <- org.data
    data <- org.data
    attr(data, "ys") <- y
    attr(data, "multivariate") <- FALSE
    attr(data, "uvarby") <- NULL
    attr(data, "subindicators") <- NULL
  }
  return(data)
}



###########


# dx <- prepare_data(data, 
#                    x = 'age',
#                    y = c('copad', 'copod'), 
#                    id = 'id',
#                    uvarby = NA, 
#                    outliers = list(remove =  TRUE, icode = c(4:6), 
#                                    limit = 5, velpower = 0.5, 
#                                    lag = 1, linearise = FALSE, 
#                                    verbose = FALSE),
#                    mvar = FALSE
#                    )
# 
# 
# head(dx)
# nrow(dx)




