

#' An internal function to set up data for \code{bsitar} model
#'
#' @param data An input data frame.
#'
#' @param xvar The predictor (typically age) variables.
#'
#' @param yvar The outcome variables. Must be a single name except when fitting a
#'   multivariate model.
#'   
#' @param sigmaxvar The predictor (typically age) variables for \code{sigma}
#'
#' @param xoffset A real number.
#'
#' @param xfuns Optional name(s) of the transformation function(s) applied to
#'   the predictor variable (typically age). Default \code{NULL}.
#'
#' @param yfuns Optional name(s) of the transformation function(s) applied to
#'   the outcome variable. Default \code{NULL}.
#'   
#' @param sigmaxfuns Optional name(s) of the transformation function(s) applied to
#'   the predictor variable (typically age)  for \code{sigma}. Default \code{NULL}.
#'  
#' @param ixfuns Optional name(s) of the inverse transformation function(s) applied to
#'   the predictor variable (typically age). Default \code{NULL}.
#'
#' @param iyfuns Optional name(s) of the inverse transformation function(s) applied to
#'   the outcome variable. Default \code{NULL}.
#'   
#' @param isigmaxfuns Optional name(s) of the inverse transformation function(s) applied to
#'   the predictor variable (typically age)  for \code{sigma}. Default \code{NULL}. 
#'   
#' @param transform A character vector to specify variables to be transformed.
#' 
#' @param itransform A character vector to specify variables to be inverse transformed.
#'   
#' @param envir A logical (default \code{TRUE})
#'
#' @return A data frame with necessary information added a attributes.
#'
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @keywords internal
#' @noRd
#'
prepare_transformations <- function(data = NULL,
                                    xvar = NULL, 
                                    yvar = NULL,
                                    sigmaxvar = NULL,
                                    xoffset = NULL,
                                    xfuns = NULL,
                                    yfuns = NULL,
                                    sigmaxfuns = NULL,
                                    ixfuns = NULL,
                                    iyfuns = NULL,
                                    isigmaxfuns = NULL,
                                    transform = c('x', 'y', 'sigma'),
                                    itransform = c('x', 'y', 'sigma'),
                                    model = NULL,
                                    envir = NULL,
                                    verbose = FALSE) {
  

  if(is.null(data) & is.null(model)) {
    stop("specify at least one of the data or model")
  } else if(!is.null(model)) {
    if(is.null(data)) {
      data <- model$data
      if(verbose) message("'data' is extracted from the 'model'")
    } 
    if(is.null(xvar)) {
      xvar <- model$model_info$xvar
      if(verbose) message("'xvar' is extracted from the 'model'")
    } 
    if(is.null(yvar)) {
      yvar <- model$model_info$yvar
      if(verbose) message("'yvar' is extracted from the 'model'")
    } 
    if(is.null(sigmaxvar)) {
      sigmaxvar <- model$model_info$sigmaxvar
      if(verbose) message("'sigmaxvar' is extracted from the 'model'")
    } 
    if(is.null(xoffset)) {
      xoffset <- model$model_info$xoffset
      if(verbose) message("'xoffset' is extracted from the 'model'")
    } 
    if(is.null(xfuns)) {
      xfuns <- model$model_info$xfuns
      if(verbose) message("'xfuns' is extracted from the 'model'")
    } 
    if(is.null(yfuns)) {
      yfuns <- model$model_info$yfuns
      if(verbose) message("'yfuns' is extracted from the 'model'")
    } 
    if(is.null(sigmaxfuns)) {
      sigmaxfuns <- model$model_info$sigmaxfuns
      if(verbose) message("'sigmaxfuns' is extracted from the 'model'")
    } 
    if(is.null(ixfuns)) {
      ixfuns <- model$model_info$ixfuns
      if(verbose) message("'ixfuns' is extracted from the 'model'")
    } 
    if(is.null(iyfuns)) {
      iyfuns <- model$model_info$iyfuns
      if(verbose) message("'iyfuns' is extracted from the 'model'")
    } 
    if(is.null(isigmaxfuns)) {
      isigmaxfuns <- model$model_info$isigmaxfuns
      if(verbose) message("'isigmaxfuns' is extracted from the 'model'")
    } 
  } # if... else if(!is.null(model)) {
  
 
  
  # the duplicate xvar name ans funs sorted here e.g., multivariate model - age
  # this because list will be with duplicate names and unique names will be left
  if(!is.null(xvar)) {
    templist <- itemplist <- list()
    for (i in 1:length(xvar)) {
      templist[[xvar[[i]]]] <- xfuns[[i]]
      itemplist[[xvar[[i]]]] <- ixfuns[[i]]
    }
    xfuns <- templist
    ixfuns <- itemplist
    if(is_emptyx(xfuns)) xfuns <- NULL
    if(is_emptyx(ixfuns)) ixfuns <- NULL
  }
  
  
  
  if(!is.null(yvar)) {
    templist <- itemplist <- list()
    for (i in 1:length(yvar)) {
      templist[[yvar[[i]]]] <- yfuns[[i]]
      itemplist[[yvar[[i]]]] <- iyfuns[[i]]
    }
    yfuns <- templist
    iyfuns <- itemplist
    if(is_emptyx(yfuns)) yfuns <- NULL
    if(is_emptyx(iyfuns)) iyfuns <- NULL
  }
  
  
  
  
  if(!is.null(sigmaxvar)) {
    templist <- itemplist <- list()
    for (i in 1:length(sigmaxvar)) {
      templist[[sigmaxvar[[i]]]] <- sigmaxfuns[[i]]
      itemplist[[sigmaxvar[[i]]]] <- isigmaxfuns[[i]]
    }
    sigmaxfuns <- templist
    isigmaxfuns <- itemplist
    if(is_emptyx(sigmaxfuns)) sigmaxfuns <- NULL
    if(is_emptyx(isigmaxfuns)) isigmaxfuns <- NULL
  }
  
  # assign xoffset names same as xvar 
  if(!is.null(xoffset)) {
    templist <- itemplist <- list()
    for (i in 1:length(xvar)) {
      templist[[xvar[[i]]]] <- xoffset[[i]]
      itemplist[[xvar[[i]]]] <- xoffset[[i]]
    }
    xoffset <- templist
    ixoffset <- itemplist
    if(is_emptyx(sigmaxfuns)) xoffset <- 0
    if(is_emptyx(isigmaxfuns)) ixoffset <- 0
  }
  
  
  
  
  # xvar         <- unlist(xvar) %>% unique()
  # yvar         <- unlist(yvar) %>% unique()
  # xfuns        <- unlist(xfuns) %>% unique()
  # yfuns        <- unlist(yfuns) %>% unique()
  # sigmaxfuns   <- unlist(sigmaxfuns) %>% unique()
  # ixfuns       <- unlist(ixfuns) %>% unique()
  # iyfuns       <- unlist(iyfuns) %>% unique()
  # isigmaxfuns  <- unlist(isigmaxfuns) %>% unique()
  
  # print(xfuns)
  # print(ixfuns)
  
  if(!is.null(xfuns) & !is.null(ixfuns)) {
    stop("specify either xfuns or ixfuns")
  }
  if(!is.null(yfuns) & !is.null(iyfuns)) {
    stop("specify either yfuns or iyfuns")
  }
  if(!is.null(sigmaxfuns) & !is.null(isigmaxfuns)) {
    stop("specify either yfuns or iyfuns")
  }
  
  # set funs to NULL if not in transform
  # if(!'x' %in% transform) xfuns <- NULL
  # if(!'y' %in% transform) yfuns <- NULL
  # if(!'sigma' %in% transform) sigmaxfuns <- NULL
  
  
  # For testing purpose, allow reverse transformation after transformation
  # But this can't be allowed for running, se set an internal flag
  allow_both_transform_itransform <- T
  
  if(!allow_both_transform_itransform) {
    if('x' %in% transform & 'x' %in% itransform) {
      stop("Either select tranformation or reverse tranformation of 'x', not both")
    } else if(!'x' %in% transform & !'x' %in% itransform) {
      if(!is.null(xvar) & (!is.null(xfuns) | !is.null(ixfuns))) {
        if(verbose) {
          message("No transformation for 'x' variable set. Did you forget to set", 
                  "'transform' or 'itransform?")
        }
      }
    }
    
    if('y' %in% transform & 'y' %in% itransform) {
      stop("Either select tranformation or reverse tranformation of 'y', not both")
    } else if(!'y' %in% transform & !'y' %in% itransform) {
      if(!is.null(yvar) & (!is.null(yfuns) | !is.null(iyfuns))) {
        if(verbose) {
          message("No transformation for 'y' variable set. Did you forget to set", 
                  "'transform' or 'itransform?")
        }
      }
    }
    
    if('sigma' %in% transform & 'sigma' %in% itransform) {
      stop("Either select tranformation or reverse tranformation of 'sigma', not both")
    } else if(!'sigma' %in% transform & !'sigma' %in% itransform) {
      if(!is.null(sigmaxvar) & (!is.null(sigmaxfuns) | !is.null(isigmaxfuns))) {
        if(verbose) {
          message("No transformation for 'sigma' variable set. Did you forget to set", 
                  "'transform' or 'itransform?")
        }
      }
    }
    
  } # if(!allow_both_transform_itransform) {
  
  
  # This already taken care above via allow_both_transform_itransform
  # if(!'x' %in% itransform) {
  
  if('x' %in% itransform) {
    if(is.null(ixfuns)) {
      if(!is.null(xfuns)) {
        ixfuns <- list()
        for (i in names(xfuns)) {
          ixfuns[[i]] <- inverse_transform(base::body(xfuns[[i]]))
        }
      } else {
        # stop("Please specify 'ixfuns'")
      }
    }
  }
  
  
  if('y' %in% itransform) {
    if(is.null(iyfuns)) {
      if(!is.null(yfuns)) {
        iyfuns <- list()
        for (i in names(yfuns)) {
          iyfuns[[i]] <- inverse_transform(base::body(yfuns[[i]]))
        }
      } else {
        # stop("Please specify 'iyfuns'")
      }
    }
  }
  
  
  if('sigma' %in% itransform) {
    if(is.null(isigmaxfuns)) {
      if(!is.null(sigmaxfuns)) {
        isigmaxfuns <- list()
        for (i in names(sigmaxfuns)) {
          isigmaxfuns[[i]] <- inverse_transform(base::body(sigmaxfuns[[i]]))
        }
      } else {
        # stop("Please specify 'isigmaxfuns'")
      }
    }
  }
  
  
  
  
  
  ## transform
  if('x' %in% transform & !is.null(xfuns)) {
    if(is.null(xfuns)) {
      stop("you requested transformation of 'x' but xfuns is 'NULL'")
    } else if(is.null(xvar)) {
      stop("you requested transformation of 'x' but xvar is 'NULL'")
    } else {
      for (i in names(xfuns)) { # don't use xvar, it might reult in duplicate exe
        evalfuns <- xfuns[[i]]
        data[[i]] <- evalfuns(data[[i]])
        evalfuns <- NULL
      }
    }
  }
  
  
  if('y' %in% transform & !is.null(yfuns)) {
    if(is.null(yfuns)) {
      stop("you requested transformation of 'y' but yfuns is 'NULL'")
    } else if(is.null(yvar)) {
      stop("you requested transformation of 'y' but yvar is 'NULL'")
    } else {
      for (i in names(yfuns)) {
        evalfuns <- yfuns[[i]]
        data[[i]] <- evalfuns(data[[i]])
        evalfuns <- NULL
      }
    }
  }
  
  if('sigma' %in% transform & !is.null(sigmaxfuns)) {
    if(is.null(sigmaxfuns)) {
      stop("you requested transformation of 'sigma' but sigmaxfuns is 'NULL'")
    } else if(is.null(sigmaxvar)) {
      stop("you requested transformation of 'sigma' but sigmaxvar is 'NULL'")
    } else {
      for (i in names(sigmaxfuns)) {
        evalfuns <- sigmaxfuns[[i]]
        data[[i]] <- evalfuns(data[[i]])
        evalfuns <- NULL
      }
    }
  }
  
  
  ## itransform -> inverse
  if('x' %in% itransform & !is.null(ixfuns)) {
    if(is.null(ixfuns)) {
      stop("you requested transformation of 'x' but ixfuns is 'NULL'")
    } else if(is.null(xvar)) {
      stop("you requested transformation of 'x' but xvar is 'NULL'")
    } else {
      for (i in names(ixfuns)) { # don't use xvar, it might reult in duplicate exe
        evalfuns <- ixfuns[[i]]
        data[[i]] <- evalfuns(data[[i]])
        evalfuns <- NULL
      }
    }
  }
  
  if('y' %in% itransform & !is.null(iyfuns)) {
    if(is.null(iyfuns)) {
      stop("you requested transformation of 'y' but iyfuns is 'NULL'")
    } else if(is.null(yvar)) {
      stop("you requested transformation of 'y' but yvar is 'NULL'")
    } else {
      for (i in names(iyfuns)) {
        evalfuns <- iyfuns[[i]]
        data[[i]] <- evalfuns(data[[i]])
        evalfuns <- NULL
      }
    }
  }
  
  if('sigma' %in% itransform & !is.null(isigmaxfuns)) {
    if(is.null(isigmaxfuns)) {
      stop("you requested transformation of 'sigma' but isigmaxfuns is 'NULL'")
    } else if(is.null(sigmaxvar)) {
      stop("you requested transformation of 'sigma' but sigmaxvar is 'NULL'")
    } else {
      for (i in names(isigmaxfuns)) {
        evalfuns <- isigmaxfuns[[i]]
        data[[i]] <- evalfuns(data[[i]])
        evalfuns <- NULL
      }
    }
  }
  
  
  return(data)
} # prepare_transformation





