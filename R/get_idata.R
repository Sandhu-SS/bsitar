


#' An internal function to imterpolate data for plotting smooth curves
#' 
#' @param model An object of class \code{bgmfit}. This is optional (default
#'   \code{NULL}) i.e., it is not neccessary to specify the model object. When
#'   \code{model} is specified, then values for \code{newdata}, \code{idVar},
#'   and \code{timeVar} are automatically taken from the \code{model}.
#'
#' @param newdata A data frame. If \code{NULL} (default), data analysed in the
#'   original model fit is used.
#'
#' @param idVar A character string to specify the group identifier. If
#'   \code{NULL} (default), \code{id} from the model fit is used.
#'
#' @param timeVar  A character string to specify the time variable. If
#'   \code{NULL} (default), \code{x} from the model fit is used.
#'
#' @param times  A numeric vector to specify the time range. Currently ignored.
#'
#' @param length.out A numeric value to specify the length of interpolation
#'   points. Default 10.
#'
#' @param xrange An integer to set the predictor range (i.e., age) when
#'   executing the interpolation via \code{ipts}. The default \code{NULL} sets
#'   the individual specific predictor range whereas code \code{xrange = 1} sets
#'   same range for all individuals within the higher order grouping variable
#'   (e.g., study). Code \code{xrange  = 2} sets the identical range
#'   dplyr::across the entire sample. Lastly, a paired numeric values can be
#'   supplied e.g., \code{xrange = c(6, 20)} will set the range between 6 and
#'   20.
#'
#' @return A data frame.
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#'
get_idata <-
  function(model = NULL,
           newdata = NULL,
           idVar = NULL,
           timeVar = NULL,
           times = NULL,
           length.out = 10,
           xrange = 1) {
    
    if (is.null(newdata)) {
      newdata <- model$data
    } else {
      newdata <- newdata
    }
    
    if(is.null(model)) {
      if (is.null(idVar)) stop("Specify model or idVar, both can not be NULL")
      if (is.null(timeVar)) stop("Specify model or timeVar, both can not be NULL")
    }

    if(!is.null(model)) {
      if (!is.null(idVar)) stop("Specify either model or idVar, not both")
      if (!is.null(timeVar)) stop("Specify either model or timeVar, not both")
    }
    
    if(!is.null(model)) {
     if(length(model$model_info$ids) > 1) {
       stop("Please specify newdata, idVar, and timeVar manullay because 
            currently value for these can not be infered from model with three 
            or more levels of hierarchy")
     }
    }
    
    
    
    if (is.null(idVar)) {
      idVar <- model$model_info$ids
    } else {
      idVar <- idVar
    }
    
    if (is.null(timeVar)) {
      timeVar <- model$model_info$xvar
    } else {
      timeVar <- timeVar
    }
    
    all_times <- TRUE
    if (is.null(xrange))
      xrange <- 1
    else
      xrange <- xrange
    times_orig <- newdata[[timeVar]]
    times_orig <- times_orig[!is.na(times_orig)]
    
    if (is.null(times) || !is.numeric(times)) {
      times <-
        seq(min(times_orig), max(times_orig), length.out = length.out)
    }
    
    # this is when no random effects and groupvar is NULL
    # therefore, an artificial group var created
    # check utils-helper function lines 60
    
    if (nlevels(newdata[[idVar]]) == 1) {
      newdata <- newdata %>%
        dplyr::distinct(newdata[[timeVar]], .keep_all = T) %>%
        dplyr::arrange(!!as.name(timeVar))
    }
    
    id <- match(newdata[[idVar]], unique(newdata[[idVar]]))
    last_time <- tapply(newdata[[timeVar]], id, max)
    first_time <- tapply(newdata[[timeVar]], id, min)
    
    newdata_nomiss <- newdata[complete.cases(newdata),]
    id_nomiss <-
      match(newdata_nomiss[[idVar]], unique(newdata_nomiss[[idVar]]))
    n <- length(unique(id_nomiss))
    
    if (xrange == 1) {
      times_to_pred <- list()
      for (i in 1:length(unique(newdata[[idVar]]))) {
        numx <- as.character(i)
        times_to_pred[[numx]] <-
          seq(first_time[i], last_time[i], length.out = length.out)
      }
    }
    
    if (xrange == 2) {
      times_to_pred <- lapply(last_time, function (t)
        if (all_times)
          times
        else
          times[times > t])
    }
    
    id_pred <- rep(seq_len(n), sapply(times_to_pred, length))
    
    right_rows <- function (data, times, ids, Q_points) {
      fids <- factor(ids, levels = unique(ids))
      if (!is.list(Q_points))
        Q_points <- split(Q_points, row(Q_points))
      ind <- mapply(findInterval, Q_points, split(times, fids))
      ind[ind < 1] <- 1
      rownams_id <- split(row.names(data), fids)
      ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
      data[c(ind),]
    }
    
    newdata_pred <-
      right_rows(newdata, newdata[[timeVar]], id, times_to_pred)
    newdata_pred[[timeVar]] <- unlist(times_to_pred)
    newdata_pred
  }
