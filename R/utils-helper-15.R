


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
#' @param keeplevels A logical in case factor variables other than \code{idVar}
#'   present
#'   
#' @param asdf
#' 
#' @param get.newdata.call
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
           xrange = 1, 
           keeplevels = FALSE, 
           asdf = FALSE,
           get.newdata.call = FALSE) {
    
    if(get.newdata.call) {
      if(is.null(newdata)) {
        stop("'newdata' can not be NULL when get.newdata.call = TRUE")
      }
    }
    
    if (is.null(newdata)) {
      newdata <- model$model_info$bgmfit.data
    } else {
      newdata <- newdata
    }
    
    
    
    if(!get.newdata.call) {
      # prepare_data2 with model = model will get all the necessary info
      newdata <- prepare_data2(data = newdata, model = model)
      newdata <- prepare_transformations(data = newdata, model = model)
    }
    
    
    if(data.table::is.data.table(newdata)) {
      setasdt <- TRUE 
      newdata <- as.data.frame(newdata)
    } else {
      setasdt <- FALSE
    }
    
    if(keeplevels) {
      is.fact <- names(newdata[, sapply(newdata, is.factor)])
      cnames  <- colnames(newdata)
    }
    
    
    if(is.null(model)) {
      if (is.null(idVar)) stop("Specify model or idVar, both can not be NULL")
      if (is.null(timeVar)) stop("Specify model or timeVar, both can't be NULL")
    }
    
    if(!is.null(model)) {
      if (!is.null(idVar)) stop("Specify either model or idVar, not both")
      if (!is.null(timeVar)) stop("Specify either model or timeVar, not both")
    }
    
    if(!is.null(model)) {
      if(length(model$model_info$idvars) > 1) {
        stop("Please specify newdata, idVar, and timeVar manullay because 
            currently value for these can not be infered from model with three 
            or more levels of hierarchy")
      }
    }
    
    `.` <- NULL;
    
    if (is.null(idVar)) {
      idVar <- model$model_info$idvars
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
    
    # This is when no random effects and groupvar is NULL
    # Therefore, an artificial group var created
    # Check utils-helper function lines 60
    
    if (nlevels(newdata[[idVar]]) == 1) {
      newdata <- newdata %>%
        dplyr::distinct(newdata[[timeVar]], .keep_all = T) %>%
        dplyr::arrange(!!as.name(timeVar))
    }
    
    id <- match(newdata[[idVar]], unique(newdata[[idVar]]))
    
    if(length( unique(newdata[[idVar]])) == 1) {
      if(length.out == 1) stop("The argument 'ipts' should be > 1")
    }
    
    last_time  <- tapply(newdata[[timeVar]], id, max)
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
      # 1.3.2025
      # Faced issue newdata[[timeVar]] not sorted 
      if(base::is.unsorted(times)) {
        times <- sort(times)
      }
      ind <- mapply(findInterval, Q_points, split(times, fids))
      ind[ind < 1] <- 1
      rownams_id <- split(row.names(data), fids)
      ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
      data[c(ind),]
    }
    
    
    # 6.03.2025
    if(is.null(length.out)) {
      newdata_pred <- newdata
    } else {
      if(length.out == 1) {
        stop("The 'ipts' should be either 'NULL' or > 1")
      } else {
        newdata_pred <- right_rows(newdata, newdata[[timeVar]], id, 
                                   times_to_pred)
      }
    }
    
    # newdata_pred <- right_rows(newdata, newdata[[timeVar]], id, times_to_pred)
    
    if(keeplevels) {
      if(length(setdiff(is.fact, idVar)) > 0) {
        newdata_pred <- newdata_pred %>% droplevels
        newdata_pred <- newdata_pred %>% 
          dplyr::select(-dplyr::all_of(setdiff(is.fact, idVar)))
        
        newdata_is.factx <- newdata %>% 
          dplyr::select(dplyr::all_of(is.fact))
        
        newdata_pred <- newdata_pred %>% 
          dplyr::left_join(., newdata_is.factx, by = idVar,
                    relationship = "many-to-many")
      }
    }
    
  
    newdata_pred[[timeVar]] <- unlist(times_to_pred)
    
    if(keeplevels) {
      newdata_pred <- newdata_pred %>% dplyr::select(dplyr::all_of(cnames))
    }
    if(setasdt) newdata_pred <- data.table::as.data.table(newdata_pred)
    if(asdf) out <- as.data.frame(newdata_pred) else out <- newdata_pred 
    # newdata_pred
   return(out)
  }
