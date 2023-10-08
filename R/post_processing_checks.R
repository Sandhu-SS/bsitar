


#' An internal function to perform checks when calling post-processing functions
#'
#' @description The \code{post_processing_checks} perform essential checks (such
#'   as the validity of model class, response etc.) during post-processing of
#'   posterior draws.
#'
#' @param model An object of class \code{bgmfit}.
#'
#' @param xcall The \code{match.call()} from the post-processing function.
#'
#' @param resp Response variable (default \code{NULL}) specified as a character
#'   string. This is needed when processing univariate-by-subgroup and
#'   multivariate model (see \code{bgmfit} function for details).
#'
#' @param deriv An integer value to specify whether to estimate distance curve
#'   (i.e., model estimated curve(s)) or the velocity curve (first derivative of
#'   the model estimated curve(s)). A value \code{0} (default) is for distance
#'   curve and  \code{1} for the velocity curve.
#'
#' @param envir Indicator to set the environment of function evaluation. The
#'   default is \code{parent.frame}.
#'
#' @return A string with the error captured or else a list with necessary
#'   information needed when executing the post-processing function
#'
#' @noRd
#'
post_processing_checks <- function(model, 
                                   xcall, 
                                   resp = NULL, 
                                   deriv = NULL,
                                   envir = parent.frame()) {
  
  if(is.null(envir)) envir <- parent.frame()
  if(is.null(deriv)) deriv <- 0
  
  if(!'bgmfit' %in% class(model)) {
    stop("The class of model object should be 'bgmfit' ")
  }
  excall_ <- c("pp_check_bgm", "loo_bgm")
  if (strsplit(deparse((xcall[1])), "\\.")[[1]][1] %in% excall_) {
    if (!is.null(as.list(xcall)[['deriv']])) {
      stop(
        "argument deriv is not allowed for the ",
        "\n ",
        " post-processing function",
        " '",
        strsplit(deparse((xcall[1])), "\\.")[[1]][1],
        "'"
      )
    }
    deriv <- ''
  }
  if (model$model_info$nys == 1 & !is.null(resp)) {
    stop(
      "You have fit a univariate model",
      " but set resp option as: ",
      resp,
      ".",
      "\n ",
      " For univariate model, the resp option should be NULL",
      "\n ",
      " (i.e., resp = NULL)"
    )
  }
  if (model$model_info$nys > 1 & is.null(resp)) {
    if (!is.na(model$model_info$univariate_by)) {
      stop(
        "You have fit a univariate-by-subset model for ",
        model$model_info$univariate_by,
        "\n ",
        " but dit not set the the resp options correctly",
        " (which is NULL at present).",
        "\n ",
        " The response options are: ",
        paste(model$model_info$ys, collapse = ", ")
      )
    }
    if (model$model_info$multivariate) {
      stop(
        "You have fit a multivariate model ",
        "\n ",
        " but dit not set the the resp options correctly",
        " (which is NULL at present).",
        "\n ",
        " The response options are: ",
        paste(model$model_info$ys, collapse = ", ")
      )
    }
  }
  if (is.null(resp)) {
    resp_ <- resp
  } else if (!is.null(resp)) {
    resp_ <- paste0(resp, "_")
  }
  
  
  # assign expose default funs 
  if(model$model_info[['expose_method']] == 'R') {
    assign(paste0(resp_, 
                  model$model_info[['namesexefuns']], 
                  '0'), 
           model$model_info$exefuns[[paste0(resp_, 
                                            model$model_info[['namesexefuns']], 
                                            '0')]], envir = envir)
    
    assign(paste0(resp_, 'getX'), 
           model$model_info$exefuns[[paste0(resp_, 'getX')]], 
           envir = envir)
    
    if(model$model_info[['select_model']] == 'sitar' |
       model$model_info[['select_model']] == 'rcs') {
      assign(paste0(resp_, 'getKnots'), 
             model$model_info$exefuns[[paste0(resp_, 'getKnots')]], 
             envir = envir)
    }
  }
  

  list(
    paste0(resp_, model$model_info[['namesexefuns']], ''),
    paste0(resp_, model$model_info[['namesexefuns']], deriv)
    ) 
  
}

