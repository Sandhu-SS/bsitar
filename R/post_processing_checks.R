


#' An internal function to perform checks when calling the post-processing
#' functions
#'  
#' @description The \code{post_processing_checks} perform essential checks 
#' (such as the validity of model class, response etc.) when post-processing 
#' functions
#'
#' @param model An object of class \code{bsitar}.
#' function.
#' @param xcall The \code{match.call()} from the post-processing function.
#' @param resp Response variable (default \code{NULL}) specified as a string 
#' character required during the post-processing of multivariate and 
#' univariate-by-subgroup model (see \code{bsitar} function for details).
#' @param deriv An integer value to specify whether to estimate distance curve
#'  (i.e., model estimated curve(s)) or velocity curve (first derivative of the 
#'  model estimated curve(s)). A 0 value (default) is for distance curve and 
#'  1 for the velocity curve.
#'
#' @return A string with the error captured or else a list with necessary  
#' information needed when executing the post-processing function 
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' #
#' post_processing_checks(model)
#' }

post_processing_checks <- function(model, xcall, resp = NULL, deriv = 0) {
  if(!'bsitar' %in% class(model)) {
    stop("The class of model object should be 'bsitar' ")
  }
  excall_ <- c("pp_check_", "loo_")
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
  list(
    paste0(resp_, model$model_info$SplineFun, ''),
    paste0(resp_, model$model_info$SplineFun, deriv)
  )
}

