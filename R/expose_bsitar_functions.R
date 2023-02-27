

#' Expose Stan function for post-processing of posterior sample obtained 
#' from the Bayesian SITAR growth curve model fit
#'
#' @param model An object of class \code{bsitar}.
#' @param scode An option argument specifying the code with the user-defined 
#' Stan function.
#' @param expose A logical (default \code{TRUE}) to indicate whether to expose
#' functions and add to the global environment. 
#'
#' @return An object of class \code{brmsfit, bsiatr} with exposed  
#' user-defined Stan functions (when \code{expose=TRUE}).
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # 
#' expose_bsitar_functions(model)
#' #
#' }
expose_bsitar_functions <- function(model, scode = NULL, expose = TRUE) {
  expose_it_ <- function(x, scode) {
    if (is.null(scode)) {
      exposecode <- stancode(x)
    } else if (!is.null(scode)) {
      exposecode <- scode
    }
    rstan::expose_stan_functions(rstan::stanc(model_code = exposecode))
  }
  if(expose) expose_it_(model, scode = scode)
  SplineFun_name <- model$model_info$SplineFun
  spfun_collect <- c(model$model_info$SplineFun,
                     paste0(model$model_info$SplineFun, "_", c("d0", 
                                                                "d1",
                                                                "d2")))
  
  nys <- model$model_info$nys
  ys <- model$model_info$ys
  if(nys > 1) {
    spfun_collect2 <- c()
    for (ysii in ys) {
      tempysi <- paste0(ysii, "_", spfun_collect)
      spfun_collect2 <- c(spfun_collect2, tempysi)
    }
    spfun_collect <- spfun_collect2
  }
  Spl_funs <- list()
  spfun_collectic <- -1
  for (spfun_collecti in spfun_collect) {
    spfun_collectic <- spfun_collectic + 1
    spfun_collecti_name <- spfun_collecti
    spfun_collecti_name <- gsub("_d0", "0", spfun_collecti_name)
    spfun_collecti_name <- gsub("_d1", "1", spfun_collecti_name)
    spfun_collecti_name <- gsub("_d2", "2", spfun_collecti_name)
    Spl_funs[[paste0(spfun_collecti_name, "")]] <-
      eval(parse(text = spfun_collecti), envir = parent.frame())
  }
  model$Spl_funs <- Spl_funs
  scode_include <- stancode(model)
  model$bmodel <- scode_include
  if (nys == 1 | nys > 1) {
    for (nys__i in 1:nys) {
      cont_ <- 0
      for (cont_i in 0:2) {
        cont_ <- cont_ + 1
        if (nys == 1) {
          gsubit <- paste0(
            "vector",
            " ",
            paste0("", "", SplineFun_name),
            "_",
            "d",
            cont_i,
            paste0(".*end of spline function", "_", ys[nys__i],
                   "d", cont_i, "")
          )
        } else if (nys > 1) {
          gsubit <-
            paste0(
              "vector",
              " ",
              paste0(ys[nys__i], "_", SplineFun_name),
              "_",
              "d",
              cont_i,
              paste0(".*end of spline function", "_", ys[nys__i],
                     "d", cont_i, "")
            )
        }
        scode_include <-
          gsub(gsubit, "", scode_include, fixed = F)
      }
    }
  }
  model$model <- scode_include
  return(model)
}

