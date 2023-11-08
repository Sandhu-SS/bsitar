

#' Expose Stan function for post-processing of posterior sample obtained 
#' from the model fit
#'
#' @param model An object of class \code{bgmfit}.
#' @param scode An option argument specifying the code with the user-defined 
#' Stan function.
#' @param expose A logical (default \code{TRUE}) to indicate whether to expose
#' functions and add to the global environment. 
#' 
#' @param select_model A string (default \code{NULL}) to indicate the model.
#'
#' @return An object of class \code{bgmfit} with exposed  
#' user-defined Stan functions (when \code{expose=TRUE}).
#' 
#' @export
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#' @examples
#' 
#' # The examples below show the use of *conditional_effects_bgm* to plot  
#' # the population average and individual-specific distance and velocity 
#' # curves.
#' 
#' # Fit Bayesian SITAR model 
#' # To avoid running the model which takes some time, model fit to the
#' # \code{berkeley_mdata} has already been saved as berkeley_mfit.rda object.
#' # Please see \code{bgm} examples.
#' 
#' model <- berkeley_mfit
#' 
#' \donttest{
#' expose_functions_bgm(model, expose = TRUE)
#' }
#' 
expose_functions_bgm <- function(model, 
                                    scode = NULL, 
                                    expose = FALSE, 
                                    select_model = NULL) {
  expose_it_ <- function(x, scode) {
    if (is.null(scode)) {
      exposecode <- brms::stancode(x)
    } else if (!is.null(scode)) {
      exposecode <- scode
    }
     rstan::expose_stan_functions(rstan::stanc(model_code = exposecode))
    # brms::expose_functions(model, vectorize = FALSE)
  }
  
  
  
  # allowed_model_names <- c('sitar', 'pb1', 'pb2', 'pb3')
  # allowed_model_names <- paste(allowed_model_names, collapse = ", " )
  # allowed_model_names <- paste0("(", allowed_model_names, ")")
  # 
  # if(!expose & is.null(select_model)) 
  #   stop('Specify atleast one of the following two arguments:',
  #        "\n ",
  #        ' expose  - a logical argument (TRUE/FALSE)',
  #        "\n ",
  #        ' select_model - a string specifying one of the models ', 
  #        allowed_model_names
  #        )
  # 
  # if( expose & !is.null(select_model)) 
  #   stop('Specify only one of the following two arguments:',
  #        "\n ",
  #        ' expose  - a logical argument (TRUE/FALSE)',
  #        "\n ",
  #        ' select_model - a string specifying one of the models ', 
  #        allowed_model_names
  #        )
  
  
  if(!expose) {
    if (is.null(model$model_info$decomp))  expose_r_from_stan <- TRUE
    if (!is.null(model$model_info$decomp)) expose_r_from_stan <- FALSE
  } else {
    expose_r_from_stan <- FALSE
  }
  
  
  if(expose) expose_it_(model, scode = scode)
  
  
  if(expose_r_from_stan) {
    for (funi in 1:length(model$model_info$funlist_r)) {
      assign(gsub("<-.*$", "", model$model_info$funlist_r[funi]),
             ept(model$model_info$funlist_r[funi]), envir = parent.frame())
    }
  }
  
  
  # SplineFun_name <- model$model_info$SplineFun
  SplineFun_name <- model$model_info[['StanFun_name']]
  spfun_collect <- c(SplineFun_name,
                     paste0(SplineFun_name, "_", 
                            c("d0", 
                              "d1",
                              "d2")
                            )
                     )
  
  
  if(expose_r_from_stan) {
    spfun_collect <- c(spfun_collect, 'getX')
    if(select_model == 'sitar' | select_model == 'rcs') {
      spfun_collect <- c(spfun_collect, 'getKnots')
    }
  }
  
  
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
  

  
  
  if(expose) {
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
  } # if(expose) {
  
  
  
  if(expose_r_from_stan) {
    Spl_funs <- list()
    spfun_collectic <- -1
    for (spfun_collecti in spfun_collect) {
      spfun_collectic <- spfun_collectic + 1
      spfun_collecti_name <- spfun_collecti
      spfun_collecti_name <- gsub("_d0", "0", spfun_collecti_name)
      spfun_collecti_name <- gsub("_d1", "1", spfun_collecti_name)
      spfun_collecti_name <- gsub("_d2", "2", spfun_collecti_name)
      getfun_ <- spfun_collecti
      # This below to change _d0 to 0 within the d2 d2 functions 
      getfun__ <- deparse(ept(getfun_))
      gsub_it <- '_d0'
      gsub_by <- "0"
      getfun__ <- gsub(gsub_it, gsub_by, getfun__, fixed = T)
      getfun__ <- paste0(getfun__, collapse =  "\n")
      Spl_funs[[paste0(spfun_collecti_name, "")]] <- 
        eval(parse(text = getfun__), envir = parent.frame())
    }
  } # if(expose_r_from_stan) {
  
  
  if(!expose & !expose_r_from_stan) Spl_funs <- NULL
  
  
  
   # This was from the wrap_function2
  # if(!is.null(select_model)) {
  #   Spl_funs <- list()
  #   # wrap_function  -> in utils-helper-3, and execute and = constructs functions'
  #   # wrap_function2 -> in utils-helper-4, and uses constructed functions
  #   # evaluatred_funs <- wrap_function(select_model)
  #   evaluatred_funs <- wrap_function2(select_model)
  #   spfun_collectic <- -1
  #   for (spfun_collecti in 1:length(evaluatred_funs)) {
  #     spfun_collectic <- spfun_collectic + 1
  #     spfun_collecti_name <- paste0(SplineFun_name, spfun_collectic)
  #     Spl_funs[[paste0(spfun_collecti_name, "")]] <- evaluatred_funs[[spfun_collecti]]
  #   }
  # }
  
  
  
  
  model$model_info[['namesexefuns']] <- SplineFun_name
  model$model_info[['exefuns']]      <- Spl_funs
  scode_include <- brms::stancode(model)
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
  model$model <- model$bmodel # scode_include
  return(model)
}

