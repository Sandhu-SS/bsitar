

#' An internal function to create data frame for post-processing
#'
#' @inheritParams  growthparameters.bgmfit 
#' 
#' @param idata_method A character string to indicate the interpolation method.
#'   Options are \emph{method 1} (specified as  \code{'m1'}, default) and
#'   \emph{method 2} (specified as \code{'m2'}). The \code{'m1'} calls an
#'   internal function \code{idatafunction} whereas \code{'m2'} calls the
#'   \code{get_idata} function for data interpolation. The \emph{method 1}
#'   (\code{'m1'}) is adapted from the the \pkg{iapvbs} and is documented here
#'   <https://rdrr.io/github/Zhiqiangcao/iapvbs/src/R/exdata.R>. The
#'   \emph{method 2} (\code{'m2'}) is adapted from the the \pkg{JMbayes} and is
#'   documented here
#'   <https://github.com/drizopoulos/JMbayes/blob/master/R/dynPred_lme.R>.
#'   If \code{NULL} (default), method \code{'m1'} is automatically set. 
#' 
#' @keywords internal
#' 
#' @return A data frame object. 
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#'
get.newdata <- function(model,
                        newdata = NULL,
                        resp = NULL,
                        numeric_cov_at = NULL,
                        aux_variables = NULL,
                        levels_id = NULL,
                        ipts = NULL,
                        xrange = NULL,
                        idata_method = NULL,
                        dummy_to_factor = NULL,
                        verbose = FALSE,
                        ...) {
  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
 
  
  if (is.null(idata_method)) {
    idata_method <- 'm2'
  }
  

  # Initiate non formalArgs()
  `:=` <- NULL
  . <- NULL;
  
  
  validate_response(model, resp)
  
  list_c <- list()
  xvar_ <- paste0('xvar', resp_rev_)
  yvar_ <- paste0('yvar', resp_rev_)
  groupvar_ <- paste0('groupvar', resp_rev_)
  xvar <- model$model_info[[xvar_]]
  yvar <- model$model_info[[yvar_]]
  hierarchical_ <- paste0('hierarchical', resp_rev_)
  
  
  
  if (is.null(levels_id)) {
    IDvar <- model$model_info[[groupvar_]]
    if (!is.null(model$model_info[[hierarchical_]])) {
      IDvar <- model$model_info[[hierarchical_]]
    } else if (is.null(model$model_info[[hierarchical_]])) {
      # if(!is.null(model$model_info[['idvars']])) {
      #   IDvar <- model$model_info[['idvars']]
      # }
    }
  } else if (!is.null(levels_id)) {
    IDvar <- levels_id
  }
  
  xfun_ <- paste0('xfun', resp_rev_)
  yfun_ <- paste0('yfun', resp_rev_)
  xfun  <- model$model_info[[xfun_]]
  yfun  <- model$model_info[[yfun_]]
  
  
  cov_       <- paste0('cov', resp_rev_)
  cov_sigma_ <- paste0('cov_sigma', resp_rev_)
  uvarby     <- model$model_info$univariate_by$by
  
  
  # When no random effects and hierarchical, IDvar <- NULL problem 02 03 2024
  #if(idata_method == 'm2') {
  if(is.null(levels_id)) {
    if(is.null(IDvar)) {
      if(!is.null(model$model_info[['idvars']])) {
        IDvar <- model$model_info[['idvars']]
      }
    }
  }
  #}
  
  
  ##################################################
  # prepare_data2 changes 
  # if (is.null(newdata)) {
  #   if(idata_method == 'm1') newdata <- model$model_info$bgmfit.data
  #   if(idata_method == 'm2') newdata <- model$data
  # } else {
  #   newdata <- newdata
  # }
  ##################################################
  
  if (is.null(newdata)) {
    if(idata_method == 'm1') newdata <- model$model_info$bgmfit.data
    if(idata_method == 'm2') newdata <- model$model_info$bgmfit.data
  } else {
    newdata <- newdata
  }
  

  if(!is.null(dummy_to_factor)) {
    if(!is.list(dummy_to_factor)) {
      stop("dummy_to_factor must be a named list as follows:",
           "\n ",
           "dummy_to_factor = list(factor.dummy = , factor.name = )",
           "\n ",
           "where factor.dummy is a vector of character strings that will" ,
           "\n ",
           "be converted to a factor variable, and ",
           "\n ",
           "factor.name is a single character string that is used to name the",
           "\n ",
           "newly created factor variable"
      )
    }
    factor.dummy <- dummy_to_factor[['factor.dummy']]
    factor.level <- dummy_to_factor[['factor.level']]
    if(is.null(dummy_to_factor[['factor.name']])) {
      factor.name <- 'factor.dummy'
    } else {
      factor.name <- dummy_to_factor[['factor.name']]
    }
    newdata <- convert_dummy_to_factor(df = newdata, 
                                       factor.dummy = factor.dummy,
                                       factor.level = factor.level,
                                       factor.name = NULL)
    
    colnames(newdata) <- gsub('factor.var', factor.name, names(newdata))
    newdata[[factor.name]] <- as.factor(newdata[[factor.name]])
  }
  
  
  # This is when no random effects and this groupvar is NULL
  # Therefore, an artificial group var created
  # see also changes made to the get_idata function lines 17
  
  # if (is.null(model$model_info$groupvar)) {
  #   name_hypothetical_id <- paste0("id", resp_rev_)
  #   model$model_info$groupvar <- name_hypothetical_id
  #   newdata[[name_hypothetical_id]] <- as.factor("tempid")
  # } else if (!is.null(model$model_info$groupvar)) {
  #   if(length(newdata[[model$model_info$groupvar]]) == 0) {
  #     # name_hypothetical_id <- paste0("hy_id", resp_rev_)
  #     if(length(IDvar) > 1) {
  #       name_hypothetical_id <- IDvar[1] 
  #     } else {
  #       name_hypothetical_id <- IDvar
  #     }
  #     model$model_info$groupvar <- name_hypothetical_id
  #     newdata[[name_hypothetical_id]] <- as.factor("tempid")
  #   }
  # }
  
  
  
  
  newdata <- check_newdata_args(model, newdata, IDvar, resp)
  
  
  ##################################################
  # prepare_data2 changes 
  # if (!is.na(model$model_info$univariate_by$by)) {
  #   if (is.symbol(model$model_info$call.bgmfit$y)) {
  #     setorgy <- deparse(model$model_info$call.bgmfit$y)
  #   } else if (is.list(model$model_info$call.bgmfit$y)) {
  #     setorgy <- unname(unlist(model$model_info$call.bgmfit$y))
  #     if (is.symbol(setorgy))
  #       setorgy <- deparse(setorgy)
  #   } else {
  #     setorgy <- model$model_info$call.bgmfit$y
  #   }
  # }
  # 
  # if (is.na(model$model_info$univariate_by$by)) {
  #   setorgy <- model$model_info$ys
  # }
  # 
  # 
  # if (idata_method == 'm1') {
  #   newdata <- prepare_data2(
  #     # newdata <- prepare_data(
  #     data = newdata,
  #     x = model$model_info$xs,
  #     y = setorgy,
  #     id = model$model_info$ids,
  #     uvarby = model$model_info$univariate_by$by,
  #     mvar = model$model_info$multivariate$mvar,
  #     xfuns = model$model_info$xfuns,
  #     yfuns = model$model_info$yfuns,
  #     outliers = model$model_info$outliers)
  # }
  ##################################################

  # prepare_data2 with model = model will get all the necessary info
  newdata <- prepare_data2(data = newdata, model = model)
  newdata <- prepare_transformations(data = newdata, model = model)
  

  
  newdata <- newdata[,!duplicated(colnames(newdata))]
  
  ##################################################
  # prepare_data2 changes 
  if (!is.na(model$model_info$univariate_by$by)) {
    # if(idata_method == 'm1') {
    #   sortbylayer <- NA
    #   newdata <- newdata %>%
    #     dplyr::mutate(sortbylayer =
    #                     forcats::fct_relevel(!!as.name(uvarby),
    #                                          (levels(
    #                                            !!as.name(uvarby)
    #                                          )))) %>%
    #     dplyr::arrange(sortbylayer) %>%
    #     dplyr::mutate(!!as.name(IDvar) := factor(!!as.name(IDvar),
    #                                              levels =
    #                                                unique(!!as.name(IDvar)))) %>%
    #     dplyr::select(-sortbylayer)
    # } # if(idata_method == 'm1') {
    ##################################################
    
    subindicatorsi <- model$model_info$subindicators[grep(resp,
                                                          model$model_info$yvars)]
    # This aux_variables used in idatafunction etc in idata_method = 1
    # otherwise subindicatorsi is not passed to fitted etc
    aux_variables <- subindicatorsi
    
    list_c[['subindicatorsi']] <- subindicatorsi
    list_c[['uvarby']]         <- model$model_info$univariate_by$by
  } # end if (!is.na(model$model_info$univariate_by$by)) {
  
  
  covars_extrcation <- function(str) {
    str <- gsub("[[:space:]]", "", str)
    for (ci in c("*", "+", ":")) {
      str <- gsub(ci, ' ', str, fixed = T)
    }
    str <- strsplit(str, " ")[[1]]
    str
  }
  
  
  cov_vars       <-  model$model_info[[cov_]]
  cov_sigma_vars <-  model$model_info[[cov_sigma_]]
  
  
  if (!is.null(cov_vars)) {
    cov_vars <- covars_extrcation(cov_vars)
  }
  if (!is.null(cov_sigma_vars)) {
    cov_sigma_vars <- covars_extrcation(cov_sigma_vars)
  }
  
  
  if(is.null(cov_vars) & is.null(cov_sigma_vars)) {
    if(!is.null(dummy_to_factor)) {
      warning("There are no covariate(s) but have specified dummy_to_factor",
              "\n ", 
              "Please check if this is an error")
    }
  }
  
  
  if(!is.null(dummy_to_factor)) {
    cov_vars <- c(cov_vars, factor.name)
  }
  
  
  # check if cov is charcater but not factor 
  checks_for_chr_fact <- c(cov_vars, cov_sigma_vars)
  checks_for_chr_fact <- unique(checks_for_chr_fact)
  for (cov_varsi in checks_for_chr_fact) {
    if(is.character(newdata[[cov_varsi]])) {
      if(!is.factor(newdata[[cov_varsi]])) {
        if(verbose) {
          message("\nVariable '", cov_varsi, "' used as a covariate in the model ",
                  "\n ",
                  " is a character but not factor. Converting it to factor.")
        }
        newdata[[cov_varsi]] <- as.factor(newdata[[cov_varsi]])
        factor_vars <- names(newdata[sapply(newdata, is.factor)])
      }
    }
  }
  
  factor_vars <- names(newdata[sapply(newdata, is.factor)])
  numeric_vars <- names(newdata[sapply(newdata, is.numeric)])
  
  # if (!is.null(cov_sigma_vars))
  #   cov_sigma_vars <- covars_extrcation(cov_sigma_vars)
  
  
  
  cov_factor_vars <- intersect(cov_vars, factor_vars)
  cov_numeric_vars <- intersect(cov_vars, numeric_vars)
  groupby_fstr <- c(cov_factor_vars)
  groupby_fistr <- c(IDvar, cov_factor_vars)
  
  cov_sigma_factor_vars <- intersect(cov_sigma_vars, factor_vars)
  cov_sigma_numeric_vars <- intersect(cov_sigma_vars, numeric_vars)
  
  if (identical(cov_factor_vars, character(0)))
    cov_factor_vars <- NULL
  if (identical(cov_numeric_vars, character(0)))
    cov_numeric_vars <- NULL
  
  if (identical(cov_sigma_factor_vars, character(0)))
    cov_sigma_factor_vars <- NULL
  if (identical(cov_sigma_numeric_vars, character(0)))
    cov_sigma_numeric_vars <- NULL
  
  # Merge here a b c covariate with sigma co variate
  # IMP: Note that groupby_fstr and groupby_fistr are stil  a b c covariate
  # This way, plot_curves and gparameters will not produce sigam cov specific
  # curves and g parameters
  
  cov_factor_vars <- c(cov_factor_vars, cov_sigma_factor_vars)
  cov_numeric_vars <- c(cov_numeric_vars, cov_sigma_numeric_vars)
  
  if (!is.na(model$model_info$univariate_by$by)) {
    if(idata_method == 'm1') groupby_fstr <- c(uvarby, groupby_fstr)
    if(idata_method == 'm1') groupby_fistr <- c(uvarby, groupby_fistr)
  }
  
  
  set_numeric_cov_at <- function(x, numeric_cov_at) {
    name_ <- deparse(substitute(x))
    if (is.null((numeric_cov_at[[name_]]))) {
      . <- mean(x, na.rm = T)
    } else if (!is.null((numeric_cov_at[[name_]]))) {
      if (numeric_cov_at[[name_]] == 'mean') {
        . <- mean(x,  na.rm = T)
      } else if (numeric_cov_at[[name_]] == 'median') {
        . <- median(x, na.rm = T)
      } else if (numeric_cov_at[[name_]] == 'min') {
        . <- min(x, na.rm = T)
      } else if (numeric_cov_at[[name_]] == 'max') {
        . <- max(x, na.rm = T)
      } else {
        . <- numeric_cov_at[[name_]]
      }
    }
    round(., 3)
  }
  
  
  
  get.data.grid <- function(data,
                            xvar,
                            yvar,
                            IDvar,
                            cov_numeric_vars,
                            numeric_cov_at,
                            aux_variables,
                            uvarby) {
    if (!is.null(IDvar))
      relocate_vars <- c(xvar, IDvar)
    if (is.null(IDvar))
      relocate_vars <- c(xvar)
    if (!is.na(uvarby))
      if(idata_method == 'm1') relocate_vars <- c(relocate_vars, uvarby)
    if(idata_method == 'm2') relocate_vars <- c(relocate_vars)
    if (!is.null(cov_numeric_vars)) {
      cov_numeric_vars__ <- cov_numeric_vars
      if (identical(cov_numeric_vars__, character(0)))
        cov_numeric_vars__ <- NULL
      if (!is.null(cov_numeric_vars__)) {
        for (cov_numeric_vars__i in cov_numeric_vars__) {
          if (!is.null(numeric_cov_at)) {
            if (!cov_numeric_vars__i %in% names(numeric_cov_at)) {
              stop(
                "You have used the argument 'numeric_cov_at' to specify the",
                "\n ",
                " value of continous covariate '",
                cov_numeric_vars__i,
                "'.",
                "\n ",
                " However, the name of this covariate is missing from the list",
                "\n ",
                " Please use argument 'numeric_cov_at' correctly as follows:",
                "\n ",
                " numeric_cov_at = list(",
                cov_numeric_vars__i,
                " = xx)"
              )
            }
          }
        }
        data <-
          data %>% dplyr::mutate_at(cov_numeric_vars__,
                                    set_numeric_cov_at,
                                    numeric_cov_at)
        
        # This is good but get.data.grid is called twice - why?
        # hence this cat("\n"... is printed twice
        # cat("Continous covariate(s) set at:\n")
        # for (cov_numeric_vars__i in cov_numeric_vars__) {
        #   cat("\n", cov_numeric_vars__i, "at",
        #       unique(data[[cov_numeric_vars__i]]))
        # }
      }
    }
    
    
    if (!is.null(yvar)) {
      if (yvar %in% colnames(data)) {
        relocate_vars <- c(yvar, relocate_vars)
      }
    }
    data %>% dplyr::relocate(dplyr::all_of(relocate_vars)) %>% data.frame()
    
  }
  
  ########
  
  
  
  i_data <-
    function(model,
             newdata,
             resp = NULL,
             cov_factor_vars = NULL,
             cov_numeric_vars = NULL,
             aux_variables = NULL,
             levels_id = NULL,
             ipts = NULL,
             xrange = NULL) {
      if (is.null(resp)) {
        resp_rev_ <- resp
      } else if (!is.null(resp)) {
        resp_rev_ <- paste0("_", resp)
      }
      xvar_ <- paste0('xvar', resp_rev_)
      yvar_ <- paste0('yvar', resp_rev_)
      groupvar_ <- paste0('groupvar', resp_rev_)
      xvar <- model$model_info[[xvar_]]
      yvar <- model$model_info[[yvar_]]
      
      hierarchical_ <- paste0('hierarchical', resp_rev_)
      if (is.null(levels_id)) {
        IDvar <- model$model_info[[groupvar_]]
        if (!is.null(model$model_info[[hierarchical_]])) {
          IDvar <- model$model_info[[hierarchical_]]
        }
      } else if (!is.null(levels_id)) {
        IDvar <- levels_id
      }
      
      uvarby <- model$model_info$univariate_by$by
      if (!is.na(uvarby))
        cov_factor_vars <- c(uvarby, cov_factor_vars)
      
      
      
      # this idatafunction i.e., 'm1'
      if (idata_method == 'm1') {
        idatafunction <- function(.x,
                                  xvar,
                                  IDvar,
                                  nmy,
                                  xrange,
                                  set_xrange,
                                  aux_var = NULL) {
          index__x <- NA
          exdata <-
            function(x,
                     id,
                     idmat,
                     nmy,
                     xrange,
                     set_xrange,
                     aux_var) {
              n <- round(nmy * diff(range(x)))
              npt <- n / diff(range(x))
              
              extage <- apply(idmat, 1, function(x1) {
                index__x <-
                  id == x1
                
                if (is.null(xrange)) {
                  id.x <- x[index__x]
                }
                if (!is.null(xrange)) {
                  if (length(xrange) == 1) {
                    if (xrange == 1 & is.null(set_xrange))
                      id.x <- x
                    if (xrange == 2 &
                        !is.null(set_xrange))
                      id.x <- set_xrange
                  }
                  if (length(xrange) == 2) {
                    id.x <- set_xrange
                  }
                }
                
                nt <- floor(npt * diff(range(id.x))) + 1
                newx <- seq(min(id.x), max(id.x), length = nt)
                
                newid <-
                  rep(x1, nt)
                extx <- data.frame(x = newx, id = newid)
                colnames(extx) <- c("x", "id")
                extx
              })
              df <-
                extage[[1]][FALSE,]
              for (dft in extage)
                df <- rbind(df, dft)
              df
            }
          
          inidnull <- FALSE
          if(is.null(.x[[IDvar]])) {
            inidnull <- TRUE
            .x[[IDvar]] <- unique(levels(newdata[[IDvar]]))[1]
          }
          
          out <- exdata(
            x = .x[[xvar]],
            id = .x[[IDvar]],
            idmat = matrix(unique(.x[[IDvar]], ncol = 1)),
            nmy = nmy,
            xrange = xrange,
            set_xrange = set_xrange,
            aux_var = aux_var
          )
          
          out <- out %>% dplyr::rename(!!IDvar := 'id') %>% data.frame()
          
          if(inidnull) out <- out %>% dplyr::select(-dplyr::all_of(IDvar))
          
          idxx <- NULL
          if (!is.null(aux_var)) {
            aux_varx <- c(aux_var, IDvar)
            newx. <- .x %>% dplyr::select(dplyr::all_of(aux_varx)) %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(IDvar))) %>%
              dplyr::mutate(idxx = dplyr::row_number()) %>%
              dplyr::ungroup()
            outx. <- out %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(IDvar))) %>%
              dplyr::mutate(idxx = dplyr::row_number()) %>%
              dplyr::ungroup()
            out <- outx. %>% dplyr::left_join(., newx.,
                                              by = c(IDvar, 'idxx')) %>%
              dplyr::select(-idxx) %>% data.frame()
          }
          
          out 
        } # end idatafunction -> m1
        
       
        
        if (!is.null(xrange)) {
          if (length(xrange) < 1 | length(xrange) > 2) {
            stop(
              "Argument xrange should be either NULL, numeric value 1 or 2",
              "\n ",
              "or else a paired values indicating the range e.g., c(6, 20)"
            )
          }
        }
        
        if (!is.null(xrange)) {
          if (length(xrange) == 1) {
            if (xrange == 1)
              set_xrange <- NULL
            if (xrange == 2)
              set_xrange <- range(newdata[[xvar]])
          }
          if (length(xrange) == 2) {
            set_xrange <- xrange
          }
        }
        
        if (is.null(xrange))
          set_xrange <- NULL
       
        
        if (is.null(model$model_info[[hierarchical_]])) {
          if (!is.null(ipts) & is.null(cov_factor_vars)) {
            newdata %>% dplyr::arrange(IDvar, xvar) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  IDvar = IDvar,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                )
              ) %>%
              dplyr::rename(!!xvar := 'x') %>%
              dplyr::mutate(!!IDvar := as.factor(eval(parse(text = IDvar)))) %>%
              dplyr::relocate(dplyr::all_of(IDvar), dplyr::all_of(xvar)) %>%
              data.frame() -> newdata
          } else if (!is.null(ipts) & !is.null(cov_factor_vars)) {
            newdata %>% dplyr::arrange(IDvar, xvar) %>%
              dplyr::group_by(dplyr::across(dplyr::all_of(cov_factor_vars))) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  IDvar = IDvar,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                )
              ) %>%
              dplyr::rename(!!xvar := 'x') %>%
              dplyr::mutate(!!IDvar := as.factor(eval(parse(text = IDvar)))) %>%
              dplyr::relocate(dplyr::all_of(IDvar), dplyr::all_of(xvar)) %>%
              data.frame() -> newdata
          }
        } # if(is.null(model$model_info[[hierarchical_]]))
        
       
        
        multiNewVar <- function(df, df2, varname) {
          df %>% dplyr::mutate(.,!!varname := df2[[varname]])
        }
        
        if (!is.null(model$model_info[[hierarchical_]])) {
          if (!is.null(ipts) & is.null(cov_factor_vars)) {
            IDvar_ <- IDvar[1]
            higher_ <- IDvar[2:length(IDvar)]
            arrange_by <- c(IDvar_, xvar)
            cov_factor_vars_by <- c(higher_, cov_factor_vars)
            newdata_o <- newdata
            newdata <-
              newdata %>% dplyr::arrange(!!as.symbol(arrange_by)) %>%
              dplyr::group_by(
                dplyr::across(dplyr::all_of(cov_factor_vars_by))) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  IDvar = IDvar_,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                ),
                .keep = F
              ) %>%
              dplyr::rename(!!xvar := 'x') %>%
              data.frame()
            
            for (i in IDvar) {
              newdata <- newdata %>% multiNewVar(df = .,
                                                 df2 = newdata,
                                                 varname = i)
            }
            
            newdata %>% dplyr::relocate(dplyr::all_of(IDvar), 
                                        dplyr::all_of(xvar)) %>%
              data.frame() -> newdata
          }
          
          if (!is.null(ipts) & !is.null(cov_factor_vars)) {
            IDvar_ <- IDvar[1]
            higher_ <- IDvar[2:length(IDvar)]
            arrange_by <- c(IDvar_, xvar)
            # cov_factor_vars_by <- c(higher_, cov_factor_vars)
            if(length(IDvar) > 1) {
              cov_factor_vars_by <- c(higher_, cov_factor_vars)
            } else {
              cov_factor_vars_by <- c(cov_factor_vars)
            }
            newdata <-
              newdata %>% dplyr::arrange(!!as.symbol(arrange_by)) %>%
              dplyr::group_by(
                dplyr::across(dplyr::all_of(cov_factor_vars_by))) %>%
              dplyr::group_modify(
                ~ idatafunction(
                  .x,
                  xvar = xvar,
                  IDvar = IDvar_,
                  nmy = ipts,
                  xrange = xrange,
                  set_xrange = set_xrange,
                  aux_var = aux_variables
                )
              ) %>%
              dplyr::rename(!!xvar := 'x') %>% data.frame()
            for (i in IDvar) {
              newdata <- newdata %>% multiNewVar(df = .,
                                                 df2 = newdata,
                                                 varname = i)
            }
            newdata %>% dplyr::relocate(dplyr::all_of(IDvar), 
                                        dplyr::all_of(xvar)) %>%
              data.frame() -> newdata
          }
        } # if(!is.null(model$model_info[[hierarchical_]])) {
        
      } # end of if(idata_method == 'm1') {
      
      
      # this is get_idata i.e., 'm2'
      if (idata_method == 'm2') {
        if (!is.null(ipts)) {
          # for 3 or more level data, idvar shoud be first of vector
          if (is.null(model$model_info[[hierarchical_]])) {
            IDvar_for_idata <- IDvar
          } else if (!is.null(model$model_info[[hierarchical_]])) {
            IDvar_for_idata <- IDvar[1]
          }
          newdata <-
            get_idata(
              model = NULL,
              newdata = newdata,
              idVar = IDvar_for_idata,
              timeVar = xvar,
              times = NULL,
              length.out = ipts,
              xrange = xrange,
              keeplevels = FALSE, 
              asdf = FALSE,
              get.newdata.call = TRUE
            )
        }
        
      } # end if(idata_method == 'm2') {
      
      
      
      if (is.null(ipts)) {
        newdata <- newdata
      }
      
      
      ##################################################
      # prepare_data2 changes 
      # if (!is.null(ipts)) {
      #   # outliers must be NULL
      #   # Because these has already been taken care of by get.newdata
      #   if (!is.na(model$model_info$univariate_by$by)) {
      #     if (is.symbol(model$model_info$call.bgmfit$y)) {
      #       setorgy <- deparse(model$model_info$call.bgmfit$y)
      #     } else if (is.list(model$model_info$call.bgmfit$y)) {
      #       setorgy <- unname(unlist(model$model_info$call.bgmfit$y))
      #       if (is.symbol(setorgy))
      #         setorgy <- deparse(setorgy)
      #     } else {
      #       setorgy <- model$model_info$call.bgmfit$y
      #     }
      #   }
      #   
      #   if (is.na(model$model_info$univariate_by$by)) {
      #     setorgy <- model$model_info$ys
      #   }
      #   
      #   if (idata_method == 'm1') {
      #     newdata <- prepare_data2(
      #       # newdata <- prepare_data(
      #       data = newdata,
      #       x = model$model_info$xs,
      #       y = setorgy,
      #       id = model$model_info$ids,
      #       uvarby = model$model_info$univariate_by$by,
      #       mvar = model$model_info$multivariate$mvar,
      #       xfuns = model$model_info$xfuns,
      #       yfuns = model$model_info$yfuns,
      #       outliers = NULL) # model$model_info$outliers
      #   } # if (idata_method == 'm1') {
      #   
      # }
      ##################################################
      
      # dont repeate prepare_data2 call 
      if (!is.null(ipts)) {
        # prepare_data2 with model = model will get all the necessary info
      #  newdata <- prepare_data2(data = newdata, model = model)
      }
      
      
      ##################################################
      # prepare_data2 changes 
      # if (!is.na(model$model_info$univariate_by$by)) {
      #   if(idata_method == 'm1') {
      #     sortbylayer <- NA
      #     unique_names <- unique(names(newdata))
      #     newdata <- newdata %>% dplyr::select(dplyr::all_of(unique_names))
      #     newdata <- newdata %>%
      #       dplyr::mutate(sortbylayer =
      #                       forcats::fct_relevel(!!as.name(uvarby),
      #                                            (levels(
      #                                              !!as.name(uvarby)
      #                                            )))) %>%
      #       dplyr::arrange(sortbylayer) %>%
      #       dplyr::mutate(!!as.name(IDvar) :=
      #                       factor(!!as.name(IDvar),
      #                              levels =
      #                                unique(!!as.name(IDvar)))) %>%
      #       dplyr::select(-sortbylayer)
      #   } # if(idata_method == 'm1') {
      # }
      ##################################################
      
      
      newdata.oo <- get.data.grid(
        data = newdata,
        xvar = xvar,
        yvar = yvar,
        IDvar = IDvar,
        cov_numeric_vars = cov_numeric_vars,
        numeric_cov_at = numeric_cov_at,
        aux_variables = aux_variables,
        uvarby = uvarby
      )
      
      
      #         newdata <- check_newdata_args(model, newdata, IDvar, resp)
      
      j_b_names <- intersect(names(newdata), names(newdata.oo))
      j_b_names__ <- c(j_b_names, cov_numeric_vars)
      j_b_names__ <- unique(j_b_names__)
      
      if(idata_method == 'm1') {
        newdata <-
          newdata %>% 
          dplyr::left_join(., 
                           newdata.oo %>%
                             dplyr::select(dplyr::all_of(j_b_names__)),
                           by = j_b_names)
        
      } else if(idata_method == 'm2') {
        if(length(IDvar) > 1) {
          name_hypothetical_id <- IDvar[1]
        } else {
          name_hypothetical_id <- IDvar
        }
        if(length(unique(newdata[[name_hypothetical_id]])) > 1) {
          newdata <-
            newdata %>% 
            dplyr::left_join(., 
                             newdata.oo %>%
                               dplyr::select(dplyr::all_of(j_b_names__)),
                             by = j_b_names)
        } else if(length(unique(newdata[[name_hypothetical_id]])) == 1) {
          newdata <- newdata
        }
      } # else if(idata_method == 'm2') {
      
      
      return(newdata)
    }
  
  
  newdata <- i_data(
    model,
    newdata,
    resp = resp,
    cov_factor_vars = cov_factor_vars,
    cov_numeric_vars = cov_numeric_vars,
    aux_variables = aux_variables,
    levels_id = levels_id,
    ipts = ipts,
    xrange = xrange
  )
  
  
  list_c[['xvar']] <- xvar
  list_c[['yvar']] <- yvar
  list_c[['IDvar']] <- IDvar
  list_c[['cov_vars']] <- cov_vars
  list_c[['cov_factor_vars']] <- cov_factor_vars
  list_c[['cov_numeric_vars']] <- cov_numeric_vars
  list_c[['groupby_fstr']] <- groupby_fstr
  list_c[['groupby_fistr']] <- groupby_fistr
  
  attr(newdata, 'list_c') <- list_c
  
  
  
  return(newdata)
} # get.newdata




