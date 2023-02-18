




get.newdata <- function(model, newdata, resp, 
                        numeric_cov_at = NULL, 
                        levels_id = NULL,
                        ipts = NULL,
                        xrange = NULL) {
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  
  validate_response(model, resp)
  
  list_c <- list()
  xvar_ <- paste0('xvar', resp_rev_)
  yvar_ <- paste0('yvar', resp_rev_)
  groupvar_ <- paste0('groupvar', resp_rev_)
  xvar <- model$model_info[[xvar_]]
  yvar <- model$model_info[[yvar_]]
  hierarchical_ <- paste0('hierarchical', resp_rev_)
 
  if(is.null(levels_id)) {
    IDvar <- model$model_info[[groupvar_]]
    if(!is.null(model$model_info[[hierarchical_]])) {
      IDvar <- model$model_info[[hierarchical_]]
    }
  } else if(!is.null(levels_id)) {
    IDvar <- levels_id
  }
  
  xfun_ <- paste0('xfun', resp_rev_)
  yfun_ <- paste0('yfun', resp_rev_)
  xfun <- model$model_info[[xfun_]]
  yfun <- model$model_info[[yfun_]]
  
  
  cov_ <- paste0('cov', resp_rev_)
  uvarby <- model$model_info$univariate_by
  
  # if (is.null(newdata) & is.na(model$model_info$univariate_by)) {
  #   newdata <- model$data
  # }
  
  
  if (is.null(newdata)) {
    newdata <- eval.parent(model$model_info$call.bsitar$data)
  } else {
    newdata = newdata
  }
  
  
  
  if(!is.null(yfun)) {
    # if(yfun == 'log') newdata[[yvar]] <- log(newdata[[yvar]])
    # if(yfun == 'sqrt') newdata[[yvar]] <- sqrt(newdata[[yvar]] )
  }
  
  
  
  if(!is.na(model$model_info$univariate_by)) {
    if (is.null(newdata)) {
      # newdata_ <- eval.parent(model$model_info$call.bsitar$data)
      newdata_ <- newdata
    } else  if (!is.null(newdata)) {
      newdata_ <- newdata
    }
    
    
    
    newdata <- model$model_info$make_bsitar_data(newdata_,
                                                 model$model_info$univariate_by,
                                                 model$model_info$multivariate,
                                                 model$model_info$ys,
                                                 model$model_info$xfuns,
                                                 model$model_info$yfuns)
  
    sortbylayer <- NA
    newdata <- newdata %>%
      dplyr::mutate(sortbylayer =
                      forcats::fct_relevel(!!as.name(uvarby),
                                           (levels(
                                             !!as.name(uvarby)
                                           )))) %>%
      dplyr::arrange(sortbylayer) %>%
      dplyr::mutate(!!as.name(IDvar) := factor(!!as.name(IDvar),
                                               levels = 
                                                 unique(!!as.name(IDvar)))) %>% 
      dplyr::select(-sortbylayer)
    
    subindicatorsi <- model$model_info$subindicators[grep(resp, model$model_info$ys)]
    list_c[['subindicatorsi']] <- subindicatorsi
    list_c[['uvarby']] <- uvarby
  }
  
  covars_extrcation <- function(str) {
    str <- gsub("[[:space:]]", "", str)
    for (ci in c("*", "+", ":")) {
      str <- gsub(ci, ' ', str, fixed = T)
    }
    str <- strsplit(str, " ")[[1]]
    str
  }
  
  factor_vars <- names(newdata[sapply(newdata, is.factor)])
  numeric_vars <- names(newdata[sapply(newdata, is.numeric)])
  cov_vars <-  model$model_info[[cov_]]
  if(!is.null(cov_vars))  cov_vars <- covars_extrcation(cov_vars)
  
  cov_factor_vars <- intersect(cov_vars, factor_vars)
  cov_numeric_vars <- intersect(cov_vars, numeric_vars)
  groupby_fstr <- c(cov_factor_vars)
  groupby_fistr <- c(IDvar, cov_factor_vars)
  
  if(identical(cov_factor_vars, character(0))) cov_factor_vars <- NULL
  if(identical(cov_numeric_vars, character(0))) cov_numeric_vars <- NULL
  
  if(!is.na(model$model_info$univariate_by)) {
    groupby_fstr <- c(uvarby, groupby_fstr)
    groupby_fistr <- c(uvarby, groupby_fistr)
  }
  
  #########
  
  
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
  
  
  #########
  
  get.data.grid <- function(data, 
                            xvar = xvar, 
                            yvar = yvar, 
                            IDvar = IDvar, 
                            cov_numeric_vars = cov_numeric_vars,
                            numeric_cov_at = numeric_cov_at,
                            uvarby = uvarby) {
    
    if(!is.null(IDvar)) relocate_vars <- c(xvar, IDvar)
    if(is.null(IDvar)) relocate_vars <- c(xvar)
    if(!is.na(uvarby)) relocate_vars <- c(relocate_vars, uvarby)
    if(!is.null(cov_numeric_vars)) {
      cov_numeric_vars__ <-
        cov_numeric_vars[!grepl(paste(xvar, collapse = "|"), cov_numeric_vars)]
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
                " Please use the argument 'numeric_cov_at' correctly as follows:",
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
                                    set_numeric_cov_at, numeric_cov_at)
        cat("Continous covariate(s) set at:")
        for (cov_numeric_vars__i in cov_numeric_vars__) {
          cat("\n", cov_numeric_vars__i, "at",
              unique(data[[cov_numeric_vars__i]]))
        }
      }
    }
    if(!is.null(yvar)) {
      if (yvar %in% colnames(data)) {
        relocate_vars <- c(yvar, relocate_vars)
      }
    }
    data %>% dplyr::relocate(all_of(relocate_vars)) %>% data.frame()
  }
  
  ########
  
  i_data <-
    function(model,
             newdata,
             resp = NULL, 
             cov_factor_vars = NULL,
             cov_numeric_vars = NULL,
             levels_id = NULL, 
             ipts = NULL,
             xrange = NULL
    ) {
      
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
      if(is.null(levels_id)) {
        IDvar <- model$model_info[[groupvar_]]
        if(!is.null(model$model_info[[hierarchical_]])) {
          IDvar <- model$model_info[[hierarchical_]]
        }
      } else if(!is.null(levels_id)) {
        IDvar <- levels_id
      }
      
      uvarby <- model$model_info$univariate_by
      if(!is.na(uvarby)) cov_factor_vars <- c(uvarby, cov_factor_vars)
      
      
      
      
      if(is.null(ipts)) newdata <- newdata
      
      
      idatafunction <- function(.x, xvar, IDvar, nmy, xrange, set_xrange) {
        exdata <- function(x, id, idmat, nmy, xrange, set_xrange) {
          n <- round(nmy * diff(range(x)))
          npt <- n / diff(range(x))
          
          extage <- apply(idmat, 1, function(x1) {
            index <-
              id == x1
            
            if(is.null(xrange)) {
              id.x <- x[index]
            }
            if(!is.null(xrange)) {
              if(length(xrange) == 1) {
                if(xrange == 1 & is.null(set_xrange)) id.x <- x
                if(xrange == 2 & !is.null(set_xrange)) id.x <- set_xrange
              }
              if(length(xrange) == 2) {
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
            extage[[1]][FALSE, ]
          for (dft in extage)
            df <- rbind(df, dft)
          df
        }
        
        exdata(.x[[xvar]], .x[[IDvar]], matrix(unique(.x[[IDvar]], ncol = 1)),
               nmy = nmy, xrange = xrange, set_xrange = set_xrange)
      }
      
      if(!is.null(xrange)) {
        if(length(xrange) < 1 | length(xrange) > 2) {
          stop("Argument xrange should be either NULL, numeric value 1 or 2",
               "\n ",
               "or else a paired values indicating the range e.g., c(6, 20)")
        }
      }
      
      if(!is.null(xrange)) {
        if(length(xrange) == 1) {
          if(xrange == 1) set_xrange <- NULL
          if(xrange == 2) set_xrange <- range(newdata[[xvar]])
        }
        if(length(xrange) == 2) {
          set_xrange <- xrange
        }
      }
      
      if(is.null(xrange)) set_xrange <- NULL
      
      if(is.null(model$model_info[[hierarchical_]])) {
        if (!is.null(ipts) & is.null(cov_factor_vars)) {
          newdata %>% dplyr::arrange(IDvar, xvar) %>% 
            dplyr::group_modify(~ idatafunction(
              .x,
              xvar = xvar,
              IDvar = IDvar,
              nmy = ipts,
              xrange = xrange, 
              set_xrange = set_xrange
            )) %>%
            dplyr::rename(!!xvar := 'x') %>%
            dplyr::mutate(!!IDvar := as.factor(eval(parse(text = IDvar)))) %>%
            dplyr::relocate(all_of(IDvar), all_of(xvar)) %>%
            data.frame() -> newdata
        } else if (!is.null(ipts) & !is.null(cov_factor_vars)) {
          newdata %>% dplyr::arrange(IDvar, xvar) %>%
            dplyr::group_by(across(all_of(cov_factor_vars))) %>%
            dplyr::group_modify(~ idatafunction(
              .x,
              xvar = xvar,
              IDvar = IDvar,
              nmy = ipts,
              xrange = xrange, 
              set_xrange = set_xrange
            )) %>%
            dplyr::rename(!!xvar := 'x') %>%
            dplyr::mutate(!!IDvar := as.factor(eval(parse(text = IDvar)))) %>%
            dplyr::relocate(all_of(IDvar), all_of(xvar)) %>%
            data.frame() -> newdata
        }
      } # if(is.null(model$model_info[[hierarchical_]]))
      
      
      multiNewVar <- function(df, df2, varname) {
        df %>% dplyr::mutate(., !!varname := df2[[varname]])
      }
     
      if(!is.null(model$model_info[[hierarchical_]])) {
        if (!is.null(ipts) & is.null(cov_factor_vars)) {
          IDvar_ <- IDvar[1]
          higher_ <- IDvar[2:length(IDvar)]
          arrange_by <- c(IDvar_, xvar)
          cov_factor_vars_by <- c(higher_, cov_factor_vars)
          newdata <- newdata %>% dplyr::arrange(!!as.symbol(arrange_by)) %>%
            dplyr::group_by(across(all_of(cov_factor_vars_by))) %>%
            dplyr::group_modify(~ idatafunction(
              .x,
              xvar = xvar,
              IDvar = IDvar_,
              nmy = ipts,  
              xrange = xrange, 
              set_xrange = set_xrange
            )) %>%
            dplyr::rename(!!xvar := 'x') %>% data.frame()
          for(i in IDvar) {
            newdata <- newdata %>% multiNewVar(df=., df2 = newdata, varname=i)
          } 
          newdata %>% dplyr::relocate(all_of(IDvar), all_of(xvar)) %>%
            data.frame() -> newdata
        } 
        
        if (!is.null(ipts) & !is.null(cov_factor_vars)) {
          IDvar_ <- IDvar[1]
          higher_ <- IDvar[2:length(IDvar)]
          arrange_by <- c(IDvar_, xvar)
          cov_factor_vars_by <- c(higher_, cov_factor_vars)
          newdata <- newdata %>% dplyr::arrange(!!as.symbol(arrange_by)) %>%
            dplyr::group_by(across(all_of(cov_factor_vars_by))) %>%
            dplyr::group_modify(~ idatafunction(
              .x,
              xvar = xvar,
              IDvar = IDvar_,
              nmy = ipts,
              xrange = xrange, 
              set_xrange = set_xrange
            )) %>%
            dplyr::rename(!!xvar := 'x') %>% data.frame()
          for(i in IDvar) {
            newdata <- newdata %>% multiNewVar(df=., df2 = newdata, varname=i)
          } 
          newdata %>% dplyr::relocate(all_of(IDvar), all_of(xvar)) %>%
            data.frame() -> newdata
        }
      } # if(!is.null(model$model_info[[hierarchical_]])) {
      
      
      
      
      # if(is.null(ipts)) newdata <- newdata
      # 
      if(!is.null(ipts) & !is.na(model$model_info$univariate_by)) {
        newdata <- model$model_info$make_bsitar_data(newdata,
                                                     model$model_info$univariate_by,
                                                     model$model_info$multivariate,
                                                     model$model_info$ys,
                                                     model$model_info$xfuns,
                                                     model$model_info$yfuns)
        sortbylayer <- NA
        newdata <- newdata %>%
          dplyr::mutate(sortbylayer =
                          forcats::fct_relevel(!!as.name(uvarby),
                                               (levels(
                                                 !!as.name(uvarby)
                                               )))) %>%
          dplyr::arrange(sortbylayer) %>%
          dplyr::mutate(!!as.name(IDvar) := factor(!!as.name(IDvar),
                                                   levels =
                                                     unique(!!as.name(IDvar)))) %>%
          dplyr::select(-sortbylayer)

      }
      
      newdata.oo <- get.data.grid(newdata, xvar, yvar, IDvar, cov_numeric_vars,
                                  numeric_cov_at, uvarby)
      
      # j_b_names <- names(newdata)
      j_b_names <- intersect(names(newdata), names(newdata.oo))
      j_b_names__ <- c(j_b_names, cov_numeric_vars)
      j_b_names__ <- unique(j_b_names__)
      
      
      
      newdata <-
        newdata %>% dplyr::left_join(., newdata.oo %>%
                                       dplyr::select(all_of(j_b_names__)),
                                     by = j_b_names)
      
      newdata 
    }
  
  
  newdata <- i_data(model, newdata, resp = resp,
                             cov_factor_vars = cov_factor_vars,
                             cov_numeric_vars = cov_numeric_vars, 
                    levels_id = levels_id,
                             ipts = ipts,
                    xrange = xrange)
  
  #########
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







get_args_ <- function(arguments, xcall) {
  `%!in%` <- Negate(`%in%`)
  f_bsitar_arg <- formals(paste0(xcall, '.bsitar'))
  nf_bsitar_arg_names <-
    intersect(names(arguments), names(f_bsitar_arg))
  arguments <-
    c(arguments, f_bsitar_arg[names(f_bsitar_arg) %!in% nf_bsitar_arg_names])
  arguments
}



get.cores <- function(cores.arg) {
  cores_ <- eval(cores.arg)
  if(!is.null(cores_)) {
    if(cores_ == "maximise") {
      max.cores <- 
        as.numeric(future::availableCores(methods = "system", omit = 0))
      if(max.cores < 1) max.cores <- 1
    } else if(cores_ == "optimize") {
      max.cores <- 
        as.numeric(future::availableCores(methods = "system", omit = 1))
      if(max.cores < 1) max.cores <- 1
    } else {
      max.cores <- eval(cores_)
    }
  } else if(is.null(cores_)) {
    max.cores <- NULL
  }
  
  if(!is.null(cores_)) {
    if(Sys.info()["sysname"] == "Windows") {
      .cores_ps <- 1
    } else {
      .cores_ps <- max.cores
    }
  } else if(is.null(cores_)) {
    .cores_ps <- 1
  }
  
  list(max.cores = max.cores, .cores_ps = .cores_ps)
}







validate_response <- function(model, resp = NULL) {
  if (model$model_info$nys == 1 & !is.null(resp)) {
    stop(
      "You have fit a univariate model",
      " but set resp option as ",
      resp,
      ".",
      "\n ",
      " The resp option should be appropriately set to NULL",
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
        " but dit not set the the resp options appropriately",
        " (which is NULL at present).",
        "\n ",
        " The response options are ",
        paste(model$model_info$ys, collapse = ", ")
      )
    }
    if (model$model_info$multivariate) {
      stop(
        "You have fit a multivariate model ",
        "\n ",
        " but dit not set the the resp options appropriately",
        " (which is NULL at present).",
        "\n ",
        " The response options are ",
        paste(model$model_info$ys, collapse = ", ")
      )
    }
  }
  if(!resp %in% model$model_info[['ys']]) {
    stop("Response should be one of the following: ", 
         paste(model$model_info[['ys']], collapse = " "),
         "\n ",
         " but you have specified: ", resp)
  }
}







make_bsitar_data <- function(data, uvarby, mvar = FALSE, 
                             ys, 
                             xfuns = NULL, yfuns = NULL) {
  
  org.data <- data
  if (!(is.na(uvarby) | uvarby == "NA")) {
    # org.ys <- ys[1]
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
      data[[l]] <- data[[ys[1]]]
    }
    unibyimat <-
      model.matrix( ~ 0 + eval(parse(text = uvarby)), data)
    subindicators <- paste0(uvarby, levels(data[[uvarby]]))
    colnames(unibyimat) <- subindicators
    ys <- levels(data[[uvarby]])
    data <- as.data.frame(cbind(data, unibyimat))
    # if (verbose) {
    #   resvcts_ <- levels(data[[uvarby]])
    #   resvcts <- paste0(resvcts_, collapse = " ")
    #   setmsgtxt <- paste0(
    #     "\n For univariate-by-subgroup model fitting for variable '",
    #     uvarby,
    #     "'",
    #     " (specified via 'univariate_by' argument)",
    #     "\n ",
    #     resvcts,
    #     " response vectors created based on the factor levels",
    #     "\n\n ",
    #     "Please check corresponding arguments list.",
    #     " E.g, df = list(4, 5) denotes that\n df = 4 is for ",
    #     resvcts_[1],
    #     ", and  df = 5 is for ",
    #     resvcts_[2],
    #     " (and similalry knots, priors, initials etc)",
    #     "\n\n ",
    #     "If it does't correspond correctly, then either reverse the list ",
    #     "arguments\n such as df = list(5, 4),",
    #     " or else reverse sort the order of factor levels"
    #   )
    #   if (displayit == 'msg') {
    #     message(setmsgtxt)
    #   } else if (displayit == 'col') {
    #     col <- setcolb
    #     cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    #   }
    # }
    attr(data, "ys") <- ys
    attr(data, "multivariate") <- FALSE
    attr(data, "uvarby") <- uvarby
    attr(data, "subindicators") <- subindicators
    data_out <- data
  } else if(mvar) { # not yet included
    for(myfunsi in 1:length(ys)) {
      mysi <- ys[[myfunsi]]
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
    data_out <- org.data
    attr(data, "ys") <- ys
    attr(data, "multivariate") <- TRUE
    attr(data, "uvarby") <- NULL
    attr(data, "subindicators") <- NULL
  } else {
    data_out <- org.data
    attr(data, "ys") <- ys
    attr(data, "multivariate") <- FALSE
    attr(data, "uvarby") <- NULL
    attr(data, "subindicators") <- NULL
  }
  return(data)
}

