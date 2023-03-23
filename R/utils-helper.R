

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
  cov_sigma_ <- paste0('cov_sigma', resp_rev_)
  uvarby <- model$model_info$univariate_by
  
  
  
  if (is.null(newdata)) {
    # newdata <- eval.parent(model$model_info$call.bsitar$data)
    newdata <- model$model_info$bsitar.data
  } else {
    newdata = newdata
  }
  
  if(!is.null(yfun)) {
    # if(yfun == 'log') newdata[[yvar]] <- log(newdata[[yvar]])
    # if(yfun == 'sqrt') newdata[[yvar]] <- sqrt(newdata[[yvar]] )
  }
  
  
 
  newdata <- prepare_data(data = newdata,
                          x = model$model_info$xs,
                          y = model$model_info$ys,  
                          id = model$model_info$ids,
                          uvarby = model$model_info$univariate_by,
                          mvar = model$model_info$multivariate,
                          xfuns = model$model_info$xfuns,
                          yfuns = model$model_info$yfuns,
                          outliers = model$model_info$outliers)
  

  if(!is.na(model$model_info$univariate_by)) {
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
  cov_sigma_vars <-  model$model_info[[cov_sigma_]]
  if(!is.null(cov_vars))  cov_vars <- covars_extrcation(cov_vars)
  
  if(!is.null(cov_sigma_vars))  cov_sigma_vars <- covars_extrcation(cov_sigma_vars)
  
  
  cov_factor_vars <- intersect(cov_vars, factor_vars)
  cov_numeric_vars <- intersect(cov_vars, numeric_vars)
  groupby_fstr <- c(cov_factor_vars)
  groupby_fistr <- c(IDvar, cov_factor_vars)
  
  cov_sigma_factor_vars <- intersect(cov_sigma_vars, factor_vars)
  cov_sigma_numeric_vars <- intersect(cov_sigma_vars, numeric_vars)
  
  if(identical(cov_factor_vars, character(0))) cov_factor_vars <- NULL
  if(identical(cov_numeric_vars, character(0))) cov_numeric_vars <- NULL
 
  if(identical(cov_sigma_factor_vars, character(0))) cov_sigma_factor_vars <- NULL
  if(identical(cov_sigma_numeric_vars, character(0))) cov_sigma_numeric_vars <- NULL
  
  
  # Merge here a b c covariate with sigma co variate
  # IMP: Note that groupby_fstr and groupby_fistr are stil  a b c covariate 
  # This way, plot_bsitar and gparameters will not produce sigam cov specific
  # curves and g parameters
  
  cov_factor_vars <- c(cov_factor_vars, cov_sigma_factor_vars)
  cov_numeric_vars <- c(cov_numeric_vars, cov_sigma_numeric_vars)
  

  # groupby_fstr <- c(groupby_fstr, groupby_fstr)
  # groupby_fistr <- c(groupby_fstr, groupby_fstr)
  
  
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
      # This cov_numeric_vars__ excludes sigma covariates
      
      # cov_numeric_vars__ <-
      #   cov_numeric_vars[!grepl(paste(xvar, collapse = "|"), cov_numeric_vars)]
     
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
        
        # This is good but get.data.grid is called twice - why? 
        # hence this cat("\n"... is printed twice
        
        # cat("Continous covariate(s) set at:\n")
        # for (cov_numeric_vars__i in cov_numeric_vars__) {
        #   cat("\n", cov_numeric_vars__i, "at",
        #       unique(data[[cov_numeric_vars__i]]))
        # }
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
      
      
      
      idatafunction <- function(.x, xvar, IDvar, nmy, xrange, set_xrange, aux_var = NULL) {
        index__x <- NA
        exdata <- function(x, id, idmat, nmy, xrange, set_xrange, aux_var) {
          n <- round(nmy * diff(range(x)))
          npt <- n / diff(range(x))
          
          extage <- apply(idmat, 1, function(x1) {
            index__x <-
              id == x1
            
            if(is.null(xrange)) {
              id.x <- x[index__x]
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
        out <- exdata(.x[[xvar]], .x[[IDvar]], matrix(unique(.x[[IDvar]], ncol = 1)),
               nmy = nmy, xrange = xrange, set_xrange = set_xrange, aux_var = aux_var)
        
        idxx <- NULL
        
        if(!is.null(aux_var)) {
          aux_varx <- c(aux_var, IDvar)
          newx. <- .x %>% dplyr::select(all_of(aux_varx)) %>% 
            dplyr::group_by(across(all_of(IDvar))) %>% 
            dplyr::mutate(idxx = dplyr::row_number()) %>% 
            dplyr::ungroup()
          outx. <- out %>% 
            dplyr::group_by(across(all_of(IDvar))) %>% 
            dplyr::mutate(idxx = dplyr::row_number()) %>% 
            dplyr::ungroup()
          out <- outx. %>% dplyr::left_join(., newx., by = c(IDvar, 'idxx')) %>% 
            dplyr::select(-idxx) %>% data.frame()
        }
        
        out 
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
          newdata_o <- newdata
          newdata <- newdata %>% dplyr::arrange(!!as.symbol(arrange_by)) %>%
            dplyr::group_by(across(all_of(cov_factor_vars_by))) %>%
            dplyr::group_modify(~ idatafunction(
              .x,
              xvar = xvar,
              IDvar = IDvar_,
              nmy = ipts,  
              xrange = xrange, 
              set_xrange = set_xrange,
              aux_var = "agem"
            ), .keep = F) %>%
            dplyr::rename(!!xvar := 'x') %>% 
            data.frame()  
     
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
      
      
      if(is.null(ipts)) {
        newdata <- newdata
      } 
      
    
      if(!is.null(ipts)) {
        # outliers must be NULL 
        # because these has already been taken care of by get.newdata
        newdata <- prepare_data(data = newdata,
                                x = model$model_info$xs,
                                y = model$model_info$ys,  
                                id = model$model_info$ids,
                                uvarby = model$model_info$univariate_by,
                                mvar = model$model_info$multivariate,
                                xfuns = model$model_info$xfuns,
                                yfuns = model$model_info$yfuns,
                                outliers = NULL) # model$model_info$outliers
      }
      
      # model$model_info$prepare_data
      if(!is.na(model$model_info$univariate_by)) { # !is.null(ipts) & 
        # newdata <- prepare_data(data = newdata,
        #                         x = model$model_info$xvar,
        #                         y = model$model_info$ys,  
        #                         id = model$model_info$groupvar_,
        #                         uvarby = model$model_info$univariate_by,
        #                         mvar = model$model_info$multivariate,
        #                         xfuns = model$model_info$xfuns,
        #                         yfuns = model$model_info$yfuns,
        #                         outliers = model$model_info$outliers)
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




########################
########################


get_args_ <- function(arguments, xcall) {
  `%!in%` <- Negate(`%in%`)
  f_bsitar_arg <- formals(paste0(xcall, '.bsitar'))
  nf_bsitar_arg_names <-
    intersect(names(arguments), names(f_bsitar_arg))
  arguments <-
    c(arguments, f_bsitar_arg[names(f_bsitar_arg) %!in% nf_bsitar_arg_names])
  arguments
}


########################
########################

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




########################
########################

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
  
  if(!is.null(resp)) {
    if(!resp %in% model$model_info[['ys']]) {
      stop("Response should be one of the following: ", 
           paste(model$model_info[['ys']], collapse = " "),
           "\n ",
           " but you have specified: ", resp)
    }
  }
  
}





########################
########################

setup_higher_priors <- function(new_prior_list) {
  o_l <- list()
  ixi = 0
  group_ <- class_ <- nlpar_ <- resp_ <- cor_check <- sd_check <- NA
  group_ <- class_ <- nlpar_ <- resp_ <- cor_check <- sd_check <- c()
  for (new_prior_listi in 1:length(new_prior_list)) {
    ixi <- ixi + 1
    pstr <- new_prior_list[[new_prior_listi]]
    if(is.null(pstr[['prior']])) prior_i <- '' else prior_i <- pstr[['prior']]
    if(is.null(pstr[['group']])) group_i <- '' else group_i <- pstr[['group']]
    if(is.null(pstr[['class']])) class_i <- '' else class_i <- pstr[['class']]
    if(is.null(pstr[['nlpar']])) nlpar_i <- '' else nlpar_i <- pstr[['nlpar']]
    if(is.null(pstr[['coef']]))  coef_i  <- '' else coef_i  <- pstr[['coef']]
    if(is.null(pstr[['resp']]))  resp_i  <- '' else resp_i  <- pstr[['resp']]
    if(is.null(pstr[['lb']]))    lb_i    <- '' else lb_i    <- pstr[['lb']]
    if(is.null(pstr[['ub']]))    ub_i    <- '' else ub_i <- pstr[['ub']]
    if(is.null(pstr[['dpar']]))  dpar_i  <- '' else dpar_i  <- pstr[['dpar']]
    if(is.null(pstr[['prior']])) prior_i <- '' else prior_i <- pstr[['prior']]
    
    group_ <- c(group_, group_i)
    class_ <- c(class_, class_i)
    nlpar_ <- c(nlpar_, nlpar_i)
    resp_ <- c(nlpar_, resp_i)
    if(class_i == 'sd') sd_check <- c(sd_check, group_i)
    if(class_i == 'cor') cor_check <- c(cor_check, group_i)
    
    if(lb_i == '' & ub_i == '' ) {
      o_l[[ixi]] <- prior_string(prior_i,
                                 group = group_i, 
                                 class = class_i, 
                                 nlpar = nlpar_i,
                                 coef = coef_i,
                                 resp = resp_i,
                                 # lb = lb_i,
                                 # ub = ub_i,
                                 dpar = dpar_i)
    } else if(lb_i != '' | ub_i != '' ) {
      o_l[[ixi]] <- prior_string(prior_i,
                                 group = group_i, 
                                 class = class_i, 
                                 nlpar = nlpar_i,
                                 coef = coef_i,
                                 resp = resp_i,
                                 # lb = lb_i,
                                 # ub = ub_i,
                                 dpar = dpar_i)
    }
    
    
  }
  o_l %>%  do.call(rbind, .)
}


########################
########################

insert_new_priors <- function(setdf_1, setdf_2) {
  index__x <- NA
  valid__x <- NA
  setdf_1 <- 
    setdf_1 %>% dplyr::mutate(index__x = interaction(class, coef, group, nlpar)) %>% 
    dplyr::mutate(order = dplyr::row_number()) %>% 
    dplyr::arrange(index__x)
  
  setdf_2 <- 
    setdf_2 %>% dplyr::mutate(index__x = interaction(class, coef, group, nlpar)) %>% 
    dplyr::mutate(order = dplyr::row_number()) %>% 
    dplyr::arrange(index__x)
  
  vi_1 <- setdf_1 %>% dplyr::mutate(valid__x = ifelse(!(class == 'sd' & coef == ""), 1, 0)) %>% 
    data.frame() %>% dplyr::filter(valid__x == 1) %>% dplyr::select(index__x) %>% unlist() %>% 
    droplevels()
  
  setdf_1 <- setdf_1 %>% dplyr::filter(!(class == 'sd' & coef == ""))
  setdf_2 <- setdf_2 %>% dplyr::filter(!(class == 'sd' & coef == ""))
  setdf_3 <- setdf_1[setdf_1$index__x %in% vi_1,]
  setdf_2 <- setdf_2 %>% dplyr::filter(!index__x %in% vi_1)
  setdf_4 <- rbind(setdf_2, setdf_3)
  # setdf_4 <- setdf_4 %>% arrange(order) %>% select(-c(index__x, order))
  setdf_4 <- setdf_4 %>% dplyr::select(-c(index__x, order))
  setdf_4
}



