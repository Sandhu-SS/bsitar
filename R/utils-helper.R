


get.newdata <- function(model, newdata, resp) {
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  list_c <- list()
  xvar_ <- paste0('xvar', resp_rev_)
  yvar_ <- paste0('yvar', resp_rev_)
  groupvar_ <- paste0('groupvar', resp_rev_)
  xvar <- model$model_info[[xvar_]]
  yvar <- model$model_info[[yvar_]]
  IDvar <- model$model_info[[groupvar_]]
  
  cov_ <- paste0('cov', resp_rev_)
  uvarby <- model$model_info$univariate_by
  
  if (is.null(newdata) & is.na(model$model_info$univariate_by)) {
    newdata <- model$data
  }
  
  if(!is.na(model$model_info$univariate_by)) {
    if (is.null(newdata)) {
      newdata_ <- eval.parent(model$model_info$call.bsitar$data)
    } else  if (!is.null(newdata)) {
      newdata_ <- newdata
    }
    newdata <- model$model_info$make_bsitar_data(newdata_,
                                                 model$model_info$univariate_by,
                                                 model$model_info$org.ycall)
    # print(head(newdata_))
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
  
  
  factor_vars <- names(newdata[sapply(newdata, is.factor)])
  numeric_vars <- names(newdata[sapply(newdata, is.numeric)])
  cov_vars <-  model$model_info[[cov_]]
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
  
  list_c[['xvar']] <- xvar
  list_c[['yvar']] <- yvar
  list_c[['IDvar']] <- IDvar
  list_c[['cov_vars']] <- cov_vars
  list_c[['cov_factor_vars']] <- cov_factor_vars
  list_c[['cov_numeric_vars']] <- cov_numeric_vars
  list_c[['groupby_fstr']] <- groupby_fstr
  list_c[['groupby_fistr']] <- groupby_fistr
  
  attr(newdata, 'list_c') <- list_c
  newdata 
}



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
      .cores_ps <- cores
    }
  } else if(is.null(cores_)) {
    .cores_ps <- 1
  }
  
  list(max.cores = max.cores, .cores_ps = .cores_ps)
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



ged.data.grid <- function(data, 
                          xvar = NULL, 
                          yvar = NULL, 
                          IDvar = NULL, 
                          cov_numeric_vars = NULL,
                          numeric_cov_at = NULL,
                          uvarby = NA) {
  
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
        data %>% dplyr::mutate_at(cov_numeric_vars__, set_numeric_cov_at, numeric_cov_at)
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



i_data <-
  function(model,
           newdata,
           newdata.oo,
           resp = NULL, 
           cov_factor_vars = NULL,
           cov_numeric_vars = NULL,
           ipts = NULL
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
    IDvar <- model$model_info[[groupvar_]] 
    uvarby <- model$model_info$univariate_by
    
    if(!is.na(uvarby)) cov_factor_vars <- c(uvarby, cov_factor_vars)
    
    
    idatafunction <- function(.x, xvar, IDvar, nmy) {
      exdata <- function(x, id, idmat, nmy = nmy) {
        n <- round(nmy * diff(range(x)))
        npt <- n / diff(range(x))
        extage <- apply(idmat, 1, function(x1) {
          index <-
            id == x1
          id.x <- x[index]
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
             nmy = nmy)
    }
    
    
    if (!is.null(ipts) & is.null(cov_factor_vars)) {
      newdata %>% dplyr::arrange(IDvar, xvar) %>% 
        dplyr::group_modify(~ idatafunction(
          .x,
          xvar = xvar,
          IDvar = IDvar,
          nmy = ipts
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
          nmy = ipts
        )) %>%
        dplyr::rename(!!xvar := 'x') %>%
        dplyr::mutate(!!IDvar := as.factor(eval(parse(text = IDvar)))) %>%
        dplyr::relocate(all_of(IDvar), all_of(xvar)) %>%
        data.frame() -> newdata
    }
    if(is.null(ipts)) newdata <- newdata
    
    if(!is.null(ipts) & !is.na(model$model_info$univariate_by)) {
      newdata <- model$model_info$make_bsitar_data(eval.parent(model$model_info$call.bsitar$data),
                                                   model$model_info$univariate_by,
                                                   model$model_info$org.ycall)
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
    
    newdata.oo <- ged.data.grid(newdata, xvar, yvar, IDvar, cov_numeric_vars,
                                numeric_cov_at, uvarby)
    
    j_b_names <- names(newdata)
    j_b_names__ <- c(j_b_names, cov_numeric_vars)
    newdata <-
      newdata %>% dplyr::left_join(., newdata.oo %>%
                                     dplyr::select(all_of(j_b_names__)),
                                   by = j_b_names)
    
    newdata 
  }




