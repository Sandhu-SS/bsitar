

# load("z.RData")

expose_optimize_fit <- function(optimize_fit, subset_list = NULL, expose_function = T) {
  optimize_fit_models <-  optimize_fit
  
  if(!is.null(subset_list)) {
    if(!is.numeric(subset_list)) stop("models must a numeric vector")
    optimize_fit_models <- optimize_fit_models[subset_list]
  } else {
    optimize_fit_models <- optimize_fit_models
  }
  
  m_list <- list()
  for (il in 1:length(optimize_fit_models)) {
    if(is.null(expose_function)) {
      m_list[[il]] <- optimize_fit_models[[il]]
    } else if(!is.null(expose_function)) {
      if(expose_function) {
        message("Exposing Stan function...", " (for model no. ", il, ")")
        m_list[[il]] <- expose_bsitar_functions(optimize_fit_models[[il]], 
                                                optimize_fit_models[[il]]$bmodel)
      } else if(!expose_function) {
        m_list[[il]] <- optimize_fit_models[[il]]
      }
    }
  }
  m_list <- m_list[!sapply(m_list, is.null)]
  m_list
}

# elist <- expose_optimize_fit(optimize_fit2$models, subset_list = c(1, 6))



plot_optimize_fit <- function(optimize_fit, subset_list = NULL, what = "plot", 
                              expose_function = F, print= T, ...) {
  
  optimize_fit_models <-  optimize_fit # optimize_fit$models
  
  dots <- list(...)
  
  for (i in names(dots)) {
    if(!i %in% formalArgs(plot_bsitar.bsitar)) stop("arguments must be be one of the following",
                                                    "\n ",
                                                    formalArgs(plot_bsitar.bsitar))
  }
  
  if(!is.null(subset_list)) {
    if(!is.numeric(subset_list)) stop("models must a numeric vector")
    optimize_fit_models <- optimize_fit_models[subset_list]
  } else {
    optimize_fit_models <- optimize_fit_models
  }
  
  
  
  m_list <- list()
  for (il in 1:length(optimize_fit_models)) {
    if(is.null(expose_function)) {
      m_list[[il]] <- optimize_fit_models[[il]]
    } else if(!is.null(expose_function)) {
      if(expose_function) {
        message("Exposing Stan function...", " (for model no. ", il, ")")
        m_list[[il]] <- expose_bsitar_functions(optimize_fit_models[[il]], 
                                                optimize_fit_models[[il]]$bmodel)
      } else if(!expose_function) {
        m_list[[il]] <- optimize_fit_models[[il]]
      }
    }
  }
  
  m_list <- m_list[!sapply(m_list, is.null)]
  
  
  
  nx <- function(.x, bx, args_) {
    message("Working on model no. ", .x)
    
    dots$model <- bx[[.x]]
    dots$... <- NULL
    
    if(is.null(what)) what <- 'plot'
    
    if(what == "plot") {
      out_ <- do.call(plot_bsitar, dots)
      title_ <- bx[[.x]]$model_info$optimization_info
      out_ <- out_ + ggplot2::labs(title = title_)
    }
    
    if(what == "gparameters") {
      out_ <- do.call(gparameters, dots)
    }
    
    if(!is.null(print)) {
      if(print) {
        print(out_)
      }
    }
    return(out_)
  }
  
  out <- purrr::map(1:length(m_list), ~nx(.x, m_list, args_))
  return(out)
}



# out_plots <- plot_optimize_fit(optimize_fit$models, subset_list = 6, resp = NULL)
# 
# 
# 
# out_plots <- plot_optimize_fit(elist, subset_list = NULL, resp = NULL)
