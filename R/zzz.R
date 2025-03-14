
.onAttach <- function(libname, pkgname) {
  
  invisible(suppressPackageStartupMessages(
    sapply(c("brms"),
           requireNamespace, quietly = TRUE)
  ))
  
  
  
  if("rethinking" %in% (.packages())){
    packageStartupMessage("Package 'rethinking' detached and unloaded as it creates conflict",
            " \nwith the current rstan version ", utils::packageVersion('rstan'))
    detach("package:rethinking", unload=TRUE) 
  }
  
}


# utils-helper-1-misc-1
# utils-helper-2-misc-2
# utils-helper-3-prepare_data2
# utils-helper-4-prepare_formula
# utils-helper-5-prepare_formula_sigma
# utils-helper-6-prepare_function
# utils-helper-7-prepare_function_nsp
# utils-helper-8-prepare_function_sigma
# utils-helper-9-set_priors_initials
# utils-helper-10-prepare_priors
# utils-helper-11-prepare_initials
# utils-helper-12-R_support_nsp
# utils-helper-13-get_hypothesis_x
# utils-helper-14-get.newdata
# utils-helper-15-get_idata
# utils-helper-16-plot_curves
# utils-helper-17-stancode_custom
# utils-helper-18-standata_custom
# utils-helper-19-prepare_transformations
# utils-helper-20-Pipe
# utils-helper-21-modelbased_growthparameters_nonS3
# utils-helper-22-modelbased_growthparameters_call

