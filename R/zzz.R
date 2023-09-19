
.onAttach <- function(libname, pkgname) {
  
  invisible(suppressPackageStartupMessages(
    sapply(c("brms"),
           requireNamespace, quietly = TRUE)
  ))
  
  packageStartupMessage("The minimum rstan version required is 2.26.2",
                        "\n ", 
                        "The latest rstan version can be installed from",
                        "\n ", 
                        "https://mc-stan.org/r-packages/")
  
  
  if("rethinking" %in% (.packages())){
    message("Package 'rethinking' detached and unloaded as it creates conflict",
            " \nwith the current rstan version ", utils::packageVersion('rstan'))
    detach("package:rethinking", unload=TRUE) 
  }
  
}

