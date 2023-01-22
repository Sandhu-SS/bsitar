
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("The minimum rstan version required is 2.26.2",
                        "\n ", 
                        "The latest rstan version can be installed from",
                        "\n ", 
                        "https://mc-stan.org/r-packages/")
}

