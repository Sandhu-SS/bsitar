

#' Prepare Stan function for fitting Bayesian SITAR growth curve model
#'
#' The \code{prepare_function}) prepares custom Stan function  which is passed
#' on to the [bsitar::bsitar] function. For univariate-by-
#' subgroup model (specified by using the \code{univariate_by}) and
#' multivariate model (specified by using the \code{multivariate}),
#' the \code{x}, \code{y}, \code{id}, \code{knots}, \code{nknots}, are
#' automatically matched with the sub-model(s).
#'
#' @param x vector of predictor (typically age in years).
#' @param y vector of outcome (i.e., repeated growth measurements).
#' @param id a factor variable identifying the groups (typically individuals).
#' @param knots vector of values for knots.
#' @param nknots an integer specifying the number of knots.
#' @param data data frame containing variables \code{x}, \code{y} and \code{id}.
#' @param internal_function_args Other internal arguments passed from the
#' [bsitar::bsitar] to the \code{prepare_formula}).
#'
#' @return An character string which later evaluated to a custom function
#'   and inserted into the Stan's functions block.
#'
#' @export
prepare_function <- function(x,
                             y,
                             id,
                             knots,
                             nknots,
                             data,
                             internal_function_args) {
  if (!is.null(internal_function_args)) {
    eout <- list2env(internal_function_args)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
  }
  

  backend <- eval(brms_arguments$backend)
  
  if (backend == 'cmdstanr' & cmdstanr::cmdstan_version() < "2.26.0") {
    stop("Please install CmdStan version 2.26 or newer.")
  }
  
  vector_X_name <- "Xp" # vector Xm to vector Xm
  
  
  if (!is.null(xfunsi[[1]][1]) & xfunsi != "NULL") {
    if (xfunsi == "log") {
      tranform_x_int <- 1
    } else if (xfunsi == "sqrt") {
      tranform_x_int <- 2
    } else if (xfunsi != "log" | xfunsi == "sqrt") {
      tranform_x_int <- 0
    }
  } else if(is.null(xfunsi[[1]][1]) | xfunsi == "NULL") {
    tranform_x_int <- 0
  }
  
  
  
  abcnames <-
    paste0(strsplit(gsub("\\+", " ", fixedsi), " ")[[1]], sep = ",")
  
  snames <- c()
  for (i in 1:(nknots - 1)) {
    if (i < (nknots - 1)) {
      name1 <- paste0("s", i, sep = ",")
    }
    else {
      name1 <- paste0("s", i, sep = "")
    }
    snames[i] <- name1
  }
  
  ept <- function(x)
    eval(parse(text = x), envir = parent.frame())
  
  # For some reasons, 'sitar' (Tim Cole) allows random only 'd' parameter
  # In fact for df > 1, it forces 'd' to be random parameter only
  if (match_sitar_d_form) {
    if (!grepl("d", fixedsi, fixed = T) &
        grepl("d", randomsi, fixed = T)) {
      abcnames <- c(abcnames, "d,")
    }
  }

  fullabcsnames <- c(abcnames, snames)
  fullabcsnames_v <- paste("vector", fullabcsnames, collapse = " ")
  
  if (grepl("b", fixedsi, fixed = T) &
      grepl("c", fixedsi, fixed = T)) {
    defineEx <- paste0("(Xm-b).*exp(c)")
  }
  if (grepl("b", fixedsi, fixed = T) &
      !grepl("c", fixedsi, fixed = T)) {
    defineEx <- paste0("(Xm-b)")
  }
  if (!grepl("b", fixedsi, fixed = T) &
      grepl("c", fixedsi, fixed = T)) {
    defineEx <- paste0("(Xm).*exp(c)")
  }
  if (!grepl("b", fixedsi, fixed = T) &
      !grepl("c", fixedsi, fixed = T)) {
    defineEx <- paste0("(Xm)")
  }
  
 
  
  add_knotinfo <- paste0(
    "\n  int N=num_elements(", vector_X_name, ");",
    paste0("\n  vector[N] Xm=", paste0(getxname, "(", vector_X_name, ")"), ";"), # vector_X_name
    paste0("\n  vector[N] X=", defineEx, ";"),
    paste0("\n  int nknots=", eval(parse(text = nknots)), ";"),
    paste0("\n  vector[nknots] knots=", paste0(getknotsname, "(", '', ")"), ";")
  #   paste0(
  #     "\n  vector[nknots] knots=",
  #     "[",
  #     paste(knots, collapse = ","),
  #     "]'"
  #   ),
  #   ";"
   )
  
  vectorA <- "\n  vector[N] A=a-(s1*min(knots));"
  
  add_knotinfo <- paste0(add_knotinfo, vectorA)
  
  
  if((backend == "rstan" & utils::packageVersion("rstan") >= "2.26.1") | 
     backend == "cmdstanr") {
    fun_body <- "
    matrix[N, nknots-1] Spl;
    matrix[nknots-1, N] rcs;
    matrix[N, nknots] Xx;
    int km1 = nknots - 1;
    int jp1;
    int j=1;
    for(ia in 1:N) {
     for(ja in 1:nknots) {
      Xx[ia,ja] = (X[ia] - knots[ja] > 0 ? X[ia] - knots[ja] : 0);
        }
    }
     Spl[,1]=X;
     while (j <= nknots - 2) {
      jp1 = j + 1;
       Spl[,jp1] = (Xx[,j]^3-(Xx[,km1]^3)*(knots[nknots]-knots[j])/
       (knots[nknots]-knots[km1]) + (Xx[,nknots]^3)*(knots[km1]-knots[j])/
       (knots[nknots]-knots[km1])) / (knots[nknots]-knots[1])^2;
        j = j + 1;
      }"
  } 
  
  if((backend == "rstan" & utils::packageVersion("rstan") < "2.26.1") & 
     backend != "cmdstanr") {
    fun_body <- "
    matrix[N, nknots-1] Spl;
    matrix[nknots-1, N] rcs;
    matrix[N, nknots] Xx;
    int km1 = nknots - 1;
    int jp1;
    int j=1;
    for(ia in 1:N) {
     for(ja in 1:nknots) {
      Xx[ia,ja] = (X[ia] - knots[ja] > 0 ? X[ia] - knots[ja] : 0);
        }
    }
     Spl[,1]=X;
     while (j <= nknots - 2) {
     for(i in 1:N) {
      jp1 = j + 1;
       Spl[i,jp1] = (pow(Xx[i,j],3)-(pow(Xx[i,km1],3))*(knots[nknots]-knots[j])/
       (knots[nknots]-knots[km1]) +
       (pow(Xx[i,nknots],3))*(knots[km1]-knots[j])/(knots[nknots]-knots[km1]))/ 
       pow((knots[nknots]-knots[1]),2);
      }
      j = j + 1;
    }
    "
  }
  
  
  name4 <- c()
  for (i in 1:(nknots - 1)) {
    name1 <- paste0("", "s", i, sep = "")
    if (i < (nknots - 1)) {
      # name2 <- paste0(' .* to_vector(Spl[,',i,"]') +")
      name2 <- paste0(' .* Spl[,', i, "] +")
    }
    else {
      # name2 <- paste0(' .* to_vector(Spl[,',i,"]') ;\n")
      name2 <- paste0(' .* Spl[,', i, "]")
    }
    name3 <- paste0(name1, name2, sep = "")
    name4[i] <- name3
  }
  name50 <- paste("", name4, collapse = " ")
  
  nameadja <- "A"
  
  ###########
  # For some reasons, 'sitar' (Tim Cole) allows random only 'd' parameter
  # In fact for df > 1, it forces 'd' to be random parameter only
  if (match_sitar_d_form) {
    if (grepl("d", randomsi, fixed = T)) {
      nameadja <- "A+(d . * Spl[,1])"
    }
  }
  
  if (!match_sitar_d_form) {
    if (grepl("d", fixedsi, fixed = T)) {
      nameadja <- "A+(d . * Spl[,1])"
    }
  }

  name5 <- paste(" (", name50, ");\n")
  
  if (grepl("c", fixedsi, fixed = T)) {
    name51 <- paste(" (", name50, ") .* exp(c) ;\n")
    name52 <- paste(" (", name50, ") .* exp(c)^2 ;\n")
  }
  if (!grepl("c", fixedsi, fixed = T)) {
    name51 <- paste(" (", name50, ") ;\n")
    name52 <- paste(" (", name50, ") ;\n")
  }
  
  returnmu <-
    paste0("return(",   paste0(nameadja, "+", gsub(";", "", name5))     , ");")
  
  # need spaces otherwise rstan 2.21 throws error: variable s1. not found
  returnmu <- gsub("\\s", "", returnmu)
  returnmu <- gsub("\\." , " \\." , returnmu, fixed = FALSE)
  returnmu <- gsub("\\*" , "\\* " , returnmu, fixed = FALSE)
  # don't create space for +
  # returnmu <- gsub("+" , " + " , returnmu, fixed = TRUE)
  
  endof_fun <-
    paste0("\n    ", returnmu, "\n  } // end of spline function", sep = " ")
  
  
  
  
  start_fun <-
    paste0("\nvector ",
           spfncname,
           "(vector ", vector_X_name, ", ",
           fullabcsnames_v,
           ") {" ,
           collapse = " ")
  rcsfun <- paste(start_fun, add_knotinfo, fun_body, endof_fun)
  
  # print(cat(rcsfun))
  # stop()
  
  ######
  
add_context_getknots_fun <- "
/* Knots
 * xoffset and Knots already transformed:
 * Returns:
 *   Knots
 */
"
  
  getknots_fun <- 
    paste0(add_context_getknots_fun,
           "\n  ",
           "vector ",
           getknotsname,
           "() {" ,
           paste0("\n  int nknots=", eval(parse(text = nknots)), ";"),
           paste0(
             "\n  vector[nknots] knots=",
             "[",
             paste(knots, collapse = ","),
             "]';"
           ),
           "\n  ",
           "return(knots);",
           "\n}  "
           , paste0("// end of ", getknotsname), 
           collapse = " ")  # %>% cat()
  
add_context_getx_fun <- "
/* Transform x variable
 * Args:
 *   Xp: x variable
 * Transformation code (tranform_x, 0 to 2) 
 * 0, no transformation, 1 log, 2 square rooot
 * Note that the xoffset  is already transformed
 * Returns:
 *   x variable with log/sqrt transformation
 */
"
  
  getx_fun <- 
    paste0(add_context_getx_fun,
           "\n  ",
           "vector ",
           getxname,
           paste0(" (vector ", vector_X_name, ") {") ,
           "\n  ",
           paste0("int N=num_elements(", vector_X_name, ");"),
           "\n  ",
           paste0("real xoffset = ", xoffset, ";"),
           paste0("\n  int tranform_x = ", eval(parse(text = tranform_x_int)), ";"),
           paste0("\n  vector[N] x;"),
           "\n",
           paste0("  if(tranform_x == 0 ) {",
                  "\n   ",
                  "x = ", vector_X_name, " - xoffset",  ";", 
                  "\n  }",
                  "\n  ",
                  "if(tranform_x == 1 ) {",
                  "\n   ",
                  "x = log(", vector_X_name, ") - xoffset;", 
                  "\n  }",
                  "\n  ",
                  "if(tranform_x == 2 ) {",
                  "\n    ",
                  "x = sqrt(", vector_X_name, ") - xoffset;", 
                  "\n  }"
           ),
           "\n  ",
           "return(x);",
           "\n}  "
           , paste0("// end of ", getxname), 
           collapse = " ")  # %>% cat()
  
  
  getx_knots_fun <- paste0(getx_fun, 
                           "\n", 
                           getknots_fun)
  
  
  
  #####
  # rcsfun <- paste(getx_knots_fun, rcsfun)
  
  
  ##########
  set_x_y_scale_factror <- function(xfunsi = NULL, yfunsi = NULL, 
                                    tranformations = "identity") {
    scale_set_comb <- tranformations
    scale_set_comb1 <- paste(scale_set_comb, scale_set_comb, sep = "_")
    scale_set_comb2 <- with(subset(expand.grid(scale_set_comb,scale_set_comb),
                                   Var1!=Var2),paste0(Var1,'_',Var2))
    scale_set_comb <- c(scale_set_comb1, scale_set_comb2)
    
    if (!is.null(xfunsi[[1]][1]) & xfunsi != "NULL") {
      if (xfunsi == "log") {
        xscale_set <- "log"
      } else if (xfunsi == "sqrt") {
        xscale_set <- "sqrt"
      } else if (xfunsi != "log" | xfunsi == "sqrt") {
        xscale_set <- "identity"
      }
    } else if(is.null(xfunsi[[1]][1]) | xfunsi == "NULL") {
      xscale_set <- "identity"
    }
    
    if (!is.null(yfunsi[[1]][1]) & yfunsi != "NULL") {
      if (yfunsi == "log") {
        yscale_set <- "log"
      } else if (yfunsi == "sqrt") {
        yscale_set <- "sqrt"
      } else if (yfunsi != "log" | yfunsi == "sqrt") {
        yscale_set <- "identity"
      }
    } else if(is.null(yfunsi[[1]][1]) | yfunsi == "NULL") {
      yscale_set <- "identity"
    }
    
    
    
    
    if (xscale_set == "identity" & yscale_set == "identity") {
      xscale_factor_str_d1 <- "rep_vector(1, N);"
      xscale_factor_str_d2 <- "rep_vector(1, N);"
      yscale_factor_str_d1 <- "rep_vector(1, N);"
      yscale_factor_str_d2 <- "rep_vector(1, N);"
    } else if (xscale_set == "log" & yscale_set == "log") {
      xscale_factor_str_d1 <- "exp(Xm + xoffset);"
      xscale_factor_str_d2 <- "exp(Xm + xoffset);"
      yscale_factor_str_d1 <- "(pred_d0);"
      yscale_factor_str_d2 <- "(pred_d0);"
    } else if (xscale_set == "sqrt" & yscale_set == "sqrt") {
      xscale_factor_str_d1 <- "(Xm + xoffset);"
      xscale_factor_str_d2 <- "(Xm + xoffset);"
      yscale_factor_str_d1 <- "(sqrt(pred_d0));"
      yscale_factor_str_d2 <- "(sqrt(pred_d0));"
    } else if(xscale_set == "log" & yscale_set == "identity") {
      xscale_factor_str_d1 <- "exp(Xm + xoffset);"
      xscale_factor_str_d2 <- "exp(Xm + xoffset);"
      yscale_factor_str_d1 <- "rep_vector(1, N);"
      yscale_factor_str_d2 <- "rep_vector(1, N);"
    } else if(xscale_set == "sqrt" & yscale_set == "identity") {
      xscale_factor_str_d1 <- "(Xm + xoffset);"
      xscale_factor_str_d2 <- "(Xm + xoffset);"
      yscale_factor_str_d1 <- "(0.5);"
      yscale_factor_str_d2 <- "(0.5);"
    } else if(xscale_set == "identity" & yscale_set == "log") {
      xscale_factor_str_d1 <- "rep_vector(1, N);"
      xscale_factor_str_d2 <- "rep_vector(1, N);"
      yscale_factor_str_d1 <- "(pred_d0);"
      yscale_factor_str_d2 <- "(pred_d0);"
    } else if(xscale_set == "sqrt" & yscale_set == "log") {
      xscale_factor_str_d1 <- "(Xm + xoffset);"
      xscale_factor_str_d2 <- "(Xm + xoffset);"
      yscale_factor_str_d1 <- "(0.5 * (pred_d0));"
      yscale_factor_str_d2 <- "(0.5 * (pred_d0));"
    } else if(xscale_set == "identity" & yscale_set == "sqrt") {
      xscale_factor_str_d1 <- "rep_vector(1, N);"
      xscale_factor_str_d2 <- "rep_vector(1, N);"
      yscale_factor_str_d1 <- "(2.0 * sqrt(pred_d0));"
      yscale_factor_str_d2 <- "(2.0 * sqrt(pred_d0));"
    } else if(xscale_set == "log" & yscale_set == "sqrt") {
      xscale_factor_str_d1 <- "exp(Xm + xoffset);"
      xscale_factor_str_d2 <- "exp(Xm + xoffset);"
      yscale_factor_str_d1 <- "(2.0 * sqrt(pred_d0));"
      yscale_factor_str_d2 <- "(2.0 * sqrt(pred_d0));"
    } 
    
    list(xscale_factor_str_d1 = xscale_factor_str_d1,
         xscale_factor_str_d2 = xscale_factor_str_d2,
         yscale_factor_str_d1 = yscale_factor_str_d1,
         yscale_factor_str_d2 = yscale_factor_str_d2)
    
  } # end set_x_y_scale_factror
  
  ##########
  
  
  
  create_internal_function <-
    function(y,
             function_str,
             fname,
             fnameout,
             spl,
             splout,
             xfunsi,
             yfunsi,
             xoffset,
             body,
             fixedsi) {
      split1 <- strsplit(function_str, gsub("\\[", "\\\\[", spl))[[1]][-1]
      split2 <- strsplit(split1, "return")[[1]][-2]
      out <- gsub(split2, body, function_str, fixed = T)
      out <- gsub(spl, splout, out, fixed = T)
      out <- gsub(fname, fnameout, out, fixed = T)
      
      if(grepl("d0", fnameout)) {
        out <- out
      } else if(grepl("d1", fnameout) | grepl("d2", fnameout)) {
        out <- gsub("return(A+", "return(0+", out, fixed = T)
      }
      
      
      
      if (grepl("c", fixedsi, fixed = T)) {
        if (grepl("d0", fnameout)) {
          out <- gsub("]));", "]));", out, fixed = T)
          out <-
            gsub(
              "end of spline function",
              paste0("end of spline function", "_", y, "d", 0),
              out,
              fixed = T
            )
        } else if (grepl("d1", fnameout)) {
          out <- gsub("]));", "]).*exp(c));", out, fixed = T)
          out <-
            gsub(
              "end of spline function",
              paste0("end of spline function", "_", y, "d", 1),
              out,
              fixed = T
            )
        } else if (grepl("d2", fnameout)) {
          out <- gsub("]));", "]).*exp(c)^2);", out, fixed = T)
          out <-
            gsub(
              "end of spline function",
              paste0("end of spline function", "_", y, "d", 2),
              out,
              fixed = T
            )
        } else {
          out <- out
        }
      }
      
      ####
      if(grepl("d0", fnameout)) {
        pattern <- "return\\(\\s*(.*?)\\s*\\);"
        result <- regmatches(out, regexec(pattern, out))
        out_unscaled <- paste0("vector[N] out_unscaled=", result[[1]][2], ";")
        
        if (yfunsi == "log") {
          out_scaled <- paste0("    vector[N] out_scaled=", "exp", "(", "out_unscaled", ")", ";")
        } else if (yfunsi == "sqrt") {
          if((backend == "rstan" & utils::packageVersion("rstan") >= "2.26.1") |
             backend == "cmdstanr") {
            out_scaled <- paste0("    vector[N] out_scaled=", "", "(", "out_unscaled", ")^2", ";")
          }
          if((backend == "rstan" & utils::packageVersion("rstan") < "2.26.1") & 
             backend != "cmdstanr") {
            out_scaled <- paste0("    vector[N] out_scaled=", "", "(", "pow(", "out_unscaled", ", 2)" , ")", ";")
          }
        } else if (yfunsi != "log" & yfunsi != "sqrt") {
          out_scaled <- paste0("    vector[N] out_scaled=", "", "(", "out_unscaled", ")", ";")
        }
        
        out <- gsub(result[[1]][2], "out_scaled", out, fixed = T)
        out_return <- paste0(out_unscaled, "\n", out_scaled)
        setxoffset <- paste0("real xoffset = ", xoffset, ";")
        out_return <- paste0(setxoffset,
                             "\n    ",
                             out_return)
        out_return_p <- paste0(out_return, "\n", "    return")
        out <- gsub("return", out_return_p, out, fixed = T)
        
      } else if (grepl("d1", fnameout) | grepl("d2", fnameout)) {
      pattern <- "return\\(\\s*(.*?)\\s*\\);"
      result <- regmatches(out, regexec(pattern, out))
      
      set_x_y_scale <- set_x_y_scale_factror(xfunsi = xfunsi, yfunsi = yfunsi, 
                                             tranformations = c("identity", "log", "sqrt"))
      
      if (grepl("d1", fnameout)) {
        xscale_factor <- set_x_y_scale[['xscale_factor_str_d1']]
        yscale_factor <- set_x_y_scale[['yscale_factor_str_d1']]
      } else if (grepl("d2", fnameout)) {
        xscale_factor <- set_x_y_scale[['xscale_factor_str_d2']]
        yscale_factor <- set_x_y_scale[['yscale_factor_str_d2']]
      }
      
      xscale_factor <- gsub(";", "", xscale_factor)
      yscale_factor <- gsub(";", "", yscale_factor)
      out_unscaled <- paste0("vector[N] out_unscaled=", result[[1]][2], ";")
      out_scaled <- paste0("    vector[N] out_scaled=", "(", "(", yscale_factor, ")", " .* ", "(", 'out_unscaled', ")", ")", " ./ ", "(", xscale_factor, ")", ";")
      out <- gsub(result[[1]][2], "out_scaled", out, fixed = T)
      out_return <- paste0(out_unscaled, "\n", out_scaled)
      addpdo <- paste0("vector[N] pred_d0=", spl_fun_ford, ";")
      setxoffset <- paste0("real xoffset = ", xoffset, ";")
      out_return <- paste0(addpdo,
                           "\n    ",
                           setxoffset,
                           "\n    ",
                           out_return)
      out_return_p <- paste0(out_return, "\n", "    return")
      out <- gsub("return", out_return_p, out, fixed = T)
      }
      
      
      # print(out_return)
      # print(cat(out))
      # stop()
      ####
      return(out)
    }
  
  
  
  
  # Create d0 function 
  
  fnameout <- paste0(spfncname, "_", "d0")
  spl <- "Spl[,1]=X;"
  splout <- spl
  spl_fun_ford <- paste0(fnameout, "(vector ", vector_X_name, ", ", fullabcsnames_v, ")")
  spl_fun_ford <- gsub("vector", "", spl_fun_ford, fixed = T)
  
  if((backend == "rstan" & utils::packageVersion("rstan") >= "2.26.1") | 
     backend == "cmdstanr") {
    body <- "
     while (j <= nknots - 2) {
      jp1 = j + 1;
      Spl[,jp1] =
        (1*Xx[,j]^3) * (1/((knots[nknots]-knots[1])^2))  -
        (1*Xx[,km1]^3) * (knots[nknots]-knots[j]) / 
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) +
        (1*Xx[,nknots]^3) * (knots[km1]-knots[j]) / 
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) ;
      j = j + 1;
    }
    "
  }
  
  if((backend == "rstan" & utils::packageVersion("rstan") < "2.26.1") & 
     backend != "cmdstanr") {
    body <- "
     while (j <= nknots - 2) {
      for(i in 1:N) {
          jp1 = j + 1;
          Spl[i,jp1] = (1*pow(Xx[i,j],3) -
          (1*pow(Xx[i,km1],3))*(knots[nknots]-knots[j]) /
          (knots[nknots]-knots[km1]) +
          (1*pow(Xx[i,nknots],3))*(knots[km1]-knots[j]) /
          (knots[nknots]-knots[km1])) /
          pow((knots[nknots]-knots[1]),2);
      }
      j = j + 1;
    }
    "
  }
  
  
  spl_d0 <- create_internal_function(
    y = y,
    function_str = rcsfun,
    fname = spfncname,
    fnameout = fnameout,
    spl = spl,
    splout = splout,
    xfunsi = xfunsi,
    yfunsi = yfunsi,
    xoffset = xoffset,
    body = body,
    fixedsi = fixedsi
  )
  
  
  
  
  
  
  # Create d1 function 
  
  fnameout <- paste0(spfncname, "_", "d1")
  spl <- "Spl[,1]=X;"
  splout <- "Spl[,1]=rep_vector(1, N);"
  
  if((backend == "rstan" & utils::packageVersion("rstan") >= "2.26.1") | 
     backend == "cmdstanr") {
    body <- "
     while (j <= nknots - 2) {
      jp1 = j + 1;
      Spl[,jp1] =
        (3*Xx[,j]^2) * (1/((knots[nknots]-knots[1])^2))  -
        (3*Xx[,km1]^2) * (knots[nknots]-knots[j]) / 
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) +
        (3*Xx[,nknots]^2) * (knots[km1]-knots[j]) / 
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) ;
      j = j + 1;
    }
    "
  }
  
  if((backend == "rstan" & utils::packageVersion("rstan") < "2.26.1") & 
     backend != "cmdstanr") {
    body <- "
     while (j <= nknots - 2) {
      for(i in 1:N) {
          jp1 = j + 1;
          Spl[i,jp1] = (3*pow(Xx[i,j],2) -
          (3*pow(Xx[i,km1],2))*(knots[nknots]-knots[j]) /
          (knots[nknots]-knots[km1]) +
          (3*pow(Xx[i,nknots],2))*(knots[km1]-knots[j]) /
          (knots[nknots]-knots[km1])) /
          pow((knots[nknots]-knots[1]),2);
      }
      j = j + 1;
    }
    "
  }
  
  
  spl_d1 <- create_internal_function(
    y = y,
    function_str = rcsfun,
    fname = spfncname,
    fnameout = fnameout,
    spl = spl,
    splout = splout,
    xfunsi = xfunsi,
    yfunsi = yfunsi,
    xoffset = xoffset,
    body = body,
    fixedsi = fixedsi
  )
  
  
  # Create d2 function
  
  fnameout <- paste0(spfncname, "_", "d2")
  spl <- "Spl[,1]=X;"
  splout <- "Spl[,1]=rep_vector(0, N);"

  if((backend == "rstan" & utils::packageVersion("rstan") >= "2.26.1") | 
     backend == "cmdstanr") {
    body <- "
     while (j <= nknots - 2) {
      jp1 = j + 1;
      Spl[,jp1] =
        (6*Xx[,j]^1) * (1/((knots[nknots]-knots[1])^2))  -
        (6*Xx[,km1]^1) * (knots[nknots]-knots[j]) / 
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) +
        (6*Xx[,nknots]^1) * (knots[km1]-knots[j]) / 
        ((knots[nknots]-knots[km1]) * (knots[nknots]-knots[1])^2) ;
      j = j + 1;
   }
    "
  }
  
  if((backend == "rstan" & utils::packageVersion("rstan") < "2.26.1") & 
     backend != "cmdstanr") {
    body <- "
     while (j <= nknots - 2) {
      for(i in 1:N) {
          jp1 = j + 1;
          Spl[i,jp1] = (6*pow(Xx[i,j],1) -
          (6*pow(Xx[i,km1],1))*(knots[nknots]-knots[j]) /
          (knots[nknots]-knots[km1]) +
          (6*pow(Xx[i,nknots],1))*(knots[km1]-knots[j]) /
          (knots[nknots]-knots[km1])) /
          pow((knots[nknots]-knots[1]),2);
      }
      j = j + 1;
    }
    "
  }
  
  
  
  spl_d2 <- create_internal_function(
    y = y,
    function_str = rcsfun,
    fname = spfncname,
    fnameout = fnameout,
    spl = spl,
    splout = splout,
    xfunsi = xfunsi,
    yfunsi = yfunsi,
    xoffset = xoffset,
    body = body,
    fixedsi = fixedsi
  )
  
  rcsfun <- paste0(rcsfun, spl_d0, spl_d1, spl_d2, sep = "\n")
  
  rcsfun <- paste(getx_knots_fun, rcsfun)
  
  
  # rcsfun <<- rcsfun
  # print(cat(rcsfun))
  # stop()
  
  return(rcsfun)
} 