

#' An internal function to set var_cov initials t0 zero#'
#' @param xscode A character string to evaluate code 
#' method. Options available are \code{'get_idata'} (default) and
#' \code{'idatafunction'}
#' @keywords internal
#' @return A prior object. 
#' @noRd
#'

set_init_gr_effects <- function(xscode,
                                xsdata,
                                full = TRUE,
                                what = 'L',
                                sd_value = 1,
                                z_value = 0,
                                L_value = 0) {
  xscode <-
    get_par_names_from_stancode(xscode, full = full, what = what)
  
  sdi_c <- c()
  parm_c <- c()
  sd_init_list_c_ <- list()
  getindexl <- 0
  for (sdi in xscode) {
    fullxscode_i <- sdi
    getparm <- tail(strsplit(sdi, split = " ")[[1]], 1)
    parm_c <- getparm # c(parm_c, getparm)
    sdi <- gsub("[", "(", sdi, fixed = T)
    sdi <- gsub("]", ")", sdi, fixed = T)
    str_d <-
      regmatches(sdi, gregexpr("(?=\\().*?(?<=\\))", sdi, perl = T))[[1]]
    # this to remove dim after par name
    parm_cm <- parm_c
    if (grepl("[", parm_cm, fixed = T)) {
      parm_cm <- gsub("[", "(", parm_cm, fixed = T)
      parm_cm <- gsub("]", ")", parm_cm, fixed = T)
      parm_cm_ <-
        regmatches(parm_cm,
                   gregexpr("(?=\\().*?(?<=\\))", parm_cm, perl = T))[[1]]
      parm_cm2 <- gsub(parm_cm_, "", parm_cm, fixed = T)
      parm_c <- parm_cm2
    } else if (!grepl("[", parm_cm, fixed = T)) {
      parm_c <- parm_c
    }
    parm_c <- parm_c
    # sdi_c <- c(sdi_c, str_d)
    sdi <- str_d
    sdi <- gsub("(", "", sdi, fixed = T)
    sdi <- gsub(")", "", sdi, fixed = T)
    str_d <- sdi
    #
    sdi <- gsub("(", "", sdi, fixed = T)
    sdi <- gsub(")", "", sdi, fixed = T)
    str_d <- gsub("[[:space:]]", "", sdi)
    # str_d_ <- strsplit(str_d, ",")[[1]]
    str_d_ <- strsplit(str_d, ",") %>% unlist()
    # print(fullxscode_i)
    
    # for student_nu distribution parameter
    if (grepl("sd_nu", parm_c, fixed = T)) {
      set_value <- 3
    }
    
    if (grepl("sd", parm_c, fixed = T)) {
      set_value <- sd_value
    }
    if (grepl("z", parm_c, fixed = T)) {
      set_value <- z_value
    }
    if (grepl("L", parm_c, fixed = T)) {
      set_value <- L_value
    }
    
    # to exclude student_nu distribution parameter
    if (grepl("sd", parm_c, fixed = T)) {
      # to exclude student_nu distribution parameter
      if (!grepl("sd_nu", parm_c, fixed = T)) {
        if (length(str_d_) == 1) {
          dim1 <- xsdata[[str_d_[1]]]
          out <- rep(set_value, dim1)
        } else if (length(str_d_) == 2) {
          dim1 <- xsdata[[str_d_[1]]]
          dim2 <- xsdata[[str_d_[2]]]
          out <- matrix(set_value, dim1, dim2)
        }
        if (is.vector(out)) {
          out <- array(out, dim = length(out))
        }
      } # if(!grepl("sd_nu", parm_c, fixed = T)) {
    } # if(grepl("sd", parm_c, fixed = T)) {
    
    if (grepl("z", parm_c, fixed = T)) {
      if (length(str_d_) == 1) {
        dim1 <- xsdata[[str_d_[1]]]
        out <- rep(set_value, dim1)
      } else if (length(str_d_) == 2) {
        dim1 <- xsdata[[str_d_[1]]]
        dim2 <- xsdata[[str_d_[2]]]
        out <- matrix(set_value, dim1, dim2)
      }
      
      if (is.vector(out)) {
        out <- array(out, dim = length(out))
      }
      
      # for cmdstanr, one dim is assigned array
      # if(grepl("array", fullxscode_i)) {
      #   out <- t(out)
      # } else if(!grepl("array", fullxscode_i)) {
      #   out <- out
      # }
      # but it seems generic problem even for rstan
      if (ncol(out) == 1) {
        out <- t(out)
      }
    } # if(grepl("sd", parm_c, fixed = T)) {
    
    if (grepl("L", parm_c, fixed = T)) {
      if (length(str_d_) == 1) {
        dim1 <- xsdata[[str_d_[1]]]
        out <- matrix(set_value, dim1, dim1)
        diag(out) <- 1
      }
      if (length(str_d_) == 2) {
        dim1 <- xsdata[[str_d_[1]]]
        dim2 <- xsdata[[str_d_[2]]]
        # out <- array(set_value, dim = c(dim2, dim1, dim1))
        outxx <- array(0, dim = c(dim2, dim1, dim1))
        for (i in 1:dim(outxx)[2]) {
          outxx[, i, i] <- 1
        }
        out <- outxx
      }
    } # if(grepl("L", parm_c, fixed = T)) {
    sd_init_list_c_[[parm_c]] <- out
  }
  sd_init_list_c_
}
