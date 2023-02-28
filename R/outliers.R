


#' Identify (and remove) outliers with abnormal growth velocity
#'
#' @description The function identifies and remove the putative outliers for y
#'   in data, Codes range from 0 (normal) to 8, where 4 and 6 are conventional
#'   outliers
#'
#' @inherit sitar::velout params description
#' @inherit sitar::zapvelout params
#' @param remove A logical (default \code{FALSE}) to indicate whether
#'   identified remove to be removed. When \code{FALSE}, outliers are set as
#'   \code{NA}. If \code{TRUE}, outliers set as \code{NA} are removed.
#' @param verbose A logical (default \code{FALSE}) to show frequency of
#'   observations with different codes (see [sitar::velout()]), and the number
#'   of observations (outliers) removed when \code{remove=TRUE}.
#'
#' @return A data frame with outliers removed when \code{remove=TRUE}
#' @export
#'
#' @examples
#' outliers(x = age, y = height, id = id, data = heights, limit=2, remove=TRUE)
#'
outliers <-
  function (x,
            y,
            id,
            data,
            icode = c(4:6),
            lag = 1,
            velpower = 0.5,
            limit = 5,
            linearise = FALSE,
            remove = FALSE,
            verbose = TRUE) {
    mcall <- match.call()
    data <- eval(mcall$data)
    
    xx_ <- deparse(substitute(x))
    yy_ <- deparse(substitute(y))
    idid_ <- deparse(substitute(id))
    # print(data)
    data <- data %>% dplyr::mutate(order = dplyr::row_number())
    
    data <- data %>% dplyr::arrange(idid_, xx_)
    
    dc <- data[, sapply(mcall[c("x", "y", "id")],
                        deparse)]
    
    colnames_data_ex <- colnames(data)
    colnames_dc_ex <- colnames(dc)
    data_ex <- data %>% dplyr::select(-colnames(dc))
    
    nrow <- nrow(data)
    dc <- na.omit(cbind(dc, count = 1:nrow))
    dc <- dc[order(dc[, 3], dc[, 1]),]
    if (linearise) {
      spline.lm <- loess(dc[, 2] ~ dc[, 1])
      dc[, 2] <- residuals(spline.lm)
    }
    dt1 <- diff(dc[, 1], lag = lag)
    vel1 <- diff(dc[, 2], lag = lag) / dt1 ^ velpower
    dt2 <- diff(dc[, 1], lag = lag * 2)
    vel2 <- diff(dc[, 2], lag = lag * 2) / dt2 ^ velpower
    dt1 <- dt1 == 0
    idlev <- as.numeric(dc[, 3])
    dt1[diff(idlev, lag = lag) != 0] <- FALSE
    vel1[diff(idlev, lag = lag) != 0] <- NA
    vel2[diff(idlev, lag = lag * 2) != 0] <- NA
    vel1 <- trunc(vel1 / mad(vel1, na.rm = TRUE))
    vel2 <- trunc(vel2 / mad(vel2, na.rm = TRUE))
    vel3 <- c(rep(NA, lag), vel2, rep(NA, lag))
    vel2 <- c(vel1, rep(NA, lag))
    vel1 <- c(rep(NA, lag), vel1)
    code <- (as.numeric(abs(vel1) >= limit) + as.numeric(abs(vel2) >=
                                                           limit)) * 2 + as.numeric(abs(vel3) >= limit)
    dt2 <- c(dt1, rep(FALSE, lag))
    dt1 <- c(rep(FALSE, lag), dt1)
    code[dt2 | dt1] <- 8
    code[dt2 & !dt1] <- 7
    t <- is.na(vel3) & !(dt1 | dt2)
    code[t] <- (as.numeric(!is.na(vel1[t]) & abs(vel1[t]) >=
                             limit) + as.numeric(!is.na(vel2[t]) &
                                                   abs(vel2[t]) >=
                                                   limit)) * 6
    dc <- cbind(dc[, c(3, 1, 2, 4)], code, vel1, vel2, vel3)
    mat <- as.data.frame(matrix(
      nrow = nrow,
      ncol = dim(dc)[2],
      dimnames = list(row.names(data), dimnames(dc)[[2]])
    ))
    attr(mat, "data") <- deparse(mcall$data)
    mat[dc$count,] <- dc
    if (is.factor(dc[, 1])) {
      mat[, 1] <- as.factor(mat[, 1])
      levels(mat[, 1]) <- levels(dc[, 1])
    }
    mat$count <- NULL
    mat$code <- factor(mat$code)
    if (verbose) {
      cat("code frequencies\n")
      print(summary(mat$code))
    }
    
    
    mat <- cbind(mat, data_ex)
    mat <- mat %>% dplyr::relocate(all_of(colnames_data_ex))
    
    if (remove) {
      zap <- mat$code %in% icode
      if (verbose) {
        cat(sum(zap), yy_, "values set missing\n")
      }
      mat[mat$code %in% icode, yy_] <- NA
      mat <- mat %>% dplyr::select(all_of(colnames_data_ex))
      mat <- mat %>% tidyr::drop_na()
    }
    mat <- mat %>% dplyr::arrange(order)
    mat <- mat %>% dplyr::select(-order)
    mat
  }

