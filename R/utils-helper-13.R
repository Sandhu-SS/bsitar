

#' An internal function to compute hypothesis from estimates
#'
#' @details Function adpated from \code{marginaleffects:::get_hypothesis}
#'   available at
#'   <https://github.com/vincentarelbundock/marginaleffects/blob/main/R/get_hypothesis.R>
#' 
#' @param x A data frame
#' @param hypothesis A character vector 
#' @param by A character or a character vector 
#' @param newdata A character vector 
#' @param draws A column name (a character)
#' 
#' @return A data frame
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
get_hypothesis_x <- function(
    x,
    hypothesis,
    by = NULL,
    newdata = NULL,
    draws = NULL) {
  
  if (is.null(hypothesis)) return(x)
  
  try(zz <- insight::check_if_installed(c("data.table"), 
                                        minimum_version = 
                                          get_package_minversion(
                                            'data.table'
                                          ), 
                                        prompt = FALSE,
                                        stop = FALSE))
  
  try(zz <- insight::check_if_installed(c("checkmate"), 
                                        minimum_version = 
                                          get_package_minversion(
                                            'checkmate'
                                          ), 
                                        prompt = FALSE,
                                        stop = FALSE))
  
  if (is.function(hypothesis)) {
    argnames <- names(formals(hypothesis))
    if (!"x" %in% argnames) stop("The `hypothesis` function must accept an `x` argument.", call. = FALSE)
    if (any(!argnames %in% c("x", "newdata", "by", "draws"))) {
      msg <- "The allowable arguments for the `hypothesis` function are: x, newdata`, by, and draws."
      stop(msg, call. = FALSE)
    }
    args <- list(x = x, newdata = newdata, by = by, draws = draws)
    args <- args[names(args) %in% argnames]
    out <- CustomDoCall(hypothesis, args)
    if (!inherits(out, "data.frame") || !"term" %in% colnames(out) || !"estimate" %in% colnames(out)) {
      msg <- "The `hypothesis` function must return a data frame with `term` and `estimate` columns."
      stop(msg, call. = FALSE)
    }
    return(out)
  }
  
  lincom <- NULL
  
  # lincom: numeric vector or matrix
  if (isTRUE(checkmate::check_numeric(hypothesis))) {
    if (isTRUE(checkmate::check_atomic_vector(hypothesis))) {
      checkmate::assert_numeric(hypothesis, len = nrow(x))
      lincom <- as.matrix(hypothesis)
    } else if (isTRUE(checkmate::check_matrix(hypothesis))) {
      lincom <- hypothesis
    }
  }
  
  # lincom: string shortcuts
  valid <- c("pairwise", "reference", "sequential", "revpairwise", "revreference", "revsequential")
  if (isTRUE(hypothesis %in% valid)) {
    if (nrow(x) > 25) {
      msg <- 'The "pairwise", "reference", and "sequential" options of the `hypotheses` argument are not supported for `marginaleffects` commands which generate more than 25 rows of results. Use the `newdata`, `by`, and/or `variables` arguments to compute a smaller set of results on which to conduct hypothesis tests.'
      insight::format_error(msg)
    }
  }
  if (isTRUE(hypothesis == "reference")) {
    lincom <- lincom_reference(x, by)
  } else if (isTRUE(hypothesis == "revreference")) {
    lincom <- lincom_revreference(x, by)
  } else if (isTRUE(hypothesis == "sequential")) {
    lincom <- lincom_sequential(x, by)
  } else if (isTRUE(hypothesis == "revsequential")) {
    lincom <- lincom_revsequential(x, by)
  } else if (isTRUE(hypothesis == "pairwise")) {
    lincom <- lincom_pairwise(x, by)
  } else if (isTRUE(hypothesis == "revpairwise")) {
    lincom <- lincom_revpairwise(x, by)
  }
  lincom <- sanitize_lincom(lincom, x)
  
  # matrix hypothesis
  if (isTRUE(checkmate::check_matrix(lincom))) {
    out <- lincom_multiply(x, lincom)
    return(out)
    
    # string hypothesis
  } else if (is.character(hypothesis)) {
    out_list <- draws_list <- list()
    lab <- attr(hypothesis, "label")
    tmp <- expand_wildcard(hypothesis, nrow(x), lab)
    hypothesis <- tmp[[1]]
    labs <- tmp[[2]]
    for (i in seq_along(hypothesis)) {
      out_list[[i]] <- eval_string_hypothesis(x, hypothesis[i], labs[i])
      draws_list[[i]] <- attr(out_list[[i]], "posterior_draws")
    }
    out <- CustomDoCall(rbind, out_list)
    attr(out, "posterior_draws") <- CustomDoCall(rbind, draws_list)
    attr(out, "label") <- if (!is.null(attr(labs, "names"))) {
      attr(labs, "names")
    } else {
      labs
    }
    return(out)
  }
  
  insight::format_error("`hypotheses` is broken. Please report this bug: https://github.com/vincentarelbundock/marginaleffects/issues.")
}


get_hypothesis_row_labels <- function(x, by = NULL) {
  lab <- grep("^term$|^by$|^group$|^value$|^contrast$|^contrast_", colnames(x), value = TRUE)
  lab <- Filter(function(z) length(unique(x[[z]])) > 1, lab)
  if (isTRUE(checkmate::check_character(by))) {
    lab <- unique(c(lab, by))
  }
  if (length(lab) == 0) {
    lab <- NULL
  } else {
    lab_df <- data.frame(x)[, lab, drop = FALSE]
    idx <- vapply(lab_df, FUN = function(x) length(unique(x)) > 1, FUN.VALUE = logical(1))
    if (sum(idx) > 0) {
      lab <- apply(lab_df[, idx, drop = FALSE], 1, paste, collapse = ", ")
    } else {
      lab <- apply(lab_df, 1, paste, collapse = ", ")
    }
  }
  
  # wrap in parentheses to avoid a-b-c-d != (a-b)-(c-d)
  if (any(grepl("-", lab))) {
    lab <- sprintf("(%s)", lab)
  }
  
  return(lab)
}


sanitize_lincom <- function(lincom, x) {
  if (isTRUE(checkmate::check_matrix(lincom))) {
    checkmate::assert_matrix(lincom, nrows = nrow(x))
    if (is.null(colnames(lincom))) {
      colnames(lincom) <- rep("custom", ncol(lincom))
    }
  }
  return(lincom)
}


lincom_revreference <- function(x, by) {
  lincom <- -1 * diag(nrow(x))
  lincom[1, ] <- 1
  lab <- get_hypothesis_row_labels(x, by = by)
  if (length(lab) == 0 || anyDuplicated(lab) > 0) {
    lab <- sprintf("Row 1 - Row %s", seq_len(ncol(lincom)))
  } else {
    lab <- sprintf("%s - %s", lab[1], lab)
  }
  colnames(lincom) <- lab
  lincom <- lincom[, 2:ncol(lincom), drop = FALSE]
  return(lincom)
}


lincom_reference <- function(x, by) {
  lincom <- diag(nrow(x))
  lincom[1, ] <- -1
  lab <- get_hypothesis_row_labels(x, by = by)
  if (length(lab) == 0 || anyDuplicated(lab) > 0) {
    lab <- sprintf("Row %s - Row 1", seq_len(ncol(lincom)))
  } else {
    lab <- sprintf("%s - %s", lab, lab[1])
  }
  colnames(lincom) <- lab
  lincom <- lincom[, 2:ncol(lincom), drop = FALSE]
  return(lincom)
}


lincom_revsequential <- function(x, by) {
  lincom <- matrix(0, nrow = nrow(x), ncol = nrow(x) - 1)
  lab <- get_hypothesis_row_labels(x, by = by)
  if (length(lab) == 0 || anyDuplicated(lab) > 0) {
    lab <- sprintf("Row %s - Row %s", seq_len(ncol(lincom)), seq_len(ncol(lincom)) + 1)
  } else {
    lab <- sprintf("%s - %s", lab[seq_len(ncol(lincom))], lab[seq_len(ncol(lincom)) + 1])
  }
  for (i in seq_len(ncol(lincom))) {
    lincom[i:(i + 1), i] <- c(1, -1)
  }
  colnames(lincom) <- lab
  return(lincom)
}


lincom_sequential <- function(x, by) {
  lincom <- matrix(0, nrow = nrow(x), ncol = nrow(x) - 1)
  lab <- get_hypothesis_row_labels(x, by = by)
  if (length(lab) == 0 || anyDuplicated(lab) > 0) {
    lab <- sprintf("Row %s - Row %s", seq_len(ncol(lincom)) + 1, seq_len(ncol(lincom)))
  } else {
    lab <- sprintf("%s - %s", lab[seq_len(ncol(lincom)) + 1], lab[seq_len(ncol(lincom))])
  }
  for (i in seq_len(ncol(lincom))) {
    lincom[i:(i + 1), i] <- c(-1, 1)
  }
  colnames(lincom) <- lab
  return(lincom)
}


lincom_revpairwise <- function(x, by) {
  lab_row <- get_hypothesis_row_labels(x, by = by)
  lab_col <- NULL
  flag <- length(lab_row) == 0 || anyDuplicated(lab_row) > 0
  mat <- list()
  for (i in seq_len(nrow(x))) {
    for (j in 2:nrow(x)) {
      if (i < j) {
        tmp <- matrix(0, nrow = nrow(x), ncol = 1)
        tmp[i, ] <- -1
        tmp[j, ] <- 1
        mat <- c(mat, list(tmp))
        if (isTRUE(flag)) {
          lab_col <- c(lab_col, sprintf("Row %s - Row %s", j, i))
        } else {
          lab_col <- c(lab_col, sprintf("%s - %s", lab_row[j], lab_row[i]))
        }
      }
    }
  }
  lincom <- CustomDoCall("cbind", mat)
  colnames(lincom) <- lab_col
  return(lincom)
}


lincom_pairwise <- function(x, by) {
  lab_row <- get_hypothesis_row_labels(x, by = by)
  lab_col <- NULL
  flag <- length(lab_row) == 0 || anyDuplicated(lab_row) > 0
  mat <- list()
  for (i in seq_len(nrow(x))) {
    for (j in 2:nrow(x)) {
      if (i < j) {
        tmp <- matrix(0, nrow = nrow(x), ncol = 1)
        tmp[j, ] <- -1
        tmp[i, ] <- 1
        mat <- c(mat, list(tmp))
        if (isTRUE(flag)) {
          lab_col <- c(lab_col, sprintf("Row %s - Row %s", i, j))
        } else {
          lab_col <- c(lab_col, sprintf("%s - %s", lab_row[i], lab_row[j]))
        }
      }
    }
  }
  lincom <- CustomDoCall("cbind", mat)
  colnames(lincom) <- lab_col
  return(lincom)
}


lincom_multiply <- function(x, lincom) {
  # bayesian
  draws <- attr(x, "posterior_draws")
  if (!is.null(draws)) {
    draws <- t(as.matrix(lincom)) %*% draws
    out <- data.table::data.table(
      term = colnames(lincom),
      tmp = apply(draws, 1, stats::median))
    data.table::setnames(out, old = "tmp", new = "estimate")
    attr(out, "posterior_draws") <- draws
    
    # frequentist
  } else {
    out <- data.table::data.table(
      term = colnames(lincom),
      tmp = as.vector(x[["estimate"]] %*% lincom))
    data.table::setnames(out, old = "tmp", new = "estimate")
  }
  
  out <- out[out$term != "1 - 1", , drop = FALSE]
  return(out)
}


eval_string_hypothesis <- function(x, hypothesis, lab) {
  # row indices: `hypotheses` includes them, but `term` does not
  if (isTRUE(grepl("\\bb\\d+\\b", hypothesis)) && !any(grepl("\\bb\\d+\\b", x[["term"]]))) {
    bmax <- regmatches(lab, gregexpr("\\bb\\d+\\b", lab))[[1]]
    bmax <- tryCatch(max(as.numeric(gsub("b", "", bmax))), error = function(e) 0)
    if (bmax > nrow(x)) {
      msg <- "%s cannot be used in `hypothesis` because the call produced just %s estimate(s). Try executing the exact same command without the `hypothesis` argument to see which estimates are available for hypothesis testing."
      msg <- sprintf(msg, paste0("b", bmax), nrow(x))
      insight::format_error(msg)
    }
    for (i in seq_len(nrow(x))) {
      tmp <- paste0("marginaleffects__", i)
      hypothesis <- gsub(paste0("b", i), tmp, hypothesis)
    }
    rowlabels <- paste0("marginaleffects__", seq_len(nrow(x)))
    
    # term names
  } else {
    if (!"term" %in% colnames(x) || anyDuplicated(x$term) > 0) {
      msg <- c(
        'To use term names in a `hypothesis` string, the same function call without `hypothesis` argument must produce a `term` column with unique row identifiers. You can use `b1`, `b2`, etc. indices instead of term names in the `hypotheses` string Ex: "b1 + b2 = 0" Alternatively, you can use the `newdata`, `variables`, or `by` arguments:',
        "",
        "mod <- lm(mpg ~ am * vs + cyl, data = mtcars)",
        'comparisons(mod, newdata = "mean", hypothesis = "b1 = b2")',
        'comparisons(mod, newdata = "mean", hypothesis = "am = vs")',
        'comparisons(mod, variables = "am", by = "cyl", hypothesis = "pairwise")')
      insight::format_error(msg)
    }
    rowlabels <- x$term
  }
  
  eval_string_function <- function(vec, hypothesis, rowlabels) {
    envir <- parent.frame()
    void <- sapply(
      seq_along(vec), function(i) {
        assign(rowlabels[i], vec[i], envir = envir)
      })
    out <- eval(parse(text = hypothesis), envir = envir)
    return(out)
  }
  
  if (!is.null(attr(lab, "names"))) {
    lab = attr(lab, "names")
  } else {
    lab = gsub("\\s+", "", lab)
  }
  
  draws <- attr(x, "posterior_draws")
  if (!is.null(draws)) {
    insight::check_if_installed("collapse", minimum_version = "1.9.0")
    tmp <- apply(
      draws,
      MARGIN = 2,
      FUN = eval_string_function,
      hypothesis = hypothesis,
      rowlabels = rowlabels)
    draws <- matrix(tmp, ncol = ncol(draws))
    out <- data.table::data.table(
      term = lab,
      tmp = collapse::dapply(draws, MARGIN = 1, FUN = collapse::fmedian))
  } else {
    out <- eval_string_function(
      x[["estimate"]],
      hypothesis = hypothesis,
      rowlabels = rowlabels)
    out <- data.table::data.table(
      term = lab,
      tmp = out)
  }
  
  data.table::setnames(out, old = "tmp", new = "estimate")
  attr(out, "posterior_draws") <- draws
  return(out)
}


expand_wildcard <- function(hyp, bmax, lab) {
  # Find all occurrences of b*
  bstar_indices <- gregexpr("b\\*", hyp)[[1]]
  if (bstar_indices[1] == -1) return(list(hyp, lab))
  bstar_count <- length(bstar_indices)
  if (bstar_count > 1) {
    insight::format_error("Error: More than one 'b*' substring found.")
  }
  
  # Replace b* with b1, b2, b3, ..., bmax
  labs <- character(bmax)
  result <- character(bmax)
  for (i in 1:bmax) {
    result[i] <- sub("b\\*", paste0("b", i), hyp)
    labs[i] <- sub("b\\*", paste0("b", i), lab)
  }
  
  return(list(result, labs))
}






#' An internal function to compute hypothesis from estimates
#'
#' @details Function adpated from \code{marginaleffects:::get_hypothesis}
#'   available at
#'   <https://github.com/vincentarelbundock/marginaleffects/blob/main/R/get_hypothesis.R>
#' 
#' @param data A data frame
#' @param from A character vector to be renamed \code{to}
#' @param to A character specifying the new names for \code{from}
#' @param to_title A character vector to be be converted to title case
#' @param to_lower A character vector to be be converted to lower case
#' @param to_upper A character vector to be be converted to upper case
#' 
#' @return A data frame
#' 
#' @author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#' 
#' @keywords internal
#' @noRd
#' 
rename_keyvars <- function(data, 
                           from, 
                           to, 
                           to_title = NULL, 
                           to_lower = NULL, 
                           to_upper = NULL) {
  
  if(length(from) != length(to)) stop("lenght mismatch")
  
  for (i in 1:length(from)) {
    if(from[i] %in% colnames(data) & !to[i] %in% colnames(data)) {
      data[[to[i]]] <- data[[from[i]]]
      data[[from[i]]] <- NULL
    }
  }
  
  if(!is.null(to_title)) {
    for (i in 1:length(to_title)) {
      if(!is.character(to_title[i])) to_title[i] <- deparse(to_title[i])
      if(any(grepl(to_title[i], colnames(data), ignore.case = T))) {
        # newname <- tools::toTitleCase(to_title[i])
        newname <- sub("(.)", "\\U\\1", to_title[i], perl = TRUE)
        colnames(data) <- gsub(to_title[i], newname, 
                               colnames(data), fixed = T)
      }
    }
  }
  
  
  if(!is.null(to_lower)) {
    for (i in 1:length(to_lower)) {
      if(!is.character(to_lower[i])) to_lower[i] <- deparse(to_lower[i])
      if(any(grepl(to_lower[i], colnames(data), ignore.case = T))) {
        newname <- base::tolower(to_lower[i])
        colnames(data) <- gsub(to_lower[i], newname, 
                               colnames(data), fixed = T)
      }
    }
  }
  
  if(!is.null(to_upper)) {
    for (i in 1:length(to_upper)) {
      if(!is.character(to_upper[i])) to_upper[i] <- deparse(to_upper[i])
      if(any(grepl(to_upper[i], colnames(data), ignore.case = T))) {
        newname <- base::toupper(to_upper[i])
        colnames(data) <- gsub(to_upper[i], newname, 
                               colnames(data), fixed = T)
      }
    }
  }
  
  return(data)
}

