


#' Plot \code{bsitar} model
#'
#' @description The \code{plot_bsitar} provides visualization of six different
#'   types of growth curves that are ploted by using the [ggplot2]. The
#'   \code{plot_bsitar} also allows users to make their own detailed plots from
#'   the data returned as a \code{data.frame}.
#'
#' @details The \code{plot_bsitar} is a generic function that allows
#'   visualization of following six curves: population average distance curve,
#'   population average velocity curve, individual-specific distance curves,
#'   individual-specific velocity curves, unadjusted individual growth curves
#'   (i.e, observed growth curves), and the adjusted individual growth curves
#'   (adjusted for the model estimated random effects). The \code{plot_bsitar}
#'   internally calls the [gparameters.bsitar] function to estimate and
#'   summaries the distance and velocity curves and to estimate growth
#'   parameters such as the age at peak growth velocity (APGV). The
#'   \code{plot_bsitar} in turn calls the [brms::fitted.brmsfit] and
#'   [brms::predict.brmsfit] functions to make inference from the posterior
#'   draws. Thus, \code{plot_bsitar} allows plotting fitted or predicted curves.
#'   See [brms::fitted.brmsfit] and [brms::predict.brmsfit] for details on these
#'   functions and the difference between fitted and predicted values.
#'
#' @param model An object of class \code{bsitar}.
#' @param opt A character string containing letter(s) corresponding to the
#'   following plotting options: 'd' for population average distance curve, 'v'
#'   for population average velocity curve, 'D' for individual-specific distance
#'   curves, 'V' for individual-specific velocity curves, 'u' for unadjusted
#'   individual-specific distance curves, and 'a' for adjusted
#'   individual-specific distance curves (adjusted for the random effects).
#'   Options 'd' and 'D' can not be specified simultaneously. Likewise, Options
#'   'v' and 'V' can not be specified simultaneously. All other combinations are
#'   allowed. For example, dvau', Dvau', dVau', DVau', or dvau'.
#' @param apv An optional logical (default \code{FALSE}) specifying whether or
#'   not to calculate and plot the age at peak velocity (APGV) when \code{opt})
#'   includes 'v' or 'V'.
#' @param bands A character string containing letter(s), or \code{NULL}
#'   (default) to indicate if CI bands to be plotted around the distance and
#'   velocity curves (and also the APGV). If \code{NULL}, no band plotted.
#'   Alternatively, user can specify a string with any one of the following or
#'   their combination(s): \code{'d'} for band around the distance curve,
#'   \code{'v} for band around the velocity curve, and \code{'p} for band around
#'   the the vertical line denoting the APGV parameter. The \code{'dvp'} will
#'   include CI bands for distance and velocity curves, and the APGV.
#' @param conf A numeric value (default \code{0.95}) to be used to compute the
#'   CI and hence the width of the \code{bands}. See [gparameters.bsitar] for
#'   further details.
#' 
#' @param trim A number (default 0) of long line segments to be excluded from
#'   plot with option 'u' or 'a'. See [sitar::plot.sitar] for details.
#' @param layout A character string defining the layout structure of the plot. A
#'   \code{'single'} (default) layout provides overlaid distance and velocity
#'   curves on a single plot when opt includes \code{'dv'}, \code{'Dv'},
#'   \code{'dV'} or \code{'DV'} options.  Similarly, when opt includes
#'   \code{'au'}, the adjusted and unadjusted curves are plotted as a single
#'   plot. When opt is a single letter (e.g., \code{'d'}. \code{'v'} \code{'D'},
#'   \code{'V'}, \code{'a'}, \code{'u'}), the \code{'single'} optiion is
#'   ignored. The alternative layout option, the \code{'facet'} uses
#'   [ggplot2::facet_wrap] to map and draw plot when \code{opt} include two or
#'   more letters.
#' @param linecolor The color of line used when layout is \code{'facet'}. The
#'   default is \code{NULL} which internally set the \code{linecolor} as
#'   \code{'grey50'}.
#' @param linecolor1 The color of first line when layout is \code{'single'}. For
#'   example, for \code{opt = 'dv'}, the color of distance line is controlled by
#'   the \code{linecolor1}. Default \code{NULL} will internally set
#'   \code{linecolor1} as \code{'orange2'}.
#' @param linecolor2 The color of second line when layout is \code{'single'}.
#'   For example, for \code{opt = 'dv'}, the color of velocity line is
#'   controlled by the \code{linecolor2}. Default \code{NULL} sets the color
#'   \code{'green4'} for \code{linecolor2}.
#' @param label.x An optional character string to label the x axis. When
#'   \code{NULL} (default), the x axis label is taken from the predictor (e.g.,
#'   age).
#' @param label.y An optional character string to label the y axis. When
#'   \code{NULL} (default), the y axis label is taken from the type of plot
#'   (e.g., distance, velocity etc.). Note that when layout option is
#'   \code{'facet'}, then y axis label is removed and instead the same label is
#'   used as a title.
#' @param legendpos An optional character string to specify the position of
#'   legends. When \code{NULL}, the legend position is set as 'bottom' for
#'   distance and velocity curves with \code{'single'} layout option.
#' @param linetype.apv An optional character string to specify the type of the
#'   vertical line drawn to mark the APGV. Default \code{NULL} sets the linetype
#'   as 'dotted'
#' @param linewidth.main An optional character string to specify the width of
#'   the the line for the distance and velocity curves. The default \code{NULL}
#'   will set it as 0.35.
#' @param linewidth.apv An optional character string to specify the width of the
#'   the vertical line drawn to mark the APGV. The default \code{NULL} will set
#'   it as 0.25.
#' @param linetype.groupby An optional argument to specify the line type for the
#'   distance and velocity curves when drawing plots for a model that includes
#'   factor covariate(s) or when visualising individual specific distance/velocity curves
#'   (default \code{NA}). Setting it to \code{NULL} will automatically sets the
#'   linetype for each factor level or individual This will also add legends for
#'   the factor level covariate or individuals whereas \code{NA} will set a
#'   'solid' line type and suppress legends. It is recommended to keep the
#'   default \code{NULL} option when plotting population average curves for when
#'   model included factor covariates because this would appropriately set the
#'   legends otherwise it is difficult to differentiate which curve belongs to
#'   which level of factor. For individual specific curves, the line type can be
#'   set to \code{NULL} when the number of individuals is small. However, when
#'   the number of individuals is large, \code{NA} is a better choice which
#'   prevents printing a large number of legends for each individual.
#'   
#' @param color.groupby An optional argument to specify the line color for
#'   distance and velocity curves when drawing plots for a model that includes
#'   factor covariate(s), or when visualising individual specific distance/velocity curves
#'   (default \code{NA}). Setting it to \code{NULL} will automatically sets the
#'   line color for each factor level or individual. This will also add legends for
#'   the factor level covariate or individuals. However, setting it as \code{NA} will set a
#'   'solid' line type and suppress legends. It is recommended to keep the
#'   default \code{NULL} option when plotting population average curves for
#'   factor covariates because this would appropriately set the
#'   legends otherwise it is difficult to differentiate which curve belongs to
#'   which level of the factor. For individual specific curves, the line color can be
#'   set to \code{NULL} when the number of individuals is small. However, when
#'   the number of individuals is large, \code{NA} is a better choice which
#'   prevents printing a large number of legends for each individual. 
#' @param band.alpha An optional numeric value to specify the transparency of
#'   the CI band(s) around the distance curve, velocity curve and the line
#'   indicating the APGV. The default \code{NULL} will set this value to 0.4.
#' @param returndata A logical (default \code{FALSE}) indicating whether to plot
#'   the data or return the data. If \code{TRUE}, the data is returned as a
#'   \code{data.frame}.
#'
#'@inheritParams  gparameters.bsitar
#'
#' @return A [ggplot2] object (default) or a \code{data.frame} when returndata
#'   is \code{TRUE}.
#'
#' @importFrom rlang .data
#' @importFrom graphics curve
#'
#' @export plot_bsitar.bsitar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Population average distance and velocity curves with default options
#' plot_bsitar(model, opt = 'dv')
#'
#' # Individual-specific distance and velocity curves with default options
#' plot_bsitar(model, opt = 'DV')
#'
#' # Population average distance and velocity curves with APGV
#' plot_bsitar(model, opt = 'dv', apv = TRUE)
#'
#' # Individual-specific distance and velocity curves with APGV
#' plot_bsitar(model, opt = 'DV', apv = TRUE)
#'
#' # Population average distance curve, velocity curve, and APGV with CI bands
#' plot_bsitar(model, opt = 'dv', apv = TRUE, bands = 'dvp')
#'
#' # Adjusted and unadjusted individual curves
#' plot_bsitar(model, opt = 'au')
#'
#' # Population average distance and velocity curves along with adjusted
#' # and unadjusted individual curves
#' plot_bsitar(model, opt = 'dvau')
#'
#' }
#' 
plot_bsitar.bsitar <- function(model,
                               opt = 'dv',
                               apv = FALSE,
                               bands = NULL,
                               conf = 0.95,
                               resp = NULL,
                               ndraws = NULL,
                               newdata = NULL,
                               summary = TRUE,
                               re_formula = NULL,
                               numeric_cov_at = NULL,
                               levels_id = NULL,
                               ipts = NULL,
                               seed = 123,
                               estimation_method = 'fitted',
                               robust = FALSE,
                               future = FALSE,
                               future_session = 'multisession',
                               cores = NULL,
                               trim = 0,
                               layout = 'single',
                               linecolor = NULL,
                               linecolor1 = NULL,
                               linecolor2 = NULL,
                               label.x = NULL,
                               label.y = NULL,
                               legendpos = NULL,
                               linetype.apv = NULL,
                               linewidth.main = NULL,
                               linewidth.apv = NULL,
                               linetype.groupby = NA,
                               color.groupby = NA,
                               band.alpha = NULL,
                               returndata = FALSE,...) {
  if (is.null(ndraws))
    ndraws  <- ndraws(model)
  else
    ndraws <- ndraws
  
  o <- post_processing_checks(model = model,
                              xcall = match.call(),
                              resp = resp)
  
  xcall <- strsplit(deparse(sys.calls()[[1]]), "\\(")[[1]][1]
  
  arguments <- get_args_(as.list(match.call())[-1], xcall)
  
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  probtitles <- probs[order(probs)] * 100
  probtitles <- paste("Q", probtitles, sep = "")
  set_names_  <- c('Estimate', 'Est.Error', probtitles)
 
  cores <- 1
  get.cores_ <- get.cores(arguments$cores)
  arguments$cores <- cores <-  get.cores_[['max.cores']] 
  .cores_ps <- get.cores_[['.cores_ps']]
  
  if (future) {
    if (future_session == 'multisession') {
      future::plan('multisession', workers = cores)
    } else if (future_session == 'multicore') {
      future::plan('multicore', workers = cores)
    }
  }
  
  newdata <- get.newdata(model, newdata = newdata, 
                         resp = resp, 
                         numeric_cov_at = numeric_cov_at,
                         levels_id = levels_id,
                         ipts = ipts)
  
  list_c <- attr(newdata, 'list_c')
  for (list_ci in names(list_c)) {
    assign(list_ci, list_c[[list_ci]])
  }
  check__ <- c('xvar', 'yvar', 'IDvar', 'cov_vars', 'cov_factor_vars', 
               'cov_numeric_vars', 'groupby_fstr', 'groupby_fistr', 'uvarby', 'subindicatorsi')
  
  for (check___ in check__) {
    if(!exists(check___)) assign(check___, NULL)
  }
  
  if(is.null(uvarby)) uvarby <- NA
  
  Xx <- xvar
  Yy <- yvar
  
  if (is.null(resp)) {
    resp_rev_ <- resp
  } else if (!is.null(resp)) {
    resp_rev_ <- paste0("_", resp)
  }
  
  
  if (is.null(bands)) {
    bands <- ''
  }
  
  if (grepl("d", opt, ignore.case = F) &
      grepl("D", opt, ignore.case = F)) {
    stop(
      "Options 'd' and 'D' can not be specified simultanously",
      "\n ",
      " Please check opt argument which is set as '",
      opt,
      "'",
      "\n ",
      " Either specify option 'd' for the population average distance curve",
      "\n ",
      "Or, else option 'D' for the individual specific distance curves",
      "\n ",
      "\n ",
      " Also note that you can combine option 'd' or 'D' with:",
      "\n ",
      "  - option 'v' (population average velocity curve) or",
      "\n ",
      "   option 'V' (individual specific velocity curves)",
      "\n ",
      " - option 'a' (adjusted curves) and 'u' (unadjusted curves)",
      "\n ",
      "\n ",
      "For example, opt = 'dvau', opt = 'dVau', opt = 'Dvau' or opt = 'DVau"
    )
  }
  
  if (grepl("v", opt, ignore.case = F) &
      grepl("V", opt, ignore.case = F)) {
    stop(
      "Options 'v' and 'V' can not be specified simultanously",
      "\n ",
      " Please check opt argument which is set as '",
      opt,
      "'",
      "\n ",
      " Either specify option 'v' for the population average velocity curve",
      "\n ",
      "Or, else option 'V' for the individual specific velocity curves",
      "\n ",
      "\n ",
      " Also note that you can combine option 'v' or 'V' with:",
      "\n ",
      "  - option 'd' (population average distance curve) or",
      "\n ",
      "   option 'D' (individual specific distance curves)",
      "\n ",
      " - option 'a' (adjusted curves) and 'u' (unadjusted curves)",
      "\n ",
      "\n ",
      "For example, opt = 'dvau', opt = 'dVau', opt = 'Dvau' or opt = 'DVau"
    )
  }
  
  if (grepl("p", bands, ignore.case = T) & summary) {
    stop(
      "To construct bands (e.g., 95%) around the parameter estimates",
      "\n ",
      " (such as APGV, PGV), they are first calculated for each",
      "\n ",
      " posterior draw and then summarised for the give conf limit.",
      "\n ",
      " Therefore,summary option must be set to FALSE"
    )
  }
  
  
 
 
  if (grepl("a", opt, ignore.case = F) |
      grepl("u", opt, ignore.case = F)) {
    
    testdata1 <- model$data %>% dplyr::select(dplyr::all_of(IDvar)) %>% droplevels() %>% 
      dplyr::mutate(groupbytest = interaction(dplyr::across(IDvar))) %>% 
      dplyr::select(groupbytest) %>% dplyr::ungroup()
    
    testdata2 <- newdata %>% dplyr::select(dplyr::all_of(IDvar)) %>% droplevels() %>% 
      dplyr::mutate(groupbytest = interaction(dplyr::across(IDvar))) %>% 
      dplyr::select(groupbytest) %>% dplyr::ungroup()
    
    
    if (!identical(testdata1, testdata2)) {
      warning(
        "Your have specified 'a' (adjusted curves) and/or 'u'",
        "\n  (unadjusted curves) in the opt argument (i.e., opt = 'au') ",
        "\n ",
        " but newdata is not identical to the data fitted.",
        "\n ",
        " Please note that adjusted and unadjusted curves will be",
        "\n ",
        " plotted using the original data fitted"
      )
    }
  }
 
  pv <- FALSE
  if (returndata & nchar(opt) > 1) {
    stop(
      "For returndata, please specify only one option at a time",
      "\n ",
      " (out of the total six optiona available, i.e., dvDVau)",
      "\n ",
      " For example, opt = 'd'"
    )
  }
  
  
  
  if (!grepl("v", opt, ignore.case = F) &
      !grepl("V", opt, ignore.case = F)) {
    apv <- arguments$apv <- FALSE
    pv <- arguments$pv <- FALSE
  }
  
  if(is.null(arguments$pv)) arguments$pv <- FALSE
  
  if(length(list(...)) != 0) arguments <- c(arguments, list(...))
  
  arguments$envir_ <- parent.frame()
  
  d. <- do.call(gparameters.bsitar, arguments)
  
  p. <- d.[['parameters']]
  probtitles <- d.[['probtitles']]
  groupby_str_d <- d.[['groupby_str_d']]
  groupby_str_v <- d.[['groupby_str_v']]
  
  
  d.[['parameters']] <- NULL
  d.[['probtitles']] <- NULL
  d.[['groupby_str_d']] <- NULL
  d.[['groupby_str_v']] <- NULL
  
  d. <- d. %>% do.call(rbind, .) %>% data.frame()
  row.names(d.) <- NULL
  
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  
  
  
  curve.d <- 'distance'
  curve.v <- 'velocity'
  
  name.apv <- "APGV"
  name.pv <- "PGV"
  
  x_minimum <- min(newdata[[Xx]])
  x_maximum <- max(newdata[[Xx]])
  x_minimum <- floor(x_minimum)
  x_maximum <- floor(x_maximum)
  
  single_plot_pair_color_dv_au <- c('black', 'red')
  
  
  
  if (nchar(opt) > 2) {
    layout <- 'facet'
  }
  
  if (is.null(linecolor)) {
    color_single <- "grey50"
  }
  
  if (is.null(linecolor1)) {
    color.d <- "orange2"
      color.adj <- "orange2"
  }
  if (is.null(linecolor2)) {
    color.v <- "green4"
      color.unadj <- "green4"
  }
  
  if (is.null(label.y)) {
    label.d <- firstup(curve.d)
    label.v <- firstup(curve.v)
    label.adj <- firstup('adjusted')
    label.unadj <- firstup('unadjusted')
  } else {
    label.d <- label.v <- label.y
    label.adj <- label.unadj <- label.y
  }
  
  if (is.null(label.x)) {
    label.x <- paste0(firstup(Xx), "")
  }
  
  if (is.null(legendpos)) {
    legendpos <- "bottom"
    legendpos.adj.unadj <- "topleft"
  } else {
    legendpos <- legendpos.adj.unadj <- legendpos
  }
  
  if (is.null(linetype.apv)) {
    linetype.apv <- 'dotted'
    linetype.pv <- 'dotted'
  } else {
    linetype.apv <- linetype.pv <- linetype.apv
  }
  
  if (is.null(linewidth.main)) {
    linewidth.main <- 0.5
  }
  
  if (is.null(linewidth.apv)) {
    linewidth.apv <- 0.5
    linewidth.pv <- 0.5
  } else {
    linewidth.apv <- linewidth.pv <- linewidth.apv
  }
  
  if (is.null(band.alpha)) {
    band.alpha <- 0.25
  }
  
  
  
  add_global_label <-
    function(pwobj,
             Xlab = NULL,
             Ylab = NULL,
             Xgap = 0.08,
             Ygap = 0.03,
             ...) {
      ylabgrob <- patchwork::plot_spacer()
      if (!is.null(Ylab)) {
        ylabgrob <- ggplot2::ggplot() +
          ggplot2::geom_text(ggplot2::aes(x = .5, y = .5),
                             label = Ylab,
                             angle = 90,
                             ...) +
          ggplot2::theme_void()
      }
      if (!is.null(Xlab)) {
        xlabgrob <- ggplot2::ggplot() +
          ggplot2::geom_text(ggplot2::aes(x = .5, y = .5), label = Xlab, ...) +
          ggplot2::theme_void()
      }
      if (!is.null(Ylab) & is.null(Xlab)) {
        return((ylabgrob + patchwork::patchworkGrob(pwobj)) +
                 patchwork::plot_layout(widths = 100 * c(Ygap, 1 - Ygap))
        )
      }
      if (is.null(Ylab) & !is.null(Xlab)) {
        return((ylabgrob + pwobj) +
                 (xlabgrob) +
                 patchwork::plot_layout(
                   heights = 100 * c(1 - Xgap, Xgap),
                   widths = c(0, 100),
                   design = "
                                   AB
                                   CC
                                   "
                 )
        )
      }
      if (!is.null(Ylab) & !is.null(Xlab)) {
        return((ylabgrob + pwobj) +
                 (xlabgrob) +
                 patchwork::plot_layout(
                   heights = 100 * c(1 - Xgap, Xgap),
                   widths = 100 * c(Ygap, 1 - Ygap),
                   design = "
                                   AB
                                   CC
                                   "
                 )
        )
      }
      return(pwobj)
    }
  
  transform.sec.axis <- function(primary, secondary, na.rm = TRUE) {
    from <- range(secondary, na.rm = na.rm)
    to   <- range(primary, na.rm = na.rm)
    zero_range <- function(x, tol = 1000 * .Machine$double.eps) {
      if (length(x) == 1) {
        return(TRUE)
      }
      if (length(x) != 2)
        stop("x must be length 1 or 2")
      if (any(is.na(x))) {
        return(NA)
      }
      if (x[1] == x[2]) {
        return(TRUE)
      }
      if (all(is.infinite(x))) {
        return(FALSE)
      }
      m <- min(abs(x))
      if (m == 0) {
        return(FALSE)
      }
      abs((x[1] - x[2]) / m) < tol
    }
    rescale.numeric_ <-
      function(x,
               to = c(0, 1),
               from = range(x, na.rm = TRUE, finite = TRUE),
               ...) {
        if (zero_range(from) || zero_range(to)) {
          return(ifelse(is.na(x), NA, mean(to)))
        }
        (x - from[1]) / diff(from) * diff(to) + to[1]
      }
    forward <- function(x) {
      rescale.numeric_(x, from = from, to = to)
    }
    reverse <- function(x) {
      rescale.numeric_(x, from = to, to = from)
    }
    list(fwd = forward, rev = reverse)
  }
  
  trimlines_ <-
    function(model,
             newdata,
             resp = NULL,
             ndraws = NULL,
             level = 0,
             trim = 0,
             ...) {
      if (is.null(ndraws))
        ndraws  <- ndraws(model)
      else
        ndraws <- ndraws
      
      o <-
        post_processing_checks(model = model,
                               xcall = match.call(),
                               resp = resp,
                               deriv = '')
      
      
      newdata.o <- newdata
      if (trim == 0)
        return(newdata)
      
      .x <- Xx
      .y <- Yy
      .id <- IDvar
      
      newdata <-
        with(newdata, newdata[order(newdata[[.id]], newdata[[.x]]), ])
      extra <- dplyr::as_tibble(diff(as.matrix(newdata[, 1:2])))
      extra[[.id]] <- newdata[[.id]][-1]
      did <- diff(as.integer(newdata[[.id]]))
      extra$dx <- extra[[.x]]
      extra[, 1:2] <- newdata[-1, 1:2] - extra[, 1:2] / 2
      extra <- extra[!did, ]
      
      if(!is.na(uvarby)) {
        extra[[subindicatorsi]] <- 1
      }
      
      if (level == 0) {
        re_formula <- NA
      } else if (level == 1) {
        re_formula <- NULL
      }
      if (estimation_method == 'fitted') {
        extra$ey <-
          fitted_.bsitar(
            model,
            resp = resp,
            newdata = extra,
            ndraws = ndraws,
            re_formula = re_formula,
            summary = TRUE
          )
      } else if (estimation_method == 'predict') {
        extra$ey <-
          predict_.bsitar(
            model,
            resp = resp,
            newdata = extra,
            ndraws = ndraws,
            re_formula = re_formula,
            summary = TRUE
          )
      }
      extra$ey <- extra$ey[, 1]
      extra <- extra %>%
        dplyr::mutate(dy = abs(extra[[.y]] - extra$ey),
                      xy = extra$dx / mad(extra$dx) + dy / mad(dy))
      outliers <- order(extra$xy, decreasing = TRUE)[1:trim]
      extra <- extra[outliers, 1:3]
      extra[[.y]] <- NA
      if(!is.na(uvarby)) {
        newdata_tt <- newdata
        common_colsnms <- intersect(colnames(newdata) , colnames(extra))
        newdata <-newdata %>% dplyr::select(all_of(common_colsnms))
      }
      newdata <- rbind(newdata, extra)
      newdata <-
        with(newdata, newdata[order(newdata[[.id]], newdata[[.x]]), ])
      
      if(!is.na(uvarby)) {
        tempotnames <- c(IDvar, Xx, Yy)
        tempot <- newdata_tt %>%  dplyr::select(-all_of(tempotnames))
        newdata <- cbind(newdata[-1, ], tempot) %>% data.frame()
      }
      
      newdata
    }
  
  
  xyadj_ <-
    function (model,
              x,
              y = NULL,
              id,
              resp = NULL,
              ndraws = NULL,
              newdata = NULL,
              abc = NULL,
              tomean = TRUE,...) {
      if (is.null(ndraws))
        ndraws  <- ndraws(model)
      else
        ndraws <- ndraws
      
     
      
      o <-
        post_processing_checks(model = model,
                               xcall = match.call(),
                               resp = resp,
                               deriv = '')
      
      newdata <- get.newdata(model, newdata = newdata,
                             resp = resp,
                             numeric_cov_at = NULL,
                             levels_id = levels_id,
                             ipts = NULL)
      
      
      
      if(!is.na(uvarby)) {
        newdata <- newdata %>%
          dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
      }
      
      xoffsetXnames <- 'xoffset'
      randomRnames <- 'random' 
      if (!is.null(resp)) randomRnames <- paste0(randomRnames, resp_rev_)
      if (!is.null(resp)) xoffsetXnames <- paste0(xoffsetXnames, resp_rev_)
      randomRnames <- model$model_info[[randomRnames]]
      xoffsetXnames <- model$model_info[[xoffsetXnames]]
      if (missing(x))
        x <- newdata[[Xx]]
      if (missing(y))
        y <- newdata[[Yy]]
      if (missing(id))
        id <- newdata[[IDvar]]
      
      if (is.null(resp)) {
        if (is.null(abc)) {
          re <- ranef(model)[[IDvar]][, 1 , ]
          abc <- re[match(id, rownames(re)), , drop = FALSE]
        }
        abc <- as.data.frame(abc)
        
      } else if (!is.null(resp)) {
        if (is.null(abc)) {
          re <- ranef(model)[[IDvar]]
          re <- re[match(id, rownames(re)) , 1, , drop = FALSE]
          re <- re[ , 1, grepl(paste0("^", resp), attr(re, "dimnames")[3][[1]])]
          abc <- re 
        }
        abc <- as.data.frame(abc)
      }
      
      colnames(abc) <- randomRnames
      
      fixef.df <- as.data.frame(fixef(model))
      fixef.df$names <- rownames(fixef.df)
      fixef.df <- fixef.df
      name_fixef.df_b <- paste0("b", "_Intercept")
      if (!is.null(resp)) name_fixef.df_b <- paste0(resp, "_", name_fixef.df_b)
      fixef.df_b <-
        subset(fixef.df, names == name_fixef.df_b)
      for (i in letters[1:4])
        if (!i %in% names(abc))
          abc[, i] <- 0
      xoffset <- xoffsetXnames
      if (!is.na(b0 <- as.numeric(fixef.df_b[1])))
        xoffset <- xoffset + b0
      x <- x - xoffset
      if (tomean) {
        x.adj <- (x - abc$b) * exp(abc$c) + xoffset
        y.adj <- y - abc$a - abc$d * x
      }
      else {
        x.adj <- x / exp(abc$c) + abc$b + xoffset
        y.adj <- y + abc$a + abc$d * x
      }
      out <- as.data.frame(as.factor(newdata[[IDvar]]))
      out <- cbind(x.adj, y.adj, out)
      colnames(out) <- c(Xx, Yy, IDvar)
      if(!is.na(uvarby)) {
        tempotnames <- c(IDvar, Xx, Yy)
        tempot <- newdata %>%  dplyr::select(-all_of(tempotnames))
        out <- cbind(out, tempot) %>% data.frame()
      }
      out
    }
  
  xyunadj_ <-
    function (model,
              x,
              y = NULL,
              id,
              resp = NULL,
              newdata = NULL,...) {
      
      o <-
        post_processing_checks(model = model,
                               xcall = match.call(),
                               resp = resp,
                               deriv = '')
      
      newdata <- get.newdata(model, newdata = newdata, resp = resp)
      
      if(!is.na(uvarby)) {
        newdata <- newdata %>%
          dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>% droplevels()
      }
      
      if (missing(x))
        x <- newdata[[Xx]]
      if (missing(y))
        y <- newdata[[Yy]]
      if (missing(id))
        id <- newdata[[IDvar]]
      out <- as.data.frame(as.factor(newdata[[IDvar]]))
      out <- cbind(x, y, out)
      colnames(out) <- c(Xx, Yy, IDvar)
      if(!is.na(uvarby)) {
        out[[uvarby]] <- resp
      }
      out
    }
  
  
  
  
  
  
  
  set_lines_colors <- function(plot, ngroups, 
                               linetype.groupby, color.groupby) {
    nrepvals <- ngroups
    
    if(is.null(linetype.groupby)) {
      linetype.groupby <- deparse(linetype.groupby)
    } else  if(is.na(linetype.groupby)) {
      linetype.groupby <- deparse(linetype.groupby)
    } else {
      linetype.groupby <- linetype.groupby
    }
    
    if(is.null(color.groupby)) {
      color.groupby <- deparse(color.groupby)
    } else  if(is.na(color.groupby)) {
      color.groupby <- deparse(color.groupby)
    } else {
      color.groupby <- color.groupby
    }
    
    # https://data.library.virginia.edu/setting-up-color-palettes-in-r/
    ggplotColors <- function(g){
      g <- g - 1
      d <- 360/g
      h <- cumsum(c(15, rep(d,g - 1)))
      O <- hcl(h = h, c = 100, l = 65)
      O <- c('black', O)
      O
    }
    
    # https://groups.google.com/g/ggplot2/c/XIcXU3KlxW0
    ggplotlines <- function(g){
      lineTypes1 <- c("solid", "22", "42", "44", "13", "1343", "73", "2262")
      # lineTypes1 <- c("solid", "solid", "solid", "13", "1343", "73", "2262")
      lineTypes2 <- apply(expand.grid(1:3, 1:3, 1:3, 1:3), 1, paste0, collapse="")
      lineTypes3 <- apply(expand.grid(1:2, 1:2, 1:2, 1:2), 1, paste0, collapse="")
      lineTypes <- c(lineTypes1, lineTypes2, lineTypes3)
      lineTypes[1:g]
    }
    
    default.set.line.groupby <- 'solid'
    default.set.color.groupby <- 'black'
      
    line.guide <- "none"
    color.guide <- "none"
    
    if(linetype.groupby == 'NA' & color.groupby == 'NA') {
      if(nrepvals > 1) {
        set.line.groupby <- rep(default.set.line.groupby, nrepvals)
        set.color.groupby <- rep(default.set.color.groupby, nrepvals)
        line.guide <- "none"
        color.guide <- "legend"
      }
    } # if(is.na(linetype.groupby) & is.na(color.groupby)) {
    
    if(linetype.groupby == 'NA' & color.groupby != 'NA') {
      set.line.groupby <- rep(default.set.line.groupby, nrepvals)
      if(nrepvals == 1) {
        if(color.groupby == 'NULL') {
          set.color.groupby <- default.set.color.groupby
        } else if(color.groupby != 'NULL') {
          set.color.groupby <- color.groupby[1]  
        }
      }
      
      if(nrepvals > 1) {
        set.line.groupby <- rep(default.set.line.groupby, nrepvals)
        if(color.groupby == 'NULL') {
          set.color.groupby <- ggplotColors(nrepvals)
        }
        
        if(color.groupby != 'NULL') {
          if(length(color.groupby) == nrepvals) {
            set.color.groupby <- color.groupby
          } else if(length(color.groupby) != nrepvals) {
            set.color.groupby <- rep(color.groupby, nrepvals)
          }
        }
        line.guide <- "none"
        color.guide <- "legend"
      }
    } # if(is.na(linetype.groupby) & !is.na(color.groupby)) {
    
    # print(linetype.groupby)
    # print(color.groupby)
    if(linetype.groupby != 'NA' & color.groupby == 'NA') {
      set.color.groupby <- rep(default.set.color.groupby, nrepvals)
      if(nrepvals == 1) {
        if(linetype.groupby == 'NULL') {
          set.line.groupby <- default.set.linetype.groupby
        } else if(linetype.groupby != 'NULL') {
          set.line.groupby <- linetype.groupby[1]  
        }
      }
      
      if(nrepvals > 1) {
        if(linetype.groupby == 'NULL') {
          set.line.groupby <- ggplotlines(nrepvals)
          if(length(set.line.groupby) < nrepvals) {
            set.line.groupby <- rep(set.line.groupby, nrepvals)
          }
        }
        
        if(linetype.groupby != 'NULL') {
          if(length(linetype.groupby) == nrepvals) {
            set.line.groupby <- linetype.groupby
          } else if(length(color.groupby) != nrepvals) {
            set.line.groupby <- rep(linetype.groupby, nrepvals)
          }
        }
        line.guide <- "legend"
        color.guide <- "none"
      }
    } # if(!is.na(linetype.groupby) & is.na(color.groupby)) {
    
    
    
    if(linetype.groupby != 'NA' & color.groupby != 'NA') {
      if(nrepvals == 1) {
        if(color.groupby == 'NULL') {
          set.color.groupby <- 'black'
        } else if(color.groupby != 'NULL') {
          set.color.groupby <- color.groupby[1]   
        }
        
        if(linetype.groupby == 'NULL') {
          set.line.groupby <- 'solid'
        } else if(linetype.groupby != 'NULL') {
          set.line.groupby <- linetype.groupby[1]   
        }
      }
      
      if(nrepvals > 1) {
        if(color.groupby == 'NULL') {
          set.color.groupby <- ggplotColors(nrepvals)
          if(length(set.color.groupby) < nrepvals) {
            set.color.groupby <- rep(set.color.groupby, nrepvals)
          }
        }
        if(linetype.groupby == 'NULL') {
          set.line.groupby <- ggplotlines(nrepvals)
          if(length(set.line.groupby) < nrepvals) {
            set.line.groupby <- rep(set.line.groupby, nrepvals)
          }
        }
        
        if(color.groupby != 'NULL') {
          if(length(color.groupby) == nrepvals) {
            set.color.groupby <- color.groupby
          } else if(length(color.groupby) != nrepvals) {
            set.color.groupby <- rep(color.groupby, nrepvals)
          }
        }
        if(linetype.groupby != 'NULL') {
          if(length(linetype.groupby) == nrepvals) {
            set.line.groupby <- linetype.groupby
          } else if(length(linetype.groupby) != nrepvals) {
            set.line.groupby <- rep(linetype.groupby, nrepvals)
          }
        }
        line.guide <- "none"
        color.guide <- "legend"
      }
    } # if(!is.na(linetype.groupby) & !is.na(set.color.groupby)) {
    
    # print(str(set.line.groupby))
    # print(str(set.color.groupby))
    suppressMessages({
      plot <- plot + 
        ggplot2::scale_linetype_manual(values=set.line.groupby, guide = line.guide) +
        ggplot2::scale_color_manual(values=set.color.groupby, guide = color.guide)
    })
    
    plot
  } # set_lines_colors
  
  
  
  set_lines_colors_ribbon <- function(plot, guideby = NULL) {
    getbuiltingg <- ggplot2::ggplot_build(plot)
    get_line_  <- getbuiltingg$data[[1]]["linetype"]
    get_color_ <- getbuiltingg$data[[1]]["colour"]
    get_fill_  <- getbuiltingg$data[[1]]["colour"]
    ngrpanels  <- getbuiltingg$data[[1]]["group"]
    get_line_  <- unique(unlist(get_line_))
    get_color_ <- unique(unlist(get_color_))
    get_fill_  <- unique(unlist(get_fill_))
    ngrpanels <- length(unique(unlist(ngrpanels)))
    
    if(length(get_line_) != ngrpanels) get_line_ <- rep(get_line_, ngrpanels)
    if(length(get_color_) != ngrpanels) get_color_ <- rep(get_color_, ngrpanels)
    if(length(get_fill_) != ngrpanels) get_fill_ <- rep(get_fill_, ngrpanels)
    
    setguide_line <- setguide_color <- setguide_fill <- 'none'
    if(is.null(guideby)) {
      setguide_line <- setguide_color <- setguide_fill <- 'none'
    } else if(guideby == 'line') {
      setguide_line <- 'legend'
    } else if(guideby == 'color') {
      setguide_color <- 'legend'
    } else if(guideby == 'fill') {
      setguide_fill <- 'legend'
    }
    
    suppressMessages({
      plot <- plot +
        ggplot2::scale_linetype_manual(values=get_line_, guide = setguide_line) +
        ggplot2::scale_color_manual(values=get_color_, guide = setguide_color) +
        ggplot2::scale_fill_manual(values=get_fill_, guide = setguide_fill)
    })
    
    plot
  }
  
  
  
  
  
  defaultW <- getOption("warn")
  options(warn = -1)
  if (grepl("d", opt, ignore.case = T) |
      grepl("v", opt, ignore.case = T)) {
    curves <- unique(d.$curve)
    if (length(curves) == 1) {
      layout <- 'facet'
    } else {
      layout <- layout
    }
    
    if (layout == 'facet')
      color.d <- color.v <- color_single
    
    if (grepl("d", opt, ignore.case = T)) {
      d.o <- d.
      index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
      dist.. <- substr(opt, index_opt, index_opt)
      if (grepl("^[[:upper:]]+$", dist..)) {
        d. <-
          d. %>% dplyr::mutate(groupby = interaction(dplyr::across(groupby_str_d)))
      } else if (!grepl("^[[:upper:]]+$", dist..)) {
        if (is.null(groupby_str_d))
          d. <- d. %>% dplyr::mutate(groupby = NA)
        if (!is.null(groupby_str_d))
          d. <-
            d. %>% dplyr::mutate(groupby =
                                   interaction(dplyr::across(groupby_str_d)))
      }
      
      
      
      
      if(is.na(d.[['groupby']][1])) {
        d.$groupby_line <- 'solid'
        d.$groupby_color <- 'black'
      } else {
        d.$groupby_line <- d.$groupby
        d.$groupby_color <- d.$groupby
      }
      
      plot.o.d <- d. %>% dplyr::filter(curve == curve.d) %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = Estimate,
            group = groupby,
            linetype = groupby_line,
            color = groupby_color
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::scale_x_continuous(breaks = seq(x_minimum, x_maximum, 1)) +
        ggplot2::labs(title = label.d) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
     
      plot.o.d <- set_lines_colors(plot.o.d, length(unique(d.[['groupby']])), 
                                   linetype.groupby = linetype.groupby, 
                                   color.groupby = color.groupby)
      
      
      

      
      if (grepl("d", bands, ignore.case = T)) {
        plot.o.d <- plot.o.d +
          ggplot2::geom_ribbon(
            data = d. %>% dplyr::filter(curve == curve.d),
            ggplot2::aes(
              ymin = .data[[paste0(probtitles[1], '')]],
              ymax = .data[[paste0(probtitles[2], '')]],
              group = groupby,
              linetype = groupby_line,
              color = groupby_color,
              fill = groupby_color
            ),
            alpha = band.alpha
          )
        plot.o.d <- set_lines_colors_ribbon(plot.o.d, guideby = 'color')
      }
      
      
      
    
      
      d. <- d.o
      if ('curve' %in% names(d.)) {
        d.out <- d. %>% dplyr::select(-curve)
      } else {
        d.out <- d.
      }
    }
    if (!grepl("d", opt, ignore.case = T)) {
      plot.o.d <- NULL
    }
    
    
    if (grepl("v", opt, ignore.case = T)) {
      d.o <- d.
      index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
      velc.. <- substr(opt, index_opt, index_opt)
      if (grepl("^[[:upper:]]+$", velc..)) {
        d. <-
          d. %>% dplyr::mutate(groupby = interaction(dplyr::across(groupby_str_v)))
      } else if (!grepl("^[[:upper:]]+$", velc..)) {
        if (is.null(groupby_str_v))
          d. <- d. %>% dplyr::mutate(groupby = NA)
        if (!is.null(groupby_str_v))
          d. <-
            d. %>% dplyr::mutate(groupby =
                                   interaction(dplyr::across(groupby_str_v)))
      }
      
      
     
      if(is.na(d.[['groupby']][1])) {
        d.$groupby_line <- 'solid'
        d.$groupby_color <- 'black'
      } else {
        d.$groupby_line <- d.$groupby
        d.$groupby_color <- d.$groupby
      }
      
      
      plot.o.v <- d. %>% dplyr::filter(curve == curve.v) %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = Estimate,
            group = groupby,
            linetype = groupby_line,
            color = groupby_color
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::scale_x_continuous(breaks = seq(x_minimum, x_maximum, 1)) +
        ggplot2::labs(title = label.v) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
      
      plot.o.v <- set_lines_colors(plot.o.v, length(unique(d.[['groupby']])), 
                                   linetype.groupby = linetype.groupby, 
                                   color.groupby = color.groupby)
      
      
      
      if (grepl("v", bands, ignore.case = T)) {
        plot.o.v <- plot.o.v +
          ggplot2::geom_ribbon(
            data = d. %>% dplyr::filter(curve == curve.v),
            ggplot2::aes(
              ymin = .data[[paste0(probtitles[1], '')]],
              ymax = .data[[paste0(probtitles[2], '')]],
              group = groupby,
              linetype = groupby_line,
              color = groupby_color,
              fill = groupby_color
            ),
            alpha = band.alpha
          )
        plot.o.v <- set_lines_colors_ribbon(plot.o.v, guideby = 'color')
      }
      
      if (pv) {
        data_hline <- p. %>% dplyr::filter(Parameter == name.pv)
        plot.o.v <- plot.o.v +
          ggplot2::geom_hline(
            data = data_hline,
            ggplot2::aes(yintercept = .data[['Estimate']]),
            linewidth = linewidth.pv,
            linetype = linetype.pv
          )
        
        if (grepl("p", bands, ignore.case = T)) {
          plot.o.v <- plot.o.v +
            ggplot2::annotate(
              "rect",
              xmin = -Inf,
              xmax = Inf,
              ymin = (data_hline[[paste0(probtitles[1], '')]]),
              ymax = (data_hline[[paste0(probtitles[2], '')]]),
              alpha = band.alpha
            )
        }
      }
      
      if (apv) {
        data_vline <- p. %>% dplyr::filter(Parameter == name.apv)
        
        plot.o.v <- plot.o.v +
          ggplot2::geom_vline(
            data = data_vline,
            ggplot2::aes(xintercept = .data[['Estimate']]),
            linewidth = linewidth.apv,
            linetype = linetype.apv
          )
        
        if (grepl("p", bands, ignore.case = T)) {
          plot.o.v <- plot.o.v +
            ggplot2::annotate(
              "rect",
              xmin = (data_vline[[paste0(probtitles[1], '')]]),
              xmax = (data_vline[[paste0(probtitles[2], '')]]),
              ymin = -Inf,
              ymax = Inf,
              alpha = band.alpha
            )
        }
      }
      d. <- d.o
      if ('curve' %in% names(d.)) {
        d.out <- d. %>% dplyr::select(-curve)
      } else {
        d.out <- d.
      }
    }
    
    if (!grepl("v", opt, ignore.case = T)) {
      plot.o.v <- NULL
    }
    
    if (length(curves) > 1 & layout == 'facet') {
      plot.o <- patchwork::wrap_plots(plot.o.d + plot.o.v) %>%
        add_global_label(
          Xlab = label.x,
          Ylab = "",
          size = 5,
          Xgap = 0.08,
          Ygap = 0.04
        ) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = legendpos,
                       legend.direction = 'horizontal')
      
    } else if (length(curves) == 1 & layout == 'facet') {
      if (!is.null(plot.o.d)) {
        plot.o <- plot.o.d +
          ggplot2::labs(x = label.x, y = label.d) +
          ggplot2::theme(plot.title = ggplot2::element_blank()) +
          ggplot2::theme(legend.position = legendpos,
                         legend.direction = 'horizontal')
      } else if (!is.null(plot.o.v)) {
        plot.o <- plot.o.v +
          ggplot2::labs(x = label.x, y = label.v) +
          ggplot2::theme(plot.title = ggplot2::element_blank()) +
          ggplot2::theme(legend.position = legendpos,
                         legend.direction = 'horizontal')
      }
    }
    
    
    if (length(curves) > 1 & layout == 'single') {
      data_d <- subset(d., curve == "distance")
      data_v <- subset(d., curve == "velocity")
      by_join_ <- c(IDvar, Xx, groupby_str_d)
      by_join_ <- unique(by_join_)
      data_dv <- dplyr::left_join(data_d, data_v, by = by_join_)
      data_dv.o <- data_dv
      if (grepl("d", opt, ignore.case = T)) {
        index_opt <- gregexpr("d", opt, ignore.case = T)[[1]]
        dist.. <- substr(opt, index_opt, index_opt)
        if (grepl("^[[:upper:]]+$", dist..)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(groupby = interaction(dplyr::across(groupby_str_d))) %>%
            dplyr::mutate(groupby.x = interaction(dplyr::across(groupby_str_d)))
        } else if (!grepl("^[[:upper:]]+$", dist..)) {
          if (is.null(groupby_str_d)) {
            data_dv <- data_dv %>% dplyr::mutate(groupby = NA) %>%
              dplyr::mutate(groupby.x = NA)
          } else if (!is.null(groupby_str_d)) {
            data_dv <- data_dv %>%
              dplyr::mutate(groupby =
                              interaction(dplyr::across(groupby_str_d))) %>%
              dplyr::mutate(groupby.x = interaction(dplyr::across(groupby_str_d)))
          }
        }
      }
      
      if (grepl("v", opt, ignore.case = T)) {
        index_opt <- gregexpr("v", opt, ignore.case = T)[[1]]
        velc.. <- substr(opt, index_opt, index_opt)
        if (grepl("^[[:upper:]]+$", velc..)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(groupby = interaction(dplyr::across(groupby_str_v))) %>%
            dplyr::mutate(groupby.y = groupby)
        } else if (!grepl("^[[:upper:]]+$", velc..)) {
          if (is.null(groupby_str_v)) {
            data_dv <- data_dv %>% dplyr::mutate(groupby = NA) %>%
              dplyr::mutate(groupby.y = NA)
          } else if (!is.null(groupby_str_v)) {
            data_dv <- data_dv %>%
              dplyr::mutate(groupby =
                              interaction(dplyr::across(groupby_str_v))) %>%
              dplyr::mutate(groupby.y = interaction(dplyr::across(groupby_str_v)))
          }
        }
      }
      
      if (grepl("^[[:upper:]]+$", dist..) &
          !grepl("^[[:upper:]]+$", velc..)) {
        if (is.null(groupby_str_v)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(groupby.x = interaction(dplyr::across(groupby_str_d)),
                          groupby.y = NA)
        } else if (!is.null(groupby_str_v)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(groupby.x = interaction(dplyr::across(groupby_str_d)),
                          groupby.y =
                            interaction(dplyr::across(groupby_str_v)))
        }
      }
      
      if (!grepl("^[[:upper:]]+$", dist..) &
          grepl("^[[:upper:]]+$", velc..)) {
        if (is.null(groupby_str_d)) {
          data_dv <-
            data_dv %>% dplyr::mutate(groupby.x = NA,
                                      groupby.y =
                                        interaction(dplyr::across(groupby_str_v)))
        } else if (!is.null(groupby_str_d)) {
          data_dv <-
            data_dv %>%
            dplyr::mutate(groupby.x = interaction(dplyr::across(groupby_str_d)),
                          groupby.y =
                            interaction(dplyr::across(groupby_str_v)))
        }
      }
      
      t.s.axis <-
        with(data_dv, transform.sec.axis(Estimate.x, Estimate.y))
      
      
      if(is.na(data_dv[['groupby.x']][1])) {
        legendlabs_mult_singel <- c('Distance', 'Velocity')
        legendlabs_mult_color <- single_plot_pair_color_dv_au
        legendlabs_mult_line <- c('solid', 'solid')
        data_dv$groupby_line.x <- legendlabs_mult_singel[1]
        data_dv$groupby_color.x <- legendlabs_mult_singel[1]
          data_dv$groupby_line.y <- legendlabs_mult_singel[2]
          data_dv$groupby_color.y <- legendlabs_mult_singel[2]
      } else {
        data_dv$groupby_line.x <- data_dv$groupby.x
        data_dv$groupby_color.x <- data_dv$groupby.x
        data_dv$groupby_line.y <- data_dv$groupby.y
        data_dv$groupby_color.y <- data_dv$groupby.y
        legendlabs_mult_mult <- unique(data_dv[['groupby.x']])
      }
      
      
      
      plot.o <- data_dv %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = Estimate.x,
            group = groupby.x,
            linetype = groupby_line.x,
            colour = groupby_color.x
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = t.s.axis$fwd(Estimate.y),
            group = groupby.y,
            linetype = groupby_line.y,
            colour = groupby_color.y
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::scale_y_continuous(sec.axis =
                                      ggplot2::sec_axis(~ t.s.axis$rev(.),
                                                        name = label.v)) +
        ggplot2::labs(x = label.x, y = label.d, color = "") +
        ggplot2::scale_x_continuous(breaks = seq(x_minimum, x_maximum, 1)) +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90))
      
      
      getbuiltingg <- ggplot2::ggplot_build(plot.o)
      get_line_  <- getbuiltingg$data[[1]]["linetype"]
      get_color_ <- getbuiltingg$data[[1]]["colour"]
      get_fill_  <- getbuiltingg$data[[1]]["colour"]
      ngrpanels  <- getbuiltingg$data[[1]]["group"]
      ngrpanels <- length(unique(unlist(ngrpanels)))
      get_line_ <- unique(unlist(get_line_))
      get_color_ <- unique(unlist(get_color_))
      if(length(get_line_) != ngrpanels) get_line_ <- rep(get_line_, ngrpanels)
      if(length(get_color_) != ngrpanels) get_color_ <- rep(get_color_, ngrpanels)
      
      # These will be carried forward for ribbon also (below)
      if(ngrpanels > 1) {
        get_line_ <- get_line_
        get_color_ <- get_color_
        legendlabs_ <- legendlabs_mult_mult
      } else if(ngrpanels == 1) {
        get_line_ <- legendlabs_mult_line
        get_color_ <- legendlabs_mult_color
        legendlabs_ <- legendlabs_mult_singel
      }
      
      plot.o <- plot.o + 
        ggplot2::scale_linetype_manual(values=get_line_, guide = 'none') +
        ggplot2::scale_color_manual(breaks=legendlabs_, values=get_color_)
      
      
      if (grepl("d", bands, ignore.case = T)) {
        plot.o <- plot.o +
          ggplot2::geom_ribbon(
            data = data_dv,
            ggplot2::aes(
              ymin = .data[[paste0(probtitles[1], '.x')]],
              ymax = .data[[paste0(probtitles[2], '.x')]],
              group = groupby.x,
             linetype = groupby_line.x,
             # color = groupby_color.x,
              fill = groupby_color.x,
            ),
            alpha = band.alpha
          )
        # Need only oncen which is done for vel plot below
        # plot.o <- plot.o +
        #   ggplot2::scale_fill_manual(values=get_color_, guide = 'none')
      }
      
      if (grepl("v", bands, ignore.case = T)) {
        plot.o <- plot.o +
          ggplot2::geom_ribbon(
            data = data_dv,
            ggplot2::aes(
              ymin = t.s.axis$fwd(.data[[paste0(probtitles[1], '.y')]]),
              ymax = t.s.axis$fwd(.data[[paste0(probtitles[2], '.y')]]),
              group = groupby.y,
              linetype = groupby_line.y,
              # color = groupby_color.y,
              fill = groupby_color.y,
            ),
            alpha = band.alpha
          )
        plot.o <- plot.o +
          ggplot2::scale_fill_manual(values=get_color_, guide = 'none')
      }
      
      
      if (pv) {
        data_hline <- p. %>% dplyr::filter(Parameter == name.pv)
        plot.o <- plot.o +
          ggplot2::geom_hline(
            data = data_hline,
            ggplot2::aes(yintercept = t.s.axis$fwd(.data[['Estimate']])),
            linewidth = linewidth.pv,
            linetype = linetype.pv
          )
        
        if (grepl("p", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::annotate(
              "rect",
              xmin = -Inf,
              xmax = Inf,
              ymin = t.s.axis$fwd(data_hline[[paste0(probtitles[1], '')]]),
              ymax = t.s.axis$fwd(data_hline[[paste0(probtitles[2], '')]]),
              alpha = band.alpha
            )
        }
      }
      
      if (apv) {
        data_vline <- p. %>% dplyr::filter(Parameter == name.apv)
        plot.o <- plot.o +
          ggplot2::geom_vline(
            data = data_vline,
            ggplot2::aes(xintercept = .data[['Estimate']]),
            linewidth = linewidth.apv,
            linetype = linetype.apv
          )
        
        if (grepl("p", bands, ignore.case = T)) {
          plot.o <- plot.o +
            ggplot2::annotate(
              "rect",
              xmin = (data_vline[[paste0(probtitles[1], '')]]),
              xmax = (data_vline[[paste0(probtitles[2], '')]]),
              ymin = -Inf,
              ymax = Inf,
              alpha = band.alpha
            )
        }
      }
      data_dv <- data_dv.o
      if ('curve' %in% names(data_dv)) {
        d.out <- data_dv %>% dplyr::select(-curve)
      } else {
        d.out <- data_dv
      }
    }
  }
  
  
  
  ############
  groupby_str_au <- groupby_fistr
  
  if (grepl("a", opt, ignore.case = T) |
      grepl("u", opt, ignore.case = T)) {
    if (!is.null(cov_vars)) {
      stop("Adjusted curves not yet supported for model with covariate(s)")
    }
    if (grepl("a", opt, ignore.case = T)) {
      
      xyadj_ed <- xyadj_(model, resp = resp)
      
      out_a_ <-
        d.out <- trimlines_(model, resp = resp, newdata = xyadj_ed, trim = trim)
      out_a_ <-
        out_a_ %>%
        dplyr::mutate(groupby = interaction(dplyr::across(groupby_str_au)))
      
      # x_minimum_a_ <- floor(min(out_a_[[Xx]]))
      # x_maximum_a_ <- ceiling(max(out_a_[[Xx]]))
      x_minimum_a_ <- x_minimum
      x_maximum_a_ <- x_maximum
      
      out_a_ <- out_a_[out_a_[[Xx]] >= x_minimum_a_ & out_a_[[Xx]] <= x_maximum_a_, ]
      
      out_a_ <- out_a_ %>% dplyr::mutate(groupby.x = groupby, groupby.y = groupby.x)
      
     
      
      if(is.na(out_a_[['groupby']][1])) {
        legendlabs_mult_singel <- c('Distance', 'Velocity')
        legendlabs_mult_color <- single_plot_pair_color_dv_au
        legendlabs_mult_line <- c('solid', 'solid')
        out_a_$groupby_line.x <- legendlabs_mult_singel[1]
        out_a_$groupby_color.x <- legendlabs_mult_singel[1]
        out_a_$groupby_line.y <- legendlabs_mult_singel[2]
        out_a_$groupby_color.y <- legendlabs_mult_singel[2]
      } else {
        out_a_$groupby_line.x <- out_a_$groupby.x
        out_a_$groupby_color.x <- out_a_$groupby.x
        out_a_$groupby_line.y <- out_a_$groupby.y
        out_a_$groupby_color.y <- out_a_$groupby.y
        legendlabs_mult_mult <- unique(out_a_[['groupby']])
      }
      
      ############# set to NULL
      # linetype.groupby <- NULL
      
      plot.o.a <- out_a_ %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = !!as.name(Yy),
            group = groupby.x,
            # linetype = groupby_line.x,
             colour = groupby_color.x
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::labs(x = label.x, y = label.d, color = "") +
        # ggplot2::scale_color_manual(values = c(color.adj)) +
        ggplot2::scale_x_continuous(breaks =
                                      seq(x_minimum_a_, x_maximum_a_, 1)) +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(y = paste0("Adjusted ", "individual curves")) +
        ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90)) +
        ggplot2::labs(title = label.adj) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      

      getbuiltingg <- ggplot2::ggplot_build(plot.o.a)
      get_line_  <- getbuiltingg$data[[1]]["linetype"]
      get_color_ <- getbuiltingg$data[[1]]["colour"]
      get_fill_  <- getbuiltingg$data[[1]]["colour"]
      ngrpanels  <- getbuiltingg$data[[1]]["group"]
      ngrpanels <- length(unique(unlist(ngrpanels)))
      get_line_ <- unique(unlist(get_line_))
      get_color_ <- unique(unlist(get_color_))
      if(length(get_line_) != ngrpanels) get_line_ <- rep(get_line_, ngrpanels)
      if(length(get_color_) != ngrpanels) get_color_ <- rep(get_color_, ngrpanels)
      
      # These will be carried forward for ribbon also (below)
      if(ngrpanels > 1) {
        get_line_ <- get_line_
        get_color_ <- get_color_
        legendlabs_ <- legendlabs_mult_mult
      } else if(ngrpanels == 1) {
        get_line_ <- legendlabs_mult_line
        get_color_ <- legendlabs_mult_color
        legendlabs_ <- legendlabs_mult_singel
      }
      
      plot.o.a <- plot.o.a +
        ggplot2::scale_linetype_manual(values=get_line_, guide = 'none') +
        ggplot2::scale_color_manual(breaks=legendlabs_, values=get_color_)
      
      
      
      if (nchar(opt) == 1) {
        plot.o.a <- plot.o.a +
          ggplot2::theme(plot.title = ggplot2::element_blank())
      } else if (nchar(opt) > 1) {
        plot.o.a <- plot.o.a +
          ggplot2::theme(axis.title.y = ggplot2::element_blank())
      }
      if (!grepl("u", opt, ignore.case = T)) {
        suppressMessages({
          # plot.o.a <- plot.o.a +
          #   ggplot2::scale_color_manual(values = c(color_single))
        })
      }
    } else if (!grepl("a", opt, ignore.case = T)) {
      plot.o.a <- NULL
    }
   
    
    if (grepl("u", opt, ignore.case = T)) {
      xyadj_ed <- xyunadj_(model, resp = resp)
      out_u_ <-
        d.out <- trimlines_(model, resp = resp, newdata = xyadj_ed, trim = trim)
      out_u_ <-
        out_u_ %>%
        dplyr::mutate(groupby = interaction(dplyr::across(groupby_str_au)))
      
      out_u_ <- out_u_ %>% dplyr::mutate(groupby.x = groupby, groupby.y = groupby.x)
      
      
      if(is.na(out_u_[['groupby']][1])) {
        legendlabs_mult_singel <- c('Distance', 'Velocity')
        legendlabs_mult_color <- single_plot_pair_color_dv_au
        legendlabs_mult_line <- c('solid', 'solid')
        out_u_$groupby_line.x <- legendlabs_mult_singel[1]
        out_u_$groupby_color.x <- legendlabs_mult_singel[1]
        out_u_$groupby_line.y <- legendlabs_mult_singel[2]
        out_u_$groupby_color.y <- legendlabs_mult_singel[2]
      } else {
        out_u_$groupby_line.x <- out_u_$groupby.x
        out_u_$groupby_color.x <- out_u_$groupby.x
        out_u_$groupby_line.y <- out_u_$groupby.y
        out_u_$groupby_color.y <- out_u_$groupby.y
        legendlabs_mult_mult <- unique(out_u_[['groupby']])
      }
      
      plot.o.u <- out_u_ %>%
        ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
        ggplot2::geom_line(
          ggplot2::aes(
            y = !!as.name(Yy),
            group = groupby,
            group = groupby.y,
            # linetype = groupby_line.y,
            colour = groupby_color.y
          ),
          linewidth = linewidth.main
        ) +
        ggplot2::labs(x = label.x, y = label.d, color = "") +
        # ggplot2::scale_color_manual(values = c(color.unadj)) +
        ggplot2::scale_x_continuous(breaks = seq(x_minimum, x_maximum, 1)) +
        jtools::theme_apa(legend.pos = legendpos) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(y = paste0("Unadjusted ", "individual curves")) +
        ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90)) +
        ggplot2::labs(title = label.unadj) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
      getbuiltingg <- ggplot2::ggplot_build(plot.o.u)
      get_line_  <- getbuiltingg$data[[1]]["linetype"]
      get_color_ <- getbuiltingg$data[[1]]["colour"]
      get_fill_  <- getbuiltingg$data[[1]]["colour"]
      ngrpanels  <- getbuiltingg$data[[1]]["group"]
      ngrpanels <- length(unique(unlist(ngrpanels)))
      get_line_ <- unique(unlist(get_line_))
      get_color_ <- unique(unlist(get_color_))
      if(length(get_line_) != ngrpanels) get_line_ <- rep(get_line_, ngrpanels)
      if(length(get_color_) != ngrpanels) get_color_ <- rep(get_color_, ngrpanels)
      
      # These will be carried forward for ribbon also (below)
      if(ngrpanels > 1) {
        get_line_ <- get_line_
        get_color_ <- get_color_
        legendlabs_ <- legendlabs_mult_mult
      } else if(ngrpanels == 1) {
        get_line_ <- legendlabs_mult_line
        get_color_ <- legendlabs_mult_color
        legendlabs_ <- legendlabs_mult_singel
      }
      
      plot.o.u <- plot.o.u +
        ggplot2::scale_linetype_manual(values=get_line_, guide = 'none') +
        ggplot2::scale_color_manual(breaks=legendlabs_, values=get_color_)
      
      
      if (!grepl("a", opt, ignore.case = T)) {
        suppressMessages({
          # plot.o.u <- plot.o.u +
          #   ggplot2::scale_color_manual(values = c(color_single))
        })
      }
      if (nchar(opt) == 1) {
        plot.o.u <- plot.o.u +
          ggplot2::theme(plot.title = ggplot2::element_blank())
      } else if (nchar(opt) > 1) {
        plot.o.u <- plot.o.u +
          ggplot2::theme(axis.title.y = ggplot2::element_blank())
      }
    } else if (!grepl("u", opt, ignore.case = T)) {
      plot.o.u <- NULL
    }
    
    
    if (grepl("a", opt, ignore.case = T) &
        !grepl("u", opt, ignore.case = T)) {
      plot.o <- plot.o.a
    } else if (!grepl("a", opt, ignore.case = T) &
               grepl("u", opt, ignore.case = T)) {
      plot.o <- plot.o.u
    } else if (grepl("a", opt, ignore.case = T) &
               grepl("u", opt, ignore.case = T)) {
      if (layout == 'facet') {
        out_a_u_ <-
          d.out <- out_a_ %>% dplyr::mutate(curve = 'Adjusted') %>%
          dplyr::bind_rows(., out_u_ %>%
                             dplyr::mutate(curve = 'Unadjusted')) %>%
          data.frame()
        
        # x_minimum_a_ <- floor(min(out_a_[[Xx]]))
        # x_maximum_a_ <- ceiling(max(out_a_[[Xx]]))
        x_minimum_a_ <- x_minimum
        x_maximum_a_ <- x_maximum
        
        out_a_u_ <- out_a_u_[out_a_u_[[Xx]] >= x_minimum_a_ & out_a_u_[[Xx]] <= x_maximum_a_, ]
        
        out_a_u_ <- out_a_u_ %>% dplyr::mutate(groupby.x = groupby, groupby.y = groupby.x)
        
        
        if(is.na(out_a_u_[['groupby']][1])) {
          legendlabs_mult_singel <- c('Distance', 'Velocity')
          legendlabs_mult_color <- single_plot_pair_color_dv_au
          legendlabs_mult_line <- c('solid', 'solid')
          out_a_u_$groupby_line.x <- legendlabs_mult_singel[1]
          out_a_u_$groupby_color.x <- legendlabs_mult_singel[1]
          out_a_u_$groupby_line.y <- legendlabs_mult_singel[2]
          out_a_u_$groupby_color.y <- legendlabs_mult_singel[2]
        } else {
          out_a_u_$groupby_line.x <- out_a_u_$groupby.x
          out_a_u_$groupby_color.x <- out_a_u_$groupby.x
          out_a_u_$groupby_line.y <- out_a_u_$groupby.y
          out_a_u_$groupby_color.y <- out_a_u_$groupby.y
          legendlabs_mult_mult <- unique(out_a_u_[['groupby']])
        }
        
        plot.o <- out_a_u_ %>%
          ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
          ggplot2::geom_line(
            ggplot2::aes(
              y = !!as.name(Yy),
              group = groupby.x,
              linetype = groupby_line.x,
              colour = groupby_color.x
            ),
            linewidth = linewidth.main
          ) +
          ggplot2::labs(x = label.x,
                        y = label.d,
                        color = "") +
          # ggplot2::scale_color_manual(values = c(color_single)) +
          ggplot2::scale_x_continuous(breaks =
                                        seq(x_minimum_a_, x_maximum_a_, 1)) +
          jtools::theme_apa(legend.pos = legendpos) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::theme(axis.title.y.right =
                           ggplot2::element_text(angle = 90)) +
          ggplot2::theme(legend.position = "none") +
          ggplot2::labs(y = paste0("Individual curves")) +
          ggplot2::facet_wrap(~ curve, scales = 'free_x')
        
        
        getbuiltingg <- ggplot2::ggplot_build(plot.o)
        get_line_  <- getbuiltingg$data[[1]]["linetype"]
        get_color_ <- getbuiltingg$data[[1]]["colour"]
        get_fill_  <- getbuiltingg$data[[1]]["colour"]
        ngrpanels  <- getbuiltingg$data[[1]]["group"]
        ngrpanels <- length(unique(unlist(ngrpanels)))
        get_line_ <- unique(unlist(get_line_))
        get_color_ <- unique(unlist(get_color_))
        if(length(get_line_) != ngrpanels) get_line_ <- rep(get_line_, ngrpanels)
        if(length(get_color_) != ngrpanels) get_color_ <- rep(get_color_, ngrpanels)
        
        # These will be carried forward for ribbon also (below)
        if(ngrpanels > 1) {
          get_line_ <- get_line_
          get_color_ <- get_color_
          legendlabs_ <- legendlabs_mult_mult
        } else if(ngrpanels == 1) {
          get_line_ <- legendlabs_mult_line
          get_color_ <- legendlabs_mult_color
          legendlabs_ <- legendlabs_mult_singel
        }
        
        plot.o <- plot.o +
          ggplot2::scale_linetype_manual(values=get_line_, guide = 'none') +
          ggplot2::scale_color_manual(breaks=legendlabs_, values=get_color_)
      
      }
      
      if (layout == 'single') {
        
        plot.o <- out_a_ %>%
          ggplot2::ggplot(., ggplot2::aes(!!as.name(Xx))) +
          ggplot2::geom_line(
            data = out_u_,
            ggplot2::aes(
              y = !!as.name(Yy),
              group = groupby,
              linetype = linetype.groupby,
              colour = label.unadj
            ),
            linewidth = linewidth.main
          ) +
          ggplot2::geom_line(
            data = out_a_,
            ggplot2::aes(
              y = !!as.name(Yy),
              group = groupby,
              linetype = linetype.groupby,
              colour = label.adj
            ),
            linewidth = linewidth.main
          ) +
          ggplot2::labs(x = label.x,
                        y = label.d,
                        color = "") +
          
          # ggplot2::scale_color_manual(values = c(color.unadj, color.adj)) +
          ggplot2::scale_color_manual(values = single_plot_pair_color_dv_au) +
          ggplot2::scale_x_continuous(breaks =
                                        seq(x_minimum_a_, x_maximum_a_, 1)) +
          jtools::theme_apa(legend.pos = legendpos.adj.unadj) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::labs(y = paste0("Individual curves")) +
          ggplot2::theme(axis.title.y.right = ggplot2::element_text(angle = 90))
      }
    }
  }
  
  if (nchar(opt) > 2) {
    if (!exists('plot.o.d'))
      plot.o.d <- NULL
    if (!exists('plot.o.v'))
      plot.o.v <- NULL
    if (!exists('plot.o.a'))
      plot.o.a <- NULL
    if (!exists('plot.o.u'))
      plot.o.u <- NULL
    
    suppressMessages({
      if (!is.null(plot.o.d)) {
        plot.o.d <- plot.o.d +
          ggplot2::scale_color_manual(values = c(color_single)) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank())
      }
      if (!is.null(plot.o.v)) {
        plot.o.v <- plot.o.v +
          ggplot2::scale_color_manual(values = c(color_single)) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank())
      }
      if (!is.null(plot.o.a)) {
        plot.o.a <- plot.o.a +
          ggplot2::scale_color_manual(values = c(color_single)) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank())
      }
      if (!is.null(plot.o.u)) {
        plot.o.u <- plot.o.u +
          ggplot2::scale_color_manual(values = c(color_single)) +
          ggplot2::theme(axis.title.x = ggplot2::element_blank())
      }
    })
    
    plot.list <- list(plot.o.d, plot.o.v, plot.o.a, plot.o.u)
    plot.list <- plot.list[lengths(plot.list) != 0]
    
    plot.o <- patchwork::wrap_plots(plot.list,
                                    ncol = 2, nrow = NULL) %>%
      add_global_label(
        Xlab = label.x,
        Ylab = "",
        size = 5,
        Xgap = 0.08,
        Ygap = 0.04
      )
    plot.o <-
      plot.o +  patchwork::plot_layout(guides = "collect")
  }
  
  if (!returndata) {
    print(plot.o)
    options(warn = defaultW)
    return(plot.o)
  } else if (returndata) {
    return(d.out)
  }
}

#' @rdname plot_bsitar.bsitar
#' @export
plot_bsitar <- function(model, ...) {
  UseMethod("plot_bsitar")
}

