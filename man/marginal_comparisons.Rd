% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal_comparisons.R
\name{marginal_comparisons.bgmfit}
\alias{marginal_comparisons.bgmfit}
\alias{marginal_comparisons}
\alias{marginal_comparison}
\title{Estimate and compare growth curves}
\usage{
\method{marginal_comparisons}{bgmfit}(
  model,
  resp = NULL,
  dpar = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  newdata = NULL,
  datagrid = NULL,
  re_formula = NA,
  allow_new_levels = FALSE,
  sample_new_levels = "gaussian",
  xrange = 1,
  digits = 2,
  numeric_cov_at = NULL,
  aux_variables = NULL,
  levels_id = NULL,
  avg_reffects = NULL,
  idata_method = NULL,
  ipts = NULL,
  seed = 123,
  future = FALSE,
  future_session = "multisession",
  future_splits = NULL,
  future_method = "future",
  future_re_expose = NULL,
  usedtplyr = FALSE,
  usecollapse = TRUE,
  cores = NULL,
  fullframe = FALSE,
  average = FALSE,
  plot = FALSE,
  showlegends = NULL,
  variables = NULL,
  deriv = 0,
  model_deriv = TRUE,
  method = "pkg",
  marginals = NULL,
  pdrawso = FALSE,
  pdrawsp = FALSE,
  pdrawsh = FALSE,
  comparison = "difference",
  type = NULL,
  by = FALSE,
  conf_level = 0.95,
  transform = NULL,
  transform_draws = NULL,
  cross = FALSE,
  wts = NULL,
  hypothesis = NULL,
  equivalence = NULL,
  eps = NULL,
  constrats_by = NULL,
  constrats_at = NULL,
  reformat = NULL,
  estimate_center = NULL,
  estimate_interval = NULL,
  dummy_to_factor = NULL,
  verbose = FALSE,
  expose_function = FALSE,
  usesavedfuns = NULL,
  clearenvfuns = NULL,
  funlist = NULL,
  xvar = NULL,
  idvar = NULL,
  itransform = NULL,
  newdata_fixed = NULL,
  envir = NULL,
  ...
)

marginal_comparisons(model, ...)

marginal_comparison(model, ...)
}
\arguments{
\item{model}{An object of class \code{bgmfit}.}

\item{resp}{A character string (default \code{NULL}) to specify the response
variable when processing posterior draws for \code{univariate_by} and
\code{multivariate} models. See \code{\link[=bsitar]{bsitar()}} for details on
\code{univariate_by} and \code{multivariate} models.}

\item{dpar}{Optional name of a predicted distributional parameter.
If specified, expected predictions of this parameters are returned.}

\item{ndraws}{A positive integer indicating the number of posterior draws to
use in estimation. If \code{NULL} (default), all draws are used.}

\item{draw_ids}{An integer specifying the specific posterior draw(s) to use
in estimation (default \code{NULL}).}

\item{newdata}{An optional data frame for estimation. If \code{NULL}
(default), \code{newdata} is retrieved from the \code{model}.}

\item{datagrid}{A data frame or named list for setting up a custom grid of
predictor values to evaluate the quantities of interest. If \code{NULL}
(default), no custom grid is used. The grid can be constructed using
\code{\link[marginaleffects:datagrid]{marginaleffects::datagrid()}}. If \code{datagrid = list()}, essential
arguments such as \code{model} and \code{newdata} are inferred
automatically from the respective arguments in the model fit.}

\item{re_formula}{Option to indicate whether or not to include
individual/group-level effects in the estimation. When \code{NA} (default),
individual-level effects are excluded, and population average growth
parameters are computed. When \code{NULL}, individual-level effects are
included in the computation, and the resulting growth parameters are
individual-specific. In both cases (\code{NA} or \code{NULL}), continuous
and factor covariates are appropriately included in the estimation.
Continuous covariates are set to their means by default (see
\code{numeric_cov_at} for details), while factor covariates remain
unaltered, allowing for the estimation of covariate-specific population
average and individual-specific growth parameters.}

\item{allow_new_levels}{A flag indicating if new levels of group-level
effects are allowed (defaults to \code{FALSE}). Only relevant if
\code{newdata} is provided.}

\item{sample_new_levels}{Indicates how to sample new levels for grouping
factors specified in \code{re_formula}. This argument is only relevant if
\code{newdata} is provided and \code{allow_new_levels} is set to
\code{TRUE}. If \code{"uncertainty"} (default), each posterior sample for a
new level is drawn from the posterior draws of a randomly chosen existing
level. Each posterior sample for a new level may be drawn from a different
existing level such that the resulting set of new posterior draws
represents the variation across existing levels. If \code{"gaussian"},
sample new levels from the (multivariate) normal distribution implied by the
group-level standard deviations and correlations. This options may be useful
for conducting Bayesian power analysis or predicting new levels in
situations where relatively few levels where observed in the old_data. If
\code{"old_levels"}, directly sample new levels from the existing levels,
where a new level is assigned all of the posterior draws of the same
(randomly chosen) existing level.}

\item{xrange}{An integer to set the predictor range (e.g., age) when
executing the interpolation via \code{ipts}. By default, \code{NULL} sets
the individual-specific predictor range. Setting \code{xrange = 1} applies
the same range for individuals within the same higher grouping variable
(e.g., study). Setting \code{xrange = 2} applies an identical range across
the entire sample. Alternatively, a numeric vector (e.g., \code{xrange =
  c(6, 20)}) can be provided to set the range within the specified values.}

\item{digits}{An integer (default \code{2}) for rounding the estimates to the
specified number of decimal places using \code{\link[base:Round]{base::round()}}.}

\item{numeric_cov_at}{An optional (named list) argument to specify the value
of continuous covariate(s). The default \code{NULL} option sets the
continuous covariate(s) to their mean. Alternatively, a named list can be
supplied to manually set these values. For example, \code{numeric_cov_at =
  list(xx = 2)} will set the continuous covariate variable 'xx' to 2. The
argument \code{numeric_cov_at} is ignored when no continuous covariates are
included in the model.}

\item{aux_variables}{An optional argument to specify the variable(s) that can
be passed to the \code{ipts} argument (see below). This is useful when
fitting location-scale models and measurement error models. If
post-processing functions throw an error such as \code{variable 'x' not
  found in either 'data' or 'data2'}, consider using \code{aux_variables}.}

\item{levels_id}{An optional argument to specify the \code{ids} for the
hierarchical model (default \code{NULL}). It is used only when the model is
applied to data with three or more levels of hierarchy. For a two-level
model, \code{levels_id} is automatically inferred from the model fit. For
models with three or more levels, \code{levels_id} is inferred from the
model fit under the assumption that hierarchy is specified from the lowest
to the uppermost level, i.e., \code{id} followed by \code{study}, where
\code{id} is nested within \code{study}. However, it is not guaranteed that
\code{levels_id} is sorted correctly, so it is better to set it manually
when fitting a model with three or more levels of hierarchy.}

\item{avg_reffects}{An optional argument (default \code{NULL}) to calculate
(marginal/average) curves and growth parameters, such as APGV and PGV. If
specified, it must be a named list indicating the \code{over} (typically a
level 1 predictor, such as age), \code{feby} (fixed effects, typically a
factor variable), and \code{reby} (typically \code{NULL}, indicating that
parameters are integrated over the random effects). For example,
\code{avg_reffects = list(feby = 'study', reby = NULL, over = 'age')}.}

\item{idata_method}{A character string to indicate the interpolation method.
The number of interpolation points is set by the \code{ipts} argument.
Available options for \code{idata_method} are \emph{method 1} (specified as
\code{'m1'}) and \emph{method 2} (specified as \code{'m2'}).
\itemize{
\item \emph{Method 1} (\code{'m1'}) is adapted from the \pkg{iapvbs} package
and is documented
\href{https://rdrr.io/github/Zhiqiangcao/iapvbs/src/R/exdata.R}{here}.
\item \emph{Method 2} (\code{'m2'}) is based on the \pkg{JMbayes} package
and is documented
\href{https://github.com/drizopoulos/JMbayes/blob/master/R/dynPred_lme.R}{here}.
The \code{'m1'} method works by internally constructing the data frame
based on the model configuration, while the \code{'m2'} method uses the
exact data frame from the model fit, accessible via \code{fit$data}. If
\code{idata_method = NULL} (default), method \code{'m2'} is automatically
selected. Note that method \code{'m1'} may fail in certain cases,
especially when the model includes covariates (particularly in
\code{univariate_by} models). In such cases, it is recommended to use
method \code{'m2'}.
}}

\item{ipts}{An integer to set the length of the predictor variable for
generating a smooth velocity curve. If \code{NULL}, the original values are
returned. If an integer (e.g., \code{ipts = 10}, default), the predictor is
interpolated. Note that these interpolations do not alter the range of the
predictor when calculating population averages and/or individual-specific
growth curves.}

\item{seed}{An integer (default \code{123}) that is passed to the estimation
method to ensure reproducibility.}

\item{future}{A logical value (default \code{FALSE}) to specify whether or
not to perform parallel computations. If set to \code{TRUE}, the
\code{\link[future.apply:future_lapply]{future.apply::future_sapply()}} function is used to summarize the posterior
draws in parallel.}

\item{future_session}{A character string specifying the session type when
\code{future = TRUE}. The \code{'multisession'} (default) option sets the
multisession environment, while the \code{'multicore'} option sets up a
multicore session. Note that \code{'multicore'} is not supported on Windows
systems. For more details, see \code{\link[future.apply:future_lapply]{future.apply::future_sapply()}}.}

\item{future_splits}{A list (default \code{NULL}) that can be an unnamed
numeric list, a logical value, or a numeric vector of length 1 or 2. It is
used to split the processing of posterior draws into smaller subsets for
parallel computation.
\itemize{
\item If passed as a list (e.g., \code{future_splits = list(1:6, 7:10)}),
each sequence of
numbers is passed to the \code{draw_ids} argument.
\item If passed as a numeric vector (e.g., \code{future_splits = c(10, 2)}),
the first element
specifies the number of draws (see \code{draw_ids}) and the second element
indicates the number of splits. The splits are created using
\code{\link[parallel:splitIndices]{parallel::splitIndices()}}.
\item If passed as a numeric vector of length 1, the first element is
internally set as the
number of draws (\code{ndraws} or \code{draw_ids}) depending on which one
is not \code{NULL}.
\item If \code{TRUE}, a numeric vector for \code{future_splits} is created
based on the number
of draws (\code{ndraws}) and the number of cores (\code{cores}).
\item If \code{FALSE}, \code{future_splits} is ignored.
The use case for \code{future_splits} is to save memory and improve
performance, especially on \code{Linux} systems when \code{future::plan()}
is set to \code{multicore}. Note: on Windows systems, R processes may not
be freed automatically when using \code{'multisession'}. In such cases, the
R processes can be interrupted using \code{\link[installr:kill_all_Rscript_s]{installr::kill_all_Rscript_s()}}.
}}

\item{future_method}{A character string (default \code{'future'}) to specify
the method for parallel computation. Options include:
\itemize{
\item \code{'future'}: Uses \code{\link[future:future]{future::future()}} along with
\code{\link[future.apply:future_lapply]{future.apply::future_lapply()}} for parallel execution.
\item \code{'foreach'}: Uses \code{\link[foreach:foreach]{foreach::foreach()}} with the
\code{'dofuture'} function from the \code{doFuture} package for parallel
execution.
}}

\item{future_re_expose}{A logical (default \code{NULL}) to indicate whether
to re-expose \code{Stan} functions when \code{future = TRUE}. This is
especially relevant when \code{\link[future:plan]{future::plan()}} is set to \code{'multisession'},
as already exposed C++ \code{Stan} functions cannot be passed across
multiple sessions.
\itemize{
\item When \code{future_re_expose = NULL} (the default), \code{future_re_expose}
is automatically set to \code{TRUE} for the \code{'multisession'} plan.
\item It is advised to explicitly set \code{future_re_expose = TRUE} for speed
gains when using parallel processing with \code{future = TRUE}.
}}

\item{usedtplyr}{A logical (default \code{FALSE}) indicating whether to use
the \pkg{dtplyr} package for summarizing the draws. This package uses
\pkg{data.table} as a back-end. It is useful when the data has a large
number of observations. For typical use cases, it does not make a
significant performance difference. The \code{usedtplyr} argument is
evaluated only when \code{method = 'custom'}.}

\item{usecollapse}{A logical (default \code{FALSE}) to indicate whether to
use the \pkg{collapse} package for summarizing the draws.}

\item{cores}{The number of cores to be used for parallel computations if
\code{future = TRUE}. On non-Windows systems, this argument can be set
globally via the \code{mc.cores} option. By default, \code{NULL}, the
number of cores is automatically determined using
\code{\link[future:re-exports]{future::availableCores()}}, and it will use the maximum number of cores
available minus one (i.e., \code{future::availableCores() - 1}).}

\item{fullframe}{A logical value indicating whether to return a
\code{fullframe} object in which \code{newdata} is bound to the summary
estimates. Note that \code{fullframe} cannot be used with \code{summary =
  FALSE}, and it is only applicable when \code{idata_method = 'm2'}. A
typical use case is when fitting a \code{univariate_by} model. This option
is mainly for internal use.}

\item{average}{A logical indicating whether to call
\code{\link[marginaleffects:comparisons]{marginaleffects::comparisons()}} (if \code{FALSE}) or
\code{\link[marginaleffects:comparisons]{marginaleffects::avg_comparisons()}} (if \code{TRUE}). Default is
\code{FALSE}.}

\item{plot}{A logical indicating whether to plot the comparisons using
\code{\link[marginaleffects:plot_comparisons]{marginaleffects::plot_comparisons()}}. Default is \code{FALSE}.}

\item{showlegends}{A logical value to specify whether to show legends
(\code{TRUE}) or not (\code{FALSE}). If \code{NULL} (default), the value of
\code{showlegends} is internally set to \code{TRUE} if \code{re_formula =
  NA}, and \code{FALSE} if \code{re_formula = NULL}.}

\item{variables}{A named list specifying the level 1 predictor, such as
\code{age} or \code{time}, used for estimating growth parameters in the
current use case. The \code{variables} list is set via the \code{esp}
argument (default value is \code{1e-6}). If \code{variables} is
\code{NULL}, the relevant information is retrieved internally from the
\code{model}. Alternatively, users can define \code{variables} as a named
list, e.g., \code{variables = list('x' = 1e-6)} where \code{'x'} is the
level 1 predictor. By default, \code{variables = list('age' = 1e-6)} in the
\pkg{marginaleffects} package, as velocity is usually computed by
differentiating the distance curve using the \code{dydx} approach. When
using this default, the argument \code{deriv} is automatically set to
\code{0} and \code{model_deriv} to \code{FALSE}. If parameters are to be
estimated based on the model's first derivative, \code{deriv} must be set
to \code{1} and \code{variables} will be defined as \code{variables =
  list('age' = 0)}. Note that if the default behavior is used (\code{deriv =
  0} and \code{variables = list('x' = 1e-6)}), additional arguments cannot be
passed to \code{variables}. In contrast, when using an alternative approach
(\code{deriv = 0} and \code{variables = list('x' = 0)}), additional options
can be passed to the \code{\link[marginaleffects:comparisons]{marginaleffects::comparisons()}} and
\code{\link[marginaleffects:comparisons]{marginaleffects::avg_comparisons()}} functions.}

\item{deriv}{A numeric value indicating whether to estimate parameters based
on the differentiation of the distance curve or the model's first
derivative.}

\item{model_deriv}{A logical value specifying whether to estimate the
velocity curve from the derivative function or by differentiating the
distance curve. Set \code{model_deriv = TRUE} for functions that require
the velocity curve, such as \code{growthparameters()} and
\code{plot_curves()}. Set it to \code{NULL} for functions that use the
distance curve (i.e., fitted values), such as \code{loo_validation()} and
\code{plot_ppc()}.}

\item{method}{A character string specifying whether to compute comparisons
using the \pkg{marginaleffects} machinery (\code{method = 'pkg'}) or via
custom functions for efficiency (\code{method = 'custom'}). Default is
\code{'pkg'}. The \code{'custom'} method is useful when testing hypotheses
and should be used cautiously.}

\item{marginals}{A \code{list}, \code{data.frame}, or \code{tibble} of
velocity returned by the \pkg{marginaleffects} functions (default
\code{NULL}). This is only evaluated when \code{method = 'custom'}. The
\code{marginals} can be the output from \pkg{marginaleffects} functions or
posterior draws from \code{marginaleffects::posterior_draws()}. The
\code{marginals} argument is primarily used for internal purposes only.}

\item{pdrawso}{A character string (default \code{FALSE}) to indicate whether
to return the original posterior draws for parameters. Options include:
\itemize{
\item \code{'return'}: returns the original posterior draws,
\item \code{'add'}: adds the original posterior draws to the outcome.
}
When \code{pdrawso = TRUE}, the default behavior is \code{pdrawso =
  'return'}. Note that the posterior draws are returned before calling
\code{\link[marginaleffects:posterior_draws]{marginaleffects::posterior_draws()}}.}

\item{pdrawsp}{A character string (default \code{FALSE}) to indicate whether
to return the posterior draws for parameters. Options include:
\itemize{
\item \code{'return'}: returns the posterior draws for parameters,
\item \code{'add'}: adds the posterior draws to the outcome.
}
When \code{pdrawsp = TRUE}, the default behavior is \code{pdrawsp =
  'return'}. The \code{pdrawsp} represent the parameter estimates for each of
the posterior samples, and the summary of these are the estimates returned.}

\item{pdrawsh}{A character string (default \code{FALSE}) to indicate whether
to return the posterior draws for parameter contrasts. Options include:
\itemize{
\item \code{'return'}: returns the posterior draws for contrasts.
}
The summary of posterior draws for parameters is the default returned
object. The \code{pdrawsh} represent the contrast estimates for each of the
posterior samples, and the summary of these are the contrast returned.}

\item{comparison}{A character string specifying the comparison type for
growth parameter estimation. Options are \code{'difference'} and
\code{'differenceavg'}. This argument sets up the internal function for
estimating parameters using \code{\link[sitar:getPeakTrough]{sitar::getPeak()}}, \code{\link[sitar:getPeakTrough]{sitar::getTakeoff()}}, and
\code{\link[sitar:getPeakTrough]{sitar::getTrough()}} functions. These options are restructured according to
the user-specified \code{hypothesis} argument.}

\item{type}{string indicates the type (scale) of the predictions used to
compute contrasts or slopes. This can differ based on the model
type, but will typically be a string such as: "response", "link", "probs",
or "zero". When an unsupported string is entered, the model-specific list of
acceptable values is returned in an error message. When \code{type} is \code{NULL}, the
first entry in the error message is used by default.}

\item{by}{Aggregate unit-level estimates (aka, marginalize, average over). Valid inputs:
\itemize{
\item \code{FALSE}: return the original unit-level estimates.
\item \code{TRUE}: aggregate estimates for each term.
\item Character vector of column names in \code{newdata} or in the data frame produced by calling the function without the \code{by} argument.
\item Data frame with a \code{by} column of group labels, and merging columns shared by \code{newdata} or the data frame produced by calling the same function without the \code{by} argument.
\item See examples below.
\item For more complex aggregations, you can use the \code{FUN} argument of the \code{hypotheses()} function. See that function's documentation and the Hypothesis Test vignettes on the \code{marginaleffects} website.
}}

\item{conf_level}{numeric value between 0 and 1. Confidence level to use to build a confidence interval.}

\item{transform}{string or function. Transformation applied to unit-level estimates and confidence intervals just before the function returns results. Functions must accept a vector and return a vector of the same length. Support string shortcuts: "exp", "ln"}

\item{transform_draws}{A function is applied to individual draws from the
posterior distribution before computing summaries (default \code{NULL}).
The argument \code{transform_draws} is derived from the
\code{\link[marginaleffects:predictions]{marginaleffects::predictions()}} function and should not be confused with
the \code{transform} argument from the deprecated
\code{\link[brms:posterior_predict.brmsfit]{brms::posterior_predict()}} function. It's important to note that for both
\code{\link[marginaleffects:predictions]{marginaleffects::predictions()}} and \code{\link[marginaleffects:predictions]{marginaleffects::avg_predictions()}},
the \code{transform_draws} argument takes precedence over the
\code{transform} argument.}

\item{cross}{\itemize{
\item \code{FALSE}: Contrasts represent the change in adjusted predictions when one predictor changes and all other variables are held constant.
\item \code{TRUE}: Contrasts represent the changes in adjusted predictions when all the predictors specified in the \code{variables} argument are manipulated simultaneously (a "cross-contrast").
}}

\item{wts}{logical, string or numeric: weights to use when computing average predictions, contrasts or slopes. These weights only affect the averaging in \verb{avg_*()} or with the \code{by} argument, and not unit-level estimates. See \code{?weighted.mean}
\itemize{
\item string: column name of the weights variable in \code{newdata}. When supplying a column name to \code{wts}, it is recommended to supply the original data (including the weights variable) explicitly to \code{newdata}.
\item numeric: vector of length equal to the number of rows in the original data or in \code{newdata} (if supplied).
\item FALSE: Equal weights.
\item TRUE: Extract weights from the fitted object with \code{insight::find_weights()} and use them when taking weighted averages of estimates. Warning: \code{newdata=datagrid()} returns a single average weight, which is equivalent to using \code{wts=FALSE}
}}

\item{hypothesis}{specify a hypothesis test or custom contrast using a number , formula, string equation, vector, matrix, or function.
\itemize{
\item Number: The null hypothesis used in the computation of Z and p (before applying \code{transform}).
\item String: Equation to specify linear or non-linear hypothesis tests. Two-tailed tests must include an equal \code{=} sign. One-tailed tests must start with \code{<} or \code{>}. If the terms in \code{coef(object)} uniquely identify estimates, they can be used in the formula. Otherwise, use \code{b1}, \code{b2}, etc. to identify the position of each parameter. The \verb{b*} wildcard can be used to test hypotheses on all estimates. When the hypothesis string represents a two-sided equation, the \code{estimate} column holds the value of the left side minus the right side of the equation. If a named vector is used, the names are used as labels in the output. Examples:
\itemize{
\item \code{hp = drat}
\item \code{hp + drat = 12}
\item \code{b1 + b2 + b3 = 0}
\item \verb{b* / b1 = 1}
\item \verb{<= 0}
\item \verb{>= -3.5}
\item \code{b1 >= 10}
}
\item Formula: \code{lhs ~ rhs | group}
\itemize{
\item \code{lhs}
\itemize{
\item \code{ratio} (null = 1)
\item \code{difference} (null = 0)
\item Leave empty for default value
}
\item \code{rhs}
\itemize{
\item \code{pairwise} and \code{revpairwise}: pairwise differences between estimates in each row.
\item \code{reference}: differences between the estimates in each row and the estimate in the first row.
\item \code{sequential}: difference between an estimate and the estimate in the next row.
\item \code{meandev}: difference between an estimate and the mean of all estimates.
\item `meanotherdev: difference between an estimate and the mean of all other estimates, excluding the current one.
\item \code{poly}: polynomial contrasts, as computed by the \code{stats::contr.poly()} function.
\item \code{helmert}: Helmert contrasts, as computed by the \code{stats::contr.helmert()} function. Contrast 2nd level to the first, 3rd to the average of the first two, and so on.
\item \code{trt_vs_ctrl}: difference between the mean of estimates (except the first) and the first estimate.
\item \code{I(fun(x))}: custom function to manipulate the vector of estimates \code{x}. The function \code{fun()} can return multiple (potentially named) estimates.
}
\item \code{group} (optional)
\itemize{
\item Column name of \code{newdata}. Conduct hypothesis tests withing subsets of the data.
}
\item Examples:
\itemize{
\item \code{~ poly}
\item \code{~ sequential | groupid}
\item \code{~ reference}
\item \code{ratio ~ pairwise}
\item \code{difference ~ pairwise | groupid}
\item \code{~ I(x - mean(x)) | groupid}
\item \verb{~ I(\\(x) c(a = x[1], b = mean(x[2:3]))) | groupid}
}
}
\item Matrix or Vector: Each column is a vector of weights. The the output is the dot product between these vectors of weights and the vector of estimates. The matrix can have column names to label the estimates.
\item Function:
\itemize{
\item Accepts an argument \code{x}: object produced by a \code{marginaleffects} function or a data frame with column \code{rowid} and \code{estimate}
\item Returns a data frame with columns \code{term} and \code{estimate} (mandatory) and \code{rowid} (optional).
\item The function can also accept optional input arguments: \code{newdata}, \code{by}, \code{draws}.
\item This function approach will not work for Bayesian models or with bootstrapping. In those cases, it is easy to use \code{get_draws()} to extract and manipulate the draws directly.
}
\item See the Examples section below and the vignette: \url{https://marginaleffects.com/chapters/hypothesis.html}
\item Warning: When calling \code{predictions()} with \code{type="invlink(link)"} (the default in some models), \code{hypothesis} is tested and p values are computed on the link scale.
}}

\item{equivalence}{Numeric vector of length 2: bounds used for the two-one-sided test (TOST) of equivalence, and for the non-inferiority and non-superiority tests. See Details section below.}

\item{eps}{NULL or numeric value which determines the step size to use when
calculating numerical derivatives: (f(x+eps)-f(x))/eps. When \code{eps} is
\code{NULL}, the step size is 0.0001 multiplied by the difference between
the maximum and minimum values of the variable with respect to which we
are taking the derivative. Changing \code{eps} may be necessary to avoid
numerical problems in certain models.}

\item{constrats_by}{A character vector (default \code{NULL}) specifying the
variables by which hypotheses should be tested (e.g., for post-draw
comparisons). These variables should be a subset of the variables in the
\code{by} argument of \code{\link[marginaleffects:comparisons]{marginaleffects::comparisons()}}.}

\item{constrats_at}{A character vector (default \code{NULL}) specifying the
values at which hypotheses should be tested. Useful for large estimates to
limit testing to fewer rows.}

\item{reformat}{A logical (default \code{TRUE}) to reformat the output from
\code{\link[marginaleffects:comparisons]{marginaleffects::comparisons()}} as a data.frame, with column names such as
\code{conf.low} as \code{Q2.5}, \code{conf.high} as \code{Q97.5}, and
dropping unnecessary columns like \code{term}, \code{contrast}, and
\code{predicted}.}

\item{estimate_center}{A character string (default \code{NULL}) specifying
how to center estimates: either \code{'mean'} or \code{'median'}. This
option sets the global options as follows:
\code{options("marginaleffects_posterior_center" = "mean")} or
\code{options("marginaleffects_posterior_center" = "median")}. These global
options are restored upon function exit using \code{\link[base:on.exit]{base::on.exit()}}.}

\item{estimate_interval}{A character string (default \code{NULL}) to specify
the type of credible intervals: \code{'eti'} for equal-tailed intervals or
\code{'hdi'} for highest density intervals. This option sets the global
options as follows: \code{options("marginaleffects_posterior_interval" =
  "eti")} or \code{options("marginaleffects_posterior_interval" = "hdi")},
and is restored on exit using \code{\link[base:on.exit]{base::on.exit()}}.}

\item{dummy_to_factor}{A named list (default \code{NULL}) to convert dummy
variables into a factor variable. The list must include the following
elements:
\itemize{
\item \code{factor.dummy}: A character vector of dummy variables to be
converted to factors.
\item \code{factor.name}: The name for the newly created factor variable
(default is \code{'factor.var'} if \code{NULL}).
\item \code{factor.level}: A vector specifying the factor levels.
If \code{NULL}, levels are taken from \code{factor.dummy}.
If \code{factor.level} is provided, its length must match
\code{factor.dummy}.
}}

\item{verbose}{A logical argument (default \code{FALSE}) to specify whether
to print information collected during the setup of the object(s).}

\item{expose_function}{A logical argument (default \code{FALSE}) to indicate
whether Stan functions should be exposed. If \code{TRUE}, any Stan
functions exposed during the model fit using \code{expose_function = TRUE}
in the \code{\link[=bsitar]{bsitar()}} function are saved and can be used in post-processing. By
default, \code{expose_function = FALSE} in post-processing functions,
except in \code{\link[=optimize_model]{optimize_model()}} where it is set to \code{NULL}. If
\code{NULL}, the setting is inherited from the original model fit. It must
be set to \code{TRUE} when adding \code{fit criteria} or \code{bayes_R2}
during model optimization.}

\item{usesavedfuns}{A logical value (default \code{NULL}) indicating whether
to use already exposed and saved Stan functions. This is typically set
automatically based on the \code{expose_functions} argument from the
\code{\link[=bsitar]{bsitar()}} call. Manual specification of \code{usesavedfuns} is rarely
needed and is intended for internal testing, as improper use can lead to
unreliable estimates.}

\item{clearenvfuns}{A logical value indicating whether to clear the exposed
Stan functions from the environment (\code{TRUE}) or not (\code{FALSE}). If
\code{NULL}, \code{clearenvfuns} is set based on the value of
\code{usesavedfuns}: \code{TRUE} if \code{usesavedfuns = TRUE}, or
\code{FALSE} if \code{usesavedfuns = FALSE}.}

\item{funlist}{A list (default \code{NULL}) specifying function names. This
is rarely needed, as required functions are typically retrieved
automatically. A use case for \code{funlist} is when \code{sigma_formula},
\code{sigma_formula_gr}, or \code{sigma_formula_gr_str} use an external
function (e.g., \code{poly(age)}). The \code{funlist} should include
function names defined in the \code{globalenv()}. For functions needing
both distance and velocity curves (e.g., \code{plot_curves(..., opt =
  'dv')}), \code{funlist} must include two functions: one for the distance
curve and one for the velocity curve.}

\item{xvar}{A character string (default \code{NULL}) specifying the
\code{'x'} variable. Rarely used because \code{xvar} is inferred
internally. A use case is when conflicting variables exist (e.g.,
\code{sigma_formula}) and user wants to set a specific variable as
\code{'x'}.}

\item{idvar}{A character string (default \code{NULL}) specifying the
\code{'id'} variable. Rarely used because \code{idvar} is inferred
internally.}

\item{itransform}{A character string (default \code{NULL}) indicating the
variables names that are reverse transformed. Options are  \code{c("x",
  "y", "sigma")}. The \code{itransform} is primarily used to get the
\code{xvar} variable at original scale i.e., \code{itransform = 'x'}. To
turn of all transformations, use \code{itransform = ""}. when
\code{itransform = NULL}, the appropriate transformation for \code{xvar} is
selected automatically. Note that when no match for \code{xvar} is found in
the \code{data,frame}, the \code{itransform} will be ignored within the
calling function, \code{'prepare_transformations()'}.}

\item{newdata_fixed}{An indicator to specify whether to check for the data
format and structure of the user provided \code{newdata} and apply needed
\code{prepare_data2} and \code{prepare_transformations}
(\code{newdata_fixed = NULL}, default), return user provided \code{newdata}
(\code{newdata = TRUE}) as it is without checking for the data format or
applying \code{prepare_data2} and \code{prepare_transformations}
(\code{newdata_fixed = 0}), check for the data format and if needed,
prepare data format using \code{prepare_data2} (\code{newdata_fixed = 1}),
or apply \code{prepare_transformations} only assuming that data format is
correct (\code{newdata_fixed = 2}). It is strongly recommended that user
either leave the \code{newdata = NULL} and \code{newdata_fixed = NULL} in
which case data used in the model fitting is automatically retrieved and
checked for the required data format and transformations, and if needed,
\code{prepare_data2} and \code{prepare_transformations} are applied
internally. The other flags provided for  \code{newdata_fixed = 0, 1, 2}
are mainly for the internal use during post-processing.}

\item{envir}{The environment used for function evaluation. The default is
\code{NULL}, which sets the environment to \code{parent.frame()}. Since
most post-processing functions rely on \pkg{brms}, it is recommended to set
\code{envir = globalenv()} or \code{envir = .GlobalEnv}, especially for
derivatives like velocity curves.}

\item{...}{Additional arguments passed to the \code{brms::fitted.brmsfit()}
and \code{brms::predict()} functions.}
}
\value{
A data frame with estimates and confidence intervals for the
specified parameters, or a list when \code{method = 'custom'} and
\code{hypothesis} is not \code{NULL}.
}
\description{
The \strong{marginal_comparisons()} function estimates and
compares growth curves such as distance and velocity. This function is a
wrapper around \code{\link[marginaleffects:comparisons]{marginaleffects::comparisons()}} and
\code{\link[marginaleffects:comparisons]{marginaleffects::avg_comparisons()}}. The \code{\link[marginaleffects:comparisons]{marginaleffects::comparisons()}}
function computes unit-level (conditional) estimates, whereas
\code{\link[marginaleffects:comparisons]{marginaleffects::avg_comparisons()}} returns average (marginal) estimates.
A detailed explanation is available \href{https://marginaleffects.com}{here}.
Note that the \pkg{marginaleffects} package is highly flexible, and users
are expected to have a strong understanding of its workings. Additionally,
since the \pkg{marginaleffects} package is evolving rapidly, results from
the current implementation should be considered experimental.
}
\examples{
\donttest{
# Fit Bayesian SITAR model 

# To avoid mode estimation which takes time, the Bayesian SITAR model fit to 
# the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
# See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.

# Note: Since no covariates are included, the 'marginal_comparisons' 
# function is being shown here as a dummy example. In practice, comparisons  
# may not make sense without relevant covariates. 

# Check and confirm whether model fit object 'berkeley_exfit' exists
berkeley_exfit <- getNsObject(berkeley_exfit)

model <- berkeley_exfit

# Call marginal_comparisons to demonstrate the function
# Note that since model has no covariate, the below example does't make sense

marginal_comparisons(model, draw_ids = 1)
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link[marginaleffects:comparisons]{marginaleffects::comparisons()}},
\code{\link[marginaleffects:comparisons]{marginaleffects::avg_comparisons()}},
\code{\link[marginaleffects:plot_comparisons]{marginaleffects::plot_comparisons()}}
}
\author{
Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
}
