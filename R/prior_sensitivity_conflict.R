

#' Flag parameters with possible prior-data conflict
#'
#' Screen a \code{prior_sensitivity} result for parameters showing notable prior
#' conflict, likelihood conflict, or both. Parameters flagged for both forms of
#' sensitivity may indicate possible prior-data conflict in the \pkg{priorsense}
#' framework, while parameters flagged primarily for prior sensitivity may
#' reflect weak likelihood informativeness.
#'
#' @param x An object returned by [prior_sensitivity()].
#' 
#' @param threshold_prior Numeric threshold used to flag prior sensitivity.
#' 
#' @param threshold_lik Numeric threshold used to flag likelihood sensitivity.
#' 
#' @param ... Ignored.
#'
#' @details
#' This helper is designed as a lightweight screening step after
#' [prior_sensitivity()]. It does not replace graphical inspection. The
#' recommended workflow is:
#' \enumerate{
#'   \item Run \code{prior_sensitivity()} to compute sensitivity.
#'   \item Use \code{prior_sensitivity_conflict()} to identify parameters of
#'         potential concern.
#'   \item Re-run \code{prior_sensitivity()} with \code{plot = "dens"},
#'         \code{plot = "ecdf"}, or \code{plot = "quantities"} for the
#'         flagged variables.
#' }
#'
#' @return A list with components:
#' 
#' \describe{
#'   \item{prior_flagged}{Character vector of variables flagged for
#'     prior sensitivity.}
#'   \item{likelihood_flagged}{Character vector of all variables flagged for
#'     likelihood sensitivity.}
#'   \item{prior_likelihood_flagged}{Character vector of variables flagged for 
#'     both prior and likelihood sensitivity}
#' }
#'
#' @seealso
#' [prior_sensitivity()]
#'
#' @examples
#' \donttest{
#' # Fit Bayesian SITAR model 
#' 
#' # To avoid model estimation which can take time, the Bayesian SITAR model fit
#' # to the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
#' # See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.
#' 
#' model <- getNsObject(berkeley_exfit)
#' 
#' ps <- prior_sensitivity(
#'   model = model,
#'   variable = c("b_a_Intercept", "sigma")
#' )
#' 
#' chk <- prior_sensitivity_conflict(ps)
#' chk$conflict_params
#' chk$prior_only_params
#' }
#'
#' @rdname prior_sensitivity_conflict
#' @export
#' 
prior_sensitivity_conflict.bgmfit <- function(x,
                                              threshold_prior = 0.05,
                                              threshold_lik = 0.05,
                                              ...) {
  
  sens <- x$sensitivity
  
  prior_flagged <- sens %>% 
    dplyr::filter(
      abs(.data[["prior"]]) > threshold_prior
    ) %>% 
    dplyr::pull(.data$variable) %>% 
    unique()
  
  likelihood_flagged <- sens %>% 
    dplyr::filter(
      abs(.data[["likelihood"]]) > threshold_lik
    ) %>% 
    dplyr::pull(.data$variable) %>% 
    unique()
  
  prior_only_params <- setdiff(prior_flagged, likelihood_flagged)
  prior_likelihood_flagged <-  intersect(prior_flagged, likelihood_flagged)
  
  if(is_emptyx(prior_flagged)) prior_flagged <- NULL
  if(is_emptyx(likelihood_flagged)) likelihood_flagged <- NULL
  if(is_emptyx(prior_likelihood_flagged)) prior_likelihood_flagged <- NULL
  
  list(
    prior_flagged = prior_flagged,
    likelihood_flagged = likelihood_flagged,
    prior_likelihood_flagged = prior_likelihood_flagged
  )
}





#' @rdname prior_sensitivity_conflict
#' @export
prior_sensitivity_conflict <- function(x, ...) {
  UseMethod("prior_sensitivity_conflict")
}


