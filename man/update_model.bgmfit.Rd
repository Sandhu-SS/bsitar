% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_model.R
\name{update_model.bgmfit}
\alias{update_model.bgmfit}
\alias{update_model}
\title{Update model}
\usage{
\method{update_model}{bgmfit}(
  model,
  newdata = NULL,
  recompile = NULL,
  expose_function = FALSE,
  verbose = FALSE,
  check_newargs = FALSE,
  envir = NULL,
  ...
)

update_model(model, ...)
}
\arguments{
\item{model}{An object of class \code{bgmfit}.}

\item{newdata}{An optional \code{data.frame} to be used when updating the
model. If \code{NULL} (default), the data used in the original model fit is
re used. Note that data-dependent default priors are not updated
automatically.}

\item{recompile}{A logical to indicate whether the Stan model should be
recompiled. When \code{NULL} (default), \strong{update_model()} tries to
figure out internally whether recompilation is required or not. Setting
\code{recompile} to \code{FALSE} will ignore Stan code changing arguments.}

\item{expose_function}{An optional logical argument to indicate whether to
expose Stan functions (default \code{FALSE}). Note that if user has already
exposed Stan functions during model fit by setting \code{expose_function =
  TRUE} in the \code{\link[=bsitar]{bsitar()}}, then those exposed functions are saved and can be
used during post processing of the posterior draws and therefore
\code{expose_function} is by default set as \code{FALSE} in all post
processing functions except \code{\link[=optimize_model]{optimize_model()}}. For \code{\link[=optimize_model]{optimize_model()}}, the
default setting is \code{expose_function = NULL}. The reason is that each
optimized model has different Stan function and therefore it need to be re
exposed and saved. The \code{expose_function = NULL} implies that the
setting for \code{expose_function} is taken from the original \code{model}
fit. Note that \code{expose_function} must be set to \code{TRUE} when
adding \code{fit criteria} and/or \code{bayes_R2} during model
optimization.}

\item{verbose}{An optional argument (logical, default \code{FALSE}) to
indicate whether to print information collected during setting up the
object(s).}

\item{check_newargs}{A logical (default \code{FALSE}) to check whether
arguments in the original \code{model} fit and the \code{update_model} are
same. When \code{check_newargs = TRUE} and arguments are same, it implies
that update is not needed and hence the original \code{model} object is
returned along with the message if \code{verbose = TRUE}.}

\item{envir}{Environment used for function evaluation. The default is
\code{NULL} which will set \code{parent.frame()} as default environment.
Note that since most of post processing functions are based on \pkg{brms},
the functions needed for evaluation should be in the \code{.GlobalEnv}.
Therefore, it is strongly recommended to set \code{ envir = globalenv()}
(or \code{envir = .GlobalEnv}). This is particularly true for the
derivatives such as velocity curve.}

\item{...}{Other arguments passed to \code{\link[brms:brm]{brms::brm()}}.}
}
\value{
An updated object of class \code{brmsfit}.
}
\description{
The \strong{update_model()} is a wrapper around the
\code{update()} function in the \pkg{brms} package which refits the model
as per the user specified updated arguments.
}
\details{
This is an adapted version of the \strong{update()} function from
available the the\pkg{brms} package.
}
\examples{

\donttest{

# Fit Bayesian SITAR model 

# To avoid mode estimation which takes time, the Bayesian SITAR model fit to 
# the 'berkeley_exdata' has been saved as an example fit ('berkeley_exfit').
# See 'bsitar' function for details on 'berkeley_exdata' and 'berkeley_exfit'.

# Check and confirm whether model fit object 'berkeley_exfit' exists
 berkeley_exfit <- getNsObject(berkeley_exfit)

model <- berkeley_exfit

# Update model
# Note that in case all arguments supplied to the update_model() call are 
# same as the original model fit (checked via check_newargs = TRUE), then  
# original model object is returned.   
# To explicitly get this information whether model is being updated or not, 
# user can set verbose = TRUE. The verbose = TRUE also useful in getting the
# information regarding what all arguments have been changed as compared to
# the original model.

model2 <- update_model(model, df = 5, check_newargs = TRUE, verbose = TRUE)

}

}
\author{
Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
}
