

#'Fit Bayesian SITAR growth curve model
#'
#'@description Fit Bayesian super imposition by translation and rotation (SITAR)
#'  model that summarizes the growth curves from early childhood through the
#'  adulthood (see @details). The frequentist version of the SITAR model can be
#'  fit by an already available R package, the *sitar*
#'  \insertCite{R-sitar}{bsitar}. Besides Bayesian estimation, the *bsitar*
#'  package greatly enhances the modelling capabilities offered earlier by the
#'  *sitar* package. For example, in addition to the univariate model fitting
#'  (i.e, modelling a single outcome as implemented in the *sitar* package), the
#'  **bsitar** allows univariate-by-subgroup and multivariate model
#'  specifications (see @details).
#'
#'@details The SITAR model is a shape-invariant nonlinear mixed effect growth
#'  curve model that fits the population average (i.e., mean average) curve to
#'  the data and then aligns each individual's growth trajectory to the
#'  population average curve by a set of three random effects (size, timing, and
#'  intensity). The concept of shape invariant model (SIM) was first described
#'  by \insertCite{Lindstrom1995}{bsitar} and later used by
#'  \insertCite{Beath2007;textual}{bsitar} to model infant growth data
#'  (birth to 2 years). The current version of the SITAR model is developed by
#'  \insertCite{Cole2010;textual}{bsitar} and has been used extensively for
#'  modelling human growth data \insertCite{@see
#'  @nembidzaneUsingSITARMethod2020; @mansukoskiLifeCourseAssociations2019;
#'  @coleFiftyYearsChild2018; @riddellClassifyingGestationalWeight2017;
#'  @Sandhu2020}{bsitar}. As mentioned earlier (see @description), the 
#'  frequentist version of the SITAR model can be fit by an already available 
#'  R package, the *sitar* \insertCite{R-sitar}{bsitar}.
#'
#'  The SITAR model specification is same in *sitar* and **bsitar** with the
#'  exception that unlike *sitar* which uses the B spline basis for the natural
#'  cubic spline design matrix (by calling the *ns* function of the *splines*
#'  package \insertCite{R-splines}{bsitar}), the *bsitar* constructs spline
#'  design matrix by using the truncated power basis approach as
#'  described by \insertCite{harrell2001regression}{bsitar} and implemented in
#'  the *rcspline.eval* function of the *Hmisc* package
#'  \insertCite{R-Hmisc}{bsitar}. Note that the **bsitar** package does not use
#'  the *rcspline.eval* but rather constructs a custom function on the fly that
#'  is included in the functions block of the *Stan* programs' and
#'  compiled (via the c++) during the model estimation.
#'
#'  Like *sitar* package, the **bsitar** fits SITAR model with (usually) up
#'  to three random effect parameters \insertCite{Cole2010}{bsitar}: the size
#'  (\code{a}), the timing (\code{b}) and the intensity (\code{c}). In addition,
#'  there is a slope parameter \code{d} that models the variability in the adult
#'  slope of the growth curve (See [sitar::sitar] for details). Please note that
#'  inclusion of \code{d} results in multicollinearity because the specification
#'  of the this \code{d} parameter involves a linear predictor term which is
#'  identical to the first term of the spline design matrix created by using 
#'  the truncated power basis approach.
#'
#'  The *bsitar* function is the main workhorse of the **bsitar** package that
#'  fits the Bayesian SITAR model. The package is a frontend to the R package
#'  *brms* \insertCite{@see @R-brms; @brms2021}{bsitar}  which can fit a wide
#'  range of hierarchical linear and nonlinear regression models including
#'  multivariate models. The *brms* itself depends on the 'Stan' software for 
#'  full Bayesian inference \insertCite{@see @teamStanReferenceManual;
#'  @gelman2015}{bsitar}. The **bsitar** package allows a wide
#'  range of prior specifications that encourage the users to apply prior
#'  distributions that actually reflect their prior knowledge about the human
#'  growth processes such as the timing of the age at peak growth velocity. The
#'  model fit to the data can evaluated by means of posterior predictive check
#'  (see [brms::pp_check()]). Furthermore, models with different priors and/or
#'  growth curves (i.e., with different \code{df} for splines) can be easily
#'  compared by using methods available in the *brms* package such as 
#'  the leave one out cross validation (see [brms::loo()]). The excellent
#'  post-processing support offered by the *brms* is further augmented by custom
#'  functions written for the **bsitar** that allows prediction and
#'  visualization of population average and individual specific growth
#'  trajectories velocity curves. Furthermore, population average and individual
#'  specific growth parameters such as age at peak growth velocity (APGV) and
#'  the peak growth velocity (PGV) can be easily computed.
#'
#'  The *bsitar* package allows three different model specifications:
#'  univariate, univariate-by-subgroup model, and multivariate. The
#'  univariate-by-subgroup approach fits two or more separate submodels for a
#'  single outcome defined by a factor variable (e.g, sex). The data are
#'  typically stacked and the factor variable is used to set-up the submodels by
#'  using the 'subset' option available in the [brms::brm()] function. The
#'  multivariate model specification allows simultaneous modelling of two or
#'  more outcomes with joint a distribution of random effects. For both
#'  univariate-by-subgroup and multivariate model fitting, the **bsitar**
#'  package allows full flexibility in specifying separate predictor (\code{x}),
#'  subject identifiers (\code{id}), degree of freedom (\code{df}) / knots
#'  (\code{knots}) as well as the prior and initial values for each submodel.
#'  Furthermore, to enhance the ease of specifying different options and to make
#'  it user-friendly, there is no need to enclose the character strings in
#'  single or double quotes. For example to specify the univariate-by-subgroup
#'  model for sex, the \code{univariate_by = sex} is same as \code{univariate_by
#'  = 'sex'} or \code{univariate_by = "sex"}. The same applies for all character
#'  string options.
#'
#'
#'@param x Predictor variable(s) in the data (typically age in years). For
#'  univariate model, the \code{x} is a single variable whereas for the
#'  univariate-by-subgroup (see \code{univariate_by}) and multivariate (see
#'  \code{multivariate}) model specifications, the \code{x} can be same for each
#'  submodel or else specified separately for each submodel. For example, for a
#'  bivariate model, the \code{x = list(x1, x2)} specifies that \code{x1} is the
#'  predictor for the first submodel and \code{x2} for the second sub model.
#'
#'@param y Response variable(s) in the data (i.e., repeated height measurements). For
#'  univariate and univariate-by-subgroup (see \code{univariate_by} argument)
#'  model specifications, \code{y} is specified as a single variable. For the
#'  univariate-by-subgroup model fitting, the outcome vectors for each sub model
#'  are created internally and named using the the factor levels. For example
#'  when fitting a univariate-by-subgroup model for sex (specified by
#'  \code{univariate_by = list(by = sex)} or simply as \code{univariate_by =
#'  sex}), the outcome vectors 'Female' and 'Male' are created automatically
#'  where 'Female' is the first level of the factor variable sex, and 'Male' is
#'  second level. For multivariate model (\code{multivariate} argument), the
#'  outcome vectors are specified as a list (e.g., \code{y = list(y1, y2}) where
#'  \code{y1} and \code{y2} are the outcomes.
#'
#'@param id A vector specifying the groups (typically individuals) assignment.
#'  For univariate-by-subgroup (see \code{univariate_by} argument) and
#'  multivariate (see \code{multivariate} argument) model specifications, the
#'  \code{id} can be same (typically) or else specified separately as a \code{id
#'  = list(id1, id2)} where \code{id1} and \code{id2} are subject identifiers.
#'
#'@param data data frame containing variables \code{x}, \code{y} and \code{id}.
#'
#'@param df degrees of freedom for natural cubic regression spline. For
#'  univariate-by-subgroup model (specified by using the \code{univariate_by}
#'  argument) and multivariate model (specified by using the \code{multivariate}
#'  argument), the \code{df} could be same \code{df = 4} or else specified
#'  separately as \code{df = list(4, 5)} where df = 4 is for the first outcome
#'  and df = 5 is for the second outcome.
#'
#'@param knots vector of values for knots (default \code{df} quantiles of
#'  \code{x} distribution). See \code{df} for specifying separate knots for
#'  univariate-by-subgroup and multivariate models.
#'
#'@param fixed character specifying a, b, c, d fixed effects. Typically
#'  specified as \code{fixed = a+b+c}. As mentioned earlier, there is no need to
#'  enclose character in quotes. I other words, \code{fixed = a+b+c},
#'  \code{fixed = 'a+b+c'}, and \code{fixed = "a+b+c"} are same. For specifying
#'  different fixed effect structures for univariate-by-subgroup and
#'  multivariate models, use list as follows: \code{fixed = list(a+b+c, a+b)}
#'  which implies that the fixed effect structure for the first outcome is
#'  \code{fixed = 'a+b+c'} and \code{fixed = 'a+b'} for the second outcome.
#'
#'@param random character specifying a, b, c, d random effects. See \code{fixed}
#'  for setting the random effects structure.
#'
#'@param xoffset optional value of offset for \code{x} allowing the origin of
#'  \code{x} to be varied (either mean, apv, or a numeric value). The default is
#'  mean.
#'
#'@param bstart optional value to set initial value for fixed effect \code{b}.
#'  Options are mean (default), apv, or a numeric value.
#'
#'@param xfun an optional argument to transform the predictor (i.e., age).
#'  Options are log and sqrt for the logarithmic and square root transformation.
#'  The default is NULL implying that no transformation is applied and the model
#'  is fit to the original scale of the predictor (e.g., years). Like other
#'  arguments, user can specify different xfun for univariate-by-subgroup
#'  (specified by using the \code{univariate_by} argument) and multivariate
#'  (specified by using the \code{multivariate} argument) models as a list i.e.,
#'  \code{xfun = list(log, sqrt)} or \code{xfun = list(NULL, sqrt)}.
#'
#'@param yfun an optional argument to transform the outcome (i.e., age). Options
#'  are log, sqrt or NULL (default). See \code{xfun} for details.
#'
#'@param bound span of \code{x} for regression spline, or a small extension of
#'  range (default 0.04). See package 'sitar' for details.
#'
#'@param terms_rhs An option argument (default \code{NULL}) to specify terms on
#'  the left hand side of the formula and after the outcome (separated by '|')
#'  that are used when fitting measurement error model. An an example, consider
#'  fitting a model with measurement error in the outcome specified as
#'  \code{bf(y | mi(sdy) ~ ..)} where \code{mi(sdy)} is passed as follows:
#'  \code{terms_rhs = mi(sdy)}. For multivariate models, each outcome can have
#'  its own measurement error varibales that are passed as a list, i.e.,
#'  \code{terms_rhs = list(mi(sdy1), mi(sdy2))}.
#'
#'@param a_formula formula for fixed effect a (default \code{~ 1}). User can
#'  specify different formula when fitting univariate-by-subgroup (specified by
#'  using the \code{univariate_by} argument) and the multivariate (specified by
#'  using the \code{multivariate} argument) models. As an example
#'  \code{a_formula = list(~1, ~1 + cov)} implies that the \code{a_formula} for
#'  the first outcome includes only an intercept whereas the \code{a_formula}
#'  for the second outcome includes an intercept plus a covariate. The covariate
#'  can be a continous variable or a factor variable (dummy variables will be
#'  created using the \code{model.matrix}). The formula can include a
#'  combination of continous and factor variables as covariate as well as their
#'  interactions.
#'
#'@param b_formula formula for fixed effect b (default \code{~ 1}). See
#'  \code{a_formula} for details.
#'
#'@param c_formula formula for fixed effect c (default \code{~ 1}). See
#'  \code{a_formula} for details.
#'
#'@param d_formula formula for fixed effect d (default \code{~ 1}). See
#'  \code{a_formula} for details.
#'
#'@param s_formula formula for fixed effect b (default \code{~ 1}). See
#'  \code{a_formula} for details.
#'
#'@param a_formula_gr formula for random effect a (default \code{~ 1}). See
#'  \code{a_formula} for details. The random effect can be specified by using
#'  the bar approach or allowing \code{group_by} argument to set group
#'  identifier and the correlation structure. For example, to specify
#'  unstructured varinace covarinace structure for random effects parameters
#'  \code{a}, \code{b}, and \code{c}, the formuale random parameter \code{a},
#'  \code{b}, and \code{c} are specified as \code{a_formula_gr = ~ 1},
#'  \code{b_formula_gr = ~ 1}, and \code{c_formula_gr = ~ 1}  then using the
#'  \code{group_by} argument to specify the group identifier and the to set the
#'  correlation structure as follows \code{group_by = list(groupvar = id, cor =
#'  un)}. Else, the bar approach can be used equivalently to set unstructured
#'  varinace covarinace structure \code{a_formula_gr = ~ (1 |c|id)},
#'  \code{b_formula_gr = ~ (1 |c|id)}, and \code{c_formula_gr = ~ (1 |c|id)}.
#'
#'@param b_formula_gr formula for random effect b (default \code{~ 1}). See
#'  \code{a_formula_gr} for details.
#'
#'@param c_formula_gr formula for random effect c (default \code{~ 1}). See
#'  \code{a_formula_gr} for details.
#'
#'@param d_formula_gr formula for random effect d (default \code{~ 1}). See
#'  \code{a_formula_gr} for details.
#'
#'@param a_formula_gr_str formula for random effect a (default \code{NULL}) when
#'  fitting model with hierarchical structure greater than two levels. See
#'  \code{a_formula} for details. The groupvar and priors specified by the
#'  argument \code{a_formula_gr} are used for the second level of hierarchy
#'  whereas groupvar and priors are manually specified for the third and beyond
#'  hierarchies by using the \code{get_priors} and \code{set_self_priors}
#'  arguments. An example of specifying formula for random effect parameter
#'  \code{a} for a three level model with repeated measures on individuals
#'  nested within the growth studies is as follows: \code{a_formula_gr_str = ~
#'  (1|a|id:study) + (1|b|istudy)}. This formulation implies (with \code{|a|}
#'  and \code{|b|}) that fully unstructured varinace covarinace structure is
#'  specified for indivuals levels as well as for study. Note that \code{|a|}
#'  and \code{|b|} need to be distict as becayse parameters can only only for
#'  correlated for at the same level of hierarchy. in other words, it is not
#'  alowed to use the letters such as \code{|a|} and \code{|a|}.
#'
#'@param b_formula_gr_str formula for random effect b (default \code{NULL}) when
#'  fitting model with hierarchical structure greater than two levels. See
#'  \code{a_formula_gr_str} for details.
#'
#'@param c_formula_gr_str formula for random effect c (default \code{NULL}) when
#'  fitting model with hierarchical structure greater than two levels. See
#'  \code{a_formula_gr_str} for details.
#'
#'@param d_formula_gr_str formula for random effect d (default \code{NULL}) when
#'  fitting model with hierarchical structure greater than two levels. See
#'  \code{a_formula_gr_str} for details.
#'  
#'@param sigma_formula formula for modelling distributional parameter sigma.
#'  (default \code{NULL}). This is only useful when including covariates(s) for
#'  sigma. The [brms::brm()] by defaults function includes an intercept for the
#'  residual standard deviation parameter (i.e,, sigma). The
#'  \code{sigma_formula} along with \code{sigma_formula_gr} and
#'  \code{sigma_formula_gr_str} arguments allows specifying hierarchical
#'  structure when modelling sigma. This set is similar to setting fixed and
#'  random effect structures for parameters \code{a}, \code{b}, and \code{c}.
#'  The \code{sigma_formula} sets up the fixed effect design matrix. It is
#'  important to note that another alternative to set up the fixed effect design
#'  matrix for distributional parameter sigma is argument \code{dpar_formula}.
#'  An advantage of \code{dpar_formula} over \code{sigma_formula} is that user
#'  can specify the linear and nonlinear formulation by using the brms'
#'  [brms::lf] and [brms::nlf] syntax. Both [brms::lf] and [brms::nlf] further
#'  allows control over centering of predictors as well as to enable or disable
#'  cell mean centering when excluding the \code{intercept} by adding \code{0}
#'  to the right-hand of model formulas. Note that \code{sigma_formula} and
#'  \code{dpar_formula} can not be specified together.
#'
#'@param sigma_formula_gr formula for setting up the random effect structure for
#'  the distributional parameter sigma (default \code{NULL}).  Note that
#'  \code{sigma_formula} can not be combined with the \code{dpar_formula}.
#'
#'@param sigma_formula_gr_str formula for setting up the random effect structure
#'  for the sigma (default \code{NULL}) when fitting model with hierarchical
#'  structure greater than two levels. The groupvar and priors specified by the
#'  argument \code{sigma_formula_gr} are used for the second level of hierarchy
#'  whereas groupvar and priors are manually specified for the third and beyond
#'  hierarchies by using the \code{get_priors} and \code{set_self_priors}
#'  arguments. An example of specifying formula for random effect parameter
#'  \code{sigma} for a three level model with repeated measures on individuals
#'  nested within the growth studies is as follows: \code{sigma_formula_gr_str = ~
#'  (1|a|id:study) + (1|b|istudy)}. This formulation implies 
#'  (with vertical bar, \code{|} ) that fully unstructured varinace covarinace 
#'  structure is specified for indivuals levels as well as for study. 
#'  See \code{a_formula_gr_str} for further details.
#'
#'@param dpar_formula formula for distributional parameter sigma (default
#'  \code{NULL}). This is only useful when modelling the sigma (i.e., residual
#'  standard deviation parameter). By default, the [brms::brm()] function includes
#'  the intercept for the residual standard deviation parameter (i.e,, sigma).
#'  The default setting for \code{dpar_formula} is NULL which implements the
#'  default behaviour of the [brms::brm()] function. Also note that
#'  \code{dpar_formula} can not be specified along with \code{sigma_formula},
#'  \code{sigma_formula_gr}, or \code{sigma_formula_gr_str}. See
#'  \code{sigma_formula} for relative advantages and disadvantages of using
#'  \code{sigma_formula} and \code{dpar_formula}
#'
#'@param autocor_formula formula for modelling autocorrelation. The default
#'  setting for \code{autocor_formula} is NULL i.e, not to model
#'  autocorrelation. Allowed options are autoregressive moving average (ARMA) of
#'  order (p, q), autoregressive (AR) of order (p) and moving average (MA) of
#'  order (q) which are specified as \code{autocor_formula = arms(p=1, q=1)},
#'  \code{autocor_formula = ar(p=1)}, and \code{autocor_formula = msq=1)}. See
#'  brms package for further details on order p and q and setting up
#'  autocorrelation structures.
#'
#'@param family response distribution and link function to be used in the model.
#'  The default is gaussian(). See [brms::brm()] function for details. For
#'  univariate-by-subgroup model (specified by using the \code{univariate_by}
#'  argument) and multivariate model (specified by using the \code{multivariate}
#'  argument), the \code{family} could be same \code{family = gaussian()} or
#'  else different such as \code{family = list(gaussian(), student()} which sets
#'  gaussian for the first outcome and student_t for the second outcome.
#'
#'@param group_arg specify group-level effects when fitting univariate models.
#'  The subptions for the \code{group_arg} are 'groupvar', 'dist', 'cor' and 'by.' The
#'  suboption groupvar specifies the subject identifier (which is typically same
#'  as \code{id}) wheres the suboption dist sets the distribution of the random
#'  effects (options are gaussian, the default and student). The correlation
#'  structure is specified by using the suboption cor. The correlation structure
#'  allowed are unstructured (default) and the diagonal. The unstructured
#'  correlation structure (cor = un) models the full varinace covarinace
#'  structure whereas the diagonal correlation structure (cor = diagonal)
#'  estimates only the variance (i.e, standard deviation) parameters
#'  (correlation parameters are set to zero). For further details, see
#'  [brms::brm()] function (\bold{Group-level terms}). Note that only the groupvar
#'  suboption of the \code{group_arg} is passed to the univariate-by-subgroup
#'  \code{univariate_by} and the multivariate (specified by using the
#'  \code{multivariate} model fittings. Lastly, the \code{group_arg} is
#'  completely ignored when user specify random formula by using the "|"
#'  option for \code{a_formula_gr}, \code{b_formula_gr}, and
#'  \code{c_formula_gr}. Also, the \code{group_arg} is redundant for
#'  \code{a_formula_gr_str}, \code{b_formula_gr_str},
#'  and \code{c_formula_gr_str}.
#'  
#'@param sigma_group_arg specify group-level effects for distributional parameter \code{sigma}.
#'  The subptions for the \code{sigma_group_arg} sames as are same \code{group_arg} i.e., 'groupvar', 'dist', 'cor' and 'by.' The
#'  suboption 'groupvar' specifies the subject identifier (which is typically same
#'  as \code{id}) whereas the suboption 'dist' sets the distribution of the random
#'  effects (options are gaussian, the default and student). The correlation
#'  structure is specified by using the suboption 'cor'. See \code{group_arg}
#'  for more details.
#'
#'@param univariate_by specify univariate-by-subgroup model fitting arguments.
#'  Suboptions include the by argument to specify the variable (which must be a
#'  factor variable) and the cor suboption to specify the correlation structure.
#'  The unstructured correlation structure (cor = un) models the full varinace
#'  covarinace structure separately for each submodel (i.e., outcome) whereas
#'  the diagonal correlation structure estimates only the variance (i.e,
#'  standard deviation) for each submodel (i.e., all outcomes). The default
#'  settings for by and cor argument are NULL and un, respectively.
#'
#'@param multivariate specify multivariate model fitting arguments. Suboptions
#'  include the mvar argument (logical, default FALSE) to specify whether or not
#'  to fit a multivariate model, cor suboption to specify the correlation
#'  structure, and rescor option (logical, default TRUE) to specify whether or
#'  not to estimate the residual correlation parameter for the outcomes. The The
#'  unstructured correlation structure (cor = un) jointly estimates all the
#'  varinace covarinace parameters between the outcomes whereas diagonal
#'  correlation structure (cor = diagonal) estimates only the variance
#'  parameters for each outcomes. Another option (cor = un_s) allows for
#'  estimating unstructured correlation structure separately for each outcome.
#'  The default setting for correlation is un.
#'
#'@param a_prior_beta Set priors on the the fixed effect a parameter. The
#'  allowed distributions are normal, student_t, cauchy, lognormal, uniform,
#'  exponential, gamma, inverse gamma. See [brms::prior()] function for details on
#'  priors. For each distibution, suboption allows for setting upper and lower
#'  bounds (default NA, i.e., lb = NA, ub = NA). For location scale based
#'  distributions which include the normal, student_t, cauchy and lognormal
#'  distributions, option autosclae (default FALSE) is provided to multiply the
#'  scale parameter as typically done in the rstanarm package. The rstanarm sets
#'  the autosclae as 2.5 whereas the brms package sets its to 1 or 2.5 depending
#'  on the standard deviation of the outcome (See [brms::prior()] function for
#'  details). The 'bsitar' package offers the flexibility of choosing the value
#'  for the autosclae. For convinience purposes, lower bound as zero is
#'  automatically set for the positive distributions (such as lognormal,
#'  exponential, gamma). For uniform distrubution, another option addrange is
#'  provided to symmetrically expand the lower and upper parameters specified as
#'  the uniforma prior. For example, uniform(a, b, addrange = 5) would be
#'  evaluated uniform(a-5, b+5). For exponential distribution, the rate
#'  parameter is evaluated as inverse. In other words, prior set as
#'  exponential(10)  would be translated to exponential(1/10). Examples of
#'  setting priors are shown below:
#' \code{a_prior_beta = normal(location = 5, scale = 1, lb = NA, ub = NA,
#' addrange = NA, autosclae = FALSE)}). This is same as as
#'  \code{a_prior_beta = normal(5, 1)})
#' \code{a_prior_beta = student_t(df = 3, location = 5, scale = 1, lb = NA,
#' ub = NA, addrange = NA, autosclae = FALSE)}) which can specified using
#'  shorthand as \code{a_prior_beta = student_t(3, 5, 1)}). For location scale
#'  based distributions, user can use specify the mean or the median of the
#'  outcome as location and the standard deviation (sd) or the median absolute
#'  deviation (mad) as scale as shown (any combination of these)
#'  \code{a_prior_beta = normal(ymean, ysd)}) and \code{a_prior_beta =
#'  normal(ymedian, ymad)}) Another option available for setting the location
#'  parameter is the linear model fit  based intercept (i.e., the lm model
#'  fitted to the outcome) as \code{a_prior_beta = normal(lm, ysd)}). For
#'  univariate-by-subgroup model (specified by using the \code{univariate_by}
#'  argument) and multivariate model (specified by using the \code{multivariate}
#'  argument), priors specified for each outcome can be same specified as a
#'  single option i.e., \code{a_prior_beta = normal(5, 1)} or can be specified
#'  different by using list as shown here \code{a_prior_beta = list(normal(5,
#'  1), normal(10, 5)} which
#'
#'@param b_prior_beta set priors on the the fixed effect b parameter. Specifying
#'  prior for the fixed effect b parameter is same as described above for the
#'  fixed effect a parameter \code{a_prior_beta} except for the fact that
#'  options ymean, ymedian, ysd and ymad are not allowed. Also, using lm as
#'  location parameter sets location as 0. (same behavior as the 'sitar'
#'  package)
#'
#'@param c_prior_beta set priors on the the fixed effect c parameter. Specifying
#'  prior for the fixed effect c parameter is exactly same as described above
#'  for the fixed effect b parameter (\code{b_prior_beta})
#'
#'@param d_prior_beta set priors on the the fixed effect c parameter. Specifying
#'  prior for the fixed effect d parameter is exactly same as described above
#'  for the fixed effect c parameter (\code{c_prior_beta})
#'
#'@param s_prior_beta set priors on the the fixed effect s parameter (i.e., the
#'  spline coefficients). The general approach to set priors for the s parameter
#'  is same as described earlier for the fixed effect a parameter. The allowed
#'  option for location and scale are lm which would set location parameter for
#'  spline coefficient based on the spline coefficients from the linear model
#'  fit to the data. The lm option for the scale parameter sets the standard
#'  deviation of spline design matrix used to fiot the linear model. For s
#'  parameter, it make sense to use only location scale based priors (i.e,
#'  normal, student_t and cauchy) or uniform priors. For uniform priors, the
#'  addrange  option can be utilized to symettrically add range to the lm based
#'  spline coefiecnt). An additional option available for the location scale
#'  based priors is sethp (logical, default set as FALSE) which, when set as
#'  TRUE, allows for setting hierarchical priors for the s parameter. In other
#'  words, instead of setting prior as s ~ normal(0, lm), the hierarchical
#'  priors are set as s ~ normal(0, hp) where hp ~ normal(0, lm). Note that the
#'  scale parameter for the hp ~ normal(0, lm) is automatically taken from the s
#'  ~ normal(0, hp). Setting sethp = TRUE impllies that the scale for spline
#'  coeficients is estimated from the data itself. The distribution of
#'  hierarchical priors is automatically matched with the prior set for the s
#'  parameter or else can be
#' set by the same sethp option. For example, \code{s_prior_beta =
#' normal(0, lm, sethp = caucy)} will be translated to s ~ normal(0, lm);
#'  hp ~ caucy(0, lm).
#'
#'@param a_cov_prior_beta set priors on the covariates for the fixed effect a
#'  parameter. The approach is same as described for the \code{a_prior_beta}
#'  except that the options ymean, ymedian, ysd and ymad are not allowed. Option
#'  lm for location parameter for normal, student_t, cauchy prior would set
#'  location based on Intercept obtained from the lm model fit. Separate priors
#'  can be specified for outcomes when fitting univariate-by-subgroup (specified
#'  by using the \code{univariate_by} argument) and the multivariate (specified
#'  by using the \code{multivariate} argument) models (see \code{a_prior_beta}).
#'
#'@param b_cov_prior_beta set priors on the covariates for the fixed effect b
#'  parameter. The approach is same as described for the \code{a_prior_beta}
#'  except that the options ymean, ymedian, ysd and ymad are not allowed. Option
#'  lm for location parameter for normal, student_t, cauchy prior would set
#'  location as 0. Separate priors can be specified for outcomes when fitting
#'  univariate-by-subgroup (specified by using the \code{univariate_by}) and the
#'  multivariate (specified by using the \code{multivariate}) (see
#'  \code{b_prior_beta}).
#'
#'@param c_cov_prior_beta set priors on the covariates for the fixed effect c
#'  parameter. The approach is same as described for the \code{a_prior_beta}
#'  except that the options ymean, ymedian, ysd and ymad are not allowed. Option
#'  lm for location parameter for normal, student_t, cauchy prior would set
#'  location as 0. Separate priors can be specified for outcomes when fitting
#'  univariate-by-subgroup (specified by using the \code{univariate_by}) and the
#'  multivariate (specified by using the \code{multivariate}) (see
#'  \code{c_prior_beta}).
#'
#'@param d_cov_prior_beta set priors on the covariates for the fixed effect d
#'  parameter. The approach is same as described for the \code{a_prior_beta}
#'  except that the options ymean, ymedian, ysd and ymad are not allowed. Option
#'  lm for location parameter for normal, student_t, cauchy prior would set
#'  location as 0. Separate priors can be specified for outcomes when fitting
#'  univariate-by-subgroup (specified by using the \code{univariate_by}) and the
#'  multivariate (specified by using the \code{multivariate}) (see
#'  \code{d_prior_beta}).
#'
#'@param s_cov_prior_beta set priors on the covariates for the fixed effect s
#'  parameter. The approach is same as described for the \code{s_prior_beta}
#'  except that the options ymean, ymedian, ysd and ymad are not allowed. Option
#'  lm for location parameter for normal, student_t, cauchy prior would set
#'  location based on Intercept obtained from the lm model fit. Separate priors
#'  can be specified for outcomes when fitting univariate-by-subgroup (specified
#'  by using the \code{univariate_by} argument) and the multivariate (specified
#'  by using the \code{multivariate} argument) models (see \code{s_prior_beta}).
#'
#'@param a_prior_sd set priors on the the random effect a parameter. The allowed
#'  distributions are normal, student_t, cauchy, lognormal, uniform,
#'  exponential, gamma, inverse gamma. For location scale based distributions
#'  which include the normal, student_t, cauchy and lognormal distributions,
#'  option autosclae (default FALSE) is provided to multiply the scale
#'  parameter. For location scale based distributions, user can use specify the
#'  standard deviation (sd) or the median absolute deviation (mad) as scale
#'  parameter. For univariate-by-subgroup model (specified by using the
#'  \code{univariate_by} argument) and multivariate model (specified by using
#'  the \code{multivariate} argument), priors specified for each outcome can be
#'  same or else can be different for each outcome (see \code{a_prior_beta} for
#'  details). The lower bound as zero is automatically set by the
#'  \code{brms::brm}.
#'
#'@param b_prior_sd set priors on the the random effect b parameter. Specifying
#'  prior for the random effect b parameter is same as described above for the
#'  random effect a parameter \code{a_prior_sd} except for the fact that
#'  optionsysd and ymad are not allowed.
#'
#'@param c_prior_sd set priors on the the random effect c parameter. Specifying
#'  prior for the random effect c parameter is same as described above for the
#'  random effect b parameter \code{b_prior_sd}.
#'
#'@param d_prior_sd set priors on the the random effect d parameter. Specifying
#'  prior for the random effect d parameter is same as described above for the
#'  random effect c parameter \code{b_prior_sd}.
#'
#'@param a_cov_prior_sd set priors on the the covariate(s) for random effect a
#'  parameter. Approach is same as described earlier for the
#'  \code{a_cov_prior_beta}. No pre-defined option (e.g., lm) is allowed to set
#'  the scale for the location scale based priors.
#'
#'@param b_cov_prior_sd set priors on the the covariate(s) for random effect b
#'  parameter. Approach is same as described earlier for the
#'  \code{a_cov_prior_beta}.
#'
#'@param c_cov_prior_sd set priors on the the covariate(s) for random effect c
#'  parameter. Approach is same as described earlier for the
#'  \code{a_cov_prior_sd}.
#'
#'@param d_cov_prior_sd set priors on the the covariate(s) for random effect d
#'  parameter. Approach is same as described earlier for the
#'  \code{a_cov_prior_sd}.
#'  
#'@param sigma_prior_beta Set priors on the the fixed effect distributional
#'  parameter \code{sigma}. The allowed distributions are normal, student_t,
#'  cauchy, lognormal, uniform, exponential, gamma, inverse gamma. See
#'  [brms::prior()] function for details on priors. See \code{a_prior_beta} for
#'  full details on how to specify priors including various subptional
#'  available.
#'
#'@param sigma_cov_prior_beta Set priors on the covariate(s) for the the fixed
#'  effect distributional parameter \code{sigma}. The allowed distributions are
#'  normal, student_t, cauchy, lognormal, uniform, exponential, gamma, inverse
#'  gamma. See [brms::prior()] function for details on priors. See
#'  \code{a_cov_prior_beta} for full details on how to specify priors including
#'  various suboptions available.
#'
#'@param sigma_prior_sd Set priors on the the random effect distributional
#'  parameter \code{sigma}. The allowed distributions are normal, student_t,
#'  cauchy, lognormal, uniform, exponential, gamma, inverse gamma. See
#'  [brms::prior()] function for details on priors. See \code{a_cov_prior_beta}
#'  for full details on how to specify priors including various subptional
#'  available.
#'
#'@param sigma_cov_prior_sd Set priors on the covariate(s) for the the random
#'  effect distributional parameter \code{sigma}. The allowed distributions are
#'  normal, student_t, cauchy, lognormal, uniform, exponential, gamma, inverse
#'  gamma. See [brms::prior()] function for details on priors. See
#'  \code{a_cov_prior_sd} for full details on how to specify priors including
#'  various suboptions available.
#'
#'@param rsd_prior_sigma set priors on the the residual standard deviation
#'  parameter sigma. This argument will only be evaluated if \code{dpar_formual}
#'  is set to NULL. For location scale based distributions, user can use specify
#'  the standard deviation (sd) or the median absolute deviation (mad) as scale
#'  parameter.
#'
#'@param dpar_prior_sigma set priors on the the distributional parameter which
#'  is sigma for location scale based family (Gaussian and student). This
#'  argument is evaluated only when \code{dpar_formual} is not set to NULL. For
#'  location scale based distributions, user can use specify the standard
#'  deviation (sd) or the median absolute deviation (mad) as scale parameter.
#'
#'@param dpar_cov_prior_sigma set priors on the the covariate(s) for the
#'  distributional parameter (i.e., sigma). The approach is same as descibed
#'  above for the \code{dpar_prior_sigma} except that options standard deviation
#'  (sd) and the median absolute deviation (mad) are not allowed to set the
#'  scale parameter for the location scale based distributions.
#'
#'@param autocor_prior_acor set priors on the the autocorrelation parameters
#'  (i.e., ar and ma parameters, see \code{autocor_formula} for details). The
#'  only allowed distribution is uniform distribution bounded between -1 and +
#'  1.
#'
#'@param gr_prior_cor set priors on the the correlations of group-level
#'  ('random') effects. The allowed distribution is lkj which has a single
#'  parameter eta to control the priors on correlation parameters (see
#'  \code{brms::prior} for details).
#'
#'@param mvr_prior_rescor set priors on the the residual correlations for
#'  multivariate model. The allowed distribution is lkj which has a single
#'  parameter eta to control the priors on correlation parameters (see
#'  \code{brms::prior} for details).
#'
#'@param init Initial values for the sampler. For \code{0}, all parameters are
#'  initialized to zero. If \code{random}, Stan will randomly generate initial
#'  values for parameters in a range specified by the \code{init_r} (see below).
#'  Another option is \code{prior} which allows setting initials based on the
#'  periors specified above. Lastly, for \code{NULL} (the default) initial
#'  value(s) for each parameter is set by the following initial argument
#'  specific for each model parameter.
#'
#'@param init_r Range for the random generatetion of initial values when
#'  \code{init} is set as \code{"random"} (see above).
#'
#'@param a_init_beta Initial values for fixed effect a parameter. Options are
#'  \code{0}, \code{random} and \code{prior}. Also, similar to the location
#'  parameter for \code{a_prior_beta}, user can specify ymean and ymedian to set
#'  initial as mean or median of the outcome. Furthermore, option \code{lm}
#'  would set initial based on the intercept of the lm model fit to the outome.
#'  Lastly, For univariate-by-subgroup model (specified by using the
#'  \code{univariate_by} argument) and multivariate model (specified by using
#'  the \code{multivariate} argument), the \code{a_init_beta} could be same
#'  \code{a_init_beta = 0} for each outcome or can be specified differently as
#'  \code{list(a_init_beta = 0, a_init_beta = lm)}.
#'
#'@param b_init_beta Initial values for fixed effect b parameter. The approach
#'  and options are same as described above for the \code{a_init_beta}. However,
#'  the option \code{lm} will set the initial as 0. An extra option for b
#'  parameter is \code{bstart} (see \code{bstart} argument for details).
#'
#'@param c_init_beta Initial values for fixed effect c parameter. The approach
#'  and options are same as described above for the \code{b_init_beta}.
#'
#'@param d_init_beta Initial values for fixed effect a parameter. The approach
#'  and options are same as described above for the \code{c_init_beta}.
#'
#'@param s_init_beta Initial values for fixed effect s parameter. The approach
#'  and options are same as described above for the \code{a_init_beta}. The
#'  option \code{lm} will set the initial based on the spline coeficients
#'  obtained from the lm model fit.
#'
#'@param a_cov_init_beta Initial values for covariate(s) for the fixed effect a
#'  parameter. The approach and options are same as described above for the
#'  \code{a_init_beta}. The option \code{lm} will set the initial based on the
#'  covariate intercepts obtained from the lm model fit.
#'
#'@param b_cov_init_beta Initial values for covariate(s) for the fixed effect b
#'  parameter. The approach and options are same as described above for the
#'  \code{a_cov_init_beta}. The option \code{lm} will set the initial as 0 for
#'  the covariate parameters.
#'
#'@param c_cov_init_beta Initial values for covariate(s) for the fixed effect c
#'  parameter. The approach and options are same as described above for the
#'  \code{b_cov_init_beta}.
#'
#'@param d_cov_init_beta Initial values for covariate(s) for the fixed effect c
#'  parameter. The approach and options are same as described above for the
#'  \code{c_cov_init_beta}.
#'
#'@param s_cov_init_beta Initial values for covariate(s) for the fixed effect b
#'  parameter. The approach and options are same as described above for the
#'  \code{a_cov_init_beta}. The option \code{lm} will set the initial for the
#'  covariate spline coeficients based on the lm model fit.
#'
#'@param a_init_sd Initial values for random effect a parameter. Options
#'  \code{0}, \code{random} and \code{prior} are same as describe earlier for
#'  \code{a_init_beta}. Additionally user can specify \code{lme_sd_a} which sets
#'  the initial based on the random intercept standard deviation obtained from
#'  the nlme::lme model fit to the data. Another option is \code{lm_sd_a} which
#'  is based on the lm model fit and essentially same as the ysd. It is
#'  mentioned here because in case nlme::lme  model fails to for some reasons
#'  and user has specified the \code{lme_sd_a}, then the option \code{lm_sd_a}
#'  will be set automatically. Additional allowed options are ysd and ymad that
#'  can be used to set standard deviation or the median absolute deviation of
#'  the outcome as initial values. See \code{a_init_beta} for specifying
#'  initials for univariate-by-subgroup (specified by using the
#'  \code{univariate_by}) and multivariate (specified by using the
#'  \code{multivariate}) models.
#'
#'@param b_init_sd Initial values for random effect a parameter. Options
#'  \code{0}, \code{random} and \code{prior} are same as describe earlier for
#'  \code{a_init_beta}. No additional option is allowed for \code{b_init_sd}.
#'
#'@param c_init_sd Initial values for random effect b parameter. The approach
#'  and options are same as described above for \code{b_init_sd}.
#'
#'@param d_init_sd Initial values for random effect d parameter. The approach
#'  and options are same as described above for \code{c_init_sd}.
#'
#'@param a_cov_init_sd Initial values for covariate(s) for the random effect a
#'  parameter. Options \code{0}, \code{random} and \code{prior} are same as
#'  described earlier for \code{a_init_beta}. No additional option is allowed
#'  for \code{a_cov_init_beta}.
#'
#'@param b_cov_init_sd Initial values for covariate(s) for the random effect b
#'  parameter. The approach and options are same as described above for the
#'  \code{a_cov_init_sd}.
#'
#'@param c_cov_init_sd Initial values for covariate(s) for the random effect c
#'  parameter. The approach and options are same as described above for the
#'  \code{b_cov_init_sd}.
#'
#'@param d_cov_init_sd Initial values for covariate(s) for the random effect d
#'  parameter. The approach and options are same as described above for the
#'  \code{c_cov_init_sd}.
#'
#'@param sigma_init_beta Initial values for fixed effect parameter sigma.
#'
#'@param sigma_cov_init_beta Initial values for covariate(s) fixed effect
#'  parameter sigma
#'
#'@param sigma_init_sd Initial values for sd of random effect parameter sigma.
#'
#'@param sigma_cov_init_sd Initial values for sd of covariate(s) for random
#'  effect parameter sigma.
#'
#'@param gr_init_cor Initial values for correlations of group-level ('random')
#'  effects. Allowed options are \code{0}, \code{random} and \code{prior}.
#'
#'@param rsd_init_sigma Initial values for residual standard deviation
#'  parameter, sigma. Allowed options are \code{0}, \code{random}, and
#'  \code{prior}. Additionally user can specify \code{lme_rsd} to sets initials
#'  based on the residual standard deviation obtained from the nlme::lme model
#'  or else \code{lm_rsd} which sets the initial value based on the lm model
#'  fit. In case user specifies the \code{lme_rsd} but for some reason model
#'  fails to converge successfully, the option \code{lme_rsd} will be set
#'  automatically. This argument is evalauted when \code{dpar_formual} is set to
#'  NULL.
#'
#'@param dpar_init_sigma Initial values for distributional parameter (i.e.,
#'  sigma). The approach and options are same as described above for the
#'  \code{rsd_init_sigma}. This argument is evalauted only when
#'  \code{dpar_formual} is not set to NULL.
#'
#'@param dpar_cov_init_sigma Initial values for covariate(s) for the
#'  distributional parameter. Allowed options are \code{0}, \code{random}, and
#'  \code{prior}.
#'
#'@param autocor_init_acor Initial values for autocorrelation parameter. see
#'  \code{autocor_formula} for details). Allowed options are \code{0},
#'  \code{random}, and \code{prior}.
#'
#'@param mvr_init_rescor Initial values for residual correlations for
#'  multivariate (specified by using the \code{multivariate}) models. Allowed
#'  options are \code{0}, \code{random}, and \code{prior}.
#'
#'@param r_init_z Initial values for standardized group level effects. These
#'  parameters are part of the central parametrisation approach adopted by the
#'  the brms package (see package brms for details).
#'
#'@param jitter_init_beta A value as proportion (between 0 and 1) to perturb the
#'  initial values for fixed effect parameters. The default is \code{NULL}
#'  implying that same initials will be used for all chains. The option which
#'  looked good during early testing is to set it as 0.1.
#'
#'@param jitter_init_sd A value as proportion (between 0 and 1) to perturb the
#'  initial values for random effect parameters. The default is \code{NULL}
#'  implying that same initials will be used for all chains. The option which
#'  looked good during early testing is to set it as 0.01.
#'
#'@param jitter_init_cor A value as proportion (between 0 and 1) to perturb the
#'  initial values for correlation parameters. The default is \code{NULL}
#'  implying that same initials will be used for all chains. The option which
#'  looked good during early testing is to set it as 0.001.
#'
#'@param prior_data An optional argument as named list to pass value for prior
#'  setting. The default is \code{NULL}. This option is particularly helpful
#'  when passing a long vector for setting priors for covariate(s) effect or a
#'  matrix. These vectors and matrices can be created in the R framework and
#'  then passed using the \code{prior_data}. For example, to pass a vector of
#'  location parameters when setting priors for covariate  with 10 dummy
#'  variables, one can create a named object prior_a_cov_beta as
#'  \code{prior_a_cov_beta = rnorm(10, 0, 1)} and then specify it as a named
#' list using \code{prior_data} as \code{prior_data = list(prior_a_cov_beta =
#' prior_a_cov_beta)} and specifying that to the \code{a_cov_prior_beta} as
#'  \code{a_cov_prior_beta = normal(prior_a_cov_beta, 100)}.
#'
#'@param init_data An optional argument as named list to pass value for setting
#'  initials. The approach is same as described above for the \code{prior_data}.
#'  As an example, create vector of initials for \code{a_cov_prior_beta} as
#' \code{init_a_cov_beta = rep(0, 10)} and then use \code{init_data =
#' list(init_a_cov_beta = init_a_cov_beta)} to pass on initials as
#'  \code{a_cov_init_beta = init_a_cov_beta}.
#'  
#'@param init_custom Set custom initials via a named list 
#'  (default \code{init}). These are are  rarely used (for example when 
#'  fitting a three level model). Note that the named list is directly passed 
#'  to the arguments without checking the dimensions.
#'
#'@param verbose An optional logical (default FALSE) argument to print step
#'  involved in preparing model formula,Stan function, priors, initials and also
#'  to report any relevant information during these processes. As an example, a
#'  user might be interested in knowing the outcomes created for the factor
#'  varibale that were used to specify the univariate-by-subgroup model. This is
#'  information then helps in matching the desired sequence of options used to
#'  pass on df, prior, initials etc.
#'
#'@param expose_function An optional logical (default TRUE) to expose Stan
#'  function used for model fitting. These functions are essential for
#'  post-processing.
#'  
#'@param get_stancode An optional logical (default \code{FALSE}) to get 
#' the stancode.
#' 
#'@param get_standata An optional logical (default \code{FALSE}) to get 
#' the standata.
#'
#'@param get_priors An optional logical (default \code{FALSE}) to get priors.
#'
#'@param get_set_priors An optional logical (default \code{FALSE}) to get 
#' priors specified by the \code{bsitar} via \code{prepare_priors}.
#'
#'@param validate_priors An optional logical (default \code{FALSE}) to
#'  validate the specified priors..
#'
#'@param set_self_priors An optional (default \code{NULL}) to specify
#'  priors manually.
#'  
#'@param set_replace_priors An optional (default \code{NULL}) to replace
#'  part of prior object.
#'  
#'@param outliers An optional (default \code{NULL}) to remove velocity
#' outliers. The argument should be a named list to pass options to the
#' [bsitar::outliers] function. See [bsitar::outliers] for details.
#'
#'@param cores Number of cores to use when executing the chains in parallel. See
#'  [brms::brm()] for details. Note that unlike [brms::brm()] which sets
#'  \code{cores=getOption("mc.cores", 1)}, the default in in \code{bsitar} is
#'  \code{cores=getOption("mc.cores", 'optimize')} which optimizes the
#'  utilization of system resources. The maximum number of cores that can be
#'  deployed is calculated as the maximum number of cores available minus 1.
#'  When the number of available cores is greater than the number chains (see
#'  \code{chains}), then number of cores is set equal to the number of chains.
#'  Another option is to set \code{cores} as \code{getOption("mc.cores",
#'  'maximise')} which set the number of cores as the maximum number of cores
#'  available from the system regardless of the number of chains specified. Note
#'  that the user can also set \code{cores} argument similar to the \code{brms}
#'  i.e., \code{getOption("mc.cores", 1)}. All these three options can be set
#'  globally as \code{options(mc.cores = x}) where x can be \code{optimize},
#'  \code{maximise} or \code{1}.
#'  Lastly, the \code{cores} can set by directly by specifying an integer e.g.,
#'  \code{cores= 4}.
#'
#'@param threads Number of threads to be used in within-chain parallelization.
#'  Note that [brms::brm()] sets this argument as
#'  \code{getOption("brms.threads", NULL)} which means that no within-chain
#'  parallelization is used by default. In contrast, to utilize the available
#'  resources from the modern computing systems, the \code{bsitar}, by default,
#'  sets \code{threads} as \code{getOption("brms.threads", 'optimize')}. The
#'  number of threads per chain is set as the maximum number of cores available
#'  minus 1. Another option is to set \code{threads} as
#'  \code{getOption("brms.threads", 'maximise')} which set the number threads
#'  per chains same as the  maximum number of cores available. User can also set
#'  the \code{threads} similar to the \code{brms} i.e.,
#'  \code{getOption("brms.threads", NULL)}. All these three options can be set
#'  globally as \code{options(brms.threads = x}) where x can be \code{optimize},
#'  \code{maximise} or \code{NULL}.
#'  Alternatively, the number of threads can BE set as \code{threads
#'  = threading(x)} where \code{X} is an integer. Other arguments that can the
#'  passed to the \code{threads} are \code{grainsize} and the \code{static}. See
#'  [brms::brm()] for further details on within-chain parallelization.
#'  
#'@param normalize Indicates whether normalization constants should be included
#'  in the Stan code (defaults to \code{TRUE}). Setting it to \code{FALSE}
#'  requires Stan version >= 2.25. If \code{FALSE}, sampling efficiency
#'  may be increased but some post processing functions such as
#'  [brms::bridge_sampler()] will not be available. Can be controlled
#'  globally for the current \R session via the `brms.normalize` option.
#' 
#'
#'@param sample_prior Indicate if draws from priors should be drawn additionally
#'  to the posterior draws. Options are \code{"no"} (the default), \code{"yes"},
#'  and \code{"only"}. Among others, these draws can be used to calculate Bayes
#'  factors for point hypotheses via [brms::hypothesis()]. Please note that
#'  improper priors are not sampled, including the default improper priors used
#'  by \code{brm}. See [brms::set_prior()] on how to set (proper) priors. Please
#'  also note that prior draws for the overall intercept are not obtained by
#'  default for technical reasons. See [brms::brmsformula()] how to obtain prior
#'  draws for the intercept. If \code{sample_prior} is set to \code{"only"},
#'  draws are drawn solely from the priors ignoring the likelihood, which allows
#'  among others to generate draws from the prior predictive distribution. In
#'  this case, all parameters must have proper priors
#'
#'@param save_model A character string or \code{NULL} (default). If not
#'  \code{NULL}, then the model's Stan code is saved via in a text file named
#'  after the string supplied in \code{save_model}.
#'
#'
#'
#'@param control A named \code{list} to control the sampler's behavior. The
#'  default are same as [brms::brm()] with the exception that the
#'  \code{max_treedepth} has been increased form 10 to 15 to allow better
#'  exploration of typically challenging posterior geometry posed by the
#'  nonlinear model. However, another control parameter, the \code{adpat_delta}
#'  which is also  often need to be increased for nonlinear model, has be set to
#'  default setting as in [brms::brm()] i.e, 0.8. This is to avoid unnecessarily
#'  increasing the sampling time. See [brms::brm()] for full details on control
#'  parameters and their default values.
#'
#'@inheritParams brms::brm
#'
#'@param ... Further arguments passed to [brms::brm()]
#'
#'@return An object of class \code{brmsfit, bsiatr}, that contains the posterior
#'  draws and other useful information about the model.
#'
#'@author Satpal Sandhu  \email{satpal.sandhu@bristol.ac.uk}
#'
#'@references
#' \insertAllCited{}
#'
#'
#'@seealso [brms::brm()] [brms::brmsformula()] [brms::prior()]
#'
#'
#' @examples
#' \dontrun{
#' # Examples below fit SITAR model to the Berkley height data obtained from
#' # 66 males and 70 females
#'
#' # First tow examples demonstrate fitting two separate univariate models for
#' # males and females and then a combined univariate-by-subgroup model. The
#' # third example show multivariate model fitting.
#'
#' # As shown below, univariate-by-subgroup model internally fits two
#' # sub-models, one for males and another for females. Advantage of fitting
#' # univariate-by-subgroup model is that posterior samples for all outcomes
#' # are contained in a single framework which can then be used for direct
#' # comparisons during the post-processing (e.g., hypothesis testing). The
#' # flexibility offered by the 'bsitar' package allows full control over
#' # the sub-models (e.g. df for spline curve, priors, initials etc.).
#' # Below we fit models with default setting with 4 chains and 2000 iter
#' #
#' # Fit 1 - Separate models for males and females with 5 df for males and
#' # 4 df for females.
#'
#' # Prepare data
#' data(heights)
#' data_males <- heights %>% filter(sex == 'Male)
#' data_females <- heights %>% filter(sex == 'Female)
#'
#' # Fit model
#' fit_males <- bsitar(x=age, y=height, id=id, data=heights, df=5)
#' fit_females <- bsitar(x=age, y=height, id=id, data=heights, df=4)
#'
#' # Generate a summary of results for males and females
#' summary(fit_males)
#' summary(fit_females)
#'
#' # Perform posterior predictive checks for males and females
#' pp_check(fit_males)
#' pp_check(fit_females)
#'
#' # plot distance and velocity curves for males and females
#' # Distance
#' plot(conditional_effects(fit_males, deriv = 0))
#' plot(conditional_effects(fit_females, deriv = 0))
#'
#' # Velocity
#' plot(conditional_effects(fit_males, deriv = 1))
#' plot(conditional_effects(fit_females, deriv = 1))
#'
#' # Fit 2 - univariate-by-subgroup model for sex (males and females) with 5 df
#' # for males and 4 df for females. Since factor variable sex is sorted
#' # alphabatically, the first reponse vector created is for females and second
#' # for males. As shown below for df, controlling any argument is as simple
#' # as enclosing it in list and separate arguments by a comma. Same approach
#' # applies for all argument including prior and initials.
#'
#' # Fit model
#' fit_male_female <- bsitar(x=age, y=height, id=id, data=heights,
#' univariate_by = sex, df=list(4,5))
#'
#' # Generate a summary of results for males and females
#' summary(fit_male_female)
#'
#' # Perform posterior predictive checks (specify response option resp = )
#' pp_check(fit_male_female, resp = 'Male')
#' pp_check(fit_male_female, resp = 'Female')
#'
#' # plot distance and velocity curves for males and females
#' # Distance
#' plot(conditional_effects(fit_male_female, deriv = 0, resp = 'Male'))
#' plot(conditional_effects(fit_male_female, deriv = 0, resp = 'Female'))
#'
#' # Velocity
#' plot(conditional_effects(fit_male_female, deriv = 1, resp = 'Male'))
#' plot(conditional_effects(fit_male_female, deriv = 1, resp = 'Female'))
#'
#'
#' # Fit 3 - multivariate model
#' # For demonstration purposes, we use the same heights data and artificially
#' # create the second outcome (height2) by rescaling the original first
#' # outcome (i.e., height). Again we use different degree of freedom (df) for
#' # each outcome. Here we also show how to use different priors and initials
#' # for some of the parameters.
#' data_heights2 <- heights %>% mutate(height2 = (height - 10) * 0.1)
#'
#' # Fit model
#' # We specified multivariate = TRUE for fitting multivariate model. By
#' # default, the cor structure will be set to un for modelling unstructured
#' # varinace covaraince with joint distribution of groop level random effects.
#' # Also, option rescor for modelling residual correlation is set to TRUE
#' # (default). These options can be modified by explicitly setting the
#' # multivariate argument as a list, e.g., multivariate = list(mvar = TRUE,
#' # cor = un, rescor = TRUE). This allows changing the cor suboptions to un_s
#' # or diagonal, and rescor to FALSE (see \code{multivariate} for details).
#'
#' # In the example shown below, we set df = 4 for the first outcome, height
#' # and df = 5 for the second outcome, height2. We set prior normal(ymean, ysd)
#' # for outcome height and cauchy(ymedian, 100) for the second outcome height2.
#' # Note data we set different autosclae values (2 for the first outcome and
#' # default FALSE for the second outcome). Also, we have set random initial
#' # for the first outcome and lm to the second outcome. Post-processing for
#' # multivariate model is same as univariate-by-subgroup model i.e., by using
#' # the resp = argument.
#'
#' # Fit model
#' fit_mutivar <- bsitar(x=age, y=list(height, height2), id=id, data=heights,
#' multivariate = TRUE,  df=list(4,5),
#' a_prior_beta = list(normal(ymean, ysd, autosclae = 2), cauchy(ymedian, 100)),
#' a_init_beta = list(random, lm))
#'
#' # Generate a summary of results for height and height2
#' summary(fit_mutivar)
#'
#' # Perform posterior predictive checks for height and height2
#' pp_check(fit_male_female, resp = 'height')
#' pp_check(fit_male_female, resp = 'height2')
#'
#' # plot distance and velocity curves for height and height2
#' # Distance
#' plot(conditional_effects(fit_male_female, deriv = 0, resp = 'height'))
#' plot(conditional_effects(fit_male_female, deriv = 0, resp = 'height2'))
#'
#' # Velocity
#' plot(conditional_effects(fit_male_female, deriv = 1, resp = 'height'))
#' plot(conditional_effects(fit_male_female, deriv = 1, resp = 'height2'))
#'
#' }
#'
#'@importFrom sitar getPeak
#'@importFrom methods formalArgs
#'@importFrom stats as.formula coef df dist filter fitted gaussian lm mad median
#'  model.matrix predict quantile rbeta sd setNames smooth.spline 
#'  
#'@importFrom stats loess na.omit residuals
#'@importFrom utils combn head installed.packages packageVersion
#'@importFrom Rdpack reprompt
#'@import brms
#'
#'@export
#'

bsitar <- function(x,
                   y,
                   id,
                   data,
                   df = 4,
                   knots = NA,
                   fixed = a + b + c,
                   random = a + b + c,
                   xoffset = mean,
                   bstart = mean,
                   xfun = NULL,
                   yfun = NULL,
                   bound = 0.04,
                   terms_rhs = NULL,
                   a_formula = ~ 1,
                   b_formula = ~ 1,
                   c_formula = ~ 1,
                   d_formula = ~ 1,
                   s_formula = ~ 1,
                   a_formula_gr = ~ 1,
                   b_formula_gr = ~ 1,
                   c_formula_gr = ~ 1,
                   d_formula_gr = ~ 1,
                   
                   a_formula_gr_str = NULL,
                   b_formula_gr_str = NULL,
                   c_formula_gr_str = NULL,
                   d_formula_gr_str = NULL,
                   
                   sigma_formula = NULL,
                   sigma_formula_gr = NULL,
                   sigma_formula_gr_str = NULL,
                   
                   dpar_formula = NULL,
                   autocor_formula = NULL,
                   family = gaussian(),
                   group_arg = list(
                     groupvar = NULL,
                     by = NULL,
                     cor = un,
                     cov = NULL,
                     dist = gaussian
                   ),
                   sigma_group_arg = list(
                     groupvar = NULL,
                     by = NULL,
                     cor = un,
                     cov = NULL,
                     dist = gaussian
                   ),
                   univariate_by = list(by = NA, cor = un),
                   multivariate = list(mvar = FALSE,
                                       cor = un,
                                       rescor = TRUE),
                   
                   a_prior_beta = normal(ymean, ysd, autoscale = 2.5),
                   b_prior_beta = normal(0, 2, autoscale = FALSE),
                   c_prior_beta = normal(0, 1, autoscale = FALSE),
                   d_prior_beta = normal(0, 1, autoscale = FALSE),
                   s_prior_beta = normal(0, lm, autoscale = 2.5),
                   a_cov_prior_beta = normal(0, 5, autoscale = FALSE),
                   b_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                   c_cov_prior_beta = normal(0, 0.1, autoscale = FALSE),
                   d_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                   s_cov_prior_beta = normal(0, 10, autoscale = FALSE),
                   a_prior_sd = normal(0, ysd, autoscale = 1),
                   b_prior_sd = normal(0, 1, autoscale = FALSE),
                   c_prior_sd = normal(0, 0.5, autoscale = FALSE),
                   d_prior_sd = normal(0, 1, autoscale = FALSE),
                   a_cov_prior_sd = normal(0, 2, autoscale = FALSE),
                   b_cov_prior_sd = normal(0, 1, autoscale = FALSE),
                   c_cov_prior_sd = normal(0, 0.15, autoscale = FALSE),
                   d_cov_prior_sd = normal(0, 0.5, autoscale = FALSE),
                   
                   sigma_prior_beta = normal(0, 1, autoscale = FALSE),
                   sigma_cov_prior_beta = normal(0, 0.5, autoscale = FALSE),
                   sigma_prior_sd = normal(0, 0.25, autoscale = FALSE),
                   sigma_cov_prior_sd = normal(0, 0.15, autoscale = FALSE),
                   
                   rsd_prior_sigma = normal(0, ysd, autoscale = FALSE),
                   dpar_prior_sigma = normal(0, ysd, autoscale = FALSE),
                   dpar_cov_prior_sigma = normal(0, 5, autoscale = FALSE),
                   autocor_prior_acor = uniform(-1, 1, autoscale = FALSE),
                   gr_prior_cor = lkj(1),
                   mvr_prior_rescor = lkj(1),
                   init = NULL,
                   init_r = NULL,
                   a_init_beta = lm,
                   b_init_beta = 0.001,
                   c_init_beta = 0.001,
                   d_init_beta = 0,
                   s_init_beta = lm,
                   a_cov_init_beta = 0,
                   b_cov_init_beta = 0,
                   c_cov_init_beta = 0,
                   d_cov_init_beta = 0,
                   s_cov_init_beta = 0,
                   a_init_sd = 1, # lme_sd_a
                   b_init_sd = 1,
                   c_init_sd = 1,
                   d_init_sd = 1,
                   a_cov_init_sd = 0.5,
                   b_cov_init_sd = 0.5,
                   c_cov_init_sd = 0.5,
                   d_cov_init_sd = 0.5,
                   
                   sigma_init_beta = 0,
                   sigma_cov_init_beta = 0,
                   sigma_init_sd = 1,
                   sigma_cov_init_sd = 1,
                   
                   gr_init_cor = 0,
                   rsd_init_sigma = 1,
                   dpar_init_sigma = 1,
                   dpar_cov_init_sigma = 1,
                   autocor_init_acor = 0.1,
                   mvr_init_rescor = 0,
                   r_init_z = 0,
                   jitter_init_beta = 0.01,
                   jitter_init_sd = NULL,
                   jitter_init_cor = NULL,
                   prior_data = NULL,
                   init_data = NULL,
                   init_custom = NULL,
                   verbose = FALSE,
                   expose_function = TRUE,
                   
                   get_stancode = FALSE,
                   get_standata = FALSE,
                   get_priors = FALSE,
                   get_set_priors = FALSE,
                   validate_priors = FALSE,
                   set_self_priors = NULL,
                   set_replace_priors = NULL,
                   outliers = NULL, 
                   
                   chains = 4,
                   iter = 2000,
                   warmup = floor(iter / 2),
                   thin = 1,
                   cores = getOption("mc.cores", "optimize"),
                   backend = getOption("brms.backend", "rstan"),
                   threads = getOption("brms.threads", "optimize"),
                   opencl = getOption("brms.opencl", NULL),
                   normalize = getOption("brms.normalize", TRUE),
                   algorithm = getOption("brms.algorithm", "sampling"),
                   control = list(adapt_delta = 0.8, max_treedepth = 15),
                   sample_prior = "no",
                   save_pars = NULL,
                   drop_unused_levels = TRUE,
                   stan_model_args = list(),
                   silent = 1,
                   seed = 123,
                   save_model = NULL,
                   fit = NA,
                   file = NULL,
                   file_refit = getOption("brms.file_refit", "never"),
                   future = getOption("future", FALSE),
                   ...) {
  mcall <- mcall_ <- match.call()
  
  # check and update if argument value taken from the global environment.
  # x,y,id and data are always taken from the data, thus excluded from checks
  
  no_default_args <- c("x", "y", "id", "data", "...")
  
  deparse_0 <- function(deparseobj) {
    deparseobj <- paste(deparse(deparseobj), collapse = "")
    deparseobj <- gsub("[[:space:]]", "", deparseobj)
    deparseobj
  }
  
  deparse_0s <- function(deparseobj) {
    deparseobj <- paste(deparse(substitute(deparseobj)), collapse = "")
    deparseobj
  }
  
  gsub_space <- function(deparseobj) {
    deparseobj <- gsub("[[:space:]]", "", deparseobj)
    deparseobj
  }
  
  
  
  terms_rhsxx <- deparse(mcall[['terms_rhs']])
  terms_rhsxx <- gsub("[[:space:]]", "", terms_rhsxx)
  terms_rhsxx <- gsub("^list\\(", "", terms_rhsxx)
  if(grepl("^list", terms_rhsxx)) {
    terms_rhsxx <- gsub(")$", "", terms_rhsxx)
  }
  terms_rhsxx <- strsplit(terms_rhsxx, ",")[[1]]
  terms_rhsxx2 <- terms_rhsxx
  terms_rhsxx <- paste(terms_rhsxx, sep = ",")
  terms_rhsxx <- gsub("\"", "", terms_rhsxx)
  if(grepl("^NULL", terms_rhsxx)) {
    terms_rhsxx <- gsub(")$", "", terms_rhsxx)
  }
  terms_rhsxx <- paste0("'", terms_rhsxx, "'")
  terms_rhsxx <- paste(terms_rhsxx, collapse = ",")
  terms_rhsxx <- paste0("list(", terms_rhsxx, ")")
  terms_rhsxx <- parse(text = terms_rhsxx)
  mcall[['terms_rhs']] <- terms_rhsxx
  
  
  
  
  
  
  xs <- ids <- dfs <- NA
  
  
  for (i in names(mcall)[-1]) {
    if (!i %in% no_default_args) {
      err. <- FALSE
      tryCatch(
        expr = {
          # checks. <- eval(mcall[[i]])
          if (is.function(eval(mcall[[i]]))) {
            checks. <- deparse_0(mcall[[i]])
          } else {
            checks. <- eval(mcall[[i]])
          }
        },
        error = function(e) {
          err. <<- TRUE
        }
      )
      if (err.) {
        mcall[[i]] <- mcall[[i]]
      } else if (!err.) {
        if (is.list(checks.)) {
          if (is.list(checks.[[1]])) {
            mcall[[i]] <- checks.[[1]]
          } else if (!is.list(checks.[[1]])) {
            if (is.list(checks.)) {
              if (is.symbol(mcall[[i]]))
                mcall[[i]] <- deparse_0(mcall[[i]]) # for set_self_priors
              mcall[[i]] <- eval(mcall[[i]])
              temp <- str2lang(deparse_0((mcall[[i]])))
              mcall[[i]] <- temp
            } else if (!is.list(checks.)) {
              mcall[[i]] <- checks.
            }
          }
        } else {
          mcall[[i]] <-  checks.
        }
      }
    } else if (i %in% no_default_args) {
      mcall[[i]] <-  mcall[[i]]
    }
  }
  
  
  
  arguments <- as.list(mcall)[-1]
  
  
  match.call.defaults <- function(...) {
    call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
    formals <- evalq(formals(), parent.frame(1))
    
    for(i in setdiff(names(formals), names(call)))
      call[i] <- list( formals[[i]] )
    
    match.call(sys.function(sys.parent()), call)
  }
  
  call.full <- match.call.defaults()
  call.full <- call.full[-length(call.full)]
  
  for (call.fulli in names(call.full)) {
    if(call.fulli != "") {
      if(call.fulli == 'family' & 
         is.language(call.full[[call.fulli]])) {
        call.full[[call.fulli]] <- deparse(call.full[[call.fulli]])
      } else if(call.fulli == 'stan_model_args')  { 
        if(length(eval(call.full[[call.fulli]])) == 0) {
          call.full[[call.fulli]] <- NULL
        } else {
          call.full[[call.fulli]] <- call.full[[call.fulli]] 
        }
      } else {
        #  call.full[[call.fulli]] <- call.full[[call.fulli]]
      }
    } else {
      # call.full[[call.fulli]] <- call.full[[call.fulli]]
    }
  }
  
  
  # combine user defined and default arguments
  `%!in%` <- Negate(`%in%`)
  f_bsitar_arg <- formals(bsitar)
  nf_bsitar_arg_names <-
    intersect(names(arguments), names(f_bsitar_arg))
  arguments <-
    c(arguments, f_bsitar_arg[names(f_bsitar_arg) %!in% nf_bsitar_arg_names])
  
  
  
  for (ip in names(arguments)) {
    if (grepl("_init_", ip)) {
      d_mcall_ <- deparse_0(mcall_[[ip]])
      if(is.symbol(mcall_[[ip]])) {
        arguments[[ip]] <- d_mcall_
      } else if(grepl("c(", d_mcall_, fixed = T)) {
        arguments[[ip]] <- gsub("c(", "list(", d_mcall_, fixed = T)
      } else if(grepl("list(", d_mcall_, fixed = T)) {
        arguments[[ip]] <- mcall_[[ip]]
      }
    }
  }
  
  # arguments2 <<- arguments
  
  remove_spaces <- c('a_formula_gr_str', 'b_formula_gr_str', 
                     'c_formula_gr_str', 'd_formula_gr_str',
                     'sigma_formula_gr_str')
  
  for (ip in remove_spaces) {
    arguments[[ip]] <-  gsub_space(arguments[[ip]] )
  }
  
  
  
  check_gr_str_form <- function(x, x__) {
    if(!is.null(x) | !is.null(x[[1]])) {
      if(!grepl("^list", x__)) {
        if(!grepl("^~", x__)) {
          stop("Argument ", deparse(substitute(x)), " should be a formula.",
               "\n ",
               " Please add '~' at the begining ")
        }
      } 
    }
  }
  
  check_gr_str_form(sigma_formula_gr_str, 
                    deparse(substitute(sigma_formula_gr_str)))
  
  
  # Separate 'brms' arguments from 'bsitar' arguments for the ease of handling
  
  brms_arguments_list <-
    c(
      'chains',
      'iter',
      'warmup',
      'thin',
      'cores',
      'backend',
      'threads',
      'opencl',
      'normalize',
      'algorithm',
      'control',
      'sample_prior',
      'save_pars',
      'drop_unused_levels',
      'stan_model_args',
      'silent',
      'seed',
      'save_model',
      'fit',
      'file',
      'file_refit',
      'future'
    )
  
  
  #######
  mc.cores_restore <- getOption("mc.cores")
  if(is.numeric(arguments$cores)) {
    options(mc.cores = arguments$cores)
    # arguments$cores <- getOption("mc.cores", arguments$cores)
  }
  
  iter <- arguments$iter
  warmup <- arguments$warmup <- eval(arguments$warmup)
  # print(arguments$iter)
  # print(arguments$warmup)
  ######
  
  
  brms_arguments <- list()
  for (brms_arguments_listi in brms_arguments_list) {
    brms_arguments[[brms_arguments_listi]] <-
      arguments[[brms_arguments_listi]]
    arguments[[brms_arguments_listi]] <- NULL
  }
  
  #### Note
  # brms_arguments <- list()
  # brms_arguments <- mget(brms_arguments_list)
  
  brms_arguments <- mget(brms_arguments_list)
  
  
  
  if (eval(brms_arguments$backend) != "rstan" &
      eval(brms_arguments$backend) != "cmdstanr") {
    stop("The backend argument must be either 'rstan' or 'cmdstanr'",
         "\n ",
         "\ Please check it which you have specified as: ", 
         eval(brms_arguments$backend))
  }
  
  
  ept <- function(x)
    eval(parse(text = x), envir = parent.frame())
  
  # Internal controls
  
  # 'd' formula control
  # 1) Match with 'sitar' package (i.e., exclude 'd' from the fixed effects)
  # 2) Or, like a,b,c parameters (i.e., include 'd' in fixed and random effects)
  # Setting default to FALSE to match scenario 2 for now
  # TODO
  match_sitar_d_form <- FALSE
  
  # Display method - either message or custom color
  # Can be moved to the arguments but not that worth
  displayit <- 'col'
  setcolh   <- 47 
  setcolb   <- 3
  
  # Quote unquoted character (e.g., sex to 'sex') for user's convineinec
  
  list_to_quoted_if_not <- function(x) {
    splitmvar <- x
    splitmvar <- gsub("\\s", "", splitmvar)
    splitmvar <- paste(splitmvar, collapse = "")
    splitmvar_w <-
      gsub("[\\(\\)]", "", regmatches(splitmvar, gregexpr("\\(.*?\\)",
                                                          splitmvar))[[1]])
    splitmvar_w2 <- strsplit(splitmvar_w, ",")[[1]]
    splitmvar_w3 <- sub("=[^=]+$", "", splitmvar_w2)
    gsubs_c_counter <- 0
    for (i in splitmvar_w3) {
      gsubs_c_counter <- gsubs_c_counter + 1
      if (gsubs_c_counter < max(length(splitmvar_w3))) {
        pattern <- paste0(i, "=", "\\s*(.*?)\\s*", ",")
      } else {
        pattern <- paste0(i, "=", "\\s*(.*?)\\s*", ")")
      }
      majors <- regmatches(splitmvar, regexec(pattern, splitmvar))
      majors2 <- majors[[1]][2]
      majors2 <- majors[[1]][2]
      if (grepl("^T$", majors2)) {
        majors2 <- gsub("^T$", "TRUE", majors2)
      }
      if (grepl("^F$", majors2)) {
        majors2 <- gsub("^F$", "FALSE", majors2)
      }
      majors2 <- gsub("\"", "", majors2)
      majors3 <- paste0("\"", majors2, "\"")
      if (gsubs_c_counter == 1) {
        splitmvar2 <- gsub(noquote(majors2), majors3, splitmvar, fixed = F)
      } else {
        splitmvar2 <- gsub(noquote(majors2), majors3, splitmvar2, fixed = F)
      }
    }
    for (i in 1:length(splitmvar_w3))
      splitmvar2 <- gsub("\"\"", "\"", splitmvar2)
    splitmvar3 <- eval(parse(text = splitmvar2))
    zzz <- splitmvar3
    for (z in names(splitmvar3)) {
      err. <- FALSE
      tryCatch(
        expr = {
          eval(parse(text = zzz[[z]]), envir = parent.frame())
        },
        error = function(e) {
          err. <<- TRUE
        }
      )
      if (!err.) {
        c_c_ <- eval(parse(text = zzz[[z]]))
        checkclass <- class(c_c_)
        if (checkclass == "NULL")
          checkclass_ <- NULL
        else
          checkclass_ <- NA
        if (is.logical(c_c_) | is.null(checkclass_))
          zzz[[z]] <- c_c_
      } else {
        zzz[[z]] <- zzz[[z]]
      }
    }
    return(zzz)
  }
  
  
  list_to_quoted_if_not_si <- function(xx) {
    xx.o <- xx
    prefix_ <- strsplit(xx, "\\(")[[1]][1]
    prefix_by <- "list"
    xx <- gsub(paste0("^", prefix_, ""), prefix_by, xx)
    if (sub("\\).*", "", sub(".*\\(", "", xx)) != "") {
      xxx <- list_to_quoted_if_not(xx)
      xxx <- gsub("\"" , "'", deparse_0(xxx))
      xxx <- gsub(paste0("^", prefix_by, ""), prefix_, xxx)
      xxx <- gsub("\\s", "", xxx)
    } else {
      xxx <- xx.o
    }
    xxx
  }
  
  
  list_to_quoted_if_not_si_lf <- function(xx) {
    xx.o <- xx
    prefix_ <- strsplit(xx, "\\(")[[1]][1]
    prefix_by <- "list"
    xx <- gsub(paste0("^", prefix_, ""), prefix_by, xx)
    if (sub("\\).*", "", sub(".*\\(", "", xx)) != "") {
      xxt <- sub("\\).*", "", sub(".*\\(", "", xx))
      xxtf <-
        strsplit(xxt, ",")[[1]][grepl("~", strsplit(xxt, ",")[[1]])]
      xxtnf <-
        strsplit(xxt, ",")[[1]][!grepl("~", strsplit(xxt, ",")[[1]])]
      xxtf <- gsub("\\s", "", xxtf)
      xxtnf <- gsub("\\s", "", xxtnf)
      xx <-
        paste0(prefix_by, "(", paste(xxtnf, collapse = ","), ")")
      xxx <- list_to_quoted_if_not(xx)
      xxx <- gsub("\"" , "'", deparse_0(xxx))
      xxx <- gsub(paste0("^", prefix_by, ""), prefix_, xxx)
      xxx <-
        gsub(paste0(prefix_, "\\("),
             paste0(prefix_, "(", xxtf, ","),
             xxx)
      xxx <- gsub("\\s", "", xxx)
    } else {
      xxx <- xx.o
    }
    xxx
  }
  
  
  # set multivariate arguments
  
  if (gsub("\\s", "",
           paste(deparse(substitute(multivariate)), collapse = "")) == "NULL" |
      gsub("\\s", "",
           paste(deparse(substitute(multivariate)), collapse = "")) == "NA" |
      gsub("\\s", "",
           paste(deparse(substitute(multivariate)), collapse = "")) == "FALSE" |
      gsub("\\s", "",
           paste(deparse(substitute(multivariate)), collapse = "")) == "F") {
    multivariate <- list()
    multivariate$mvar <- FALSE
  } else if (gsub("\\s", "",
                  paste(deparse(substitute(multivariate)),
                        collapse = "")) == "TRUE" |
             gsub("\\s", "",
                  paste(deparse(substitute(multivariate)),
                        collapse = "")) == "T") {
    multivariate <- list()
    multivariate$mvar <- TRUE
  } else if (!grepl("^list", gsub("\\s", "", paste(deparse(
    substitute(multivariate)
  ), collapse = ""))) &
  !is.null(gsub("\\s", "", paste(deparse(
    substitute(multivariate)
  ), collapse = "")))) {
    if (is.symbol(substitute(multivariate))) {
      multivariate <-
        gsub("\\s", "", paste(deparse(substitute(multivariate)), collapse = ""))
      if (multivariate == "T")
        multivariate <- eval(parse(text = multivariate))
      multivariate <- multivariate
      multivariate <- as.list(multivariate)
      names(multivariate) <- 'mvar'
    } else if (is.character(substitute(multivariate))) {
      multivariate <- multivariate
      multivariate <- as.list(multivariate)
      names(multivariate) <- 'mvar'
    }
  }
  if (grepl("^list", gsub("\\s", "", paste(deparse(
    substitute(multivariate)
  ), collapse = ""))) &
  length(strsplit(gsub("\\s", "", paste(
    deparse(substitute(multivariate)), collapse = ""
  )), ",")[[1]]) == 1) {
    if (!is.null(gsub("\\s", "", paste(deparse(
      substitute(multivariate)
    ), collapse = "")))) {
      if (is.language(substitute(multivariate))) {
        ttt <-
          gsub("\\s", "", paste(deparse(substitute(
            multivariate
          )), collapse = ""))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("mvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        multivariate <- list_to_quoted_if_not(ttt)
        
      } else if (grepl("^list", multivariate)) {
        ttt <- deparse_0(as.name(substitute(multivariate)))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("mvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        multivariate <- list_to_quoted_if_not(ttt)
        for (multivariatei in 1:length(multivariate)) {
          if (!is.null(multivariate[[multivariatei]])) {
            multivariate[[multivariatei]] <-
              gsub("'", "", multivariate[[multivariatei]])
          }
        }
      } else {
        if (!is.null(multivariate)) {
          if (!grepl("^list", multivariate)) {
            multivariate <- multivariate
            multivariate <- as.list(multivariate)
            names(multivariate) <- 'mvar'
          }
        } else if (is.null(multivariate)) {
          multivariate <- as.list(multivariate)
        }
      }
    }
  }
  if (length(strsplit(gsub("\\s", "", paste(
    deparse(substitute(multivariate)), collapse = ""
  )), ",")[[1]]) > 1) {
    ttt <-
      gsub("\\s", "", paste(deparse(substitute(multivariate)), collapse = ""))
    temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
    temp <- gsub("\\s", "", temp)
    if (!grepl("^mvar=", temp[1])) {
      temp[1] <- paste0("mvar=", temp[1])
    }
    temp <- paste(temp, collapse = ",")
    temp <- paste0("list(", temp, ")")
    multivariate <- list_to_quoted_if_not(temp)
  }
  
  
  # Set univariate_by arguments (i.e., univariate-by=subgroup model fitting)
  
  if (gsub("\\s", "",
           paste(deparse(substitute(univariate_by)), 
                 collapse = "")) == "NULL" |
      gsub("\\s", "",
           paste(deparse(substitute(univariate_by)), 
                 collapse = "")) == "NA" |
      gsub("\\s", "",
           paste(deparse(substitute(univariate_by)), 
                 collapse = "")) == "FALSE" |
      gsub("\\s", "",
           paste(deparse(substitute(univariate_by)), 
                 collapse = "")) == "F") {
    univariate_by <- list()
    univariate_by$by <- NA
  } else if (!grepl("^list", gsub("\\s", "", paste(deparse(
    substitute(univariate_by)
  ), collapse = ""))) &
  !is.null(gsub("\\s", "", paste(deparse(
    substitute(univariate_by)
  ), collapse = "")))) {
    if (is.symbol(substitute(univariate_by))) {
      univariate_by <-
        gsub("\\s", "", paste(deparse(substitute(
          univariate_by
        )), collapse = ""))
      univariate_by <- univariate_by
      univariate_by <- as.list(univariate_by)
      names(univariate_by) <- 'by'
    } else if (is.character(substitute(univariate_by))) {
      univariate_by <- univariate_by
      univariate_by <- as.list(univariate_by)
      names(univariate_by) <- 'by'
    }
  }
  if (grepl("^list", gsub("\\s", "", paste(deparse(
    substitute(univariate_by)
  ), collapse = ""))) &
  length(strsplit(gsub("\\s", "", paste(
    deparse(substitute(univariate_by)), collapse = ""
  )), ",")[[1]]) == 1) {
    if (!is.null(gsub("\\s", "", paste(deparse(
      substitute(univariate_by)
    ), collapse = "")))) {
      if (is.language(substitute(univariate_by))) {
        ttt <-
          gsub("\\s", "", paste(deparse(substitute(
            univariate_by
          )), collapse = ""))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("by=", strsplit(temp, "=")[[1]]), ttt)
        }
        univariate_by <- list_to_quoted_if_not(ttt)
        
      } else if (grepl("^list", univariate_by)) {
        ttt <- deparse_0(as.name(substitute(univariate_by)))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("by=", strsplit(temp, "=")[[1]]), ttt)
        }
        univariate_by <- list_to_quoted_if_not(ttt)
        for (univariate_byi in 1:length(univariate_by)) {
          if (!is.null(univariate_by[[univariate_byi]])) {
            univariate_by[[univariate_byi]] <-
              gsub("'", "", univariate_by[[univariate_byi]])
          }
        }
      } else {
        if (!is.null(univariate_by)) {
          if (!grepl("^list", univariate_by)) {
            univariate_by <- univariate_by
            univariate_by <- as.list(univariate_by)
            names(univariate_by) <- 'by'
          }
        } else if (is.null(univariate_by)) {
          univariate_by <- as.list(univariate_by)
        }
      }
    }
  }
  if (length(strsplit(gsub("\\s", "", paste(
    deparse(substitute(univariate_by)), collapse = ""
  )), ",")[[1]]) > 1) {
    ttt <-
      gsub("\\s", "", paste(deparse(substitute(univariate_by)), collapse = ""))
    temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
    temp <- gsub("\\s", "", temp)
    if (!grepl("^by=", temp[1])) {
      temp[1] <- paste0("by=", temp[1])
    }
    temp <- paste(temp, collapse = ",")
    temp <- paste0("list(", temp, ")")
    univariate_by <- list_to_quoted_if_not(temp)
  }
  
  
  
  # Set group_arg arguments (for univariate model fitting)
  
  if (!paste(deparse(substitute(group_arg)), collapse = "") == "NULL"  &
      !any(grepl("^list", gsub("\\s", "", paste(
        deparse(substitute(group_arg)), collapse = ""
      )))) &
      any(gsub("\\s", "", paste(deparse(
        substitute(group_arg)
      ), collapse = "")) == "NULL")) {
    group_arg <- list()
    group_arg$groupvar <- NULL
  } else if (!any(grepl("^list", gsub("\\s", "", paste(
    deparse(substitute(group_arg)), collapse = ""
  )))) &
  any(gsub("\\s", "", paste(deparse(
    substitute(group_arg)
  ), collapse = "")) != "NULL")) {
    if (paste(deparse(substitute(group_arg)), collapse = "") == "T" |
        paste(deparse(substitute(group_arg)), collapse = "") == "TRUE" |
        paste(deparse(substitute(group_arg)), collapse = "") == "F" |
        paste(deparse(substitute(group_arg)), collapse = "") == "FALSE" |
        paste(deparse(substitute(group_arg)), collapse = "") == "NA") {
      stop("group_arg should be either NULL or a character",
           " denoting the group idetifier")
    }
    if (is.symbol(substitute(group_arg))) {
      group_arg <-
        gsub("\\s", "", paste(deparse(substitute(group_arg)), collapse = ""))
      group_arg <- group_arg
      group_arg <- as.list(group_arg)
      names(group_arg) <- 'groupvar'
    } else if (is.character(substitute(group_arg))) {
      group_arg <- group_arg
      group_arg <- as.list(group_arg)
      names(group_arg) <- 'groupvar'
    }
  }
  if (any(grepl("^list", gsub("\\s", "",
                              paste(
                                deparse(substitute(group_arg)),
                                collapse = ""
                              )))) &
      length(strsplit(gsub("\\s", "",
                           paste(
                             deparse(substitute(group_arg)),
                             collapse = ""
                           )), ",")[[1]]) == 1) {
    if (!is.null(gsub("\\s", "", paste(deparse(
      substitute(group_arg)
    ),
    collapse = "")))) {
      if (is.language(substitute(group_arg))) {
        ttt <- gsub("\\s", "", paste(deparse(substitute(group_arg)),
                                     collapse = ""))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "T" |
            temp == "TRUE" |
            temp == "F" |
            temp == "FALSE") {
          stop(
            "group_arg should be either NULL or a character",
            " denoting the group idetifier"
          )
        }
        if (temp == "") {
          stop("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("groupvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        group_arg <- list_to_quoted_if_not(ttt)
      } else if (grepl("^list", group_arg)) {
        ttt <- deparse_0(as.name(substitute(group_arg)))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("groupvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        group_arg <- list_to_quoted_if_not(ttt)
        for (group_argi in 1:length(group_arg)) {
          if (!is.null(group_arg[[group_argi]])) {
            group_arg[[group_argi]] <- gsub("'", "", group_arg[[group_argi]])
          }
        }
      } else {
        if (!is.null(group_arg)) {
          if (!grepl("^list", gsub("\\s", "",
                                   paste(
                                     deparse(substitute(group_arg)),
                                     collapse = ""
                                   )))) {
            group_arg <- group_arg
            group_arg <- as.list(group_arg)
            names(group_arg) <- 'groupvar'
          }
        } else if (is.null(group_arg)) {
          group_arg <- as.list(group_arg)
        }
      }
    } else if (is.null(group_arg)) {
      group_arg <- list()
      group_arg$groupvar <- NULL
    }
  }
  if (any(grepl("^list", gsub("\\s", "",
                              paste(
                                deparse(substitute(group_arg)),
                                collapse = ""
                              )))) &
      length(strsplit(gsub("\\s", "",
                           paste(
                             deparse(substitute(group_arg)),
                             collapse = ""
                           )), ",")[[1]]) > 1) {
    ttt <-
      gsub("\\s", "", paste(deparse(substitute(group_arg)), collapse = ""))
    temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
    temp <- gsub("\\s", "", temp)
    if (!grepl("^groupvar=", temp[1])) {
      temp[1] <- paste0("groupvar=", temp[1])
    }
    temp <- paste(temp, collapse = ",")
    temp <- paste0("list(", temp, ")")
    group_arg <- list_to_quoted_if_not(temp)
  }
  if (length(group_arg) == 0) {
    group_arg <- list()
    group_arg$groupvar <- NULL
  }
  if (!is.null(group_arg$groupvar) &
      !is.character(group_arg$groupvar)) {
    stop("group_arg should be either NULL or a character",
         " denoting the group idetifier")
  }
  
  
  
  
  
  
  
  
  # Set sigma_group_arg arguments (for univariate model fitting)
  
  if (!paste(deparse(substitute(sigma_group_arg)), collapse = "") == "NULL"  &
      !any(grepl("^list", gsub("\\s", "", paste(
        deparse(substitute(sigma_group_arg)), collapse = ""
      )))) &
      any(gsub("\\s", "", paste(deparse(
        substitute(sigma_group_arg)
      ), collapse = "")) == "NULL")) {
    sigma_group_arg <- list()
    sigma_group_arg$groupvar <- NULL
  } else if (!any(grepl("^list", gsub("\\s", "", paste(
    deparse(substitute(sigma_group_arg)), collapse = ""
  )))) &
  any(gsub("\\s", "", paste(deparse(
    substitute(sigma_group_arg)
  ), collapse = "")) != "NULL")) {
    if (paste(deparse(substitute(sigma_group_arg)), collapse = "") == "T" |
        paste(deparse(substitute(sigma_group_arg)), collapse = "") == "TRUE" |
        paste(deparse(substitute(sigma_group_arg)), collapse = "") == "F" |
        paste(deparse(substitute(sigma_group_arg)), collapse = "") == "FALSE" |
        paste(deparse(substitute(sigma_group_arg)), collapse = "") == "NA") {
      stop("sigma_group_arg should be either NULL or a character",
           " denoting the group idetifier")
    }
    if (is.symbol(substitute(sigma_group_arg))) {
      sigma_group_arg <-
        gsub("\\s", "", paste(deparse(substitute(sigma_group_arg)), collapse = ""))
      sigma_group_arg <- sigma_group_arg
      sigma_group_arg <- as.list(sigma_group_arg)
      names(sigma_group_arg) <- 'groupvar'
    } else if (is.character(substitute(sigma_group_arg))) {
      sigma_group_arg <- sigma_group_arg
      sigma_group_arg <- as.list(sigma_group_arg)
      names(sigma_group_arg) <- 'groupvar'
    }
  }
  if (any(grepl("^list", gsub("\\s", "",
                              paste(
                                deparse(substitute(sigma_group_arg)),
                                collapse = ""
                              )))) &
      length(strsplit(gsub("\\s", "",
                           paste(
                             deparse(substitute(sigma_group_arg)),
                             collapse = ""
                           )), ",")[[1]]) == 1) {
    if (!is.null(gsub("\\s", "", paste(deparse(
      substitute(sigma_group_arg)
    ),
    collapse = "")))) {
      if (is.language(substitute(sigma_group_arg))) {
        ttt <- gsub("\\s", "", paste(deparse(substitute(sigma_group_arg)),
                                     collapse = ""))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "T" |
            temp == "TRUE" |
            temp == "F" |
            temp == "FALSE") {
          stop(
            "sigma_group_arg should be either NULL or a character",
            " denoting the group idetifier"
          )
        }
        if (temp == "") {
          stop("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("groupvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        sigma_group_arg <- list_to_quoted_if_not(ttt)
      } else if (grepl("^list", sigma_group_arg)) {
        ttt <- deparse_0(as.name(substitute(sigma_group_arg)))
        temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
        if (temp == "") {
          stop("empty list")
        }
        if (length(strsplit(temp, "=")[[1]]) == 1) {
          ttt <- gsub(strsplit(temp, "=")[[1]],
                      paste0("groupvar=", strsplit(temp, "=")[[1]]),
                      ttt)
        }
        sigma_group_arg <- list_to_quoted_if_not(ttt)
        for (sigma_group_argi in 1:length(sigma_group_arg)) {
          if (!is.null(sigma_group_arg[[sigma_group_argi]])) {
            sigma_group_arg[[sigma_group_argi]] <- gsub("'", "", sigma_group_arg[[sigma_group_argi]])
          }
        }
      } else {
        if (!is.null(sigma_group_arg)) {
          if (!grepl("^list", gsub("\\s", "",
                                   paste(
                                     deparse(substitute(sigma_group_arg)),
                                     collapse = ""
                                   )))) {
            sigma_group_arg <- sigma_group_arg
            sigma_group_arg <- as.list(sigma_group_arg)
            names(sigma_group_arg) <- 'groupvar'
          }
        } else if (is.null(sigma_group_arg)) {
          sigma_group_arg <- as.list(sigma_group_arg)
        }
      }
    } else if (is.null(sigma_group_arg)) {
      sigma_group_arg <- list()
      sigma_group_arg$groupvar <- NULL
    }
  }
  if (any(grepl("^list", gsub("\\s", "",
                              paste(
                                deparse(substitute(sigma_group_arg)),
                                collapse = ""
                              )))) &
      length(strsplit(gsub("\\s", "",
                           paste(
                             deparse(substitute(sigma_group_arg)),
                             collapse = ""
                           )), ",")[[1]]) > 1) {
    ttt <-
      gsub("\\s", "", paste(deparse(substitute(sigma_group_arg)), collapse = ""))
    temp <- sub("\\).*", "", sub(".*\\(", "", ttt))
    temp <- gsub("\\s", "", temp)
    if (!grepl("^groupvar=", temp[1])) {
      temp[1] <- paste0("groupvar=", temp[1])
    }
    temp <- paste(temp, collapse = ",")
    temp <- paste0("list(", temp, ")")
    sigma_group_arg <- list_to_quoted_if_not(temp)
  }
  if (length(sigma_group_arg) == 0) {
    sigma_group_arg <- list()
    sigma_group_arg$groupvar <- NULL
  }
  if (!is.null(sigma_group_arg$groupvar) &
      !is.character(sigma_group_arg$groupvar)) {
    stop("sigma_group_arg should be either NULL or a character",
         " denoting the group idetifier")
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # If not already specified by the user, add default values to the
  # univariate_by, multivariate, and group_arg arguments
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    univariate_by$by <- gsub("\\s", "", univariate_by$by)
  }
  if (identical(univariate_by$by, character(0))) {
    univariate_by$by <- NA
  }
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    if (univariate_by$by == "" |
        univariate_by$by == FALSE | is.null(univariate_by$by)) {
      univariate_by$by <- NA
    }
    if (univariate_by$by == TRUE) {
      stop(
        "For univeriate-by-subgroup model fitting (via univariate_by argument)",
        "\n ",
        "argument 'by' should be a variable name, '', NULL, or FALSE"
      )
    }
  }
  
  
  
  if (multivariate$mvar &
      !(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    stop(
      "You have set multivariate as TRUE and also specified ",
      "\n ",
      " univeriate-by-subgroup model (see univariate_by argument)",
      "\n ",
      " Please specify either multivariate or univariate_by argument"
    )
  }
  
  
  if (is.symbol(arguments[["y"]]) |
      is.character(arguments[["y"]])) {
    nys <- length(arguments[["y"]])
  } else {
    nys <- length(arguments[["y"]]) - 1
  }
  
  
  if (multivariate$mvar & nys == 1) {
    stop(
      "You have set multivariate as TRUE but provided only one outcome ",
      "\n ",
      " Please set y as list or vector of multiple outcomes such as ",
      "\n ",
      " list(outcome1, outcome2) or y = c(outcome1, outcome2)"
    )
  }
  
  if (!multivariate$mvar & nys > 1) {
    stop(
      "You have set multivariate as FALSE but provided more than one outcome",
      "\n ",
      " Please set y as a symbol / list / vector of single outcome such as",
      "\n ",
      " y = outcome, y = list(outcome1) or y = c(outcome1)"
    )
  }
  
  if (!(is.na(univariate_by$by) |
        univariate_by$by == "NA") & nys > 1) {
    stop(
      "You have specified univariate_by model for ",
      univariate_by$by,
      "for which ",
      "\n ",
      " only one outcome varibale should be specified but have provided ",
      nys,
      " outcomes",
      "\n ",
      " Please set y as a symbol / list / vector of single outcome ",
      "\n ",
      " such as y = outcome, y = list(outcome1) or y = c(outcome1)"
    )
  }
  
  if (multivariate$mvar) {
    if (is.null(multivariate$cor))
      multivariate$cor <- "un"
    if (is.null(multivariate$rescor))
      multivariate$rescor <- TRUE
  }
  if (!multivariate$mvar) {
    if (is.null(multivariate$cor))
      multivariate$cor <- "un"
    if (is.null(multivariate$rescor))
      multivariate$rescor <- TRUE
  }
  
  
  
  if (is.na(univariate_by$by)) {
    if (is.null(univariate_by$cor))
      univariate_by$cor <- "un"
  }
  if (!is.na(univariate_by$by)) {
    if (is.null(univariate_by$cor))
      univariate_by$cor <- "un"
  }
  
  
  
  if (is.null(group_arg$groupvar))
    group_arg$groupvar <- NULL
  if (is.null(group_arg$by))
    group_arg$by <- NULL
  if (is.null(group_arg$cor))
    group_arg$cor <- "un"
  if (is.null(group_arg$dist))
    group_arg$dist <- "gaussian"
  
  
  
  
  
  if (is.null(sigma_group_arg$groupvar))
    sigma_group_arg$groupvar <- NULL
  if (is.null(sigma_group_arg$by))
    sigma_group_arg$by <- NULL
  if (is.null(sigma_group_arg$cor))
    sigma_group_arg$cor <- "un"
  if (is.null(sigma_group_arg$dist))
    sigma_group_arg$dist <- "gaussian"
  
  
  
  
  
  
  multivariate$verbose <-
    univariate_by$verbose <- group_arg$verbose <- verbose
  
  
  sigma_group_arg$verbose <- verbose
  
  
  # Temporary placeholder for the number of outcomes when fitting
  # univariate-by-subgroup model
  
  if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
    temp_ <- univariate_by$by
    if (!temp_ %in% colnames(data)) {
      stop(
        paste(
          "\nvariable",
          temp_,
          "used for setting univariate_by-univariate submodels is missing"
        )
      )
    }
    if (!is.factor(data[[temp_]])) {
      stop(temp_, "should be a factor variable")
    }
    nlevtemp_ <- nlevels(data[[temp_]])
    nys <- nlevtemp_
  }
  
  
  
  
  # Perform checks and set-up the 'to convert arguments' 
  
  to_list_if_not <- function(.x, nys, arguments, ...) {
    if (nys == 1) {
      if (!is.symbol(arguments[[.x]]) & !is.character(arguments[[.x]])) {
        arguments[[.x]] <- deparse_0(arguments[[.x]])
      } else {
        arguments[[.x]] <- arguments[[.x]]
      }
      if (is.symbol(arguments[[.x]]) &
          !is.character(arguments[[.x]])) {
        arguments[[.x]] <- deparse_0(arguments[[.x]])
      } else {
        arguments[[.x]] <- arguments[[.x]]
      }
      if (!is.character(.x)) {
        .xx <- eval(parse(text = .x))
      } else {
        .xx <- .x
      }
    }
    if (nys > 1) {
      .xx <- .x
    }
    if (!is.character(.xx) & !is.list(.xx)) {
      .xx <- deparse_0(.xx)
    } else {
      .xx <- .xx
    }
    . <- lapply(.xx, function(x)
      if (is.list(x))
        x <- x
      else
        x <- list(x)[[1]])
    assign(.x, ., envir = parent.frame())
  }
  
  eval_c_list_args <- function(.x, nys, arguments, ...) {
    if (is.language(arguments[[.x]]) &
        (strsplit(deparse_0(arguments[[.x]]), "\\(")[[1]][1] !=
         "c" &
         strsplit(deparse_0(arguments[[.x]]), "\\(")[[1]][1] != "list")) {
      arguments[[.x]] <- deparse(arguments[[.x]])
    } else {
      arguments[[.x]] <- (arguments[[.x]])
    }
    if (is.logical(arguments[[.x]])) {
      arguments[[.x]] <- deparse_0(arguments[[.x]])
    } else {
      arguments[[.x]] <- arguments[[.x]]
    }
    
    if (is.numeric(arguments[[.x]])) {
      arguments[[.x]] <- deparse_0(arguments[[.x]])
    } else {
      arguments[[.x]] <- (arguments[[.x]])
    }
    .xo <- .x
    .x <- arguments[[.x]]
    fun_ <- function(.x) {
      if (!is.character(.x))
        .x <- deparse_0(.x)
      else
        .x <- .x
      .x <- gsub("\\s", "", .x)
    }
    if (is.symbol(arguments[[.xo]]))
      .x <- deparse_0(.x)
    else
      .x <- .x
    if (is.symbol(arguments[[.xo]]))
      args_s <- mapply(fun_, .x)
    if (!is.symbol(arguments[[.xo]]))
      args_s <- mapply(fun_, .x)[-1]
    if (is.character(arguments[[.xo]]))
      args_s <- mapply(fun_, .x)
    if (length(args_s) > 1)
      args_s <- mapply(fun_, .x)[-1]
    attr(args_s, "names") <- NULL
    if (length(args_s) < nys)
      args_s <- rep(args_s, nys)
    if (length(args_s) > nys)
      args_s <- args_s[1:nys]
    assign(paste0(.xo, "s"), args_s, envir = parent.frame())
  }
  
  
  getArgNames <-
    function(value)
      formalArgs(deparse_0(substitute(value)[[1]]))
  
  convert_to_list <- getArgNames(bsitar())
  
  
  
  for (ip in convert_to_list) {
    if (grepl("_init_", ip)) {
      err. <- FALSE
      tryCatch(
        expr = {
          out <- suppressWarnings(ept(ip))
        },
        error = function(e) {
          err. <<- TRUE
        }
      )
      if (!err.) {
        if (length(out) > 1 & !is.list(out)) {
          stop(
            "Initials specified as vector [e.g, c(1, 2)] but must be a list, ",
            "\n ",
            " Note, initials can also be specified by using a single character",
            "\n ",
            " such as 0, random, or an object defined in the init_data",
            "\n ",
            " please check the following init arg: ",
            ip
          )
        }
      }
    }
  }
  
  
  # Convert arguments to the required format for setting sub-options for
  # the univariate-by-subgroup and multivariate model fitting
  
  single_args <- c(
    "data",
    "group_arg",
    "sigma_group_arg",
    "univariate_by",
    "multivariate",
    "prior_data",
    "init_data",
    "init_custom",
    "jitter_init_beta",
    "jitter_init_sd",
    "jitter_init_cor",
    "expose_function",
    "verbose",
    "normalize",
    "seed",
    "brms_arguments",
    "get_stancode",
    "get_standata",
    "get_priors",
    "get_set_priors",
    "validate_priors",
    "set_self_priors",
    "set_replace_priors",
    "outliers",
    "..."
  )
  
  
  
  for (i in convert_to_list) {
    if (!i %in% single_args) {
      to_list_if_not(i, nys, arguments)
    }
  }
  
  for (i in convert_to_list) {
    if (!i %in% single_args) {
      eval_c_list_args(i, nys, arguments)
    }
  }
  
  
  less_args <- extra_args <- c()
  outcomes_l <- paste0(" (", paste(ys, collapse = ", "), ")")
  for (i in convert_to_list) {
    .x <- i
    if (is.call(arguments[[.x]]))
      nl <- length(arguments[[.x]]) - 1
    if (!is.call(arguments[[.x]]))
      nl <- length(arguments[[.x]])
    if (nl > nys)
      extra_args <- c(extra_args, .x)
    if (nl < nys)
      less_args <- c(less_args, .x)
  }
  
  
  # Prepare data for 'bsitar'
  
  if (verbose) {
    setmsgtxt <- paste0("\n Preparing data")
    if (displayit == 'msg') {
      message(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolh
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  
  # org.ycall <- ys[1]
  
  if(is.list(xfuns) & length(xfuns) == 0) {
    xfuns <- rep('NULL', length(ys))
  }
  if(is.list(yfuns) & length(yfuns) == 0) {
    yfuns <- rep('NULL', length(ys))
  }
  
  # This data.org.in as specified will be saved in model_info
  
  if(!is.null(outliers)) {
    if(is.null(outliers$remove))    outliers$remove <- TRUE
    if(is.null(outliers$icode))     outliers$icode <- c(4,5,6)
    if(is.null(outliers$limit))     outliers$limit <- 5
    if(is.null(outliers$velpower))  outliers$velpower <- 0.5
    if(is.null(outliers$lag))       outliers$lag <- 1
    if(is.null(outliers$linearise)) outliers$linearise <- FALSE
    if(is.null(outliers$verbose))   outliers$verbose <- FALSE
  }
  
  data.org.in <- data
  
  data <- prepare_data(data = data,
                       x = xs, 
                       y = ys, 
                       id = ids,
                       uvarby = univariate_by$by, 
                       mvar = multivariate$mvar,
                       xfuns = xfuns, 
                       yfuns = yfuns,
                       outliers = outliers)
  
  ys <- attr(data, "ys")
  subindicators <- attr(data, "subindicators")
  
  
  if (univariate_by$verbose) {
    resvcts_ <- levels(data[[univariate_by$by]])
    resvcts <- paste0(resvcts_, collapse = " ")
    setmsgtxt <- paste0(
      "\n For univariate-by-subgroup model fitting for variable '",
      uvarby,
      "'",
      " (specified via 'univariate_by' argument)",
      "\n ",
      resvcts,
      " response vectors created based on the factor levels",
      "\n\n ",
      "Please check corresponding arguments list.",
      " E.g, df = list(4, 5) denotes that\n df = 4 is for ",
      resvcts_[1],
      ", and  df = 5 is for ",
      resvcts_[2],
      " (and similalry knots, priors, initials etc)",
      "\n\n ",
      "If it does't correspond correctly, then either reverse the list ",
      "arguments\n such as df = list(5, 4),",
      " or else reverse sort the order of factor levels"
    )
    if (displayit == 'msg') {
      message(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolb
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  
  
  # yfuns <<- yfuns
  # ys <<- ys
  
  # if(multivariate$mvar) {
  #   if(is.list(yfuns) & length(yfuns) == 0) {
  #     yfuns <- rep('NULL', length(ys))
  #   }
  #   for(myfunsi in 1:length(ys)) {
  #     mysi <- ys[[myfunsi]]
  #     myfunsi <- yfuns[[myfunsi]]
  #     if(grepl('.Primitive', myfunsi, fixed = T) & 
  #        grepl('log', myfunsi, fixed = T)) {
  #       myfunsi <- 'log'
  #     }
  #     if(grepl('.Primitive', myfunsi, fixed = T) & 
  #        grepl('sqrt', myfunsi, fixed = T)) {
  #       myfunsi <- 'sqrt'
  #     }
  #     if(myfunsi == 'log') data[[mysi]] <- log(data[[mysi]])
  #     if(myfunsi == 'sqrt') data[[mysi]] <- sqrt(data[[mysi]])
  #   }
  # }
  
  
  
  
  # print(head(data))
  # stop()
  
  # print(terms_rhs)
  # print(terms_rhss)
  
  
  # Initiate loop over outcome(s) 
  # First, create empty lists, vector etc. to collect these elements
  
  dataout <- priorlist <- NULL
  
  bflist <- list()
  bflist <- initialslist <- initialslist_s <- initsilist <- bflist
  blanketinitslist <- prior_stanvarlist <- auxillary_stanvarlist <- bflist
  
  funlist <- c()
  xoffsetvaluelist <- xoffsetnamelist <- knotsvaluelist <- funlist
  knotsnamelist <- spfun_collect <- xfunvaluelist <- xfunnamelist <- funlist
  yfunvaluelist <- yfunnamelist <- yyfunvaluelist <- yyfunnamelist <- funlist
  xxfunvaluelist <- xxfunnamelist <- fixedvaluelist <- fixednamelist <- funlist
  randomvaluelist <- randomnamelist <- groupvarvaluelist <- funlist
  yvarvaluelist <- ynamelist <- covvaluelist <- covnamelist <- funlist
  groupvarnamelist <- xvarvaluelist <- xnamelist <- funlist
  hierarchicalvarnamelist <- hierarchicalvarvaluelist <- funlist
  
  sigmacovnamelist <- sigmacovvaluelist <- funlist
  
  sigma_groupvarnamelist <- sigma_groupvarvaluelist <- funlist
  sigma_hierarchicalvarnamelist <- sigma_hierarchicalvarvaluelist <- funlist
  
  
  
  
  
  # Start loop over the outcome(s)
  
  for (ii in 1:length(ys)) {
    if (nys > 1)
      resp <- ys[ii]
    else
      resp <- ""
    subindicatorsi <- subindicators[ii]
    
    for (i in convert_to_list) {
      if (!i %in% single_args) {
        assign(paste0(i, "s", "i"), eval(parse(text = paste0(i, "s")))[ii])
      }
    }
    
    if (is.null(group_arg$groupvar))
      group_arg$groupvar <- idsi
    
    
    if (is.null(sigma_group_arg$groupvar))
      sigma_group_arg$groupvar <- idsi
    
    
    if (!is.numeric(ept(dfsi)) & !is.numeric(ept(knotssi))) {
      stop("either df or knots must be specified")
    }
    if (is.numeric(ept(dfsi)) & is.numeric(ept(knotssi))) {
      stop("both df and knots specified. Specify one of them\n")
    }
    
    if (!grepl("a", fixedsi, fixed = T) &
        grepl("a", randomsi, fixed = T)) {
      stop(
        "Parameter 'a' is missing in the fixed effects part of the model ",
        "\n ",
        " but specified in the random effects part of the model ",
        "\n ",
        " Either include 'a' in the fixed effects too or else ",
        "\n ",
        " remove it from the random effect part of the model"
      )
    }
    if (!grepl("b", fixedsi, fixed = T) &
        grepl("b", randomsi, fixed = T)) {
      stop(
        "Parameter 'b' is missing in the fixed effects part of the model ",
        "\n ",
        " but specified in the random effects part of the model ",
        "\n ",
        " Either include 'b' in the fixed effects too or else ",
        "\n ",
        " remove it from the random effect part of the model"
      )
    }
    if (!grepl("c", fixedsi, fixed = T) &
        grepl("c", randomsi, fixed = T)) {
      stop(
        "Parameter 'c' is missing in the fixed effects part of the model ",
        "\n ",
        " but specified in the random effects part of the model ",
        "\n ",
        " Either include 'c' in the fixed effects too or else ",
        "\n ",
        " remove it from the random effect part of the model"
      )
    }
    
    # For some reasons, 'sitar' (Tim Cole) allows the 'd' parameter to be 
    # random only. In fact for df > 1, it forces 'd' to be random parameter only
    
    if (!match_sitar_d_form) {
      if (!grepl("d", fixedsi, fixed = T) &
          grepl("d", randomsi, fixed = T)) {
        stop(
          "Parameter 'd' is missing in the fixed effects part of the model ",
          "\n ",
          " but specified in the random effects part of the model ",
          "\n ",
          " Either include 'd' in the fixed effects too or else ",
          "\n ",
          " remove it from the random effect part of the model"
        )
      }
    }
    
    # Covariate not allowed when matching to sitar 'd' form
    
    if (match_sitar_d_form) {
      if ((grepl("d", fixedsi, fixed = T) |
           grepl("d", randomsi, fixed = T)) &
          (!grepl("^~1$", d_formulasi) |
           !grepl("^~1$", d_formula_grsi))) {
        stop(
          "Parameter 'd' is missing in the fixed effects part of the model ",
          "\n ",
          " and is specified only in the random effects part of the model ",
          "\n ",
          " (This is to match with the 'sitar' package's formulation)",
          "\n ",
          " For this formulation (i.e., 'd' is missing in the fixed effects)",
          "\n ",
          " covariate(s) are not allowed"
        )
      }
    }
    
    
    # Add missing parameters to the dpar_formula
    
    if (!is.null(dpar_formulasi)) {
      if (grepl("^1$", dpar_formulasi)) {
        dpar_formulasi <- paste0("lf(", "sigma", "~", dpar_formulasi, ")")
      } else if (grepl("^~1", dpar_formulasi)) { # if (grepl("^~1$", dpar_formulasi)) {
        dpar_formulasi <- paste0("lf(", "sigma", dpar_formulasi, ")")
      } else if (grepl("^sigma~1", dpar_formulasi)) { # f (grepl("^sigma~1$", dpar_formulasi)) {
        dpar_formulasi <- paste0("lf(", "", dpar_formulasi, ")")
      } else {
        dpar_formulasi <- dpar_formulasi
      }
      if (grepl("lf\\(", dpar_formulasi) |
          grepl("nlf\\(", dpar_formulasi)) {
        if (grepl("^lf\\(", dpar_formulasi) &
            !grepl("nlf\\(", dpar_formulasi)) {
          lf_list <- c('flist',
                       'dpar',
                       'resp',
                       'center',
                       'cmc',
                       'sparse',
                       'decomp') #
        } else if (!grepl("^lf\\(", dpar_formulasi) &
                   grepl("^nlf\\(", dpar_formulasi)) {
          lf_list <- c('flist', 'dpar', 'resp', 'loop ') 
        }
        lf_list_c <- c()
        for (lf_listi in lf_list) {
          if (!grepl(lf_listi, dpar_formulasi)) {
            if (lf_listi == 'center') {
              o. <- paste0(lf_listi, "=", 'TRUE')
            } else if (lf_listi == 'cmc') {
              o. <- paste0(lf_listi, "=", 'TRUE')
            } else if (lf_listi == 'resp') {
              if (nys > 1) {
                o. <- paste0(lf_listi, "=", paste0("'", ysi, "'"))
                # o. <- paste0(lf_listi, "=", 'NULL')
              } else {
                o. <- paste0(lf_listi, "=", 'NULL')
              }
            } else {
              o. <- paste0(lf_listi, "=", 'NULL')
            }
            lf_list_c <- c(lf_list_c, o.)
          }
        }
        lf_list_c <- paste(lf_list_c, collapse = ",")
        if (lf_list_c != "")
          lf_list_c <- paste0(",", lf_list_c)
        dpar_formulasi <- gsub(")$", lf_list_c, dpar_formulasi)
        dpar_formulasi <- paste0(dpar_formulasi, ")")
      }
    }
    
    # print(dpar_formulasi)
    
    
    # Check for higher level model and update level 2 random formula
    
    f_checks_gr_gr_str <- function(a, b) {
      gr_st_id <- sub(".*\\|", "", a) 
      a_ <- paste0("'", deparse_0(substitute(a)), "'")
      b_ <- paste0("'", deparse_0(substitute(b)), "'")
      b_out <- NULL
      if(is.null(b[[1]])) {
        if(grepl(":", gr_st_id, fixed = T) | grepl("/", gr_st_id, fixed = T)) {
          stop("Models beyound two levels of hierarchy are not supported yet",
               "\n ",
               "An alternative to argument ", a_, " is to use ",
               "\n ",
               "argument ", b_, " to directly pass on the ", 
               "\n ",
               "random formual to the brms and then either accept",
               "\n ",
               "default priors placed by the brms for those varinace covarinace", 
               "\n ",
               "or else use get_prios to place priors manually and the pass them",
               "\n ",
               "to the bsitar by using argument 'set_self_priors'"
          )
        }
      } else if(!is.null(b[[1]])) {
        b_out <- b
      }
      b_out
    }
    
    a_fcgs_out <- f_checks_gr_gr_str(a_formula_grsi, a_formula_gr_strsi)
    b_fcgs_out <- f_checks_gr_gr_str(b_formula_grsi, b_formula_gr_strsi)
    c_fcgs_out <- f_checks_gr_gr_str(c_formula_grsi, c_formula_gr_strsi)
    d_fcgs_out <- f_checks_gr_gr_str(d_formula_grsi, d_formula_gr_strsi)
    
    
    
    # First, if sigma_formula_gr_strsi not NULL but sigma_formula_grsi NULL
    # Then set sigma_formula_grsi to ~1 because then only first part of the 
    # sigma_formula_gr_strsi (i.e., before first + ) will be copied to the 
    # sigma_formula_grsi
    
    
    if(sigma_formula_gr_strsi != 'NULL') {
      if(sigma_formula_grsi == 'NULL') {
        sigma_formula_grsi <- "~1"
      } else {
        sigma_formula_grsi <- sigma_formula_grsi
      }
    } else {
      sigma_formula_grsi <- sigma_formula_grsi
    }
    
    
    if(sigma_formula_gr_strsi != 'NULL') {
      if(!grepl("^~", sigma_formula_gr_strsi)) {
        sigma_formula_gr_strsi <- paste0("~", sigma_formula_gr_strsi)
      }
    }
    if(is.null(sigma_formula_gr_strsi[[1]])) {
      sigma_formula_gr_strsi <- 'NULL'
    }
    
    if(sigma_formula_grsi != 'NULL') {
      if(!grepl("^~", sigma_formula_grsi)) {
        sigma_formula_grsi <- paste0("~", sigma_formula_grsi)
      }
    }
    if(is.null(sigma_formula_grsi[[1]])) {
      sigma_formula_grsi <- 'NULL'
    }
    
    
    # print(sigma_formula_grsi)
    # print(sigma_formula_gr_strsi)
    # stop()
    
    
    sigma_fcgs_out <- f_checks_gr_gr_str(sigma_formula_grsi, sigma_formula_gr_strsi)
    
    if(!is.null(a_fcgs_out)) {
      if(a_formula_grsi == "~1" & !is.null(a_formula_gr_strsi[[1]])) {
        a_formula_grsi <- strsplit(a_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    if(!is.null(b_fcgs_out)) {
      if(b_formula_grsi == "~1" & !is.null(b_formula_gr_strsi[[1]])) {
        b_formula_grsi <- strsplit(b_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    if(!is.null(c_fcgs_out)) {
      if(c_formula_grsi == "~1" & !is.null(c_formula_gr_strsi[[1]])) {
        c_formula_grsi <- strsplit(c_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    if(!is.null(d_fcgs_out)) {
      if(!is.null(d_formula_grsi[[1]]) & !is.null(d_formula_gr_strsi[[1]])) {
        if(d_formula_grsi == "~1" & !is.null(d_formula_gr_strsi[[1]])) {
          d_formula_grsi <- strsplit(d_formula_gr_strsi, "+(", fixed = T)[[1]][1]
        }
      }
    }
    
    if(!is.null(sigma_fcgs_out) & sigma_fcgs_out != 'NULL') {
      if(sigma_formula_grsi == "~1" & !is.null(sigma_formula_gr_strsi[[1]])) {
        sigma_formula_grsi <- strsplit(sigma_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    a_formula_grsi <- gsub("[()]", "", a_formula_grsi)
    b_formula_grsi <- gsub("[()]", "", b_formula_grsi)
    c_formula_grsi <- gsub("[()]", "", c_formula_grsi)
    if(!is.null(d_formula_grsi)) d_formula_grsi <- gsub("[()]", "", d_formula_grsi)
    
    sigma_formula_grsi <- gsub("[()]", "", sigma_formula_grsi)
    
    set_higher_levels <- TRUE
    if(is.null(a_fcgs_out) & 
       is.null(b_fcgs_out) & 
       is.null(c_fcgs_out) & 
       is.null(d_fcgs_out)) {
      set_higher_levels <- FALSE
    }
    
    sigma_set_higher_levels <- TRUE
    if(is.null(sigma_fcgs_out) | sigma_fcgs_out == 'NULL') {
      sigma_set_higher_levels <- FALSE
    }
    
    
    # Add intercept ~ 1 if missing
    
    check_formuals <-
      c(
        "a_formulasi",
        "b_formulasi",
        "c_formulasi",
        "d_formulasi",
        "s_formulasi",
        "a_formula_grsi",
        "b_formula_grsi",
        "c_formula_grsi",
        "d_formula_grsi",
        "sigma_formulasi",
        "sigma_formula_grsi"
      )
    
    for (check_formualsi in check_formuals) {
      if (!grepl("~1", ept(check_formualsi)) &
          !grepl("~0", ept(check_formualsi))) {
        check_formualsi_with1 <-
          gsub("^~", "~1+", ept(check_formualsi), fixed = F)
        assign(check_formualsi, check_formualsi_with1)
      }
    }
    
    
    if (is.null(sigma_formula_gr_strsi[[1]][1]) |
        sigma_formula_gr_strsi == "NULL") {
      sigma_formula_gr_strsi <- NULL
    }
    
    if (is.null(sigma_formula_grsi[[1]][1]) |
        sigma_formula_grsi == "NULL") {
      sigma_formula_grsi <- NULL
    }
    
    
    
    
    if (is.null(dpar_formulasi[[1]][1]) |
        dpar_formulasi == "NULL") {
      dpar_formulasi <- NULL
    }
    
    if (is.null(autocor_formulasi[[1]][1]) |
        autocor_formulasi == "NULL") {
      autocor_formi <- NULL
    } else {
      autocor_formi <- autocor_formulasi
    }
    
    
    if (is.null(familysi[[1]][1]) |
        familysi == "NULL") {
      familysi <- NULL
    }
    
    if (!is.null(familysi)) {
      familysi <- list_to_quoted_if_not_si(familysi)
    }
    
    # lf edited 
    if (!is.null(dpar_formulasi)) {
      if (grepl("^lf\\(", dpar_formulasi) |
          grepl("^nlf\\(", dpar_formulasi)) {
        #  dpar_formulasi <- list_to_quoted_if_not_si_lf(dpar_formulasi)
      } else {
        dpar_formulasi <- dpar_formulasi
      }
    }
    
    # print(dpar_formulasi)
    
    
    N_J_all <- length(unique(data[[idsi]]))
    
    if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
      sortbylayer <- NA
      data <- data %>%
        dplyr::mutate(sortbylayer =
                        forcats::fct_relevel(!!as.name(univariate_by$by),
                                             (levels(
                                               !!as.name(univariate_by$by)
                                             )))) %>%
        dplyr::arrange(sortbylayer) %>%
        dplyr::mutate(!!as.name(idsi) := factor(!!as.name(idsi),
                                                levels = 
                                                  unique(!!as.name(idsi)))) %>% 
        dplyr::select(-sortbylayer)
      
      
      datai <- data %>%
        dplyr::filter(eval(parse(text = subindicatorsi)) == 1) %>%
        droplevels()
      if (!subindicatorsi %in% colnames(datai)) {
        stop("variable ", subindicatorsi, " not in the dataframe")
      }
      if (!xsi %in% colnames(datai))
        stop("variable ", xsi, " not in the dataframe")
      if (!idsi %in% colnames(datai))
        stop("variable ", idsi, " not in the dataframe")
    }
    
    if ((is.na(univariate_by$by) | univariate_by$by == "NA")) {
      datai <- data %>%
        droplevels()
      if (!ysi %in% colnames(datai))
        stop("variable ", ysi, " not in the dataframe")
      if (!xsi %in% colnames(datai))
        stop("variable ", xsi, " not in the dataframe")
      if (!idsi %in% colnames(datai))
        stop("variable ", idsi, " not in the dataframe")
    }
    
    
    
    
    if (!is.null(xfunsi[[1]][1]) & xfunsi != "NULL") {
      if (xfunsi != "log" & xfunsi != "sqrt") {
        stop("only log and sqrt options allowed for xfun argument")
      }
    }
    
    
    
    
    if (!is.null(xfunsi[[1]][1]) & xfunsi != "NULL") {
      if (xfunsi == "log") {
        datai[[xsi]] <- log(datai[[xsi]])
      } else if (xfunsi == "sqrt") {
        datai[[xsi]] <- sqrt(datai[[xsi]])
      } else {
        stop("only log and sqrt options allowed for xfun argument")
      }
    }
    
    
    # for multivariate, yvar is transformed before initiating loop over ys
    # Changed now, yvar transformations for univariate, univariate_by, 
    # multivariate, done within the prepare_data
    
    # if(!multivariate$mvar) {
    #   if (!is.null(yfunsi[[1]][1]) & yfunsi != "NULL") {
    #     if (yfunsi == "log") {
    #       datai[[ysi]] <- log(datai[[ysi]])
    #     } else if (yfunsi == "sqrt") {
    #       datai[[ysi]] <- sqrt(datai[[ysi]])
    #     } else {
    #       stop("only log or sqrt tranformation allowed for yfun argument")
    #     }
    #   }
    # }
    # 
    # 
    
    
    
    
    
    
    
    
    
    
    
    
    
    if (!is.null(xfunsi[[1]][1]) & xfunsi != "NULL") {
      xfunvalue <- xfunsi
    } else {
      xfunvalue <- NULL
    }
    
    if (!is.null(yfunsi[[1]][1]) & yfunsi != "NULL") {
      yfunvalue <- yfunsi
    } else {
      yfunvalue <- NULL
    }
    
    
    if (nys == 1) {
      xfun_name <- "xfun"
      yfun_name <- "yfun"
      xxfun_name <- "xvar_xfun"
      yyfun_name <- "yvar_yfun"
    } else if (nys > 1) {
      xfun_name <- paste0("xfun", "_", ysi)
      yfun_name <- paste0("yfun", "_", ysi)
      xxfun_name <- paste0("xvar_xfun", "_", ysi)
      yyfun_name <- paste0("yvar_yfun", "_", ysi)
    }
    
    
    xfunvaluelist[[ii]] <- xfunvalue
    xfunnamelist[[ii]] <- xfun_name
    
    yfunvaluelist[[ii]] <- yfunvalue
    yfunnamelist[[ii]] <- yfun_name
    
    if (!is.null(xfunsi[[1]][1]) & xfunsi != "NULL") {
      xxfunvaluelist[[ii]] <- paste0(xfunsi, "(", xsi, ")")
    } else {
      xxfunvaluelist[[ii]] <- NULL
    }
    
    
    if (!is.null(yfunsi[[1]][1]) & yfunsi != "NULL") {
      yyfunvaluelist[[ii]] <- paste0(yfunsi, "(", ysi, ")")
    } else {
      yyfunvaluelist[[ii]] <- NULL
    }
    
    xxfunnamelist[[ii]] <- xxfun_name
    yyfunnamelist[[ii]] <- xxfun_name
    
    
    gkn <- function(x, df, bounds) {
      c(min(x) - bounds * (max(x) - min(x)),
        quantile(x, (1:(df - 1)) / df),
        max(x) +
          bounds * (max(x) - min(x)))
    }
    
    if (is.numeric(ept(knotssi))) {
      knots <- ept(knotssi)
    }
    if (is.numeric(ept(dfsi))) {
      knots <- (unname(gkn(datai[[xsi]], ept(dfsi), ept(boundsi))))
    }
    
    
    ########### For some reasons, 'sitar' (Tim Cole) allows random only 'd'
    ########### parameter In fact for df > 1, it forces 'd' to be random
    ########### parameter only
    
    if (match_sitar_d_form) {
      if (length(knots) > 2) {
        itemp <- strsplit(gsub("\\+", " ", fixedsi), " ")[[1]]
        itemp <- itemp[!grepl("d", itemp)]
        fixedsi <- paste(itemp, collapse = "+")
      }
    }
    
    
    ###########
    
    nabci <- length(strsplit(gsub("\\+", " ", fixedsi), " ")[[1]])
    nabcrei <-
      length(strsplit(gsub("\\+", " ", randomsi), " ")[[1]])
    
    
    # define function to evaluate offset and bstart arguments
    
    make_spline_matrix <- function(x, knots) {
      X <- x
      N <- length(X)
      nk <- length(knots)
      basis_evals <- matrix(0, N, nk - 1)
      basis_evals[, 1] <- X
      basis_evals[, 1] <- X
      Xx <- matrix(0, N, nk)
      km1 <- nk - 1
      j <- 1
      knot1 <- knots[1]
      knotnk <- knots[nk]
      knotnk1 <- knots[nk - 1]
      kd <- (knotnk - knot1) ^ (2)
      for (ia in 1:N) {
        for (ja in 1:nk) {
          Xx[ia, ja] <- ifelse(X[ia] - knots[ja] > 0, X[ia] - knots[ja], 0)
        }
      }
      while (j <= nk - 2) {
        jp1 <- j + 1
        basis_evals[, jp1] <-
          (
            Xx[, j] ^ 3 - (Xx[, km1] ^ 3) * (knots[nk] - knots[j]) /
              (knots[nk] - knots[km1]) + (Xx[, nk] ^ 3) *
              (knots[km1] - knots[j]) / (knots[nk] - knots[km1])
          ) /
          (knots[nk] - knots[1]) ^ 2
        j <- j + 1
      }
      return(basis_evals)
    }
    
    
    
    eval_xoffset_bstart_args <-
      function(x, y, knots, data, eval_arg, xfunsi) {
        if (eval_arg == "mean") {
          eval_arg.o <- mean(data[[x]])
        } else if (eval_arg == "apv") {
          mat_s <- make_spline_matrix(data[[x]], knots)
          lmform <- as.formula(paste0(y, "~1+", "mat_s"))
          lmfit <- lm(lmform, data = data)
          eval_arg.o <- sitar::getPeak(data[[x]],
                                       predict(smooth.spline(data[[x]],
                                                             fitted(lmfit)),
                                               data[[x]], deriv = 1)$y)[1]
        } else {
          eval_arg.o <- ept(eval_arg)
        }
        return(as.numeric(eval_arg.o))
      }
    
    xoffset <-
      eval_xoffset_bstart_args(xsi, ysi, knots, datai, xoffsetsi, xfunsi)
    bstart <-
      eval_xoffset_bstart_args(xsi, ysi, knots, datai, bstartsi, xfunsi)
    bstart <- bstart - xoffset
    
    
    xoffset <- round(xoffset, 4)
    datai[[xsi]] <- datai[[xsi]] - xoffset
    knots <- knots - xoffset
    knots <- round(knots, 6)
    nknots <- length(knots)
    df <- length(knots) - 1
    
    mat_s <- make_spline_matrix(datai[[xsi]], knots)
    
    
    SplineFun_name <- "SplineFun"
    getX_name <- "getX"
    getKnots_name <- "getKnots"
    if (nys > 1) {
      spfncname <- paste0(ysi, "_", SplineFun_name)
      getxname <- paste0(ysi, "_", getX_name)
      getknotsname <- paste0(ysi, "_", getKnots_name)
    } else if (nys == 1) {
      spfncname <- SplineFun_name
      getxname <- getX_name
      getknotsname <- getKnots_name
    }
    
    
    spfun_collect <-
      c(spfun_collect, c(spfncname, paste0(spfncname, "_", c("d1",
                                                             "d2"))))
    
    internal_function_args_names <-
      c(
        "fixedsi",
        "randomsi",
        "spfncname",
        "getxname",
        "getknotsname",
        'match_sitar_d_form',
        'xfunsi',
        'yfunsi',
        'xoffset',
        'brms_arguments',
        "verbose"
      )
    
    internal_function_args <- list()
    internal_function_args <- mget(internal_function_args_names)
    
    if (verbose) {
      if (ii == 1) {
        setmsgtxt <- paste0("\n Preparing function")
        if (displayit == 'msg') {
          message(setmsgtxt)
        } else if (displayit == 'col') {
          col <- setcolh
          cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
        }
      }
    }
    
    
    
    funlist[ii] <-
      prepare_function(
        x = xsi,
        y = ysi,
        id = idsi,
        knots = knots,
        nknots = nknots,
        data = datai,
        internal_function_args = internal_function_args
      )
    
    
    
    
    
    internal_formula_args_names <-
      c(
        "a_formulasi",
        "b_formulasi",
        "c_formulasi",
        "d_formulasi",
        "s_formulasi",
        "a_formula_grsi",
        "b_formula_grsi",
        "c_formula_grsi",
        "d_formula_grsi",
        
        "terms_rhssi",
        
        "sigma_formulasi",
        "sigma_formula_grsi",
        
        "dpar_formulasi",
        "autocor_formi",
        "subindicatorsi",
        "fixedsi",
        "randomsi",
        "univariate_by",
        "multivariate",
        "group_arg",
        "sigma_group_arg",
        "df",
        "mat_s",
        "spfncname",
        'nys',
        'ysi',
        'familysi',
        'xfunsi',
        'xoffset',
        'match_sitar_d_form',
        
        "a_formula_gr_strsi",
        "b_formula_gr_strsi",
        "c_formula_gr_strsi",
        "d_formula_gr_strsi",
        "sigma_formula_gr_strsi",
        "set_higher_levels",
        "sigma_set_higher_levels",
        "verbose"
      )
    
    
    
    internal_formula_args <- list()
    internal_formula_args <- mget(internal_formula_args_names)
    
    if (verbose) {
      if (ii == 1) {
        setmsgtxt <- paste0("\n Preparing formula")
        if (displayit == 'msg') {
          message(setmsgtxt)
        } else if (displayit == 'col') {
          col <- setcolh
          cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
        }
      }
    }
    
    formula_bf <-
      prepare_formula(
        x = xsi,
        y = ysi,
        id = idsi,
        knots = knots,
        nknots = nknots,
        data = datai,
        internal_formula_args = internal_formula_args
      )
    
    
    
    list_out <- attr(formula_bf, "list_out")
    
    attributes(formula_bf) <- NULL
    
    eout <- list2env(list_out)
    for (eoutii in names(eout)) {
      assign(eoutii, eout[[eoutii]])
    }
    
    
    group_arg$groupvar <- group_arg_groupvar
    multivariate$rescor <- multivariate_rescor
    univariate_by$by <- univariate_by_by
    covariates_ <- covariates_
    covariates_sigma_ <- covariates_sigma_
    set_higher_levels <- set_higher_levels
    
    sigma_set_higher_levels <- sigma_set_higher_levels
    
    
    sigma_group_arg$groupvar <- sigma_group_arg_groupvar
    
    
    lm_val_list <-
      names(eout)[grep(pattern = "^lm_|^lme_", names(eout))]
    lm_val_list <- sort(lm_val_list)
    
    lm_val_list_not <-
      names(eout)[!names(eout) %in%
                    names(eout)[grep(pattern = "^lm_|^lme_", names(eout))]]
    lm_val_list_not <- sort(lm_val_list_not)
    
    
    
    cov_list_names <- ls()[grepl(pattern = "_cov", ls())]
    cov_list_names <-
      cov_list_names[!grepl(pattern = "_init_", cov_list_names)]
    cov_list_names <-
      cov_list_names[!grepl(pattern = "_prior_", cov_list_names)]
    cov_list_names <-
      cov_list_names[!grepl(pattern = "^lm_", cov_list_names)]
    cov_list_names <- sort(cov_list_names)
    
    bflist[[ii]] <- formula_bf
    
    
    
    ####
    ymean <- mean(datai[[ysi]])
    ymedian <- median(datai[[ysi]])
    ysd <- sd(datai[[ysi]])
    ymad <- mad(datai[[ysi]])
    xsd <- sd(datai[[xsi]])
    
    # TODO
    # acov_sd etc setting numeric but later can be worked out to infer from
    # some meaningful way - until then, these will not be used in anywhere
    # and are included here just as placeholders
    lm_a_cov_sd <- 10
    lm_b_cov_sd <- 1
    lm_c_cov_sd <- 0.5
    
    
    prior_data_internal_names <-
      c(
        lm_val_list,
        "ymean",
        "ymedian",
        "ysd",
        "ymad",
        "xsd",
        "lm_a_cov_sd",
        "lm_b_cov_sd",
        "lm_c_cov_sd",
        "bstart"
      )
    
    
    prior_args_internal_names <-
      c(
        lm_val_list_not,
        cov_list_names,
        "a_formulasi",
        "b_formulasi",
        "c_formulasi",
        "d_formulasi",
        "s_formulasi",
        "a_formula_grsi",
        "fixedsi",
        "b_formula_grsi",
        "c_formula_grsi",
        "d_formula_grsi",
        "sigma_formulasi",
        "sigma_formula_grsi",
        "sigma_formula_gr_strsi",
        "autocor_formi",
        "randomsi",
        "nabci",
        "nabcrei",
        "univariate_by",
        "multivariate",
        "group_arg",
        "sigma_group_arg",
        "initsi",
        "df",
        "idsi",
        "ys",
        "resp",
        "ii",
        "nys",
        "N_J_all",
        "dpar_formulasi",
        "normalize",
        "seed",
        'match_sitar_d_form',
        "verbose"
      )
    
    
    prior_data_internal <- list()
    prior_data_internal <- mget(prior_data_internal_names)
    
    
    prior_args_internal <- list()
    prior_args_internal <- mget(prior_args_internal_names)
    
    
    
    if (!is.null(prior_data[[1]])) {
      get_common_names_lists <-
        intersect(names(prior_data_internal), names(prior_data))
      ttt_n1 <- paste(names(prior_data_internal), collapse = ", ")
      ttt_nn2 <- paste(get_common_names_lists, collapse = ", ")
      if (!identical(get_common_names_lists, character(0))) {
        stop(
          "Names in prior_data list should not be following reserved names:",
          "\n ",
          ttt_n1,
          "\n ",
          " Please change the following name(s) ",
          "\n ",
          ttt_nn2
        )
      }
    }
    
    
    init_data_internal <- prior_data_internal
    init_args_internal <- prior_args_internal
    
    
    init_arguments <-
      list(
        a_init_beta = a_init_betasi,
        b_init_beta = b_init_betasi,
        c_init_beta = c_init_betasi,
        d_init_beta = d_init_betasi,
        s_init_beta = s_init_betasi,
        a_cov_init_beta = a_cov_init_betasi,
        b_cov_init_beta = b_cov_init_betasi,
        c_cov_init_beta = c_cov_init_betasi,
        d_cov_init_beta = d_cov_init_betasi,
        s_cov_init_beta = s_cov_init_betasi,
        a_init_sd = a_init_sdsi,
        b_init_sd = b_init_sdsi,
        c_init_sd = c_init_sdsi,
        d_init_sd = d_init_sdsi,
        a_cov_init_sd = a_cov_init_sdsi,
        b_cov_init_sd = b_cov_init_sdsi,
        c_cov_init_sd = c_cov_init_sdsi,
        d_cov_init_sd = d_cov_init_sdsi,
        
        sigma_init_beta = sigma_init_betasi,
        sigma_cov_init_beta = sigma_cov_init_betasi,
        sigma_init_sd = sigma_init_sdsi,
        sigma_cov_init_sd = sigma_cov_init_sdsi,
        
        rsd_init_sigma = rsd_init_sigmasi,
        dpar_init_sigma = dpar_init_sigmasi,
        dpar_cov_init_sigma = dpar_cov_init_sigmasi,
        autocor_init_acor = autocor_init_acorsi,
        gr_init_cor = gr_init_corsi,
        mvr_init_rescor = mvr_init_rescorsi,
        r_init_z = r_init_zsi
      )
    
    
    if (verbose) {
      if (ii == 1) {
        setmsgtxt <- paste0("\n Preparing priors and initials")
        if (displayit == 'msg') {
          message(setmsgtxt)
        } else if (displayit == 'col') {
          col <- setcolh
          cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
        }
      }
    }
    
    bpriors <-
      set_priors_initials(
        a_prior_beta = a_prior_betasi,
        b_prior_beta = b_prior_betasi,
        c_prior_beta = c_prior_betasi,
        d_prior_beta = d_prior_betasi,
        s_prior_beta = s_prior_betasi,
        a_cov_prior_beta = a_cov_prior_betasi,
        b_cov_prior_beta = b_cov_prior_betasi,
        c_cov_prior_beta = c_cov_prior_betasi,
        d_cov_prior_beta = d_cov_prior_betasi,
        s_cov_prior_beta = s_cov_prior_betasi,
        a_prior_sd = a_prior_sdsi,
        b_prior_sd = b_prior_sdsi,
        c_prior_sd = c_prior_sdsi,
        d_prior_sd = d_prior_sdsi,
        a_cov_prior_sd = a_cov_prior_sdsi,
        b_cov_prior_sd = b_cov_prior_sdsi,
        c_cov_prior_sd = c_cov_prior_sdsi,
        d_cov_prior_sd = d_cov_prior_sdsi,
        gr_prior_cor = gr_prior_corsi,
        sigma_prior_beta = sigma_prior_betasi, 
        sigma_cov_prior_beta = sigma_cov_prior_betasi,
        sigma_prior_sd = sigma_prior_sdsi,
        sigma_cov_prior_sd = sigma_cov_prior_sdsi,
        rsd_prior_sigma = rsd_prior_sigmasi,
        dpar_prior_sigma = dpar_prior_sigmasi,
        dpar_cov_prior_sigma = dpar_cov_prior_sigmasi,
        autocor_prior_acor = autocor_prior_acorsi,
        mvr_prior_rescor = mvr_prior_rescorsi,
        prior_data = prior_data,
        prior_data_internal = prior_data_internal,
        prior_args_internal = prior_args_internal,
        init_arguments = init_arguments,
        init_data = init_data,
        init_data_internal = init_data_internal,
        init_args_internal = init_args_internal
      )
    
    
    priorlist <- rbind(priorlist, bpriors)
    
    stanvar_priors <- attr(bpriors, "stanvars")
    
    stanvar_priors_names <- names(stanvar_priors)
    
    prior_stanvarlist[[ii]] <- stanvar_priors
    
    initials <- attr(bpriors, "initials")
    
    scode_auxillary <- attr(bpriors, "scode_auxillary")
    auxillary_stanvarlist[[ii]] <- scode_auxillary
    
    initialslist[[ii]] <- initials
    initialslist_s[[ii]] <- initsi
    
    if (initsi == "random") {
      blanketinits <- "random"
    } else if (initsi == "0") {
      blanketinits <- "0"
    } else {
      blanketinits <- "no"
    }
    
    blanketinitslist[[ii]] <- blanketinits
    
    
    if (nys == 1) {
      xoffset_name <- "xoffset"
      knots_name <- "knots"
      fixed_name <- "fixed"
      random_name <- "random"
      groupvar_name <- "groupvar"
      sigma_groupvar_name <- "sigma_groupvar"
      hierarchical_name <- "hierarchical"
      sigma_hierarchical_name <- "sigma_hierarchical"
      xvar_name <- "xvar"
      yvar_name <- "yvar"
      cov_name <- "cov"
      cov_name_sigma <- "cov_sigma"
    } else if (nys > 1) {
      xoffset_name <- paste0("xoffset", "_", ysi)
      knots_name <- paste0("knots", "_", ysi)
      fixed_name <- paste0("fixed", "_", ysi)
      random_name <- paste0("random", "_", ysi)
      groupvar_name <- paste0("groupvar", "_", ysi)
      sigma_groupvar_name <- paste0("sigma_groupvar", "_", ysi)
      hierarchical_name <- paste0("hierarchical", "_", ysi)
      sigma_hierarchical_name <- paste0("sigma_hierarchical", "_", ysi)
      xvar_name <- paste0("xvar", "_", ysi)
      yvar_name <- paste0("yvar", "_", ysi)
      cov_name <- paste0("cov", "_", ysi)
      cov_name_sigma <- paste0("cov_sigma", "_", ysi)
    }
    
    
    xoffsetnamelist[[ii]] <- xoffset_name
    xoffsetvaluelist[[ii]] <- xoffset
    
    knotsnamelist[[ii]] <- knots_name
    knotsvaluelist[[ii]] <- knots
    
    fixednamelist[[ii]] <- fixed_name
    fixedvaluelist[[ii]] <-
      strsplit(gsub("\\+", " ", fixedsi), " ")[[1]]
    
    randomnamelist[[ii]] <- random_name
    randomvaluelist[[ii]] <-
      strsplit(gsub("\\+", " ", randomsi), " ")[[1]]
    
    
    groupvarnamelist[[ii]] <- groupvar_name
    groupvarvaluelist[[ii]] <- group_arg_groupvar
    
    
    sigma_groupvarnamelist[[ii]] <- sigma_groupvar_name
    sigma_groupvarvaluelist[[ii]] <- sigma_group_arg_groupvar
    
    hierarchicalvarnamelist[[ii]] <- hierarchical_name
    hierarchicalvarvaluelist[[ii]] <- hierarchical_gr_names
    
    sigma_hierarchicalvarnamelist[[ii]] <- sigma_hierarchical_name
    sigma_hierarchicalvarvaluelist[[ii]] <- sigma_hierarchical_gr_names
    
    xnamelist[[ii]] <- xvar_name
    xvarvaluelist[[ii]] <- xsi
    
    ynamelist[[ii]] <- yvar_name
    yvarvaluelist[[ii]] <- ysi
    
    covnamelist[[ii]] <- cov_name
    covvaluelist[[ii]] <- covariates_
    
    sigmacovnamelist[[ii]] <- cov_name_sigma
    sigmacovvaluelist[[ii]] <- covariates_sigma_
    
    
    
    # Restore x var for data
    # Imp: Note that only xvar is reverse transformed and not the yvar
    # This is because the xvar is again transformed in the Stan function block
    # Whereas the yvar is passed to bsitar as such 
    
    
    if (!is.null(xfunsi[[1]][1]) & xfunsi != "NULL") {
      if (xfunsi == "log") {
        datai[[xsi]] <- exp(datai[[xsi]] + xoffset)
      } else if (xfunsi == "sqrt") {
        datai[[xsi]] <- (datai[[xsi]] + xoffset) ^ 2
      }
    } else if (is.null(xfunsi[[1]][1]) | xfunsi == "NULL") {
      datai[[xsi]] <- (datai[[xsi]] + xoffset)
    }
    
    
    
    
    if (!(is.na(univariate_by$by) | univariate_by$by == "NA"))
      dataout <- rbind(dataout, datai)
    else
      dataout <- datai
    
    # print(head(dataout))
    # stop()
    # "kkkk"
    
    if (!(is.na(univariate_by$by) | univariate_by$by == "NA"))
      uvarbyTF <- TRUE
    else
      uvarbyTF <- FALSE
  }  # end of loop over outcome(s)
  
  
  
  if (verbose) {
    if (multivariate$mvar) {
      setmsgtxt <-
        paste0(
          "\n Combining formula, function, priors and initials",
          "\n ",
          "for multivariate model fitting"
        )
    }
    if (!is.na(univariate_by$by)) {
      setmsgtxt <-
        paste0(
          "\n Combining formula, function, priors and initials",
          "\n ",
          "for univariate-by-subgroup model fitting"
        )
    }
    if (displayit == 'msg') {
      message(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolh
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  
  
  # Now collect and combine elemets for univariate-by-sybgroup and multivariate
  # models
  
  brmsdata <- dataout
  brmspriors <- priorlist
  
  bflist_c_list <- list()
  bflist_c <- c()
  for (il in 1:length(bflist)) {
    bflist_c_list[[il]] <- ept(bflist[[il]])
    bflist_c <- c(bflist_c, paste0("bflist_c_list[[", il, "]]"))
  }
  
  bformula <- ept(paste(bflist_c, collapse = "+"))
  
  
  if (nys > 1) {
    if (!(is.na(univariate_by$by) | univariate_by$by == "NA")) {
      bformula <- bformula + set_rescor(FALSE)
    }
    if (multivariate$mvar && multivariate$rescor) {
      bformula <- bformula + set_rescor(TRUE)
    }
    if (multivariate$mvar && !multivariate$rescor) {
      bformula <- bformula + set_rescor(FALSE)
    }
  }
  
  
  bstanvars <-
    stanvar(scode = paste(funlist, collapse = "\n"), block = "function")
  
  prior_stanvarlistlist <- c()
  for (i in 1:nys) {
    prior_stanvarlistlist[i] <- paste0("prior_stanvarlist[[", i, "]]")
  }
  
  bstanvars <-
    bstanvars + eval(parse(text = paste(prior_stanvarlistlist, collapse = "+")))
  
  
  if (length(auxillary_stanvarlist) != 0) {
    auxillary_stanvarlistlist <- c()
    for (i in 1:nys) {
      auxillary_stanvarlistlist[i] <-
        paste0("auxillary_stanvarlist[[", i, "]]")
    }
    bstanvars <-
      bstanvars + eval(parse(text = paste(
        auxillary_stanvarlistlist, collapse = "+"
      )))
  }
  
  
  if (is.list(initialslist) & length(initialslist) == 0) {
    brmsinits <- NULL
  } else if (is.list(initialslist) & length(initialslist) > 0) {
    clistlist <- c()
    for (i in 1:length(initialslist)) {
      clistlist <- c(clistlist, ept(paste0("initialslist[[", i, "]]")))
    }
    brmsinits <- clistlist
  }
  
  
  if (!is.null(brmsinits)) {
    if (multivariate$mvar & multivariate$cor == "un") {
      # combine sd
      c_it <- "sd_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      temppp <- unlist(unname(temppp))
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      brmsinits[[keys[1]]] <- temppp %>%
        unname()
      # combine cor
      c_it <- "L_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      t_names <- l_comb <- d_comb <- c()
      for (lnamei in names(temppp)) {
        t <- temppp[[lnamei]]
        l <- t[lower.tri(t)]
        names(l) <-
          apply(combn(colnames(t), 2), 2, paste, collapse = "_")
        d_comb <- c(d_comb, ncol(t))
        l_comb <- c(l_comb, l)
        t_names <- c(t_names, colnames(t))
      }
      #
      create_cor_mat <- function(n, cor = NULL) {
        n_elements <- n
        m <- diag(n_elements)
        m_upper <- m_lower <- matrix(0, n_elements, n_elements)
        nc <- n_elements * (n_elements - 1) / 2
        if (is.null(cor)) {
          x <- rep(0, nc)
        } else {
          x <- cor
          if (length(x) != nc) {
            stop("length of correlation vector must be ",
                 nc,
                 "\n, ",
                 ", but found ",
                 length(x))
          }
        }
        m_lower[lower.tri(m_lower, diag = FALSE)] <- x
        m_upper <- t(m_lower)
        M <- m_lower + m + m_upper
        M
      }
      #
      tt_names <- apply(combn(t_names, 2), 2, paste, collapse = "_")
      tt_dims <- sum(d_comb)
      tt_nc <- (tt_dims * (tt_dims - 1) / 2)
      tt_12 <- create_cor_mat(tt_dims, rep(0, tt_nc))
      colnames(tt_12) <- rownames(tt_12) <- t_names
      tt_ll <- tt_12[lower.tri(tt_12)]
      names(tt_ll) <-
        apply(combn(colnames(tt_12), 2), 2, paste, collapse = "_")
      tt_ll[names(l_comb)] <- l_comb
      tt_ll[!names(tt_ll) %in% names(l_comb)] <- 0
      brmsinits[[keys[1]]] <- create_cor_mat(tt_dims, tt_ll)
      # combine std z
      c_it <- "z_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      brmsinits[[keys[1]]] <- do.call(rbind, temppp)
    } else if (multivariate$mvar &
               (multivariate$cor == "un" | multivariate$cor ==
                "un_s") &
               !any(grepl("^L_", names(brmsinits))))
    {
      # combine sd
      c_it <- "sd_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      temppp <- unlist(unname(temppp))
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      for (sdi in 1:length(temppp)) {
        brmsinits[[paste0(c_it, sdi)]] <- temppp[sdi] %>%
          unname()
      }
    }  # else if(!any(grepl('^L_', names(brmsinits)))) {
    
    
    # keep only one Lrescor
    if (multivariate$mvar & multivariate$rescor) {
      c_it <- "Lrescor_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      brmsinits[["Lrescor"]] <- temppp[[1]]
    }
    
    if ((multivariate$mvar &
         multivariate$cor == "diagonal") |
        (!is.na(univariate_by$by) &
         univariate_by$cor == "diagonal") |
        group_arg$cor == "diagonal" |
        sigma_group_arg$cor == "diagonal") {
      # combine sd
      c_it <- "sd_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      temppp <- unlist(unname(temppp))
      brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      xxx <- temppp
      ilc <- list()
      ilc_c <- 0
      for (nysi_ in 1:nys) {
        for (il in letters[1:4]) {
          ilc_c <- ilc_c + 1
          na <- paste0("^", il, nysi_)
          nb <- paste0("^", il, "cov", nysi_)
          nanb <- paste0(na, "|", nb)
          if (length(xxx[grepl(nanb, names(xxx))]) > 0) {
            ilc[[ilc_c]] <- xxx[grepl(nanb, names(xxx))]
          }
        }
      }
      ilc <- ilc[lengths(ilc) != 0]
      names(ilc) <- paste0("sd_", 1:length(ilc))
      
      for (sdi in names(ilc)) {
        brmsinits[[sdi]] <- ilc[[sdi]]
      }
      
      
      c_it <- "z_"
      brmsinits_names <- names(brmsinits)
      keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
      temppp <- brmsinits[names(brmsinits) %in% keys]
      
      for (zi in 1:length(ilc)) {
        brmsinits[[paste0(c_it, zi)]] <- temppp[[zi]]
      }
    }
  }  # if(!is.null(brmsinits)) {
  
  
  # For multivariate, it makes sense to keep initials for betas only otherwise
  # dimensional mismatch
  if (!is.null(brmsinits) & length(initialslist) != nys) {
    if (multivariate$mvar & multivariate$cor == "un") {
      c_it_names <- c("sd_", "L_", "z_", "Lrescor")
      for (c_it in c_it_names) {
        brmsinits_names <- names(brmsinits)
        keys <- brmsinits_names[grepl(c_it, brmsinits_names)]
        brmsinits <- brmsinits[!names(brmsinits) %in% keys]
      }
    }
  }
  
  
  
  if (all(sapply("random", grepl, initialslist_s))) {
    brmsinits <- "random"
    brmsinits_r <- ept(init_rsi)
    brmsinits_ <- NULL
  } else if (all(sapply("0", grepl, initialslist_s))) {
    brmsinits <- "0"
    brmsinits_r <- ept(init_rsi)
    brmsinits_ <- NULL
  } else {
    brmsinits <- brmsinits
    brmsinits_r <- NULL
    brmsinits_ <- ""
  }
  
  
  # Strip off initials attributes now
  for (inm in names(brmsinits)) {
    if (is.matrix(brmsinits[[inm]])) {
      colnames(brmsinits[[inm]]) <- rownames(brmsinits[[inm]]) <- NULL
      t__ <- brmsinits[[inm]]
      if (!is.null(attr(t__, "dimnames")))
        attr(t__, "dimnames") <- NULL
      brmsinits[[inm]] <- t__
    }
    if (is.vector(brmsinits[[inm]])) {
      t__ <- brmsinits[[inm]]
      if (!is.null(attr(t__, "names")))
        attr(t__, "names") <- NULL
      brmsinits[[inm]] <- t__
    }
  }
  
  
  
  if (!is.null(brmsinits_)) {
    eval_inits_fun <-
      function(inits,
               jitter_init_beta,
               jitter_init_sd,
               jitter_init_cor,
               digits) {
        if (is.character(jitter_init_beta))
          jitter_init_beta <- ept(jitter_init_beta)
        if (is.character(jitter_init_sd))
          jitter_init_sd <- ept(jitter_init_sd)
        if (is.character(jitter_init_cor))
          jitter_init_cor <- ept(jitter_init_cor)
        
        if (!is.null(jitter_init_beta) &
            !is.numeric(jitter_init_beta)) {
          stop("Argument jitter_init_beta should be NULL or a numeric value")
        }
        if (!is.null(jitter_init_sd) &
            !is.numeric(jitter_init_sd)) {
          stop("Argument jitter_init_sd should be NULL or a numeric value")
        }
        if (!is.null(jitter_init_cor) &
            !is.numeric(jitter_init_cor)) {
          stop("Argument jitter_init_cor should be NULL or a numeric value")
        }
        
        jitter_x <- function(x, a, digits) {
          x <- unname(x)
          col <- c()
          for (i in 1:length(x)) {
            amount <- abs(x[i]) * a
            col <- c(col, jitter(x[i], factor = 1, amount = amount))
          }
          col <- round(col, digits)
          col
        }
        
        jitter_mat <- function(x, a, digits) {
          mat_out <- x
          x <- x[lower.tri(x)]
          col <- c()
          for (i in 1:length(x)) {
            amount <- abs(x[i]) * a
            col <- c(col, jitter(x[i], factor = 1, amount = amount))
          }
          col <- round(col, digits)
          col <- ifelse(col > 1, 1, col)
          col <- ifelse(col < -1, 1, col)
          mat_out[lower.tri(mat_out)] <-
            mat_out[upper.tri(mat_out)] <- col
          return(mat_out)
        }
        
        eval_inits <- c()
        for (i_init in names(inits)) {
          if (grepl("^b_", i_init)) {
            if (!is.null(jitter_init_beta)) {
              values_i <-
                jitter_x(inits[[i_init]], jitter_init_beta, digits = digits)
            } else if (is.null(jitter_init_beta)) {
              values_i <- inits[[i_init]]
              values_i <- round(values_i, digits)
            }
            eval_inits[[i_init]] <- values_i
          } else if (grepl("^sd_", i_init)) {
            if (!is.null(jitter_init_sd)) {
              values_i <- jitter_x(inits[[i_init]], jitter_init_sd, digits)
              values_i <- abs(values_i)
              values_i <-
                ifelse(values_i <= 0, values_i + 0.01, values_i)
            } else if (is.null(jitter_init_sd)) {
              values_i <- inits[[i_init]]
              values_i <- abs(round(values_i, digits))
              values_i <-
                ifelse(values_i <= 0, values_i + 0.01, values_i)
            }
            eval_inits[[i_init]] <- values_i
          } else if (grepl("^L_", i_init)) {
            if (!is.null(jitter_init_cor)) {
              values_i <- jitter_mat(inits[[i_init]], jitter_init_cor, digits)
            } else if (is.null(jitter_init_cor)) {
              values_i <- inits[[i_init]]
              values_i <- round(values_i, digits)
            }
            eval_inits[[i_init]] <- values_i
          } else {
            eval_inits[[i_init]] <- inits[[i_init]]
          }
        }  # for(i_init in names(inits)) {
        eval_inits
        return(eval_inits)
      }
    
    brmsinits <- lapply(1:brms_arguments$chains, function(id) {
      eval_inits_fun(
        inits = brmsinits,
        jitter_init_beta = jitter_init_beta,
        jitter_init_sd = jitter_init_sd,
        jitter_init_cor = jitter_init_cor,
        digits = 4
      )
    })
  }
  
  
  
  # Set brm arguments and fit model
  
  setup_brms_args <-
    function(formula,
             prior,
             stanvars,
             data,
             init_set,
             init_str,
             init_r,
             seed,
             verbose,
             setarguments,
             brmsdots) {
      
      exc_args <- c("formula", "prior", "stanvars", "init", "data")
      if (eval(setarguments$backend) == "rstan")
        exc_args <- c(exc_args, "stan_model_args")
      for (exc_argsi in exc_args) {
        if (exc_argsi %in% names(setarguments))
          setarguments[[exc_argsi]] <- NULL
      }
      setarguments$formula <- formula
      setarguments$prior <- prior
      setarguments$stanvars <- stanvars
      setarguments$data <- data
      
      if (eval(setarguments$backend) == "cmdstanr") {
        if (all(sapply("0", grepl, init_str))) {
          setarguments$init <- 0
          custom_init <- FALSE
        } else if (all(sapply("random", grepl, init_str))) {
          setarguments$init <- NULL
          custom_init <- FALSE
        } else {
          setarguments$init <- init_set
          custom_init <- TRUE
        }
        if (!custom_init & !is.null(init_r)) {
          setarguments$init <- init_r
        }
      }
      
      if (eval(setarguments$backend) == "rstan") {
        if (all(sapply("0", grepl, init_str))) {
          setarguments$init <- "0"
          custom_init <- FALSE
        } else if (all(sapply("random", grepl, init_str))) {
          setarguments$init <- "random"
          custom_init <- FALSE
        } else {
          setarguments$init <- init_set
          custom_init <- TRUE
        }
        if (!custom_init & !is.null(init_r)) {
          setarguments$init_r <- init_r
        }
      }
      
      if (eval(setarguments$backend) == "rstan" | 
          eval(setarguments$backend) == "cmdstanr") {
        if (is.null(eval(setarguments$control))) {
          setarguments$control <- list(adapt_delta = 0.8, max_treedepth = 15)
        }
        if (is.na(eval(setarguments$seed))) {
          setarguments$seed <- seed
        }
        
        cores_ <- eval(setarguments$cores)
        threads_ <- eval(setarguments$threads)
        
        if(cores_ == "maximise") {
          max.cores <- 
            as.numeric(future::availableCores(methods = "system", omit = 0))
          if(max.cores < 1) max.cores <- 1
        } else if(cores_ == "optimize") {
          max.cores <- 
            as.numeric(future::availableCores(methods = "system", omit = 1))
          if(max.cores < 1) max.cores <- 1
          if(max.cores > eval(setarguments$chains)) {
            max.cores <- eval(setarguments$chains)
          }
        } else if(!is.null(getOption('mc.cores')) &
                  cores_ != "maximise" &
                  cores_ != "optimize") {
          max.cores <- getOption('mc.cores')
        } else {
          max.cores <- eval(setarguments$cores)
        }
        setarguments$cores <-  max.cores
        
        
        
        if(!is.list(threads_)) {
          if( is.character(threads_) & threads_ == "maximise") {
            max.threads <- 
              as.numeric(future::availableCores(methods = "system", omit = 0))
            if(max.threads < 1) max.threads <- 1
          } else if( is.character(threads_) & threads_ == "optimize") {
            max.threads <- 
              as.numeric(future::availableCores(methods = "system", omit = 1))
            if(max.threads < 1) max.threads <- 1
            max.threads <- floor(max.threads /  eval(setarguments$chains))
          } else if(!is.null(getOption('brms.threads')) &
                    (is.character(threads_) & threads_ != "maximise") &
                    (is.character(threads_) & threads_ != "optimize")) {
            max.threads <- getOption('brms.threads')
          } else if(is.null(getOption('brms.threads')) &
                    (is.character(threads_) & threads_ != "maximise") &
                    (is.character(threads_) & threads_ != "optimize")) {
            max.threads <- getOption('brms.threads')
          } else {
            max.threads <- eval(setarguments$cores)
          }
          setarguments$threads <-  threading(max.threads)
        }
      } 
      
      
      if (eval(setarguments$backend) == "cmdstanr") {
        if (is.list(eval(setarguments$stan_model_args)) &
            eval(length(setarguments$stan_model_args)) == 0) {
          setarguments$stan_model_args <- list(stanc_options = list("O1"))
        }
      }
      
      if (eval(setarguments$backend) == "rstan" & 
          packageVersion("rstan") < "2.26.1") {
        # placeholder, will be updated later when rstan > 2.26.1 on CRAN
        # brms takes care of threads options for rstan by the version
        setarguments$threads <- setarguments$threads 
      }
      
      if (length(brmsdots) > 0) {
        setarguments <- c(setarguments, brmsdots)
      }
      return(setarguments)
    }
  
  
  
  if (verbose) {
    setmsgtxt <- paste0("\n Setting-up brms arguments")
    if (displayit == 'msg') {
      message(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolh
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  
  
  
  
  brmsdots_ <- list(...)
  
  brm_args <-
    setup_brms_args(
      formula = bformula,
      prior = brmspriors,
      stanvars = bstanvars,
      data = brmsdata,
      init = brmsinits,
      init_str = initialslist_s,
      init_r = brmsinits_r,
      seed,
      verbose,
      setarguments = brms_arguments,
      brmsdots = brmsdots_
    )
  
  
  if (verbose) {
    setmsgtxt <- paste0("\n Fitting model")
    if (displayit == 'msg') {
      message(setmsgtxt)
    } else if (displayit == 'col') {
      col <- setcolh
      cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
    }
  }
  
  cat("\n")
  
  
  ######################
  # not using insert_new_priors
  
  insert_new_priors <- function(setdf_1, setdf_2) {
    index <- row_number <- valid <- NA
    setdf_1 <- 
      setdf_1 %>% dplyr::mutate(index = interaction(class, coef, group, nlpar)) %>% 
      dplyr::mutate(order = row_number()) %>% 
      dplyr::arrange(index)
    setdf_2 <- 
      setdf_2 %>% dplyr::mutate(index = interaction(class, coef, group, nlpar)) %>% 
      dplyr::mutate(order = row_number()) %>% 
      dplyr::arrange(index)
    vi_1 <- setdf_1 %>% dplyr::mutate(valid = ifelse(!(class == 'sd' & coef == ""), 1, 0)) %>% 
      data.frame() %>% dplyr::filter(valid == 1) %>% dplyr::select(index) %>% unlist() %>% 
      droplevels()
    setdf_1 <- setdf_1 %>% dplyr::filter(!(class == 'sd' & coef == ""))
    setdf_2 <- setdf_2 %>% dplyr::filter(!(class == 'sd' & coef == ""))
    setdf_3 <- setdf_1[setdf_1$index %in% vi_1,]
    setdf_2 <- setdf_2 %>% dplyr::filter(!index %in% vi_1)
    setdf_4 <- rbind(setdf_2, setdf_3)
    setdf_4 <- setdf_4 %>% dplyr::select(-c(index, order))
    setdf_4
  }
  
  
  insert_new_priors <- function(setdf_1, setdf_2) {
    cc <- zz <- list()
    for (i in 1:nrow(setdf_1)) {
      getx <- setdf_1[i,]
      zz[[i]] <- setdf_2 %>% dplyr::mutate(prior = dplyr::if_else(class == getx[['class']] &
                                                                    coef == getx[['coef']] &
                                                                    group == getx[['group']] &
                                                                    resp == getx[['resp']] &
                                                                    dpar == getx[['dpar']] &
                                                                    nlpar == getx[['nlpar']] ,
                                                                  getx$prior,
                                                                  setdf_2$prior)) %>% 
        dplyr::filter(class == getx[['class']] &
                        coef == getx[['coef']] &
                        group == getx[['group']] &
                        resp == getx[['resp']] &
                        dpar == getx[['dpar']] &
                        nlpar == getx[['nlpar']])
      
      cc[[i]] <- setdf_2 %>% dplyr::mutate(prior = dplyr::if_else(class != getx[['class']] &
                                                                    coef != getx[['coef']] &
                                                                    group != getx[['group']] &
                                                                    resp != getx[['resp']] &
                                                                    dpar != getx[['dpar']] &
                                                                    nlpar != getx[['nlpar']] ,
                                                                  setdf_2$prior,
                                                                  getx$prior)) %>% 
        
        dplyr::filter(class != getx[['class']] & 
                        coef != getx[['coef']] &
                        group != getx[['group']] &
                        resp != getx[['resp']] &
                        dpar != getx[['dpar']] &
                        nlpar != getx[['nlpar']])
      
    }
    
    p1 <- cc %>% do.call(rbind, .)
    p2 <- zz %>% do.call(rbind, .)
    p1p2 <- rbind(p1, p2)
    p1p2
  }
  
  ######################
  
  
  if(set_higher_levels) {
    brmspriors_sdcor <- brmspriors %>% 
      dplyr::filter(class == 'sd' | class == 'cor')
    brmspriors_sdcor_gr <- brmspriors_sdcor$group
    
    brmsfit_sdcor <- do.call(get_prior, brm_args) %>% 
      dplyr::filter(class == 'sd' | class == 'cor')
    
    brmsfit_sdcor_prior_gr <- brmsfit_sdcor %>% 
      dplyr::filter(!group %in%  brmspriors_sdcor_gr)
    
    brmspriors_brmsfit_sdcor <- brmspriors %>% 
      dplyr::bind_rows(., brmsfit_sdcor_prior_gr)
    
    brmspriors <- brmspriors_brmsfit_sdcor
  }
  
  
  
  # if(sigma_set_higher_levels) {
  #   brmspriors_sdcor <- brmspriors %>% 
  #     dplyr::filter(class == 'sd' | class == 'cor')
  #   brmspriors_sdcor_gr <- brmspriors_sdcor$group
  #   
  #   brmsfit_sdcor <- do.call(get_prior, brm_args) %>% 
  #     dplyr::filter(class == 'sd' | class == 'cor')
  #   
  #   brmsfit_sdcor_prior_gr <- brmsfit_sdcor %>% 
  #     dplyr::filter(!group %in%  brmspriors_sdcor_gr)
  #   
  #   brmspriors_brmsfit_sdcor <- brmspriors %>% 
  #     dplyr::bind_rows(., brmsfit_sdcor_prior_gr)
  #   
  #   brmspriors <- brmspriors_brmsfit_sdcor
  # }
  
  
  
  brm_args$prior <- brmspriors
  
  
  
  if(!is.null(set_self_priors) & !is.null(set_replace_priors)) {
    stop("Amongst 'set_self_priors' and 'set_replace_priors' arguments,",
         "\n ",
         " only one can be specified at a time")
  }
  
  if(get_priors & get_set_priors & validate_priors & 
     get_stancode & get_standata) {
    stop("Amongst 'get_priors' 'get_set_priors', 'validate_priors' ",
         "\n ",
         "'get_stancode' and 'get_standata' ",
         "\n ",
         " arguments, only one can be set to TRUE at a time")
  }
  
  
  
  exe_model_fit <- TRUE
  if(get_stancode |
     get_standata |
     get_priors |
     get_set_priors |
     validate_priors) {
    exe_model_fit <- FALSE
  }
  
  
  
  # if(exe_model_fit) {
  #   if(!is.null(set_self_priors)) {
  #     brm_args$prior <- set_self_priors
  #   } else if(!is.null(set_replace_priors)) {
  #     # brm_args$prior <- insert_new_priors(set_replace_priors, brmspriors)
  #     brm_args$prior <- brmspriors %>% 
  #       dplyr::filter(source == 'user') %>% 
  #       dplyr::bind_rows(., set_replace_priors)
  #   } else if(is.null(set_self_priors) & is.null(set_replace_priors)) {
  #     brm_args$prior <- brmspriors
  #   }
  
  
  
  
  # IMP - brms does not allow different lb conditions for sd parsm (e.e, all to be NA)
  # Because prior function automatically sets lb 0 for positive priors such as exponentials
  # the following is need
  lbbb_ <- ubbb_ <- NULL
  tempprior_hold <- brmspriors # brm_args$prior 
  setpriornamesorder <- colnames(tempprior_hold)
  tempprior_hold$lbbb_ <- tempprior_hold$lb
  tempprior_hold$ubbb_ <- tempprior_hold$ub
  tempprior_hold$lb <- tempprior_hold$ub <- NULL
  tempprior_hold <- tempprior_hold %>% dplyr::mutate(lbbb_ = dplyr::if_else(class == 'sd', NA, lbbb_))
  tempprior_hold <- tempprior_hold %>% dplyr::mutate(ubbb_ = dplyr::if_else(class == 'sd', NA, ubbb_))
  tempprior_hold$lb <- tempprior_hold$lbbb_
  tempprior_hold$ub <- tempprior_hold$ubbb_
  tempprior_hold$lbbb_ <- tempprior_hold$ubbb_ <- NULL
  tempprior_hold <- tempprior_hold %>% dplyr::relocate(dplyr::all_of(setpriornamesorder))
  # brm_args$prior <- tempprior_hold
  brmspriors <-   tempprior_hold
  
  brm_args$prior <- brmspriors
  
  
  
  if(!exe_model_fit) {
    if(get_priors) {
      options(mc.cores = mc.cores_restore)
      return(do.call(get_prior, brm_args))
    } else if(get_standata) {
      options(mc.cores = mc.cores_restore)
      return(do.call(make_standata, brm_args))
    } else if(get_stancode) {
      options(mc.cores = mc.cores_restore)
      return(do.call(make_stancode, brm_args))
    } else if(get_set_priors) {
      options(mc.cores = mc.cores_restore)
      return(brm_args$prior)
    } else if(validate_priors) {
      options(mc.cores = mc.cores_restore)
      return(do.call(validate_prior, brm_args))
    }
  } 
  
  
  
  if(exe_model_fit) {
    if(!is.null(set_self_priors)) {
      brm_args$prior <- set_self_priors
    } else if(!is.null(set_replace_priors)) {
      # brm_args$prior <- insert_new_priors(set_replace_priors, brmspriors)
      brm_args$prior <- brmspriors %>%
        dplyr::filter(source == 'user') %>%
        dplyr::bind_rows(., set_replace_priors)
    } else if(is.null(set_self_priors) & is.null(set_replace_priors)) {
      brm_args$prior <- brmspriors
    }
    
    
    
    # txx <<- brm_args$prior
    # stop()
    
    # If initials are 0 or random, then set custom init to NULL
    
    if(brm_args$backend == "rstan") {
      if(length(brm_args$init) == 1) {
        if(brm_args$init == "0") {
          init_custom <- NULL
        } else if(brm_args$init == "random") {
          init_custom <- NULL
        } else {
          init_custom <- init_custom
        }
      } else {
        init_custom <- init_custom
      }
    }
    
    
    
    if(brm_args$backend == "cmdstanr") {
      if(is.null(brm_args$init)) {
        init_custom <- NULL
      } else if(length(brm_args$init) == 1) {
        if(brm_args$init == "0") {
          init_custom <- NULL
        } else if(brm_args$init == "random") {
          init_custom <- NULL
        } else if(brm_args$init == 0) {
          init_custom <- NULL
        } else {
          init_custom <- init_custom
        }
      } else {
        init_custom <- init_custom
      }
    }
    
    
    
    if(!is.null(init_custom)) {
      init_fun <- function(chain_id = 1) init_custom
      if(!is.list(init_custom[[1]])) {
        init_custom <- 
          lapply(1:brm_args$chains, function(id) init_fun(chain_id = id))
      } else if(is.list(init_custom[[1]]) & length(init_custom) == 1) {
        init_custom <- rep(init_custom, brm_args$chains)
      } else {
        if(length(init_custom) != length(brm_args$init)) {
          stop("Custom initials specified via 'init_custom' argument must",
               "\n ", 
               " be a single named list (e.g., custom_init = list(x= 2, xx=5)) ",
               "\n ", 
               " or else a list of list matching the number of chains")
        }
      }
      new_init_append <- list()
      init_old <- brm_args$init
      init_append <- init_custom
      for (ilen in 1:length(init_old)) {
        new_init_append[[ilen]] <- c(init_old[[ilen]], init_append[[ilen]])
      }
      brm_args$init <- new_init_append
    } 
    
    
    brmsfit <- do.call(brm, brm_args)
    
    
    # Add model info for post-processing
    
    model_info <- list()
    
    for (i in 1:length(xoffsetnamelist)) {
      model_info[[xoffsetnamelist[[i]]]] <- xoffsetvaluelist[[i]]
    }
    
    for (i in 1:length(knotsnamelist)) {
      model_info[[knotsnamelist[[i]]]] <- knotsvaluelist[[i]]
    }
    
    for (i in 1:length(fixednamelist)) {
      model_info[[fixednamelist[[i]]]] <- fixedvaluelist[[i]]
    }
    
    for (i in 1:length(randomnamelist)) {
      model_info[[randomnamelist[[i]]]] <- randomvaluelist[[i]]
    }
    
    for (i in 1:length(xfunnamelist)) {
      model_info[[xfunnamelist[[i]]]] <- xfunvaluelist[[i]]
    }
    
    for (i in 1:length(yfunnamelist)) {
      model_info[[yfunnamelist[[i]]]] <- yfunvaluelist[[i]]
    }
    
    for (i in 1:length(xxfunnamelist)) {
      model_info[[xxfunnamelist[[i]]]] <- xxfunvaluelist[[i]]
    }
    
    for (i in 1:length(yyfunnamelist)) {
      model_info[[yyfunnamelist[[i]]]] <- yyfunvaluelist[[i]]
    }
    
    for (i in 1:length(groupvarnamelist)) {
      model_info[[groupvarnamelist[[i]]]] <- groupvarvaluelist[[i]]
    }
    
    for (i in 1:length(hierarchicalvarnamelist)) {
      model_info[[hierarchicalvarnamelist[[i]]]] <- hierarchicalvarvaluelist[[i]]
    }
    
    for (i in 1:length(xnamelist)) {
      model_info[[xnamelist[[i]]]] <- xvarvaluelist[[i]]
    }
    
    for (i in 1:length(ynamelist)) {
      model_info[[ynamelist[[i]]]] <- yvarvaluelist[[i]]
    }
    
    for (i in 1:length(covnamelist)) {
      model_info[[covnamelist[[i]]]] <- covvaluelist[[i]]
    }
    
    for (i in 1:length(sigmacovnamelist)) {
      model_info[[sigmacovnamelist[[i]]]] <- sigmacovvaluelist[[i]]
    }
    
    if(!is.na(univariate_by$by)) {
      model_info[['subindicators']] <- subindicators
    } 
    
    model_info[[SplineFun_name]] <- SplineFun_name
    model_info[['multivariate']] <- multivariate$mvar
    model_info[['univariate_by']] <- univariate_by$by
    model_info[['nys']] <- nys
    model_info[['ys']] <- ys
    
    model_info[['xs']] <- xs
    model_info[['ids']] <- ids
    
    model_info[['dfs']] <- dfs
    
    
    model_info[['xfuns']] <- xfuns
    model_info[['yfuns']] <- yfuns
    
    model_info[['outliers']] <- outliers
    
    model_info[['bsitar.data']] <- data.org.in
    
    model_info[['call.full.bsitar']] <- call.full
    
    model_info[['call.bsitar']] <- mcall_
    
    brmsfit$model_info <- model_info
    
    # Expose Stan function
    
    if (expose_function) {
      if (verbose) {
        setmsgtxt <-
          paste0("\n Exposing Stan functions for post-processing\n")
        if (displayit == 'msg') {
          message(setmsgtxt)
        } else if (displayit == 'col') {
          col <- setcolh
          cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
        }
      }
      
      if (!verbose) {
        setmsgtxt <-
          paste0("\n Exposing Stan functions for post-processing..\n")
        message(setmsgtxt)
      }
      
      brmsfit <- expose_bsitar_functions(brmsfit)
    }
    
    if (!expose_function) {
      brmsfit$model <- brmsfit$bmodel <- stancode(brmsfit)
    }
    
    if (verbose) {
      setmsgtxt <- paste0("\nModel Fitting complete")
      if (displayit == 'msg') {
        message(setmsgtxt)
      } else if (displayit == 'col') {
        col <- setcolh
        cat(paste0("\033[0;", col, "m", setmsgtxt, "\033[0m", "\n"))
      }
    }
    
    attr(brmsfit, 'class') <- c(attr(brmsfit, 'class'), 'bsitar')
    options(mc.cores = mc.cores_restore)
    return(brmsfit)
  } # exe_model_fit
}
