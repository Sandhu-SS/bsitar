

#'Fit Bayesian SITAR growth curve model
#'
#'@description Fit Bayesian super imposition by translation and rotation (SITAR)
#'  model that summarizes the growth curves from early childhood through the
#'  adulthood (see @details). The frequentist version of the SITAR model can be
#'  fit by an already available R package, the *sitar*
#'  \insertCite{R-sitar}{bsitar}. Besides Bayesian estimation of the SITAR
#'  model, the *bsitar* package greatly enhances the modelling capabilities
#'  offered by the *sitar* package. For example, in addition to the univariate
#'  model fitting (i.e, modelling a single response as implemented in the
#'  *sitar* package), the **bsitar** allows univariate-by-subgroup and
#'  multivariate model specifications (see @details).
#'
#'@details The SITAR is a shape-invariant nonlinear mixed effect growth curve
#'  model that fits a population average (i.e., mean average) curve to the data
#'  and aligns each individual's growth trajectory to the population average
#'  curve via a set of three random effects (size, timing, and intensity). The
#'  concept of shape invariant model (SIM) was first described by
#'  \insertCite{Lindstrom1995}{bsitar} and later used by
#'  \insertCite{Beath2007;textual}{bsitar} to model infant growth data (birth to
#'  2 years). The current version of the SITAR model is developed by
#'  \insertCite{Cole2010;textual}{bsitar} and has been used extensively for
#'  modelling human growth data \insertCite{@see
#'  @nembidzaneUsingSITARMethod2020; @mansukoskiLifeCourseAssociations2019;
#'  @coleFiftyYearsChild2018; @riddellClassifyingGestationalWeight2017;
#'  @Sandhu2020}{bsitar}. As mentioned earlier (see @description), the
#'  frequentist version of the SITAR model can be fit by an already available
#'  R package, the *sitar* \insertCite{R-sitar}{bsitar}.
#'
#'  The model specification is same in both *sitar* and **bsitar** with the
#'  exception that unlike *sitar* which uses the B spline basis for the natural
#'  cubic spline design matrix (by calling the *ns* function of the *splines*
#'  package \insertCite{R-splines}{bsitar}), the *bsitar* constructs spline
#'  design matrix by using the truncated power basis approach as
#'  described by \insertCite{harrell2001regression}{bsitar}, and implemented in
#'  the *rcspline.eval* function of the *Hmisc* package
#'  \insertCite{R-Hmisc}{bsitar}. Note that the **bsitar** package does not use
#'  the *rcspline.eval* but rather constructs a custom function on the fly that
#'  is included in the functions block of the *Stan* programs' and thus
#'  compiled (via the c++) during the model estimation.
#'
#'  Like *sitar* package, the **bsitar** fits SITAR model with (usually) up to
#'  three random effect parameters \insertCite{Cole2010}{bsitar}: the size
#'  (\code{a}), the timing (\code{b}) and the intensity (\code{c}). In addition,
#'  there is a slope parameter \code{d} that models the variability in the adult
#'  slope of the growth curve (See [sitar::sitar] for details). Please note that
#'  inclusion of \code{d} results in multicollinearity because inclusion of this
#'  this \code{d} parameter involves a linear predictor term which is identical
#'  to the first term of the spline design matrix created by the truncated power
#'  basis approach.
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
#'@param x Predictor variable in the data (typically age in years). For
#'  univariate model, the \code{x} is a single variable whereas for the
#'  univariate-by-subgroup (see \code{univariate_by}) and multivariate (see
#'  \code{multivariate}) model specifications, the \code{x} can be same for each
#'  sub model or else specified separately for each sub model. For example, for
#'  a bivariate model, the \code{x = list(x1, x2)} specifies that \code{x1} is
#'  the predictor for the first sub model and \code{x2} for the second sub
#'  model.
#'
#'@param y Response variable in the data (e.g., repeated height
#'  measurements). For univariate and univariate-by-subgroup (see
#'  \code{univariate_by} below) model specifications, \code{y} is specified as a
#'  single variable. For the univariate-by-subgroup model fitting, the outcome
#'  vectors for each sub model are created internally and named using the the
#'  factor levels. For example when fitting a univariate-by-subgroup model for
#'  sex specified as \code{univariate_by = list(by = sex)} or simply as
#'  \code{univariate_by = sex}, the outcome vectors 'Female' and 'Male' are
#'  created automatically  where 'Female' is the first level of the factor
#'  variable sex, and 'Male' is the second level. For multivariate model (see
#'  \code{multivariate} below), the response vectors are specified as a list
#'  (e.g., \code{y = list(y1, y2}) where \code{y1} and \code{y2} are the
#'  responses. Note that for \code{multivariate} model fitting, data are not
#'  stacked but rather response vectors are variables in the data frame.
#'  
#'@param id A vector specifying a unique group identifier for each individual
#'  (typically \code{id} denoting the individual) For univariate-by-subgroup
#'  (see \code{univariate_by} below) and multivariate (see \code{multivariate}
#'  below) model specifications, the \code{id} can be same (typically) or
#'  specified separately for each sub model like \code{id = list(id1, id2)}
#'  where \code{id1} and \code{id2} are individual identifiers.
#'
#'@param data Data frame containing variables \code{x}, \code{y} and \code{id}.
#'
#'@param df Set degrees of freedom for natural cubic regression spline. For
#'  univariate-by-subgroup model (see \code{univariate_by} below) and
#'  multivariate model (see \code{multivariate} below), the \code{df} can be
#'  same \code{df = 4} for each sub model or else specified separately as
#'  \code{df = list(4, 5)} where df is 4 is for the first model and 5 for the
#'  second model.
#'
#'@param knots A vector of knots used for constructing the spline design matrix.
#'  (default \code{df} quantiles of \code{x} distribution). See \code{df} for
#'  specifying separate knots for univariate-by-subgroup and multivariate
#'  models.
#'
#'@param fixed A character string to specify the fixed effects structure.
#'  Typically specified as \code{fixed = a+b+c}. As mentioned earlier, there is
#'  no need to enclose character in quotes. In other words, \code{fixed =
#'  a+b+c}, \code{fixed = 'a+b+c'}, and \code{fixed = "a+b+c"} are same. For
#'  specifying different fixed effect structures for univariate-by-subgroup and
#'  multivariate models, use list as follows: \code{fixed = list(a+b+c, a+b)}
#'  which implies that the fixed effect structure for the first sub model is
#'  \code{fixed = 'a+b+c'} and \code{fixed = 'a+b'} for the second sub model.
#'
#'@param random A character string to specify the random effects structure. The 
#'  approach is same as described above (see \code{fixed}) for setting the 
#'  fixed effects structure.
#'  
#'@param select_model A symbol or a character string specifying the model to be
#'  fitted. Allowed models are SITAR (\code{sitar}), Preece-Baines model 1
#'  (\code{pb1}), Preece-Baines model 2 (\code{pb2}) and Preece-Baines model 3
#'  (\code{pb3}). The option \code{sitar} fits the default three parameter SITAR
#'  model i.e, \code{a+b+c} (see @details). To fit four parameter SITAR model
#'  formulation that includes an additional parameter \code{d} (i.e.,
#'  \code{a+b+c+d}), use \code{select_model = sitar4}. Note that the option
#'  should be specified as either as \code{select_model = sitar4fr} or
#'  \code{select_model = sitar4r} that controls whether or not to include the
#'  parameter \code{d} in the fixed andom effect structure of the model (see
#'  @details). The \code{select_model = sitar4fr} indicates that the parameter
#'  \code{d} is allowed to be included both in the fixed and the random effects
#'  by \code{fixed} and \code{random} arguments. On the contrary,
#'  \code{select_model = sitar4r} implies that parameter \code{d} will be
#'  included only in the random effects part of the model and removed
#'  from the fixed effects structure.
#'
#'@param xoffset An optional character string or a numeric value to set up the
#'  offset for \code{x}. The \code{offset} allows the origin of \code{x} to be
#'  varied by a finite value. The options are \code{mean(x)}, \code{min(x)},
#'  \code{max(x)}, \code{age at peak velocity i.e., apv}, or a numeric values.
#'  These options are specified as \code{xoffset = mean}, \code{xoffset = min},
#'  \code{xoffset = max}, or \code{xoffset = apv} or a numeric values such as
#'  \code{xoffset = 12.5}. The default is \code{xoffset = mean}.
#'
#'@param bstart An optional character string or a numeric value to set up the
#'  initial value for the fixed effect parameter \code{b}. Options are same as
#'  described above for the \code{xoffset}. The default is to use the same value
#'  that has been set up for the \code{xoffset} i.e.,  \code{bstart = xoffset}.
#'  
#'@param pgv An optional numeric value to set up the initial value for the 
#'  fixed effect parameter \code{c}.
#'  
#'@param apgv An optional numeric value to set up the initial value for the 
#'  fixed effect parameter \code{b}. This is an alternative the \code{xoffset}.
#'  
#'@param xfun An optional character string to set up the transformation of the
#'  predictor variable, \code{x}. Options are \code{log} and \code{sqrt} for the
#'  logarithmic and square root transformation, respectively. The default is
#'  \code{NULL} indicating that no transformation is applied i.e., the model is
#'  fit to the original scale of the predictor variable, \code{x}. Like other
#'  arguments, user can specify different xfun for univariate-by-subgroup (see
#'  \code{univariate_by}) and multivariate (see \code{multivariate}) models as a
#'  list i.e., \code{xfun = list(log, sqrt)} or \code{xfun = list(NULL, sqrt)}.
#'
#'@param yfun An optional character string to set up the transformation of the
#'  response variable, \code{y} (default \code{NULL}). Options are \code{log}
#'  and \code{sqrt}. See \code{xfun} for details.
#'
#'@param bound An optional numeric value to extend the span of \code{x} by a
#'  small range (default 0.04). See package 'sitar' for details.
#'
#'@param terms_rhs An optional argument (default \code{NULL}) to specify terms
#'  on the right side of the response variable (separated by '|') but before the
#'  formula. The \code{terms_rhs} is used to fit measurement error model. As an
#'  example, consider fitting a model with measurement error in the response
#'  specified as \code{bf(y | mi(sdy) ~ ..)} where \code{mi(sdy)} is passed to
#'  the [brms::brmsformula()] as \code{terms_rhs = mi(sdy)}. For multivariate
#'  model, each outcome can have its own measurement error variable passed as a
#'  list, i.e., \code{terms_rhs = list(mi(sdy1), mi(sdy2))}.
#'
#'@param a_formula Formula for the fixed effect parameter, \code{a} (default
#'  \code{~ 1}). User can specify different formula when fitting
#'  univariate-by-subgroup (see \code{univariate_by}) and the multivariate (see
#'  \code{multivariate}) models. As an example \code{a_formula = list(~1, ~1 +
#'  cov)} implies that the \code{a_formula} for the first outcome includes only
#'  an intercept whereas the \code{a_formula} for the second outcome includes an
#'  intercept plus and covariate \code{cov}. The covariate can be a continuous
#'  variable or a factor variable (dummy variables will be created using the
#'  \code{model.matrix}). The formula can include a combination of continuous
#'  and factor variables and their interactions.
#'
#'@param b_formula Formula for the fixed effect parameter, \code{b} (default
#'  \code{~ 1}). See \code{a_formula} for details.
#'
#'@param c_formula Formula for the fixed effect parameter, \code{c} (default
#'  \code{~ 1}). See \code{a_formula} for details.
#'
#'@param d_formula Formula for the fixed effect parameter, \code{d} (default
#'  \code{~ 1}). See \code{a_formula} for details.
#'  
#'@param e_formula Formula for the fixed effect parameter, \code{e} (default
#'  \code{~ 1}). See \code{a_formula} for details.
#'  
#'@param f_formula Formula for the fixed effect parameter, \code{f} (default
#'  \code{~ 1}). See \code{a_formula} for details.
#'  
#'@param g_formula Formula for the fixed effect parameter, \code{g} (default
#'  \code{~ 1}). See \code{a_formula} for details.
#'  
#'@param h_formula Formula for the fixed effect parameter, \code{h} (default
#'  \code{~ 1}). See \code{a_formula} for details.
#'  
#'@param i_formula Formula for the fixed effect parameter, \code{i} (default
#'  \code{~ 1}). See \code{a_formula} for details.
#' 
#'@param s_formula Formula for the fixed effect parameter, \code{s} (default
#'  \code{~ 1}) i.e., spline coefficients. See \code{a_formula} for details.
#'
#'@param a_formula_gr Formula for the random effect parameter, \code{a} (default
#'  \code{~ 1}). User can set up the group identifier and the correlation
#'  structure for random effects by using the \code{group_by} argument or the
#'  vertical bar approach. As an example of random effects \code{a}, \code{b},
#'  and \code{c} specified as \code{a_formula_gr = ~1}, \code{b_formula_gr = ~1}
#'  and \code{c_formula_gr = ~1} the group identifier \code{id} and unstructured
#'  correlation structure  is set up via the \code{group_by} argument as
#'  follows: \code{group_by = list(groupvar = id, cor = un)}. The  vertical bar
#'  approach can be used equivalently as \code{a_formula_gr = ~ (1 |i|id)},
#'  \code{b_formula_gr = ~ (1|i|id)}, and \code{c_formula_gr = ~ (1 |i|id)}
#'  where i is just a placeholder.
#'
#'@param b_formula_gr Formula for the random effect parameter, \code{b} (default
#'  \code{~ 1}). See \code{a_formula_gr} for details.
#'
#'@param c_formula_gr Formula for the random effect parameter, \code{c} (default
#'  \code{~ 1}). See \code{a_formula_gr} for details.
#'
#'@param d_formula_gr Formula for the random effect parameter, \code{d} (default
#'  \code{~ 1}). See \code{a_formula_gr} for details.
#'  
#'@param e_formula_gr Formula for the random effect parameter, \code{e} (default
#'  \code{~ 1}). See \code{a_formula_gr} for details.
#'  
#'@param f_formula_gr Formula for the random effect parameter, \code{f} (default
#'  \code{~ 1}). See \code{a_formula_gr} for details.
#'  
#'@param g_formula_gr Formula for the random effect parameter, \code{g} (default
#'  \code{~ 1}). See \code{a_formula_gr} for details.
#'  
#'@param h_formula_gr Formula for the random effect parameter, \code{h} (default
#'  \code{~ 1}). See \code{a_formula_gr} for details.
#'  
#'@param i_formula_gr Formula for the random effect parameter, \code{i} (default
#'  \code{~ 1}). See \code{a_formula_gr} for details.
#'
#'@param a_formula_gr_str Formula for the random effect parameter, \code{a}
#'  (default \code{NULL}) when fitting model with hierarchical structure greater
#'  than two levels (e.g., three level data with repeated measurements (level 1)
#'  on individuals (level 2) nested further within the studies (level 3) in
#'  which individuals participated). See \code{a_formula} for details on
#'  specifying the group identifier and the correlation structure for the random
#'  effects via the \code{group_by} argument and the vertical bar approach. When
#'  using the \code{a_formula_gr_str} argument, only the vertical bar can be
#'  used to set up the group level terms. An example of specifying formula for
#'  random effect parameter \code{a} for a three level model
#'  is as follows: \code{a_formula_gr_str = ~
#'  (1|i|id:study) + (1|i2|study)}. This formulation implies (with \code{|i|}
#'  and \code{|i2|}) that fully unstructured varinace covarinace structure is
#'  specified individuals and study levels. Note that \code{|i|} and \code{|i2|}
#'  need to be distinct as because parameters at different hierarchy levels are
#'  not allowed to be correlated
#'  
#'@param b_formula_gr_str Formula for the random effect parameter, \code{b}
#'  (default \code{NULL}) when fitting model with hierarchical structure greater
#'  than two levels. See \code{a_formula_gr_str} for details.  
#'
#'@param c_formula_gr_str Formula for the random effect parameter, \code{c}
#'  (default \code{NULL}) when fitting model with hierarchical structure greater
#'  than two levels. See \code{a_formula_gr_str} for details.
#'
#'@param d_formula_gr_str Formula for the random effect parameter, \code{d}
#'  (default \code{NULL}) when fitting model with hierarchical structure greater
#'  than two levels. See \code{a_formula_gr_str} for details.
#'  
#'@param e_formula_gr_str Formula for the random effect parameter, \code{e}
#'  (default \code{NULL}) when fitting model with hierarchical structure greater
#'  than two levels. See \code{a_formula_gr_str} for details.
#'  
#'@param f_formula_gr_str Formula for the random effect parameter, \code{f}
#'  (default \code{NULL}) when fitting model with hierarchical structure greater
#'  than two levels. See \code{a_formula_gr_str} for details.
#'  
#'@param g_formula_gr_str Formula for the random effect parameter, \code{g}
#'  (default \code{NULL}) when fitting model with hierarchical structure greater
#'  than two levels. See \code{a_formula_gr_str} for details.
#'  
#'@param h_formula_gr_str Formula for the random effect parameter, \code{h}
#'  (default \code{NULL}) when fitting model with hierarchical structure greater
#'  than two levels. See \code{a_formula_gr_str} for details.
#'  
#'@param i_formula_gr_str Formula for the random effect parameter, \code{i}
#'  (default \code{NULL}) when fitting model with hierarchical structure greater
#'  than two levels. See \code{a_formula_gr_str} for details.
#'
#'@param sigma_formula Formula for the fixed effect parameter, \code{sigma}
#'  (default \code{NULL}) i.e., distributional parameter. The
#'  \code{sigma_formula} is only useful when including covariates(s) for sigma.
#'  The [brms::brm()] by default models the \code{sigma} (i.e., the residual
#'  standard deviation parameter) parameter at the link scale. The
#'  \code{sigma_formula} along with \code{sigma_formula_gr} and
#'  \code{sigma_formula_gr_str} arguments allows specifying hierarchical
#'  structure when modelling sigma. This set up is similar to setting fixed and
#'  random effect parameters such as \code{a}, \code{b}, and \code{c}. The
#'  \code{sigma_formula} sets up the fixed effect design matrix. It is important
#'  to note that another alternative to set up the fixed effect design matrix
#'  for distributional parameter \code{sigma} is argument \code{dpar_formula}.
#'  An advantage of \code{dpar_formula} over \code{sigma_formula} is that user
#'  can specify the linear and nonlinear formulation by using the brms'
#'  [brms::lf] and [brms::nlf] syntax. Both [brms::lf] and [brms::nlf] further
#'  allows control over centering of predictors as well as enabling / disabling
#'  cell mean centering when excluding \code{intercept} via \code{0 + }
#'  formulation. Note that \code{sigma_formula} and \code{dpar_formula} can not
#'  be specified together.
#'
#'@param sigma_formula_gr Formula for the random effect parameter, \code{sigma}
#'  (default \code{NULL}).  
#'
#'@param sigma_formula_gr_str Formula for the random effect parameter, 
#'  \code{sigma} (default \code{NULL}) when fitting model with hierarchical 
#'  structure greater than two levels. See \code{a_formula_gr_str} for details.
#'
#'@param dpar_formula Formula for the distributional parameter \code{sigma}
#'  (default \code{NULL}). This (\code{dpar_formula}) is only useful when
#'  including covariates(s) for sigma. The [brms::brm()] by default models the
#'  \code{sigma} (i.e., the residual standard deviation parameter) parameter at
#'  the link scale. Also note that \code{dpar_formula} can not be specified
#'  along with \code{sigma_formula}, \code{sigma_formula_gr}, or
#'  \code{sigma_formula_gr_str}. See \code{sigma_formula} for relative
#'  advantages and disadvantages of using \code{sigma_formula} and
#'  \code{dpar_formula}.
#'
#'@param autocor_formula Formula for modelling the autocorrelation or residuals.
#'  (default \code{NULL}). Allowed options are: autoregressive moving average
#'  (ARMA) of (ARMA) order \code{p}, \code{q}, autoregressive (AR) of order
#'  \code{p} and moving average (MA) of order \code{q}, and unstructured
#'  (\code{unstr}) over time across individuals that are specified as
#'  \code{autocor_formula = arms(p=1, q=1)}, \code{autocor_formula = ar(p=1)},
#'  \code{autocor_formula = msq=1)} and \code{autocor_formula = unstr(time,
#'  id))} See brms package for further details.
#'
#'@param family Family distribution (default \code{gaussian}) and link function
#'  \ (default \code{identity}) for response variable. See [brms::brm()] for
#'  available distribution and link function and how to specify them. For
#'  univariate-by-subgroup (see \code{univariate_by} ) and multivariate (see
#'  \code{multivariate}) models, the \code{family} could be same \code{family =
#'  gaussian()} or different for each response such as \code{family =
#'  list(gaussian(), student()} which sets gaussian distribution for the first
#'  response variable and the student_t distribution for the second response
#'  variable.
#'
#'@param group_arg Specify named list as group-level sub-arguments for random
#'  effects. These include \code{groupvar}, \code{dist}, \code{cor}, and
#'  \code{by}. The \code{groupvar} specifies the subject identifier. In case
#'  \code{groupvar = NULL} (default), the \code{groupvar} is taken from the
#'  \code{id}). The default \code{dist} is \code{gaussian} and the \code{by} is
#'  \code{NULL}. The default \code{cor} is \code{un} for all three model
#'  settings, i.e., \code{univariate}, \code{univariate_by_subgroup} and
#'  \code{multivariate}. The alternative correlation structure when fitting
#'  \code{univariate} and \code{univariate_by_subgroup} models is
#'  \code{diagonal}. The \code{cor = un} models the full unstructured varinace
#'  covarinace structure whereas \code{cor = diagonal} specifies the diagonal
#'  correlation structure that  estimates only the variance (i.e, standard
#'  deviation) parameters and the covariance (i.e., correlation) parameters are
#'  set to zero. For further details, see brms::brm()] (\bold{Group-level
#'  terms}). Note that only the groupvar suboption of the \code{group_arg} is
#'  passed to the univariate-by-subgroup \code{univariate_by} and the
#'  multivariate (specified by using the \code{multivariate} model fittings
#'  Lastly, the \code{group_arg} is completely ignored when user specify random
#'  formula using the vertical bar ("|") approach. Also, the \code{group_arg} is
#'  ignored for \code{a_formula_gr_str}, \code{b_formula_gr_str}, and
#'  \code{c_formula_gr_str}.
#'  
#'@param sigma_group_arg Specify named list as group-level sub-arguments for 
#'  random effects for the distributional parameter \code{sigma}. The arguments 
#'  are same as \code{group_arg}. See above \code{group_arg} for details.
#'
#'@param univariate_by Specify sub-arguments for univariate-by-subgroup model
#'  (default  \code{NULL}). These include the \code{by} and the
#'  \code{cor} arguments. The \code{by} specifies the variable (which must be a
#'  factor variable) used to set up the sub-models. The \code{cor} specifies the
#'  correlation structure. The \code{un} (default) option i.e., \code{cor = un}
#'  sets up the unstructured correlation structure for each submodel. The
#'  alternative diagonal correlation structure (specified as \code{cor =
#'  diagonal}) estimates only the variance (i.e, standard deviation) for each
#'  submodel whereas covariance (i.e., correlation) parameters are set as zero.
#'
#'@param multivariate Specify sub-arguments for multivariate model fitting
#'  (default \code{NULL}). These include the \code{mvar} argument (logical,
#'  default \code{FALSE}) to specify whether to set up a multivariate model,
#'  \code{cor} to set up the random effect correlation structure, and the
#'  \code{rescor} (logical, default \code{TRUE}) to specify whether or not to
#'  estimate the residual correlation parameter across responses. The \code{cor
#'  = un} (default) sets the correlation as unstructured implying that random
#'  effects across responses are drawn for a joint multivariate normal
#'  distribution with shared variance covarinace parameters. The diagonal
#'  correlation structure \code{cor = diagonal} estimates only the variance
#'  parameters for each response wheras all correlation parameters are set to
#'  \code{NULL} Another option \code{cor = un_s} allows for estimating
#'  unstructured variance covariance parameters separately for each response.
#'
#'@param a_prior_beta Specify priors for the fixed effect parameter, \code{a}.
#'  The allowed distributions  \code{normal},  \code{student_t}, \code{cauchy},
#'  \code{lognormal},  \code{uniform}, \code{exponential}, \code{gamma} and
#'  \code{inv_gamma} (inverse gamma). See [brms::prior()] for details. For each
#'  distribution, sub options \code{lb} and \code{ub} are used to set the upper
#'  and lower bounds (default \code{NA} for both \code{lb} and \code{ub}). For
#'  location scale based distributions (\code{normal}, \code{student_t},
#'  \code{cauchy}, and \code{lognormal}), option \code{autosclae} (default
#'  \code{FALSE}) allows autoscaling of the the scale parameter by a numeric
#'  value. While \code{rstanarm} package earlier used to set it as 2.5 (recently
#'  the authors changed this behavior to \code{FALSE}), the \code{brms} package
#'  sets its to 1 or 2.5 depending on the standard deviation of the response
#'  variable (See [brms::prior()]). The \code{bsitar} package offers the
#'  flexibility of choosing the scaling value as any real number (e.g.,
#'  \code{autosclae = FALSE}, \code{autosclae = 2} or \code{autosclae = 2.5}).
#'  For strictly positive distributions (\code{exponential}, \code{gamma} and
#'  \code{inv_gamma}), the lower bound is automatically set to zero (i.e.,
#'  \code{lb = 0}) For uniform distribution, an option \code{addrange} is
#'  available to symmetrically widen the range between the lower and upper
#'  limits of prior. For example, prior \code{uniform(a, b, addrange = 5)}
#'  implies that the lower and upper limits will be evaluated as
#'  \code{uniform(a-5, b+5)}. For exponential distribution, the rate parameter
#'  is evaluated as inverse. In other words, prior set as \code{exponential(10)}
#'  is translated to \code{exponential(1.0 / 10.0)}. Also, note that the user
#'  need not to specify each option explicitly because the missing information
#'  is added automatically. For example, the \code{normal} prior specified as
#' \code{a_prior_beta = normal(location = 5, scale = 1, lb = NA, ub = NA,
#' addrange = NA, autosclae = FALSE)}) is same as
#'  \code{a_prior_beta = normal(5, 1)}). Lastly, the location parameter for
#'  location scale based distributions can use specified as mean or the median
#'  of the response variable. Similalry, the scale can be set as the standard
#'  deviation (sd) or the median absolute deviation (mad) of the response
#'  variable. As an example, \code{normal} prior can be specified as
#'  \code{a_prior_beta = normal(ymean,ysd)}, \code{a_prior_beta = normal(ymean,
#'   ysd)} or \code{a_prior_beta = normal(ymedian, ymad)}. Another option
#'  available for setting the location is to use the coefficients from the
#'  simple linear model applied to the response (e.g., \code{y ~ age, data
#'  = data}) as follows: \code{a_prior_beta = normal(lm, ysd)}). This is true
#'  even when model has covariates i.e., \code{y ~ age + cov, data = data}. Note
#'  that options \code{ymean}, \code{ymedian}, \code{ysd}, \code{ymad} and 
#'  \code{ymad} are available only for the fixed effect parameter, \code{a} and 
#'  not for any other parameter. For univariate-by-subgroup ( see
#'  \code{univariate_by}) and multivariate ( see \code{multivariate}) models,
#'  priors specified for each response can be same (specified as a single option
#'  i.e., \code{a_prior_beta = normal(5,1)}) or separate for each response
#'  specified as a list such as \code{a_prior_beta = list(normal(5,1),
#'  normal(10, 5)}).
#'  
#'@param b_prior_beta Specify priors for the fixed effect parameter, \code{b}.
#'  See \code{a_prior_beta} for details. 
#'
#'@param c_prior_beta Specify priors for the fixed effect parameter, \code{c}.
#'  See \code{a_prior_beta} for details. 
#'
#'@param d_prior_beta Specify priors for the fixed effect parameter, \code{d}.
#'  See \code{a_prior_beta} for details. 
#'  
#'@param e_prior_beta Specify priors for the fixed effect parameter, \code{e}.
#'  See \code{a_prior_beta} for details. 
#'  
#'@param f_prior_beta Specify priors for the fixed effect parameter, \code{f}.
#'  See \code{a_prior_beta} for details. 
#'  
#'@param g_prior_beta Specify priors for the fixed effect parameter, \code{g}.
#'  See \code{a_prior_beta} for details. 
#'  
#'@param h_prior_beta Specify priors for the fixed effect parameter, \code{h}.
#'  See \code{a_prior_beta} for details. 
#'  
#'@param i_prior_beta Specify priors for the fixed effect parameter, \code{i}.
#'  See \code{a_prior_beta} for details. 
#'
#'@param s_prior_beta  Specify priors for the fixed effect parameter, \code{s}
#'  (i.e., the spline coefficients). The general approach for setting priors for
#'  parameter \code{s} is same as described earlier for the fixed effect
#'  parameters (see \code{a_prior_beta}). The allowed option for location and
#'  scale is \code{lm} which sets location parameter based on the spline
#'  coefficient obatined from the simple linear model fit to the data. The
#'  \code{lm} option for scale parameter sets the standard deviation of the
#'  spline design matrix used to fit the simple linear model. For parameter
#'  \code{s}, it make sense to use only the location scale based prior
#'  distributions (e.g, \code{normal}, \code{student_t}, and \code{cauchy}) or
#'  the \code{uniform} distribution based priors. For \code{uniform} priors, the
#'  \code{addrange} option can be utilized to symmetrically add range to the
#'  \code{lm} based spline coefficients). An additional option available for the
#'  location scale based priors is \code{sethp} (logical, default set as
#'  \code{FALSE}) which, when set as \code{TRUE}, allows for setting
#'  hierarchical priors for the \code{s} parameter. In other words, instead of
#'  setting prior \code{s ~ normal(0, lm)} the hierarchical priors are set as
#'  \code{s ~ normal(0, hp)}\ where \code{hp ~ normal(0, lm)}. Note that the
#'  scale parameter for the \code{hp ~ normal(0, lm)} is automatically taken
#'  from the code{s ~ normal(0, hp)} Setting  \code{sethp = TRUE} implies that
#'  the scale for spline coefficients is estimated from the data itself. The
#'  distribution of hierarchical priors is automatically matched with the prior
#'  set for the \code{s} parameter or else can be set by the same sethp option
#'  set by the same \code{sethp} option. For example, \code{s_prior_beta =
#'  normal(0, lm, sethp = caucy)} will be translated to \code{s ~ normal(0, lm)},
#'  \code{hp ~ caucy(0, lm)}.
#'
#'@param a_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{a}. The approach is same as described for the
#'  \code{a_prior_beta} except that the options \code{ymean}, \code{ymedian},
#'  \code{ysd}, and \code{ymad} are not allowed. The Option \code{lm} for the
#'  location parameter sets \code{Intercept} coefficient obtained from the lm
#'  model fit. Note that options \code{ymean}, \code{ymedian}, \code{ysd},
#'  \code{ymad} and \code{ymad} for covariate coeeficients are are available
#'  only for the fixed effect parameter, \code{a} and not for any other
#'  parameter. Separate priors an be specified for responses when fitting
#'  univariate-by-subgroup (see \code{univariate_by} argument) and the
#'  multivariate (see \code{a_prior_beta}).
#'
#'@param b_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{b}. The approach is same as described for the
#'  the fixed effect parameter, \code{a} (see \code{a_cov_prior_beta})
#'
#'@param c_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{c}. The approach is same as described for the
#'  the fixed effect parameter, \code{a} (see \code{a_cov_prior_beta}).
#'
#'@param d_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{d}. The approach is same as described for the
#'  the fixed effect parameter, \code{a} (see \code{a_cov_prior_beta}).
#'  
#'@param e_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{e}. The approach is same as described for the
#'  the fixed effect parameter, \code{a} (see \code{a_cov_prior_beta}).
#'  
#'@param f_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{f}. The approach is same as described for the
#'  the fixed effect parameter, \code{a} (see \code{a_cov_prior_beta}).
#'  
#'@param g_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{g}. The approach is same as described for the
#'  the fixed effect parameter, \code{a} (see \code{a_cov_prior_beta}).
#'  
#'@param h_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{h}. The approach is same as described for the
#'  the fixed effect parameter, \code{a} (see \code{a_cov_prior_beta}).
#'  
#'@param i_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{i}. The approach is same as described for the
#'  the fixed effect parameter, \code{a} (see \code{a_cov_prior_beta}).
#'
#'@param s_cov_prior_beta Specify priors for the covariate(s) included for the
#'  fixed effect parameter, \code{s}. The approach is same as described for the
#'  \code{s_prior_beta}.
#'
#'@param a_prior_sd Specify priors  for the random effect parameter, \code{a}.
#'  Note that prior is on the standard deviation (which is the square root of
#'  the variance) and not on the variance. The approach is same as earlier for
#'  the fixed effect parameter, \code{a} (See \code{a_prior_beta}) with the
#'  exception that location parameter is always zero. As described earlier (see
#'  \code{a_prior_beta}), priors for the univariate-by-subgroup (see
#'  \code{univariate_by} argument) and multivariate (see \code{multivariate})
#'  models can be same for each response variable or separate for each response.
#'  The lower bound as zero is automatically set by the \code{brms::brm}.
#'
#'@param b_prior_sd  Specify priors  for the random effect parameter, \code{b}.
#'  See \code{a_prior_sd} for details.
#'
#'@param c_prior_sd set Specify priors  for the random effect parameter,
#'  \code{c}. See \code{a_prior_sd} for details.
#'
#'@param d_prior_sd set Specify priors  for the random effect parameter,
#'  \code{d}. See \code{a_prior_sd} for details.
#'
#'@param e_prior_sd sset Specify priors  for the random effect parameter,
#'  \code{e}. See \code{a_prior_sd} for details.
#'
#'@param f_prior_sd sset Specify priors  for the random effect parameter,
#'  \code{f}. See \code{a_prior_sd} for details.
#'
#'@param g_prior_sd sset Specify priors  for the random effect parameter,
#'  \code{g}. See \code{a_prior_sd} for details.
#'  
#'@param h_prior_sd sset Specify priors  for the random effect parameter,
#'  \code{h}. See \code{a_prior_sd} for details.
#'
#'@param i_prior_sd sset Specify priors  for the random effect parameter,
#'  \code{i}. See \code{a_prior_sd} for details.
#'
#'@param a_cov_prior_sd Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{a}. The approach is same as described earlier
#'  for the \code{a_cov_prior_beta} except that no pre-defined option (e.g.,
#'  \code{lm}) is allowed to set the scale parameter for the location scale
#'  based priors.
#'
#'@param b_cov_prior_sd Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{b}. The approach is same as described above
#'  for the \code{a_cov_prior_sd}.
#'
#'@param c_cov_prior_sd Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{c}. The approach is same as described above
#'  for the \code{a_cov_prior_sd}.
#'
#'@param d_cov_prior_sd Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{d}. The approach is same as described above
#'  for the \code{a_cov_prior_sd}.
#'
#'@param e_cov_prior_sd Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{e}. The approach is same as described above
#'  for the \code{a_cov_prior_sd}.
#'
#'@param f_cov_prior_sd Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{f}. The approach is same as described above
#'  for the \code{a_cov_prior_sd}.
#'
#'@param g_cov_prior_sd Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{g}. The approach is same as described above
#'  for the \code{a_cov_prior_sd}.
#'
#'@param h_cov_prior_sd Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{h}. The approach is same as described above
#'  for the \code{a_cov_prior_sd}.
#'
#'@param i_cov_prior_sd Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{i}. The approach is same as described above
#'  for the \code{a_cov_prior_sd}.
#'
#'@param a_prior_sd_str Specify priors for the random effect parameter, \code{a}
#'  when fitting a model to the data with hierarchy level 3 and beyond. The
#'  approach of setting up the priors is same as described earlier (see
#'  \code{a_prior_sd}).
#'
#'@param b_prior_sd_str Specify priors for the random effect parameter, \code{b}
#'  when fitting a model to the data with hierarchy level 3 and beyond. The
#'  approach of setting up the priors is same as described earlier (see
#'  \code{a_prior_sd}).
#'
#'@param c_prior_sd_str Specify priors for the random effect parameter, \code{c}
#'  when fitting a model to the data with hierarchy level 3 and beyond. The
#'  approach of setting up the priors is same as described earlier (see
#'  \code{a_prior_sd}).
#'
#'@param d_prior_sd_str Specify priors for the random effect parameter, \code{d}
#'  when fitting a model to the data with hierarchy level 3 and beyond. The
#'  approach of setting up the priors is same as described earlier (see
#'  \code{a_prior_sd}). 
#'  
#'@param e_prior_sd_str Specify priors for the random effect parameter, \code{e}
#'  when fitting a model to the data with hierarchy level 3 and beyond. The
#'  approach of setting up the priors is same as described earlier (see
#'  \code{a_prior_sd}).
#'  
#'@param f_prior_sd_str Specify priors for the random effect parameter, \code{f}
#'  when fitting a model to the data with hierarchy level 3 and beyond. The
#'  approach of setting up the priors is same as described earlier (see
#'  \code{a_prior_sd}).
#'  
#'@param g_prior_sd_str Specify priors for the random effect parameter, \code{g}
#'  when fitting a model to the data with hierarchy level 3 and beyond. The
#'  approach of setting up the priors is same as described earlier (see
#'  \code{a_prior_sd}).
#'  
#'@param h_prior_sd_str Specify priors for the random effect parameter, \code{h}
#'  when fitting a model to the data with hierarchy level 3 and beyond. The
#'  approach of setting up the priors is same as described earlier (see
#'  \code{a_prior_sd}).
#'  
#'@param i_prior_sd_str Specify priors for the random effect parameter, \code{i}
#'  when fitting a model to the data with hierarchy level 3 and beyond. The
#'  approach of setting up the priors is same as described earlier (see
#'  \code{a_prior_sd}).
#'
#'@param a_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{a} when fitting a model to the data with
#'  hierarchy level 3 and beyond. The approach is same as described above for
#'  the \code{a_cov_prior_sd}.
#'
#'@param b_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{b} when fitting a model to the data with
#'  hierarchy level 3 and beyond. The approach is same as described above for
#'  the \code{a_cov_prior_sd}.
#'
#'@param c_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{c} when fitting a model to the data with
#'  hierarchy level 3 and beyond. The approach is same as described above for
#'  the \code{a_cov_prior_sd}. 
#'
#'@param d_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{d} when fitting a model to the data with
#'  hierarchy level 3 and beyond. The approach is same as described above for
#'  the \code{a_cov_prior_sd}.
#'  
#'@param e_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{e} when fitting a model to the data with
#'  hierarchy level 3 and beyond. The approach is same as described above for
#'  the \code{a_cov_prior_sd}.
#'  
#'@param f_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{f} when fitting a model to the data with
#'  hierarchy level 3 and beyond. The approach is same as described above for
#'  the \code{a_cov_prior_sd}.  
#'  
#'@param g_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{g} when fitting a model to the data with
#'  hierarchy level 3 and beyond. The approach is same as described above for
#'  the \code{a_cov_prior_sd}.
#'  
#'@param h_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{h} when fitting a model to the data with
#'  hierarchy level 3 and beyond. The approach is same as described above for
#'  the \code{a_cov_prior_sd}.
#'  
#'@param i_cov_prior_sd_str Specify priors for the covariate(s) included in the
#'  random effect parameter, \code{i} when fitting a model to the data with
#'  hierarchy level 3 and beyond. The approach is same as described above for
#'  the \code{a_cov_prior_sd}.
#'  
#'@param sigma_prior_beta Specify priors for the fixed effect distributional
#'  parameter, \code{sigma}. The approach is same as described earlier for the 
#'  fixed effect parameter, \code{a}. See \code{a_prior_beta} for details.
#'
#'@param sigma_cov_prior_beta Specify priors for the covariate(s) included in
#'  the fixed effect distributional parameter, \code{sigma}. The approach is
#'  same as described earlier for the covariate(s) included in the fixed effect
#'  parameter, \code{a}. See \code{a_cov_prior_beta} for details.
#'
#'@param sigma_prior_sd Specify priors for the random effect distributional
#'  parameter, \code{sigma}. The approach of setting up the priors is same as
#'  described earlier the random effect parameter \code{a} (see
#'  \code{a_prior_sd}).
#'
#'@param sigma_cov_prior_sd Specify priors for the covariate(s) included in the
#'  random effect distributional parameter, \code{sigma}. The approach of
#'  setting up the priors is same as described earlier the covariate(s) included
#'  in the random effect parameter \code{a} (see \code{a_cov_prior_sd}).
#'  
#'@param sigma_prior_sd_str Specify priors for the the random effect
#'  distributional parameter, \code{sigma} when fitting a model to the data with
#'  hierarchy level 3 and beyond. The approach is same as described above for
#'  the random effect parameter, \code{a} (See \code{a_prior_sd_str}).
#'  
#'@param sigma_cov_prior_sd_str Specify priors for the covariate(s) included in
#'  the random effect distributional parameter, \code{sigma} when fitting a
#'  model to the data with hierarchy level 3 and beyond. The approach is same as
#'  described above for the random effect parameter, \code{a} (See
#'  \code{a_cov_prior_sd_str}).
#'
#'@param rsd_prior_sigma Specify priors for the residual standard deviation
#'  parameter \code{sigma}. This argument will only be evaluated if
#'  \code{dpar_formual} is set to \code{NULL}. For location scale based
#'  distributions, user can use specify the standard deviation (sd) or the
#'  median absolute deviation (mad) as scale parameter.
#'
#'@param dpar_prior_sigma Specify priors for the residual standard deviation
#'  parameter \code{sigma}. The argument is evaluated only when
#'  \code{dpar_formual} is not set to \code{NULL}.
#'
#'@param dpar_cov_prior_sigma Specify priors for the covariate(s) included in
#'  the residual standard deviation parameter \code{sigma}. The argument is
#'  evaluated only when \code{dpar_formual} is not set to \code{NULL}.
#'
#'@param autocor_prior_acor Specify priors for the the autocorrelation
#'  parameters (i.e., \code{ar} and \code{ma} parameters, see
#'  \code{autocor_formula} for details). The only allowed distribution is
#'  uniform distribution bounded between -1 and + 1. For the recently added 
#'  unstructured residual correlation, the allowed prior is \code{LKJ}. For 
#'  this unstructured residual correlation structure, a separate argument  
#'  \code{autocor_prior_unstr_acor} is included (see below). 
#'  
#' @param autocor_prior_unstr_acor Specify priors for the unstructured
#' residual autocorrelation structure. The allowed prior distribution is 
#' \code{LKJ}. See \code{gr_prior_cor} for details on \code{LKJ} prior.
#'
#'@param gr_prior_cor Specify priors for the the correlations of group-level
#'  random effects. The allowed distribution is \code{LKJ} which has a
#'  single parameter \code{eta} that sets priors on correlation parameters (see
#'  \code{brms::prior} for details).'
#'  
#'@param gr_prior_cor_str Specify priors for the the correlations of group-level
#'  random effects when fitting a model to the data with hierarchy level 3 and
#'  beyond. The approach is same as described above for the correlations of
#'  group-level random effects (See \code{gr_prior_cor}).
#'  
#'@param sigma_prior_cor Specify priors for the correlations of random effects
#'  for the distributional parameter \code{sigma}. The allowed distribution is
#'  is \code{LKJ} (see \code{gr_prior_cor}). Note that currently
#'  \code{brms::brm()} does not allow for setting separate \code{LKJ} priors for
#'  the distribution and group level random effects that share the same group
#'  (because  \code{brms::brm()} does not assign group for sigma). Therefore,
#'  either create a copy of group identifier and use that but then this will not
#'  allow correlation parameter across group random effects and sigma. Another
#'  hack, which is used currently, is to remove group argument i.e., \code{group
#'  = ""} when evaluating the \code{sigma_prior_cor}. See the relevant portion
#'  of code \code{if(sigma_dpar == "sigma") group <- ""} in function
#'  \code{set_priors_initials (line 2095)}.
#'  
#'@param sigma_prior_cor_str Specify priors for the the correlations of random
#'  effects for distributional parameter \code{sigma} when fitting a model to
#'  the data with hierarchy level 3 and beyond. The approach is same as
#'  described above for the correlations of distributional parameter random
#'  effects (See \code{sigma_prior_cor}).
#'
#'@param mvr_prior_rescor Specify priors for the residual correlation parameter
#'  when fitting a multivariate model. The allowed distribution is \code{LKJ} 
#'  (see \code{gr_prior_cor}).
#'
#'@param init Specify initial values for the sampler. For \code{0}, all
#'  parameters are initialized to zero. If \code{random}, Stan will randomly
#'  generate initial values for each parameter in a range specified by the
#'  \code{init_r} (see below). Another option is \code{prior} which allows
#'  setting initials based on the prior specified for each parameter. Lastly,
#'  \code{NULL} option (the default) will let all the following init arguments
#'  to be evaluated and set initial for each parameter as specified by the
#'  corresponding initial argument.
#'
#'@param init_r Set range for the random generation of initial values. This
#'  argument is evaluated only when the \code{init} is set as \code{"random"} 
#'  (see above).
#'
#'@param a_init_beta Specify initial values for the fixed effect parameter,
#'  \code{a}. Options available are \code{0}, \code{random} and \code{prior}. In
#'  addition, user can specify \code{ymean} and \code{ymedian} to set initial as
#'  the mean or median of the response variable, or option \code{lm} that sets
#'  initials obtained from the simple linear model fitted to the data (similar
#'  to the location parameter, see \code{a_prior_beta}). Note that these options
#'  are available only for the fixed effect parameter \code{a} and not for other
#'  parameters described below. Lastly, For univariate-by-subgroup model (see
#'  \code{univariate_by}) and multivariate (see \code{multivariate}) models, the
#'  initials can be same (e.g., \code{a_init_beta = 0}) for each response
#'  variable different for each response (e.g., \code{list(a_init_beta = 0,
#'  a_init_beta = lm)}).
#'
#'@param b_init_beta Specify initial values for the fixed effect parameter,
#'  \code{b}. See \code{a_init_beta} for details. 
#'
#'@param c_init_beta Specify initial values for the fixed effect parameter,
#'  \code{c}. See \code{a_init_beta} for details. 
#'
#'@param d_init_beta Specify initial values for the fixed effect parameter,
#'  \code{d}. See \code{a_init_beta} for details.
#'  
#'@param e_init_beta Specify initial values for the fixed effect parameter,
#'  \code{e}. See \code{a_init_beta} for details. 
#'  
#'@param f_init_beta Specify initial values for the fixed effect parameter,
#'  \code{f}. See \code{a_init_beta} for details.
#'  
#'@param g_init_beta Specify initial values for the fixed effect parameter,
#'  \code{g}. See \code{a_init_beta} for details.
#'  
#'@param h_init_beta Specify initial values for the fixed effect parameter,
#'  \code{h}. See \code{a_init_beta} for details.
#'  
#'@param i_init_beta Specify initial values for the fixed effect parameter,
#'  \code{i}. See \code{a_init_beta} for details.
#'
#'@param s_init_beta  Specify initial values for the fixed effect parameter,
#'  \code{s} i.e., spline coefficients. Options available are \code{0},
#'  \code{random} and \code{prior}, and \code{lm}. See \code{a_init_beta} for
#'  details.
#'
#'@param a_cov_init_beta Specify initial values for covariate(s) included in the
#'  fixed effect parameter, \code{a}. parameter. Options available are \code{0},
#'  \code{random}, \code{prior} and \code{lm}. See \code{a_init_beta} for
#'  details. The \code{lm} is available only for the \code{a_cov_init_beta} and
#'  not for the covariate(s) for other fixed effect parameters.
#'
#'@param b_cov_init_beta Specify initial values for covariate(s) included in the
#'  fixed effect parameter, \code{b}. See \code{a_cov_init_beta} for details.
#'
#'@param c_cov_init_beta Specify initial values for covariate(s) included in the
#'  fixed effect parameter, \code{c}. See \code{a_cov_init_beta} for details.
#'
#'@param d_cov_init_beta Specify initial values for covariate(s) included in the
#'  fixed effect parameter, \code{d}. See \code{a_cov_init_beta} for details.
#'  
#'@param e_cov_init_beta Specify initial values for covariate(s) included in the
#'  fixed effect parameter, \code{e}. See \code{a_cov_init_beta} for details.
#'  
#'@param f_cov_init_beta Specify initial values for covariate(s) included in the
#'  fixed effect parameter, \code{f}. See \code{a_cov_init_beta} for details.
#'  
#'@param g_cov_init_beta Specify initial values for covariate(s) included in the
#'  fixed effect parameter, \code{g}. See \code{a_cov_init_beta} for details.
#'  
#'@param h_cov_init_beta Specify initial values for covariate(s) included in the
#'  fixed effect parameter, \code{h}. See \code{a_cov_init_beta} for details.
#'  
#'@param i_cov_init_beta Specify initial values for covariate(s) included in the
#'  fixed effect parameter, \code{i}. See \code{a_cov_init_beta} for details.
#'
#'@param s_cov_init_beta Specify initial values for covariate(s) included in the
#'  fixed effect parameter, \code{s} (spline coefficients). See
#'  \code{a_cov_init_beta} for details. The option \code{lm} will set the
#'  initial obtained from the simple linear model fit to the data.
#'
#'@param a_init_sd Specify initial values for the random effect parameter,
#'  \code{a}. Options available are \code{0}, \code{random} and \code{prior} as
#'  described above for setting initials for the fixed effect parameters.
#'  Additional options available are \code{ysd}, \code{ymad}, \code{lme_sd_a},
#'  and \code{lm_sd_a}. The \code{ysd} and \code{ymad} options set standard
#'  deviation (sd) variable and the median absolute deviation (mad) of the the
#'  response as initial value. The  \code{lme_sd_a} sets initial value based on
#'  the standard deviation of intercept from the linear mixed model
#'  (\code{nlme::lme()}) applied to the data. The initial value set by the
#'  \code{lm_sd_a} is square root of the residual variance from the simple
#'  linear model applied to the data. Note that in case \code{nlme::lme()} fails
#'  to converge for some reasons, then option \code{lm_sd_a} is set
#'  automatically. Also note that additional options are available only for the
#'  andom effect parameter \code{a}. See \code{a_init_beta} for further details
#'  like setting up the initials for the univariate-by-subgroup
#'  (\code{univariate_by}) and multivariate (\code{multivariate}) models.
#'
#'@param b_init_sd Specify initial values for the random effect parameter,
#'  \code{b}. See \code{a_init_sd} for details.
#'
#'@param c_init_sd Specify initial values for the random effect parameter,
#'  \code{c}. See \code{a_init_sd} for details.
#'
#'@param d_init_sd Specify initial values for the random effect parameter,
#'  \code{d}. See \code{a_init_sd} for details.
#'  
#'@param e_init_sd Specify initial values for the random effect parameter,
#'  \code{e}. See \code{a_init_sd} for details.
#'
#'@param f_init_sd Specify initial values for the random effect parameter,
#'  \code{f}. See \code{a_init_sd} for details.
#'  
#'@param g_init_sd Specify initial values for the random effect parameter,
#'  \code{g}. See \code{a_init_sd} for details.
#'  
#'@param h_init_sd Specify initial values for the random effect parameter,
#'  \code{h}. See \code{a_init_sd} for details.
#'  
#'@param i_init_sd Specify initial values for the random effect parameter,
#'  \code{i}. See \code{a_init_sd} for details.
#'
#'@param a_cov_init_sd Specify initial values for the covariate(s) included in
#'  the random effect parameter, \code{a}. Options available are \code{0},
#'  \code{random} and \code{prior}. See \code{a_cov_init_beta} for further
#'  details on setting up the initials for the univariate-by-subgroup
#'  (\code{univariate_by}) and multivariate (\code{multivariate}) models.
#'
#'@param b_cov_init_sd Specify initial values for the covariate(s) included in
#'  the random effect parameter, \code{b}. See \code{a_cov_init_sd} for details.
#'
#'@param c_cov_init_sd Specify initial values for the covariate(s) included in
#'  the random effect parameter, \code{c}. See \code{a_cov_init_sd} for details.
#'
#'@param d_cov_init_sd Specify initial values for the covariate(s) included in
#'  the random effect parameter, \code{d}. See \code{a_cov_init_sd} for details.
#'
#'@param e_cov_init_sd Specify initial values for the covariate(s) included in
#'  the random effect parameter, \code{e}. See \code{a_cov_init_sd} for details.
#'
#'@param f_cov_init_sd Specify initial values for the covariate(s) included in
#'  the random effect parameter, \code{f}. See \code{a_cov_init_sd} for details.
#'
#'@param g_cov_init_sd Specify initial values for the covariate(s) included in
#'  the random effect parameter, \code{g}. See \code{a_cov_init_sd} for details.
#'  
#'@param h_cov_init_sd Specify initial values for the covariate(s) included in
#'  the random effect parameter, \code{h}. See \code{a_cov_init_sd} for details.
#'  
#'@param i_cov_init_sd Specify initial values for the covariate(s) included in
#'  the random effect parameter, \code{i}. See \code{a_cov_init_sd} for details.
#'
#'@param sigma_init_beta Specify initial values for the fixed effect
#'  distributional parameter, \code{sigma}. The approach is same as described
#'  earlier for the fixed effect parameter \code{a}.See \code{a_init_beta} for
#'  details.
#'
#'@param sigma_cov_init_beta Specify initial values for the covariate(s)
#'  included in the fixed effect distributional parameter, \code{sigma}. See
#'  \code{a_cov_init_beta} for details.
#'
#'@param sigma_init_sd Specify initial values for the random effect
#'  distributional parameter, \code{sigma}. The approach is same as described
#'  earlier for the random effect parameter \code{a}.See \code{a_init_sd} for
#'  details.
#'
#'@param sigma_cov_init_sd Specify initial values for the covariate(s) included
#'  in the random effect distributional parameter, \code{sigma}. See
#'  \code{a_cov_init_sd} for details.
#'
#'@param gr_init_cor Specify initial values for correlations of group-level
#'  random effects parameters. Allowed options are \code{0}, \code{random} and
#'  \code{prior}.
#'
#'@param sigma_init_cor Specify initial values for correlations of
#'  distributional random effects parameter \code{sigma}. Allowed options are
#'  \code{0}, \code{random} and \code{prior}.
#'
#'@param rsd_init_sigma Specify initial values for residual standard deviation
#'  parameter, \code{sigma}. Options available are \code{0}, \code{random} and
#'  \code{prior} as described above for setting initials for the fixed effect
#'  parameters. Additional options available are \code{lme_rsd}, and
#'  \code{lm_rsd}. The \code{lme_rsd} sets initial value based on the standard
#'  deviation of residuals obatined from the linear mixed model
#'  (\code{nlme::lme()}) applied to the data. The initial value set by the
#'  \code{lm_rsd} is square root of the residual variance from the simple linear
#'  model applied to the data. Note that in case \code{nlme::lme()} fails to
#'  converge for some reasons, then option \code{lm_sd_a} is set automatically.
#'  The argument \code{rsd_init_sigma} is evaluated when \code{dpar_formual} is
#'  set to \code{NULL}.
#'
#'@param dpar_init_sigma Specify initial values for the distributional parameter
#'  \code{sigma}. The approach and options available are same as described above
#'  for the \code{rsd_init_sigma}. This argument is evaluated only when
#'  \code{dpar_formual} is not set to NULL.
#'
#'@param dpar_cov_init_sigma Specify initial values for the covariate(s)
#'  included in the distributional parameter, \code{sigma}. Allowed options are
#'  \code{0}, \code{random}, and \code{prior}.
#'
#'@param autocor_init_acor Specify initial values for autocorrelation parameter.
#'  see \code{autocor_formula} for details). Allowed options are \code{0},
#'  \code{random}, and \code{prior}.
#'  
#'@param autocor_init_unstr_acor Specify initial values for unstructured
#'  residual autocorrelation parameter. Allowed options are \code{0},
#'  \code{random}, and \code{prior}. Note that the procedure to set initials
#'  using \code{autocor_init_acor} is identical to the \code{gr_init_cor}.
#'
#'@param mvr_init_rescor Specify initial values for the residual correlations
#'  for multivariate (\code{multivariate}) model. Allowed options are \code{0},
#'  \code{random}, and \code{prior}.
#'
#'@param r_init_z Specify initial values for the standardized group level
#'  effects. These parameters are part of the central parameterisation  approach
#'  adopted by the the brms package (see [brms::brm()] for details).
#'  
#'@param vcov_init_0 A logical (default \code{TRUE}) to set initials for  
#'  variance covariance (i.e, standard deviation and correlation) parameters as
#'  zero. This allows for setting initials for the fixed effects and zero for 
#'  random effects. 
#'
#'@param jitter_init_beta A value as proportion (between 0 and 1) used to
#'  perturb the initials for fixed effect parameters. The default is \code{NULL}
#'  indicating that same initials are used across for all chains. An option of 
#'  setting \code{jitter_init_beta = 0.1} looked good during early testing. 
#'
#'@param jitter_init_sd A value as proportion (between 0 and 1) used to perturb
#'  the initials for standard deviation of random effect parameters. The default
#'  is \code{NULL} indicating that same initials are used across for all chains.
#'  An option of setting \code{jitter_init_beta = 0.01} looked good during early
#'  testing.
#'
#'@param jitter_init_cor A value as proportion (between 0 and 1) used to perturb
#'  the initials for correlation of random effect parameters. The default is
#'  \code{NULL} indicating that same initials are used across for all chains. An
#'  option of setting \code{jitter_init_beta = 0.001} looked good during early
#'  testing.
#'
#'@param prior_data An optional argument (a named list) pass values that can be
#'  used in the prior arguments The default is \code{NULL}. This option is
#'  particularly helpful when passing a long vector or a matrix for setting
#'  priors for covariate(s). These vectors and matrices can be created in the R
#'  framework and then passed using the \code{prior_data}. For example, to pass
#'  a vector of location parameters when setting priors for covariates  with 10
#'  dummy variables, one can create a named object prior_a_cov_location
#'  (\code{prior_a_cov_location = rnorm(10, 0, 1)}) and prior_a_cov_scale
#'  (\code{prior_a_cov_scale = rnorm(5, 0, 1)})  and then pass it as a named
#'  list \code{prior_data = list(prior_a_cov_location = prior_a_cov_location,
#'  prior_a_cov_scale = prior_a_cov_scale,)} and specifying that to the
#'  \code{a_cov_prior_beta} as \code{a_cov_prior_beta =
#'  normal(prior_a_cov_location, prior_a_cov_scale)}.
#'
#'@param init_data An optional argument (a named list) pass values that can be
#'  used in the initial arguments. The approach is exact same as described above
#'  for the \code{prior_data}.
#'  
#'@param init_custom Set a custom initials object (a named list (default
#'  \code{init}). Note that the named list is directly passed to the arguments
#'  without checking the dimensions.
#'
#'@param verbose An optional logical (default \code{FALSE}) argument to print
#'  important information during the steps involved in preparing model formula,
#'  Stan function, priors, initials. As an example, a user might be interested
#'  in knowing the responses created for the factor variable that is used to
#'  specify the univariate-by-subgroup model. This information then helps in
#'  matching the desired sequence of options used to pass on to arguments such
#'  as df, prior, initials etc.
#'
#'@param expose_function An optional logical (default TRUE) to expose Stan
#'  function used for model fitting. These functions are essential for
#'  post-processing.
#'  
#'@param get_stancode An optional logical (default \code{FALSE}) to get the
#'  stancode.
#'
#'@param get_standata An optional logical (default \code{FALSE}) to get the
#'  standata.
#'
#'@param get_formula An optional logical (default \code{FALSE}) to get formula.
#'
#'@param get_stanvars An optional logical (default \code{FALSE}) to get
#'  stanvars.
#'
#'@param get_priors An optional logical (default \code{FALSE}) to get priors.
#'
#'@param get_set_priors An optional logical (default \code{FALSE}) to get priors
#'  specified by the \code{bsitar} via \code{prepare_priors}.
#'
#'@param validate_priors An optional logical (default \code{FALSE}) to
#'  validate the specified priors..
#'  
#'@param get_set_init An optional logical (default \code{FALSE}) to get 
#' initials specified by the \code{bsitar} via \code{prepare_initials}.
#'
#'@param set_self_priors An optional (default \code{NULL}) to specify
#'  priors manually.
#'  
#'@param set_replace_priors An optional (default \code{NULL}) to replace
#'  part of prior object.
#' 
#'@param set_same_priors_hierarchy An optional (default \code{NULL}) to replace
#' part of prior object.
#'  
#'@param outliers An optional (default \code{NULL}) to remove velocity
#' outliers. The argument should be a named list to pass on to the
#' [bsitar::outliers] function. See [bsitar::outliers] for details.
#'
#'@param cores Number of cores to be used when executing the chains in parallel.
#'  See [brms::brm()] for details. Note that unlike [brms::brm()] which sets
#'  \code{cores=getOption("mc.cores", 1)}, the default in \code{bsitar} is
#'  \code{cores=getOption("mc.cores", 'optimize')} which optimizes the
#'  utilization of system resources. The maximum number of cores that can be
#'  deployed is calculated as the maximum number of available cores minus 1.
#'  When the number of available cores is greater than the number of chains (see
#'  \code{chains}), then number of cores is set equal to the number of chains.
#'  Another option is to set \code{cores} as \code{getOption("mc.cores",
#'  'maximise')} which sets the number of cores as the maximum number of cores
#'  available from the system regardless of the number of chains specified. Note
#'  that the user can also set \code{cores} argument similar to the
#'  [brms::brm()] i.e., \code{getOption("mc.cores", 1)}. All these three options
#'  can be set globally as \code{options(mc.cores = x}) where x can be
#'  \code{optimize}, \code{maximise} or \code{1}. Lastly, the \code{cores} can
#'  set by directly by specifying an integer e.g., \code{cores= 4}.
#'  
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
#'  Alternatively, the number of threads can be set as \code{threads
#'  = threading(x)} where \code{X} is an integer. Other arguments that can the
#'  passed to the \code{threads} are \code{grainsize} and the \code{static}. See
#'  [brms::brm()] for further details on within-chain parallelization.
#'
#'@param normalize Indicates whether normalization constants should be included
#'  in the Stan code (default \code{TRUE}). Setting it to \code{FALSE} requires
#'  Stan version >= 2.25. If \code{FALSE}, sampling efficiency may be increased
#'  but some post processing functions such as [brms::bridge_sampler()] will not
#'  be available. This option can be controlled globally via the
#'  \code{brms.normalize} option.
#' 
#'
#'@param sample_prior Indicates whether to draw sample from priors in addition
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
#'  this case, all parameters must have proper priors.
#'
#'@param save_model A character string or \code{NULL} (default). If not
#'  \code{NULL}, then the model's Stan code is saved via in a text file named
#'  after the string supplied in \code{save_model}.
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
#'
#' @param unused An optional formula which contains variables that are unused 
#' in the model but should still be stored in the model's data frame. 
#' This can be useful, for example, if those variables are required for 
#' post-processing the model.
#' 
#' @param chains Number of Markov chains (defaults to 4).
#' 
#' @param iter Number of total iterations per chain, including warmup (defaults
#'   2000)
#' 
#' @param warmup A positive integer specifying number of warmup (aka burnin)
#'   iterations. This also specifies the number of iterations used for stepsize
#'   adaptation, so warmup draws should not be used for inference. The number
#'   of warmup should not be larger than \code{iter} and the default is
#'   \code{iter/2}.
#'   
#' @param thin Thinning rate. Must be a positive integer. Set \code{thin > 1} to
#'   save memory and computation time if \code{iter} is large.
#'   
#' @param backend Character string naming the package to use as the backend for
#'   fitting the Stan model. Options are \code{"rstan"} (the default) or
#'   \code{"cmdstanr"}. Can be set globally for the current \R session via the
#'   \code{"brms.backend"} option (see \code{\link{options}}). Details on the
#'   \pkg{rstan} and \pkg{cmdstanr} packages are available at
#'   \url{https://mc-stan.org/rstan/} and \url{https://mc-stan.org/cmdstanr/},
#'   respectively. Additionally a \code{"mock"} backend is available to make
#'   testing \pkg{brms} and packages that depend on it easier. The \code{"mock"}
#'   backend does not actually do any fitting, it only checks the generated Stan
#'   code for correctness and then returns whatever is passed in an additional
#'   \code{mock_fit} argument as the result of the fit.
#'   
#' @param opencl The platform and device IDs of the OpenCL device to use for
#'   fitting using GPU support. If you don't know the IDs of your OpenCL device,
#'   \code{c(0,0)} is most likely what you need. For more details, see
#'   \code{\link{opencl}}. Can be set globally for the current \R session via
#'   the \code{"brms.opencl"} option.
#'
#' @param algorithm Character string naming the estimation approach to use.
#'   Options are \code{"sampling"} for MCMC (the default), \code{"meanfield"}
#'   for variational inference with independent normal distributions,
#'   \code{"fullrank"} for variational inference with a multivariate normal
#'   distribution, or \code{"fixed_param"} for sampling from fixed parameter
#'   values. Can be set globally for the current \R session via the
#'   \code{"brms.algorithm"} option (see \code{\link{options}}).
#'   
#' @param save_pars An object generated by \code{\link{save_pars}} controlling
#'   which parameters should be saved in the model. The argument has no
#'   impact on the model fitting itself.
#'   
#' @param drop_unused_levels Should unused factors levels in the data be
#'   dropped? Defaults to \code{TRUE}.
#'   
#' @param stan_model_args A \code{list} of further arguments passed to
#'   \code{\link[rstan:stan_model]{rstan::stan_model}} for \code{backend =
#'   "rstan"} or to \code{cmdstanr::cmdstan_model} for \code{backend =
#'   "cmdstanr"}, which allows to change how models are compiled.
#'   
#' @param silent Verbosity level between \code{0} and \code{2}. If \code{1} (the
#'   default), most of the informational messages of compiler and sampler are
#'   suppressed. If \code{2}, even more messages are suppressed. The actual
#'   sampling progress is still printed. Set \code{refresh = 0} to turn this off
#'   as well. If using \code{backend = "rstan"} you can also set
#'   \code{open_progress = FALSE} to prevent opening additional progress bars.
#'   
#' @param seed The seed for random number generation to make results
#'   reproducible. If \code{NA} (the default), \pkg{Stan} will set the seed
#'   randomly.
#'   
#' @param fit An instance of S3 class \code{brmsfit} derived from a previous
#'   fit; defaults to \code{NA}. If \code{fit} is of class \code{brmsfit}, the
#'   compiled model associated with the fitted result is re-used and all
#'   arguments modifying the model code or data are ignored. It is not
#'   recommended to use this argument directly, but to call the
#'   \code{\link[brms:update.brmsfit]{update}} method, instead.
#'   
#' @param file Either \code{NULL} or a character string. In the latter case, the
#'   fitted model object is saved via \code{\link{saveRDS}} in a file named
#'   after the string supplied in \code{file}. The \code{.rds} extension is
#'   added automatically. If the file already exists, \code{brm} will load and
#'   return the saved model object instead of refitting the model.
#'   Unless you specify the \code{file_refit} argument as well, the existing
#'   files won't be overwritten, you have to manually remove the file in order
#'   to refit and save the model under an existing file name. The file name
#'   is stored in the \code{brmsfit} object for later usage.
#'   
#' @param file_refit Modifies when the fit stored via the \code{file} argument
#'   is re-used. Can be set globally for the current \R session via the
#'   \code{"brms.file_refit"} option (see \code{\link{options}}).
#'   For \code{"never"} (default) the fit is always loaded if it
#'   exists and fitting is skipped. For \code{"always"} the model is always
#'   refitted. If set to \code{"on_change"}, brms will
#'   refit the model if model, data or algorithm as passed to Stan differ from
#'   what is stored in the file. This also covers changes in priors,
#'   \code{sample_prior}, \code{stanvars}, covariance structure, etc. If you
#'   believe there was a false positive, you can use
#'   \code{\link{brmsfit_needs_refit}} to see why refit is deemed necessary.
#'   Refit will not be triggered for changes in additional parameters of the fit
#'   (e.g., initial values, number of iterations, control arguments, ...). A
#'   known limitation is that a refit will be triggered if within-chain
#'   parallelization is switched on/off.
#'   
#' @param future Logical; If \code{TRUE}, the \pkg{\link[future:future]{future}}
#'   package is used for parallel execution of the chains and argument
#'   \code{cores} will be ignored. Can be set globally for the current \R
#'   session via the \code{"future"} option. The execution type is controlled
#'   via \code{\link[future:plan]{plan}} (see the examples section below).
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
#'@importFrom stats loess na.omit residuals complete.cases deriv formula update
#' 
#'@importFrom utils combn head installed.packages packageVersion tail
#'@importFrom Rdpack reprompt
#'@import brms
#'
#'@export
#'
#'
bsitar <- function(x,
                   y,
                   id,
                   data,
                   df = 4,
                   knots = NA,
                   fixed = a + b + c + d + e + f,
                   random = a + b + c + d + e + f,
                   select_model = 'sitar',
                   xoffset = mean,
                   bstart = mean,
                   apgv = 13,
                   pgv = 4,
                   xfun = NULL,
                   yfun = NULL,
                   bound = 0.04,
                   terms_rhs = NULL,
                   
                   a_formula = ~ 1,
                   b_formula = ~ 1,
                   c_formula = ~ 1,
                   d_formula = ~ 1,
                   e_formula = ~ 1,
                   f_formula = ~ 1,
                   g_formula = ~ 1,
                   h_formula = ~ 1,
                   i_formula = ~ 1,
                   
                   s_formula = ~ 1,
                   
                   a_formula_gr = ~ 1,
                   b_formula_gr = ~ 1,
                   c_formula_gr = ~ 1,
                   d_formula_gr = ~ 1,
                   e_formula_gr = ~ 1,
                   f_formula_gr = ~ 1,
                   g_formula_gr = ~ 1,
                   h_formula_gr = ~ 1,
                   i_formula_gr = ~ 1,
                   
                   a_formula_gr_str = NULL,
                   b_formula_gr_str = NULL,
                   c_formula_gr_str = NULL,
                   d_formula_gr_str = NULL,
                   e_formula_gr_str = NULL,
                   f_formula_gr_str = NULL,
                   g_formula_gr_str = NULL,
                   h_formula_gr_str = NULL,
                   i_formula_gr_str = NULL,
                   
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
                   
                   a_prior_beta = normal(lm, ysd, autoscale = FALSE),
                   b_prior_beta = normal(0, 2, autoscale = FALSE),
                   c_prior_beta = normal(0, 0.5, autoscale = FALSE),
                   d_prior_beta = normal(0, 1, autoscale = FALSE),
                   e_prior_beta = normal(0, 1, autoscale = FALSE),
                   f_prior_beta = normal(0, 1, autoscale = FALSE),
                   g_prior_beta = normal(0, 1, autoscale = FALSE),
                   h_prior_beta = normal(0, 1, autoscale = FALSE),
                   i_prior_beta = normal(0, 1, autoscale = FALSE),
                   
                   s_prior_beta = normal(lm, lm, autoscale = FALSE),
                   
                   a_cov_prior_beta = normal(0, 5, autoscale = FALSE),
                   b_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                   c_cov_prior_beta = normal(0, 0.1, autoscale = FALSE),
                   d_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                   e_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                   f_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                   g_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                   h_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                   i_cov_prior_beta = normal(0, 1, autoscale = FALSE),
                   
                   s_cov_prior_beta = normal(0, 10, autoscale = FALSE),
                   
                   a_prior_sd = normal(0, ysd, autoscale = 1),
                   b_prior_sd = normal(0, 1, autoscale = FALSE),
                   c_prior_sd = normal(0, 0.25, autoscale = FALSE),
                   d_prior_sd = normal(0, 1, autoscale = FALSE),
                   e_prior_sd = normal(0, 1, autoscale = FALSE),
                   f_prior_sd = normal(0, 1, autoscale = FALSE),
                   g_prior_sd = normal(0, 1, autoscale = FALSE),
                   h_prior_sd = normal(0, 1, autoscale = FALSE),
                   i_prior_sd = normal(0, 1, autoscale = FALSE),
                   
                   a_cov_prior_sd = normal(0, 2, autoscale = FALSE),
                   b_cov_prior_sd = normal(0, 1, autoscale = FALSE),
                   c_cov_prior_sd = normal(0, 0.15, autoscale = FALSE),
                   d_cov_prior_sd = normal(0, 0.5, autoscale = FALSE),
                   e_cov_prior_sd = normal(0, 0.5, autoscale = FALSE),
                   f_cov_prior_sd = normal(0, 0.5, autoscale = FALSE),
                   g_cov_prior_sd = normal(0, 0.5, autoscale = FALSE),
                   h_cov_prior_sd = normal(0, 0.5, autoscale = FALSE),
                   i_cov_prior_sd = normal(0, 0.5, autoscale = FALSE),
                   
                   a_prior_sd_str = NULL,
                   b_prior_sd_str = NULL,
                   c_prior_sd_str = NULL,
                   d_prior_sd_str = NULL,
                   e_prior_sd_str = NULL,
                   f_prior_sd_str = NULL,
                   g_prior_sd_str = NULL,
                   h_prior_sd_str = NULL,
                   i_prior_sd_str = NULL,
                   
                   a_cov_prior_sd_str = NULL,
                   b_cov_prior_sd_str = NULL,
                   c_cov_prior_sd_str = NULL,
                   d_cov_prior_sd_str = NULL,
                   e_cov_prior_sd_str = NULL,
                   f_cov_prior_sd_str = NULL,
                   g_cov_prior_sd_str = NULL,
                   h_cov_prior_sd_str = NULL,
                   i_cov_prior_sd_str = NULL,
                   
                   sigma_prior_beta = normal(0, 1, autoscale = FALSE),
                   sigma_cov_prior_beta = normal(0, 0.5, autoscale = FALSE),
                   sigma_prior_sd = normal(0, 0.25, autoscale = FALSE),
                   sigma_cov_prior_sd = normal(0, 0.15, autoscale = FALSE),
                   
                   sigma_prior_sd_str = NULL,
                   sigma_cov_prior_sd_str = NULL,
                   
                   rsd_prior_sigma = normal(0, ysd, autoscale = FALSE),
                   dpar_prior_sigma = normal(0, ysd, autoscale = FALSE),
                   dpar_cov_prior_sigma = normal(0, 5, autoscale = FALSE),
                   autocor_prior_acor = uniform(-1, 1, autoscale = FALSE),
                   autocor_prior_unstr_acor = lkj(1),
                   
                   gr_prior_cor = lkj(1),
                   gr_prior_cor_str = lkj(1),
                   sigma_prior_cor = lkj(1),
                   sigma_prior_cor_str = lkj(1),
                   
                   mvr_prior_rescor = lkj(1),
                   init = NULL,
                   init_r = NULL,
                   a_init_beta = lm,
                   b_init_beta = 0,
                   c_init_beta = 0,
                   d_init_beta = 0,
                   e_init_beta = 0,
                   f_init_beta = 0,
                   g_init_beta = 0,
                   h_init_beta = 0,
                   i_init_beta = 0,
                   
                   s_init_beta = lm,
                   
                   a_cov_init_beta = 0,
                   b_cov_init_beta = 0,
                   c_cov_init_beta = 0,
                   d_cov_init_beta = 0,
                   e_cov_init_beta = 0,
                   f_cov_init_beta = 0,
                   g_cov_init_beta = 0,
                   h_cov_init_beta = 0,
                   i_cov_init_beta = 0,
                   
                   s_cov_init_beta = lm,
                   
                   a_init_sd = 1,
                   b_init_sd = 1,
                   c_init_sd = 1,
                   d_init_sd = 1,
                   e_init_sd = 1,
                   f_init_sd = 1,
                   g_init_sd = 1,
                   h_init_sd = 1,
                   i_init_sd = 1,
                   
                   a_cov_init_sd = 1,
                   b_cov_init_sd = 1,
                   c_cov_init_sd = 1,
                   d_cov_init_sd = 1,
                   e_cov_init_sd = 1,
                   f_cov_init_sd = 1,
                   g_cov_init_sd = 1,
                   h_cov_init_sd = 1,
                   i_cov_init_sd = 1,
                   
                   sigma_init_beta = 0.01,
                   sigma_cov_init_beta = 0.01,
                   sigma_init_sd = 1,
                   sigma_cov_init_sd = 1,
                   
                   gr_init_cor = 0,
                   sigma_init_cor = 0,
                   rsd_init_sigma = 1,
                   dpar_init_sigma = 1,
                   dpar_cov_init_sigma = 1,
                   autocor_init_acor = 0.5,
                   autocor_init_unstr_acor = 0,
                   mvr_init_rescor = 0,
                   r_init_z = 0,
                   vcov_init_0 = TRUE,
                   jitter_init_beta = NULL,
                   jitter_init_sd = NULL,
                   jitter_init_cor = NULL,
                   prior_data = NULL,
                   init_data = NULL,
                   init_custom = NULL,
                   verbose = FALSE,
                   expose_function = FALSE,
                   
                   get_stancode = FALSE,
                   get_standata = FALSE,
                   get_formula = FALSE,
                   get_stanvars = FALSE,
                   get_priors = FALSE,
                   get_set_priors = FALSE,
                   validate_priors = FALSE,
                   get_set_init = FALSE,
                   
                   set_self_priors = NULL,
                   set_replace_priors = NULL,
                   set_same_priors_hierarchy = FALSE,
                   outliers = NULL, 
                   unused = NULL,
                   
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
  
  
  # check and set alias argument for formuale 
  dots_allias <- list(...)
  collect_dot_names <- c()
  for (ia in letters[1:10]) {
    set_name_dot <- paste0(ia, ".", 'formula')
    set_name_uns <- paste0(ia, "_", 'formula')
    collect_dot_names <- c(collect_dot_names, set_name_dot)
    if (set_name_dot %in% names(dots_allias)) {
      # if (missing(a_formula)) {
      if (eval(bquote(missing(.(set_name_uns)))) ) { 
        mcall[[set_name_uns]] <- dots_allias[[set_name_dot]]
        # else if (!missing(a_formula)) {
      } else if (!eval(bquote(missing(.(set_name_uns)))) ) { 
        err_msg <- paste0("both '", set_name_uns, "' and '" , 
                          set_name_dot, "' found, ignoring '",set_name_dot, "'")
        if(verbose) warning(err_msg)
      }
    }
  }
  
  # Rememeber to remove these from the brms args also : brmsdots_  line 5552
  for (collect_dot_namesi in collect_dot_names) {
    if(!is.null(mcall[[collect_dot_namesi]])) 
      mcall[[collect_dot_namesi]] <- NULL
  }
  rm(dots_allias)
  
  
  mcall <- mcall_ <- mcall
  
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
  
  
  # Problem with rethinking occurs during the expose_bsitar_function
  if("rethinking" %in% (.packages())){
    message("Package 'rethinking' detached and unloaded as it creates conflict",
            " \nwith the rstan version ", utils::packageVersion('rstan'))
    detach("package:rethinking", unload=TRUE) 
  }
 
  
  
  # Some checks 
  if(utils::packageVersion('rstan') < 2.26) {
    if(expose_function) stop("Argument 'expose_function' not allowed ",
                             "for this rstan version ",
                             utils::packageVersion('rstan'))
  }
  
  
  # This all done to not let init = random to evaluate to random effect str 
  # i.e, "a+b+c..."
  
  temp_init_call_in <- mcall$init
  if(is.null(temp_init_call_in)) temp_init_call_c <- temp_init_call_in
  if(is.symbol(temp_init_call_in) | is.numeric(temp_init_call_in)) {
    if(!is.character(temp_init_call_in)) {
      temp_init_call_c <- deparse(temp_init_call_in)
    } else {
      temp_init_call_c <- temp_init_call_in
    }
  } else if(is.character(temp_init_call_in)) {
    temp_init_call_c <- temp_init_call_in
  }
  
  if(is.language(temp_init_call_in)) {
    temp_init_call <- deparse(temp_init_call_in)
    temp_init_call <- gsub("[[:space:]]", "", temp_init_call)
    temp_init_call <- regmatches(temp_init_call, 
                                 gregexpr("(?<=\\().*?(?=\\))", 
                                          temp_init_call, perl=T))[[1]]
    # this check if inits = list(xx = xx etc)
    if(length(temp_init_call) != 0) {
      temp_init_call <- strsplit(temp_init_call, ",")[[1]]
      temp_init_call_c <- c()
      for (temp_init_calli in temp_init_call) {
        if(!grepl("\"", temp_init_calli)) {
          temp_init_call2 <- deparse(temp_init_calli)
        } else {
          temp_init_call2 <- temp_init_calli
        }
        temp_init_call_c <- c(temp_init_call_c, temp_init_call2)
      }
      temp_init_call_c <- gsub("[[:space:]]", "", temp_init_call_c)
      temp_init_call_c <- paste0("list(", 
                                paste(temp_init_call_c, collapse = ",") , ")")
      temp_init_call_c <- str2lang(temp_init_call_c)
    } else { # if(length(temp_init_call) != 0) {
      temp_init_call_c <- mcall$init
    }
  } # if(is.language(temp_init_call_in)) {
  
  mcall$init <- temp_init_call_c
  
  
  
  
  xs <- ids <- dfs <- NA
  
  
  for (i in names(mcall)[-1]) {
    # don't let family argument also to be evaluated
    no_default_args_plus_family <- c(no_default_args, "family")
    if (!i %in% no_default_args_plus_family) {
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
      if(length(checks.) == 0) err. <- TRUE # This one line added for update
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
    } else if (i %in% no_default_args_plus_family) {
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
  
  
  
  # e_formulasi <- NULL
  # 
  # e_formula_grsi <- NULL
  # e_formula_gr_strsi <- NULL
  # e_prior_betasi <- NULL
  # e_init_betasi <- NULL
  
  # argumentsx <<- arguments
  # 
  # for (agsxi in letters[1:26]) {
  #   if(is.null(paste0(agsxi, "_", "formula"))) {
  #     arguments[[paste0(agsxi, "_", "formula_gr")]] <- NULL
  #     arguments[[paste0(agsxi, "_", "formula_gr_str")]] <- NULL
  #     arguments[[paste0(agsxi, "_", "prior_beta")]] <- NULL
  #     arguments[[paste0(agsxi, "_", "init_beta")]] <- NULL
  #   }
  # }
  
  
  
  if(is.character(arguments$select_model)) {
    select_model <- arguments$select_model
  } else if(is.symbol(arguments$select_model)) {
    select_model <- deparse(arguments$select_model)
  } else if(!is.character(arguments$select_model) |
            !is.symbol(arguments$select_model)
            ) {
    stop("The 'select_model' must be a symbol or single character string")
  }
  
  
  
  
  # Better below, control match_sitar_d_form from select_model arg
  if(select_model == 'sitar') {
    # 'd' formula control
    # 1) Match with 'sitar' package (i.e., exclude 'd' from the fixed effects)
    # 2) Or, include 'd' in fixed and random effects
    # Setting default to FALSE to match scenario 2 for now
    # TODO
    # match_sitar_d_form <- FALSE
  }
  
  if(select_model == 'sitar3') select_model <- 'sitar'
  
  if(grepl('sitar4', select_model)) {
    if(select_model == 'sitar4fr') match_sitar_d_form <- FALSE
    if(select_model == 'sitar4r')  match_sitar_d_form <- TRUE
    if(select_model == 'sitar4')   match_sitar_d_form <- FALSE # default
    sitar_nparms <- 4
    select_model <- 'sitar'
  } else if(select_model == 'sitar') {
    sitar_nparms <- 3
    select_model <- 'sitar'
    match_sitar_d_form <- FALSE
  } else {
    match_sitar_d_form <- FALSE
  }
  
  
  
  
  allowed_model_names   <- c('sitar', 'sitar3', 'sitar4', 'sitar4fr', 'sitar4r',
                             'pb1', 'pb2', 'pb3')
  allowed_model_names_  <- paste(allowed_model_names, collapse = ", " )
  allowed_model_names__ <- paste0("(", allowed_model_names_, ")")

  if(!select_model %in% allowed_model_names) {
    stop("Currently supported models (via 'select_model' argument) are:",
         "\n ",
         " ", paste(paste0("'", allowed_model_names, "'"), collapse = ", ")
         )
  }
  
  
  
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
  
  remove_spaces <- c('a_formula_gr_str', 'b_formula_gr_str', 
                     'c_formula_gr_str', 'd_formula_gr_str',
                     'e_formula_gr_str', 'f_formula_gr_str', 
                     'sigma_formula_gr_str')
  
  for (ip in remove_spaces) {
    arguments[[ip]] <-  gsub_space(arguments[[ip]] )
  }
  
  
  
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
  iter <-  arguments$iter
  warmup <-  arguments$warmup <- eval(arguments$warmup)
 
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
        gsub("\\s", "", 
             paste(deparse(substitute(sigma_group_arg)), collapse = ""))
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
            sigma_group_arg[[sigma_group_argi]] <- 
              gsub("'", "", sigma_group_arg[[sigma_group_argi]])
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
      gsub("\\s", "", paste(deparse(substitute(sigma_group_arg)), 
                            collapse = ""))
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
    "get_formula",
    "get_stanvars",
    "get_priors",
    "get_set_priors",
    "validate_priors",
    "get_set_init",
    "set_self_priors",
    "set_replace_priors",
    "set_same_priors_hierarchy",
    "outliers",
    "select_model",
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
  
  if (!is.na(univariate_by$by) & univariate_by$verbose) {
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
  
  
  
  
  # First assign NULL to avoid global vars issue in Package
  set_env <- environment()
  for (agsxi in letters[1:26]) {
    assign(paste0(agsxi, "", "" , "") , NULL, 
           envir = set_env)
      assign(paste0(agsxi, "_", "formula" , "si") , NULL, 
             envir = set_env)
      assign(paste0(agsxi, "_", "formula_gr" , "si") , NULL,
             envir = set_env)
      assign(paste0(agsxi, "_", "formula_gr_str" , "si") , NULL,
             envir = set_env)
      assign(paste0(agsxi, "_", "prior_beta" , "si") , NULL,
             envir = set_env)
      assign(paste0(agsxi, "_", "cov_prior_beta" , "si") , NULL, 
             envir = set_env)
      assign(paste0(agsxi, "_", "prior_sd" , "si") , NULL,
             envir = set_env)
      assign(paste0(agsxi, "_", "cov_prior_sd" , "si") , NULL,
             envir = set_env)
      assign(paste0(agsxi, "_", "init_beta", "si" ) , NULL,
             envir = set_env)
      assign(paste0(agsxi, "_", "cov_init_beta" , "si") , NULL,
             envir = set_env)
      assign(paste0(agsxi, "_", "init_sd", "si" ) , NULL,
             envir = set_env)
      assign(paste0(agsxi, "_", "cov_init_sd" , "si") , NULL, 
             envir = set_env)
  }
  
  
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
  
  funlist_r <- funlist_rnamelist <- funlist_rvaluelist <- list()
  # funlist_r <- list()
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
    
    

    # Set to NULL for those not included
    for (agsxi in letters[1:26]) {
      if(is.null(arguments[[paste0(agsxi, "_", "formula" , "")]])) {
        assign(paste0(agsxi, "_", "formula" , "si") , NULL)
        assign(paste0(agsxi, "_", "formula_gr" , "si") , NULL)
        assign(paste0(agsxi, "_", "formula_gr_str" , "si") , NULL)
        assign(paste0(agsxi, "_", "prior_beta" , "si") , NULL)
        assign(paste0(agsxi, "_", "cov_prior_beta" , "si") , NULL)
        assign(paste0(agsxi, "_", "prior_sd" , "si") , NULL)
        assign(paste0(agsxi, "_", "cov_prior_sd" , "si") , NULL)
        assign(paste0(agsxi, "_", "init_beta", "si" ) , NULL)
        assign(paste0(agsxi, "_", "cov_init_beta" , "si") , NULL)
        assign(paste0(agsxi, "_", "init_sd", "si" ) , NULL)
        assign(paste0(agsxi, "_", "cov_init_sd" , "si") , NULL)
      }
    }
    

    #################
    validate_fixed_random_parms <- function(fixedsi, randomsi, 
                                            allowed_parm_letters, 
                                            select_model) {
      
      parm_letters_fixed <- strsplit(gsub("\\+", " ", fixedsi), " ")[[1]]
      parm_letters_fixed <- sort(parm_letters_fixed)
      parm_letters_fixed <- parm_letters_fixed[1:length(allowed_parm_letters)]
      parm_letters_fixed <- parm_letters_fixed[!is.na(parm_letters_fixed)]
      
      parm_letters_random <- strsplit(gsub("\\+", " ", randomsi), " ")[[1]]
      parm_letters_random <- sort(parm_letters_random)
      parm_letters_random <- parm_letters_random[1:length(allowed_parm_letters)]

      if(select_model == 'pb1' | 
         select_model == 'pb2' | 
         select_model == 'pb3') {
        if(length(parm_letters_fixed) != length(allowed_parm_letters))
          stop("For model '", select_model, "'", ", 
               the number of parameters must be ",
               length(allowed_parm_letters),
               " \n ", 
               "(parameters ", 
               paste(paste0("'", allowed_parm_letters, "'"), collapse = " "),
               ")"
          )
      }
      
      if(select_model == 'sitar') {
        if(length(parm_letters_fixed) > length(allowed_parm_letters))
          stop("For model '", select_model, "'", ", 
               the maximum number of parameters is ",
               length(allowed_parm_letters),
               " \n ", 
               "(parameters ", 
               paste(paste0("'", allowed_parm_letters, "'"), collapse = " "),
               ")"
          )
      }
      
      get_parm_letters <- parm_letters_fixed
      sub_parm_letters_fixed <- intersect(allowed_parm_letters, 
                                          get_parm_letters)
      sub_parm_letters <- sub_parm_letters_fixed
      inv_parm_letters <- get_parm_letters[!get_parm_letters %in% 
                                             sub_parm_letters]
      inv_parm_letters <- sort(inv_parm_letters)
      if(length(inv_parm_letters) > 0) {
        see_what_formual <- paste0(" Please see and correct ", "'", 
                                   'fixed', "'", " argument")
        not_allowed_parsm <- paste(paste0("'", inv_parm_letters, "'"), 
                                   collapse = " ")
        msg_1 <- paste0("Parameter ", not_allowed_parsm, 
                        " not allowed for ", "'", select_model, "'", " model" )
        msg_2 <- paste0(" Allowed parameters are ", 
                        paste(paste0("'", 
                                     allowed_parm_letters, "'"), 
                              collapse = " "))
        stop(msg_1, "\n ", msg_2, " \n ", see_what_formual)
      }
      
      get_parm_letters <- parm_letters_random
      sub_parm_letters_random <- intersect(allowed_parm_letters, 
                                           get_parm_letters)
      sub_parm_letters <- sub_parm_letters_random
      inv_parm_letters <- get_parm_letters[!get_parm_letters %in% 
                                             sub_parm_letters]
      inv_parm_letters <- sort(inv_parm_letters)
      if(length(inv_parm_letters) > 0) {
        see_what_formual <- paste0(" Please see and correct ", "'", 
                                   'random', "'", " argument")
        not_allowed_parsm <- paste(paste0("'", inv_parm_letters, "'"), 
                                   collapse = " ")
        msg_1 <- paste0("Parameter ", 
                        not_allowed_parsm, " not allowed for ", "'", 
                        select_model, "'", " model" )
        msg_2 <- paste0(" Allowed parameters are ", 
                        paste(paste0("'", 
                                     allowed_parm_letters, "'"), 
                              collapse = " "))
        stop(msg_1, "\n ", msg_2, " \n ", see_what_formual)
      }
      
      # if no error, then check if fixed parm is present for each random
      
      sub_parm_letters_fixed_random <- intersect(parm_letters_fixed, 
                                                 parm_letters_random)
      inv_parm_letters_fixed_random <- 
        parm_letters_random[!parm_letters_random %in% parm_letters_fixed]
      inv_parm_letters_fixed_random <- sort(inv_parm_letters_fixed_random)
      if(length(inv_parm_letters_fixed_random) > 0) {
        not_allowed_parsm <- paste(paste0("'", 
                                          inv_parm_letters_fixed_random, "'"), 
                                   collapse = " ")
        stop(
          "Parameter ", not_allowed_parsm , " included in the random part of ",
          "\n ",
          " the model but missing from the fixed effects.",
          "\n ",
          " Please check and correct 'fixed' / 'random' arguments"
        )
      }
      
      # if all checks pass, then reassign fixed and random 
      sub_parm_letters_fixed <- sort(sub_parm_letters_fixed)
      sub_parm_letters_random <- sort(sub_parm_letters_random)
      
      out_fixed <- paste(sub_parm_letters_fixed, collapse = "+")
      out_random <- paste(sub_parm_letters_random, collapse = "+")
      list(fixed = out_fixed, random = out_random)
    } # validate_fixed_random_parms
    
    
    
    
    
    
    # Model (select_model) specifc fixedsi and randomsi
    
    allowed_parm_letters <- NULL
    if(select_model == 'sitar') allowed_parm_letters <- letters[1:sitar_nparms]
    if(select_model == 'pb1')   allowed_parm_letters <- letters[1:5]
    if(select_model == 'pb2')   allowed_parm_letters <- letters[1:6]
    if(select_model == 'pb3')   allowed_parm_letters <- letters[1:6]
    
    # fixedsi <- "a+c"
    # randomsi <- "a+b+d+e"
    
    fixedsi_randomsi <- validate_fixed_random_parms(fixedsi, randomsi,
                                                    allowed_parm_letters, 
                                                    select_model)
    fixedsi <- fixedsi_randomsi[['fixed']]
    randomsi <- fixedsi_randomsi[['random']]
    
    
    # Covariate not allowed when matching to sitar 'd' form
    
    if(select_model == 'sitar') {
      
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
      
      
      if (match_sitar_d_form) {
        if ((grepl("d", fixedsi, fixed = T) |
             grepl("d", randomsi, fixed = T)) &
            (!grepl("^~1$", d_formulasi) |
             !grepl("^~1$", d_formula_grsi))) {
          stop(
            "Parameter 'd' is missing in the fixed effects part of the model ",
            "\n ",
            " but specified in the random effects part of the model ",
            "\n ",
            " (This is to match with the 'sitar' package's formulation)",
            "\n ",
            " For this formulation (i.e., 'd' is missing in the fixed effects)",
            "\n ",
            " covariate(s) are not allowed"
          )
        }
      }
    } # if(select_model == 'sitar') {
    
    
    
    
    # Add missing parameters to the dpar_formula
    
    if (!is.null(dpar_formulasi)) {
      if (grepl("^1$", dpar_formulasi)) {
        dpar_formulasi <- paste0("lf(", "sigma", "~", dpar_formulasi, ")")
      } else if (grepl("^~1", dpar_formulasi)) { 
        dpar_formulasi <- paste0("lf(", "sigma", dpar_formulasi, ")")
      } else if (grepl("^sigma~1", dpar_formulasi)) { 
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
    

    
    # Check for higher level model and update level 2 random formula
    
    f_checks_gr_gr_str <- function(a, b) {
      if(!is.null(a)) {
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
                 "default priors placed by brms for those varinace covarinace", 
                 "\n ",
                 "or else use get_prios to place priors manually and the pass ",
                 "\n ",
                 "to the bsitar by using argument 'set_self_priors'"
            )
          }
        } else if(!is.null(b[[1]])) {
          b_out <- b
        }
        out <- b_out
      } # if(!is.null(a)) {
      if(is.null(a)) {
        out <- NULL
      }
      out
    } # f_checks_gr_gr_str
    
    
    # First, if a,b,c,d or e not NULL but sigma_formula_grsi NULL
    # Then set to ~1 because then only first part of the 
    # (i.e., before first + ) will be copied to the 
    # _grsi
    
    test_gr_sr_str_function <- function(x_grsi, x_gr_strsi) {
      if(!is.null(x_grsi)) {
        if(x_gr_strsi != 'NULL') {
          if(x_grsi == 'NULL') {
            x_grsi <- "~1"
          } else if(x_grsi != "~1") {
            if(verbose) {
              message("Argument '", 
                      substitute(x_grsi), "' changed from '", 
                      x_grsi , "' to  '~1'.")
              message("Instead of '", 
                      substitute(x_grsi), " = ", x_grsi, "', 
                    the covariates are now specified as '", 
                      substitute(x_gr_strsi), " = ", x_grsi, "'")
            }
            x_grsi <- "~1"
            
          } else {
            x_grsi <- x_grsi
          }
        } 
        out <- x_grsi
      } # if(!is.null(x_grsi)) {
      if(is.null(x_grsi)) {
        out <- NULL
      }
      out
    } # end test_gr_sr_str_function
    
    
    a_formula_grsi <- 
      test_gr_sr_str_function(a_formula_grsi, a_formula_gr_strsi)
    b_formula_grsi <- 
      test_gr_sr_str_function(b_formula_grsi, b_formula_gr_strsi)
    c_formula_grsi <- 
      test_gr_sr_str_function(c_formula_grsi, c_formula_gr_strsi)
    d_formula_grsi <- 
      test_gr_sr_str_function(d_formula_grsi, d_formula_gr_strsi)
    e_formula_grsi <- 
      test_gr_sr_str_function(e_formula_grsi, e_formula_gr_strsi)
    f_formula_grsi <- 
      test_gr_sr_str_function(f_formula_grsi, f_formula_gr_strsi)
    
    
    a_fcgs_out <- f_checks_gr_gr_str(a_formula_grsi, a_formula_gr_strsi)
    b_fcgs_out <- f_checks_gr_gr_str(b_formula_grsi, b_formula_gr_strsi)
    c_fcgs_out <- f_checks_gr_gr_str(c_formula_grsi, c_formula_gr_strsi)
    d_fcgs_out <- f_checks_gr_gr_str(d_formula_grsi, d_formula_gr_strsi)
    e_fcgs_out <- f_checks_gr_gr_str(e_formula_grsi, e_formula_gr_strsi)
    f_fcgs_out <- f_checks_gr_gr_str(f_formula_grsi, f_formula_gr_strsi)
    
    
    
    # First, if sigma_formula_gr_strsi not NULL but sigma_formula_grsi NULL
    # Then set sigma_formula_grsi to ~1 because then only first part of the 
    # sigma_formula_gr_strsi (i.e., before first + ) will be copied to the 
    # sigma_formula_grsi
    
    
    # when no a, b, c, d, or e random effect, then sigma_formula_gr or 
    # sigma_formula_gr_str are not allowed
   
    sigma_formula_grsi_NULL <- sigma_formula_gr_strsi_NULL <- FALSE
    if (is.null(sigma_formula_grsi[[1]][1]) |
        sigma_formula_grsi == "NULL") {
      sigma_formula_grsi_NULL <- TRUE
    }
    if (is.null(sigma_formula_gr_strsi[[1]][1]) |
        sigma_formula_gr_strsi == "NULL") {
      sigma_formula_gr_strsi_NULL <- TRUE
    }
    
    if(randomsi == "") {
      if (!sigma_formula_grsi_NULL |
          !sigma_formula_gr_strsi_NULL) {
        stop("Random effect for parameter 'sigma' are not allowed",
             " \n ",
             " if no group level random effect is specified i.e., random = ''",
             " \n ",
             " Therefore, 
             please set argument sigma_formula_gr/sigma_formula_gr_str to NULL")
      }
    }
    
  
    
    
    sigma_formula_grsi <- test_gr_sr_str_function(sigma_formula_grsi, 
                                                  sigma_formula_gr_strsi)
   
    
   
    
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
    

    
    sigma_fcgs_out <- f_checks_gr_gr_str(sigma_formula_grsi, 
                                         sigma_formula_gr_strsi)
    
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
          d_formula_grsi <- strsplit(d_formula_gr_strsi, 
                                     "+(", fixed = T)[[1]][1]
        }
      }
    }
    
    if(!is.null(e_fcgs_out)) {
      if(e_formula_grsi == "~1" & !is.null(e_formula_gr_strsi[[1]])) {
        e_formula_grsi <- strsplit(e_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    
    if(!is.null(f_fcgs_out)) {
      if(f_formula_grsi == "~1" & !is.null(f_formula_gr_strsi[[1]])) {
        f_formula_grsi <- strsplit(f_formula_gr_strsi, "+(", fixed = T)[[1]][1]
      }
    }
    
    
    
    if(!is.null(sigma_fcgs_out) & sigma_fcgs_out != 'NULL') {
      if(sigma_formula_grsi == "~1" & !is.null(sigma_formula_gr_strsi[[1]])) {
        sigma_formula_grsi <- strsplit(sigma_formula_gr_strsi, 
                                       "+(", fixed = T)[[1]][1]
      }
    }
    
    a_formula_grsi <- gsub("[()]", "", a_formula_grsi)
    b_formula_grsi <- gsub("[()]", "", b_formula_grsi)
    c_formula_grsi <- gsub("[()]", "", c_formula_grsi)
    if(!is.null(d_formula_grsi)) d_formula_grsi <- gsub("[()]", "", 
                                                        d_formula_grsi)
    e_formula_grsi <- gsub("[()]", "", e_formula_grsi)
    f_formula_grsi <- gsub("[()]", "", f_formula_grsi)
    
   
    
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
        "e_formulasi",
        "f_formulasi",
        "s_formulasi",
        "a_formula_grsi",
        "b_formula_grsi",
        "c_formula_grsi",
        "d_formula_grsi",
        "e_formula_grsi",
        "f_formula_grsi",
        "sigma_formulasi",
        "sigma_formula_grsi"
      )
    
    for (check_formualsi in check_formuals) {
      if(!is.null(ept(check_formualsi)) & length(ept(check_formualsi)) !=0 ) {
        if (!grepl("~1", ept(check_formualsi)) &
            !grepl("~0", ept(check_formualsi))) {
          check_formualsi_with1 <-
            gsub("^~", "~1+", ept(check_formualsi), fixed = F)
          # new added on 25 8 2023 to a_formula_grsi etc. with _str format
          # BUT dont let it for sigma_formulasi and sigma_formula_grsi
          if(!grepl("^~", ept(check_formualsi))) {
            if(!grepl("^sigma", check_formualsi))
              check_formualsi_with1 <- paste0("~", check_formualsi_with1)
          }
          assign(check_formualsi, check_formualsi_with1)
        }
      } # if(!is.null(ept(check_formualsi))) {
      if(is.null(ept(check_formualsi)) | length(ept(check_formualsi)) ==0 ) {
        assign(check_formualsi, NULL)
      }
    } # for (check_formualsi in check_formuals) {
    
    
    

    
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
      autocor_formulasi <- gsub("\"", "", autocor_formulasi)
      if(!grepl("^~", autocor_formulasi)) {
        stop('autocor_formula argument should be a formula. E.g.,',
             "\n ",
             " autocor_formula = ~arms(p=1,q=1)",
             "\n ", 
             " It seems you forgot to add '~' before the autocor structure")
      }
      autocor_formi <- autocor_formulasi
    } # if (is.null(autocor_formulasi[[1]][1]) |
    

    if(!is.null(autocor_formi)) {
      tempunstx <- autocor_formi # '~unstr(time=visit, patient)'
      tempunstx <- gsub("[[:space:]]", "", tempunstx)
      if(grepl("unstr(", tempunstx, fixed = T)) {
        tempunstx_1 <- regmatches(tempunstx, gregexpr("(?<=\\().*?(?=\\))", 
                                                      tempunstx, perl=T))[[1]]
        tempunstx_2 <- strsplit(tempunstx_1, ",")[[1]][1]
        if(grepl("time=", tempunstx_2, fixed = T)) {
          tempunstx_3 <- sub(".*time=", "", tempunstx_2) 
        } else if(!grepl("time=", tempunstx_2, fixed = T)) {
          tempunstx_3 <- tempunstx_2
        }
        cortimeNlags_var <- tempunstx_3
      } # if(grepl("unstr(", tempunstx, fixed = T)) {
      
      if(!grepl("unstr(", tempunstx, fixed = T)) {
        cortimeNlags_var <- NULL
      }
    } # if(!is.null(autocor_formi)) {
      
    
    if(is.null(autocor_formi)) {
      cortimeNlags_var <- NULL
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
    
    
    if(!is.null(cortimeNlags_var)) {
      if(!is.factor(datai[[cortimeNlags_var]])) {
        datai[[cortimeNlags_var]] <- as.factor(datai[[cortimeNlags_var]])
        datai[[cortimeNlags_var]] <- droplevels(datai[[cortimeNlags_var]])
      } else {
        datai[[cortimeNlags_var]] <-  datai[[cortimeNlags_var]]
      }
      cortimeNlags <- nlevels(unique(datai[[tempunstx_3]]))
    } else if(is.null(cortimeNlags_var)) {
      cortimeNlags <- NULL
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
    
    if(select_model == "sitar") {
      if (match_sitar_d_form) {
        if (length(knots) > 2) {
          itemp <- strsplit(gsub("\\+", " ", fixedsi), " ")[[1]]
          itemp <- itemp[!grepl("d", itemp)]
          fixedsi <- paste(itemp, collapse = "+")
        }
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
        } else if (eval_arg == "min") {
          eval_arg.o <- min(data[[x]])
        } else if (eval_arg == "max") {
          eval_arg.o <- max(data[[x]])
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
    
    # bstartx <<- bstart
    # dataix <<- datai
    
    mat_s <- make_spline_matrix(datai[[xsi]], knots)
    
    
    # SplineFun_name <- "SplineFun" 
    SplineFun_name <- "DefFun" 
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
        'select_model',
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
    
  
    
    
    get_s_r_funs <- 
      prepare_function(
      x = xsi,
      y = ysi,
      id = idsi,
      knots = knots,
      nknots = nknots,
      data = datai,
      internal_function_args = internal_function_args
    )
    
    funlist[ii] <- get_s_r_funs[['rcsfun']]
    funlist_r[[ii]] <- get_s_r_funs[['r_funs']]
    
    
    #funlist_rx <<- funlist_r
    # stop()

    
    
    internal_formula_args_names <-
      c(
        "a_formulasi",
        "b_formulasi",
        "c_formulasi",
        "d_formulasi",
        "e_formulasi",
        "f_formulasi",
        "s_formulasi",
        "a_formula_grsi",
        "b_formula_grsi",
        "c_formula_grsi",
        "d_formula_grsi",
        "e_formula_grsi",
        "f_formula_grsi",
        
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
        "e_formula_gr_strsi",
        "f_formula_gr_strsi",
        "sigma_formula_gr_strsi",
        "set_higher_levels",
        "sigma_set_higher_levels",
        "select_model",
        "verbose",
        "unusedsi"
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
    
    # list_outx <<- list_out
    
    
    attributes(formula_bf) <- NULL
    
    temp_stancode <- temp_standata <- NULL
    
    
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
    
    
    sigma_group_arg$groupvar <- sigma_arg_groupvar
    
    
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
    ymean   <- mean(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
    ymedian <- median(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
    if(select_model == 'sitar') {
      ymax  <- max(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
      ymin  <- min(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
      ymaxs <- NULL
    } else if(select_model != 'sitar') {
      ymax <- round(max(predict(loess_fitx)), 2)
      ymaxs <- round(ymax * 0.95, 2)
    }
    
    ysd     <- sd(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
    ymad    <- mad(datai[[ysi]], na.rm = TRUE) %>% round(., 2)
    xsd     <- sd(datai[[xsi]], na.rm = TRUE) %>% round(., 2)
    
    
    
    loess_fit <- paste0("loess(", ysi, "~", xsi, ",", 'datai', ")")
    loess_fitx <- eval(parse(text = loess_fit))
    
    
    
    
    
    
    
    if (!is.null(pgvsi[[1]][1]) & pgvsi != "NULL") {
      setpgv <- eval(parse(text = pgvsi))
      cstart <- log(setpgv) / 5.0
      dstart <- log(setpgv)
    } else {
      cstart <- 0.01
      dstart <- 0.01
    }
    
    if (!is.null(apgvsi[[1]][1]) & apgvsi != "NULL") {
      setapgv <- eval(parse(text = apgvsi))
      estart <- setapgv
    } else {
      estart <- 13
    }
    
    
    cstart <- round(cstart, 2)
    dstart <- round(dstart, 2)
    estart <- round(estart, 1)
    
    
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
        "ymax",
        "ymaxs",
        "ymin",
        "ysd",
        "ymad",
        "xsd",
        "lm_a_cov_sd",
        "lm_b_cov_sd",
        "lm_c_cov_sd",
        "bstart",
        "cstart",
        "dstart",
        "estart"
      )
    
    
    prior_args_internal_names <-
      c(
        lm_val_list_not,
        cov_list_names,
        "a_formulasi",
        "b_formulasi",
        "c_formulasi",
        "d_formulasi",
        "e_formulasi",
        "f_formulasi",
        "s_formulasi",
        "a_formula_grsi",
        "fixedsi",
        "b_formula_grsi",
        "c_formula_grsi",
        "d_formula_grsi",
        "e_formula_grsi",
        "f_formula_grsi",
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
        'cortimeNlags_var',
        'cortimeNlags',
        "select_model",
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
        e_init_beta = e_init_betasi,
        f_init_beta = f_init_betasi,
        s_init_beta = s_init_betasi,
        a_cov_init_beta = a_cov_init_betasi,
        b_cov_init_beta = b_cov_init_betasi,
        c_cov_init_beta = c_cov_init_betasi,
        d_cov_init_beta = d_cov_init_betasi,
        e_cov_init_beta = e_cov_init_betasi,
        f_cov_init_beta = f_cov_init_betasi,
        s_cov_init_beta = s_cov_init_betasi,
        a_init_sd = a_init_sdsi,
        b_init_sd = b_init_sdsi,
        c_init_sd = c_init_sdsi,
        d_init_sd = d_init_sdsi,
        e_init_sd = e_init_sdsi,
        f_init_sd = f_init_sdsi,
        a_cov_init_sd = a_cov_init_sdsi,
        b_cov_init_sd = b_cov_init_sdsi,
        c_cov_init_sd = c_cov_init_sdsi,
        d_cov_init_sd = d_cov_init_sdsi,
        e_cov_init_sd = e_cov_init_sdsi,
        f_cov_init_sd = f_cov_init_sdsi,
        
        sigma_init_beta = sigma_init_betasi,
        sigma_cov_init_beta = sigma_cov_init_betasi,
        sigma_init_sd = sigma_init_sdsi,
        sigma_cov_init_sd = sigma_cov_init_sdsi,
        
        rsd_init_sigma = rsd_init_sigmasi,
        dpar_init_sigma = dpar_init_sigmasi,
        dpar_cov_init_sigma = dpar_cov_init_sigmasi,
        autocor_init_acor = autocor_init_acorsi,
        autocor_init_unstr_acor = autocor_init_unstr_acorsi,
        gr_init_cor = gr_init_corsi,
        sigma_init_cor = sigma_init_corsi,
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
    
    vcov_init_0e <- eval(parse(text =  "vcov_init_0si" ))
    vcov_init_0e <- eval(parse(text =  vcov_init_0e ))
    
  

    
    
    #######################################
    
    set_priors_initials_agrs <- list()
    

    set_priors_initials_agrs $ a_prior_beta <- a_prior_betasi
    set_priors_initials_agrs $ b_prior_beta <- b_prior_betasi
    set_priors_initials_agrs $ c_prior_beta <- c_prior_betasi
    set_priors_initials_agrs $ d_prior_beta <- d_prior_betasi
    set_priors_initials_agrs $ e_prior_beta <- e_prior_betasi
    set_priors_initials_agrs $ f_prior_beta <- f_prior_betasi
    
    set_priors_initials_agrs $ s_prior_beta <- s_prior_betasi
    
    set_priors_initials_agrs $ a_cov_prior_beta <- a_cov_prior_betasi
    set_priors_initials_agrs $ b_cov_prior_beta <- b_cov_prior_betasi
    set_priors_initials_agrs $ c_cov_prior_beta <- c_cov_prior_betasi
    set_priors_initials_agrs $ d_cov_prior_beta <- d_cov_prior_betasi
    set_priors_initials_agrs $ e_cov_prior_beta <- e_cov_prior_betasi
    set_priors_initials_agrs $ f_cov_prior_beta <- f_cov_prior_betasi
    
    set_priors_initials_agrs $ s_cov_prior_beta <- s_cov_prior_betasi
    
    set_priors_initials_agrs $ a_prior_sd <- a_prior_sdsi
    set_priors_initials_agrs $ b_prior_sd <- b_prior_sdsi
    set_priors_initials_agrs $ c_prior_sd <- c_prior_sdsi
    set_priors_initials_agrs $ d_prior_sd <- d_prior_sdsi
    set_priors_initials_agrs $ e_prior_sd <- e_prior_sdsi
    set_priors_initials_agrs $ f_prior_sd <- f_prior_sdsi
    
    set_priors_initials_agrs $ a_cov_prior_sd <- a_cov_prior_sdsi
    set_priors_initials_agrs $ b_cov_prior_sd <- b_cov_prior_sdsi
    set_priors_initials_agrs $ c_cov_prior_sd <- c_cov_prior_sdsi
    set_priors_initials_agrs $ d_cov_prior_sd <- d_cov_prior_sdsi
    set_priors_initials_agrs $ e_cov_prior_sd <- e_cov_prior_sdsi
    set_priors_initials_agrs $ f_cov_prior_sd <- f_cov_prior_sdsi
    
    set_priors_initials_agrs $ gr_prior_cor         <- gr_prior_corsi
    set_priors_initials_agrs $ sigma_prior_cor      <- sigma_prior_corsi
    set_priors_initials_agrs $ sigma_prior_beta     <- sigma_prior_betasi
    set_priors_initials_agrs $ sigma_cov_prior_beta <- sigma_cov_prior_betasi
    
    set_priors_initials_agrs $ sigma_prior_sd      <- sigma_prior_sdsi
    set_priors_initials_agrs $ sigma_cov_prior_sd  <- sigma_cov_prior_sdsi
    set_priors_initials_agrs $ rsd_prior_sigma     <- rsd_prior_sigmasi
    set_priors_initials_agrs $ dpar_prior_sigma    <- dpar_prior_sigmasi
    
  
    set_priors_initials_agrs $ dpar_cov_prior_sigma     <- 
      dpar_cov_prior_sigmasi
    set_priors_initials_agrs $ autocor_prior_acor       <- autocor_prior_acorsi
    set_priors_initials_agrs $ autocor_prior_unstr_acor <- 
      autocor_prior_unstr_acorsi
    set_priors_initials_agrs $ mvr_prior_rescor         <- mvr_prior_rescorsi
    set_priors_initials_agrs $ prior_data               <- prior_data
    set_priors_initials_agrs $ prior_data_internal      <- prior_data_internal
    set_priors_initials_agrs $ prior_args_internal      <- prior_args_internal
    set_priors_initials_agrs $ init_arguments           <- init_arguments
    set_priors_initials_agrs $ init_data                <- init_data
    set_priors_initials_agrs $ init_data_internal       <- init_data_internal
    set_priors_initials_agrs $ init_args_internal       <- init_args_internal
    set_priors_initials_agrs $ temp_stancode            <- temp_stancode
    set_priors_initials_agrs $ temp_standata            <- temp_standata

    set_priors_initials_agrs $ custom_order_prior_str   <- ""
    
    
    bpriors <- do.call(set_priors_initials, set_priors_initials_agrs)
    
    stanvar_priors <- attr(bpriors, "stanvars")
    
    initials <- attr(bpriors, "initials")
    

    
    
    
    
    #############################################################
    
    # check and add hierarchical prior (for 3 level and more)
    
    #############################################################    
    
    ##########################
    # First, sd
    
    set_class_what <- 'sd'
    
    set_org_priors_initials_agrs_what <- set_priors_initials_agrs
    
    set_randomsi_higher_levsl <- strsplit(gsub("\\+", " ", randomsi), " ")[[1]]
    
    check_sigma_str <- 
      eval(parse(text = paste0('sigma', "covcoefnames_gr_str")))
    
    if(!is.null(check_sigma_str)) {
      set_randomsi_higher_levsl <- c(set_randomsi_higher_levsl, 'sigma')
    }
    
    evaluate_higher_level_sd_priors <- function(set_nlpar_, 
                                                set_class,
                                                set_prior,
                                                set_cov_prior,
                                                org_priors_initials_agrs,
                                                set_env,
                                                ...) {
      
      custom_order_prior_str <- c(paste0(set_nlpar_, "_prior_sd"),
                                  paste0(set_nlpar_, "_cov_prior_sd"))
     
      
      eval_what <- eval(parse(text = paste0(set_nlpar_, 
                                            "covcoefnames_gr_str")), 
                        envir = set_env_what)
      if(!is.null(eval_what)) {
        gr_str_id   <- eval(parse(text = paste0(set_nlpar_, 
                                                "covcoefnames_gr_str_id")), 
                            envir = set_env_what)
        gr_str_coef <- eval(parse(text = paste0(set_nlpar_, 
                                                "covcoefnames_gr_str")), 
                            envir = set_env_what)
        gr_str_form <- eval(parse(text = paste0(set_nlpar_, 
                                                "covcoefnames_gr_str_form")), 
                            envir = set_env_what)
        gr_str_ncov <- eval(parse(text = paste0(set_nlpar_, 
                                                "ncov_gr_str")), 
                            envir = set_env_what)
        temp_gr_str_priors <- list()
        temp_gr_str_stanvars <- c()
        temp_gr_str_inits <- c()
        set_priors_initials_agrs_str <- org_priors_initials_agrs 
        # this for adding _prior_cor 
        
       counter_start_from_one_for_prior <- 0
        for (istrx in 2:length(eval_what)) {
          counter_start_from_one_for_prior <- 
            counter_start_from_one_for_prior + 1
          if(set_nlpar_ == 'sigma') {
            assign('sigma_arg_groupvar', gr_str_id[[istrx]], 
                   envir = set_env_what)
          } else {
            assign('group_arg_groupvar', gr_str_id[[istrx]], 
                   envir = set_env_what)
          }
          
          assign( paste0(set_nlpar_, "_formula_grsi"), 
                  gr_str_form[[istrx]], envir = set_env_what)
          assign( paste0(set_nlpar_, "covcoefnames_gr"), 
                  gr_str_coef[[istrx]], envir = set_env_what)
          assign( paste0(set_nlpar_, "ncov_gr"), 
                  gr_str_coef[[istrx]], envir = set_env_what)
          
          prior_args_internal_str <- list()
          prior_args_internal_str <- mget(prior_args_internal_names, 
                                          envir = set_env_what)
          set_priors_initials_agrs_str $ prior_args_internal <- 
            prior_args_internal_str

          set_priors_initials_agrs_str $ custom_order_prior_str <- 
            custom_order_prior_str
       
          set_priors_initials_agrs_str [[paste0(set_nlpar_, 
                                                "_prior_sd")]]  <- 
            set_prior[counter_start_from_one_for_prior]
          set_priors_initials_agrs_str [[paste0(set_nlpar_, 
                                                "_cov_prior_sd")]] <- 
            set_cov_prior[counter_start_from_one_for_prior]

          bpriors_str <- do.call(set_priors_initials, 
                                 set_priors_initials_agrs_str, 
                                 envir = set_env_what)

          stanvars_str <- attr(bpriors_str, "stanvars")
          initials_str <- attr(bpriors_str, "initials")
          temp_gr_str_stanvars <- c(temp_gr_str_stanvars, stanvars_str)
          temp_gr_str_priors[[istrx]] <- bpriors_str
        }
        temp_gr_str_priors <- temp_gr_str_priors %>% do.call(rbind, .)
        out <- list(temp_gr_str_priors = temp_gr_str_priors,
                    temp_gr_str_stanvars = temp_gr_str_stanvars,
                    temp_gr_str_inits = temp_gr_str_inits)
       # bpriors <-  rbind(bpriors, temp_gr_str_priors)
      }
      out
    } # evaluate_higher_level_sd_priors
    
    
    temp_gr_str_priors_sd <- list()
    temp_gr_str_stanvars_sd <-  temp_gr_str_inits_sd <- c()
    for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
      set_nlpar_what <- set_randomsi_higher_levsli
      set_env_what   <- environment()
      n_higher_str   <- length(eval(parse(text = paste0(set_nlpar_what,
                                                        "covcoefnames_gr_str")),
                                    envir = set_env_what))
      n_higher_str   <- n_higher_str - 1
      # Evaluate hierarchy structure only if levels are 3 or more 
      # This because for the second level, the first argument of _str already 
      # evaluated with _sd priors
      
      if(n_higher_str > 0) {
        # assign _prior
        set_assign_prior_what <- '_prior'
        check_prior_ifp <- 
          extract_prior_str_lv(ept(paste0(set_nlpar_what, 
                                          paste0(set_assign_prior_what, 
                                                 "_", set_class_what, 
                                                 "_strsi"))))
        check_prior_ifp_true_false <- FALSE
        if(length(check_prior_ifp) == 1) {
          if(!is.null(ept(paste0(set_nlpar_what, 
                                 paste0(set_assign_prior_what, "_", 
                                        set_class_what, "_strsi") ))[[1]][1])) {
            check_prior_ifp_true_false <- TRUE 
          } else {
            check_prior_ifp_true_false <- FALSE
          }
          if(ept(paste0(set_nlpar_what, 
                        paste0(set_assign_prior_what, "_", 
                               set_class_what, "_strsi") )) != "NULL") {
            check_prior_ifp_true_false <- TRUE
          } else {
            check_prior_ifp_true_false <- FALSE
          }
        } else if(length(check_prior_ifp) > 1) {
          check_prior_ifp_true_false <- TRUE
        }
        if(check_prior_ifp_true_false) {
          assign(paste0(set_nlpar_what, paste0(set_assign_prior_what, "_",
                                               set_class_what, "si")),  
                 check_prior_ifp)
        }
        set_prior_what <- ept(paste0(set_nlpar_what, 
                                     paste0(set_assign_prior_what, "_", 
                                            set_class_what, "si") ))
        
        
        paste_message <- paste("Length of prior elements for random effect ",
             "'", set_nlpar_what, "'",
             " \n",
             "  specified by using the argument ", 
             "'", paste0(set_nlpar_what, 
                         paste0(set_assign_prior_what, "_", 
                                set_class_what, "_str")), "'",  " ",
             " \n",
             "  should be one or same as the levels of hierarchy minus one.",
             " \n",
             "  (minus one because the prior for the second level of hierarchy",
             " \n",
             "  is taken from the ", 
             "'", paste0(set_nlpar_what, paste0(set_assign_prior_what, 
                                                "_", set_class_what, "")), "'"
        )
        
        if(length(set_prior_what) > 1 & 
           length(set_prior_what) != n_higher_str) {
          stop(paste_message)
        } else if(length(set_prior_what) == 1) {
          set_prior_what <- rep(set_prior_what, n_higher_str)
        }
        paste_message <- NULL
        # end assign _prior
        
        
        
        
        # assign _cov_prior
        set_assign_prior_what <- '_cov_prior'
        check_prior_ifp <- 
          extract_prior_str_lv(ept(paste0(set_nlpar_what, 
                                          paste0(set_assign_prior_what, 
                                                 "_", set_class_what, 
                                                 "_strsi"))))
        check_prior_ifp_true_false <- FALSE
        if(length(check_prior_ifp) == 1) {
          if(!is.null(ept(paste0(set_nlpar_what, 
                                 paste0(set_assign_prior_what, "_", 
                                        set_class_what, "_strsi") ))[[1]][1])) {
            check_prior_ifp_true_false <- TRUE 
          } else {
            check_prior_ifp_true_false <- FALSE
          }
          if(ept(paste0(set_nlpar_what, 
                        paste0(set_assign_prior_what, "_", 
                               set_class_what, "_strsi") )) != "NULL") {
            check_prior_ifp_true_false <- TRUE
          } else {
            check_prior_ifp_true_false <- FALSE
          }
        } else if(length(check_prior_ifp) > 1) {
          check_prior_ifp_true_false <- TRUE
        }
        if(check_prior_ifp_true_false) {
          assign(paste0(set_nlpar_what, paste0(set_assign_prior_what, "_",
                                               set_class_what, "si")),  
                 check_prior_ifp)
        }
        set_cov_prior_what <- ept(paste0(set_nlpar_what, 
                                         paste0(set_assign_prior_what, "_", 
                                                set_class_what, "si") ))
        
        if(length(set_cov_prior_what) > 1 & 
           length(set_cov_prior_what) != n_higher_str) {
          stop("Length of prior elements for random effect parameter ",
               "'", set_nlpar_what, "'",
               " \n",
               "  specified by using the argument ", 
               "'", paste0(set_nlpar_what, 
                           paste0(set_assign_prior_what, "_", 
                                  set_class_what, "_str")), "'",  " ",
               " \n",
               "  should be one or same as the levels of hierarchy minus one.",
               " \n",
               "  (minus one because prior for the first level of hierarchy",
               " \n",
               "  is taken from the ", 
               "'", paste0(set_nlpar_what, 
                           paste0(set_assign_prior_what, "_", 
                                  set_class_what, "")), "'"
          )
        } else if(length(set_prior_what) == 1) {
          set_cov_prior_what <- rep(set_cov_prior_what, n_higher_str)
        }
        # end assign _cov_prior
      
       
        out2 <- evaluate_higher_level_sd_priors(set_nlpar_ = set_nlpar_what, 
                                        set_class  = set_class_what,
                                        set_prior = set_prior_what,
                                        set_cov_prior = set_cov_prior_what,
                                        set_env = set_env_what,
                                        org_priors_initials_agrs = 
                                          set_org_priors_initials_agrs_what)
      
      temp_gr_str_priors_sd[[set_randomsi_higher_levsli]] <- 
        out2 $ temp_gr_str_priors
      temp_gr_str_stanvars_sd <- 
        c(temp_gr_str_stanvars_sd, out2 $ temp_gr_str_stanvars)
      temp_gr_str_inits_sd <- 
        c(temp_gr_str_inits_sd,    out2 $ temp_gr_str_inits)
      } # end if(n_higher_str > 0) {
    } # end for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
    
    
    
    
     # Add temp_gr_str_priors_sd to bpriors above 
    higher_level_priors <- temp_gr_str_priors_sd %>% do.call(rbind, .)
    bpriors             <- rbind(bpriors, higher_level_priors)
    
    # Add temp_gr_str_stanvars_sd to stanvar_priors to above 
    if(length(temp_gr_str_stanvars_sd) > 0) {
      stanvar_priors_c <- temp_gr_str_stanvars_sd_c <- c()
      for (i in 1:length(stanvar_priors)) {
        stanvar_priors_c <- c(stanvar_priors_c, stanvar_priors[i])
      }
      for (i in 1:length(temp_gr_str_stanvars_sd)) {
        temp_gr_str_stanvars_sd_c <- c(temp_gr_str_stanvars_sd_c, 
                                       temp_gr_str_stanvars_sd[i])
      }
      stanvar_priors <- c(stanvar_priors_c, temp_gr_str_stanvars_sd_c)
    } # if(length(temp_gr_str_stanvars_sd_c) > 0) {
    # Add temp_gr_str_inits_sd to initials to above 
    # But this is not good because it sets initials only for the outset level
    # Instead, when 3 or more levels, set init = '0' or  init = 'random' or else 
    # vcov_init_0 = TRUE
    
    if(length(temp_gr_str_inits_sd) > 0) {
      initials_c <- temp_gr_str_inits_sd_c <- c()
      for (i in 1:length(initials)) {
        initials_c <- c(initials_c, initials[i])
      }
      for (i in 1:length(temp_gr_str_inits_sd)) {
        temp_gr_str_inits_sd_c <- c(temp_gr_str_inits_sd_c, 
                                    temp_gr_str_inits_sd[i])
      }
      initials <- c(initials_c, temp_gr_str_inits_sd) 
    } # if(length(temp_gr_str_inits_sd) > 0) {
    
    
    
    
    
    
    
    
    ##########################
    # Now, cor priors    
    
    # Adding cor priors is tricky because of complex |x| formulations possible
    
    set_class_what <- 'cor'
    
    set_org_priors_initials_agrs_what <- set_priors_initials_agrs
    
    set_randomsi_higher_levsl <- 'gr'
    
    check_sigma_str <- eval(parse(text = paste0('sigma', 
                                                "covcoefnames_gr_str")))
    
    if(!is.null(check_sigma_str)) {
      set_randomsi_higher_levsl <- c(set_randomsi_higher_levsl, 'sigma')
    }
    
    
    evaluate_higher_level_corr_priors <- function(set_nlpar_, 
                                                set_class,
                                                set_prior,
                                                id_higher_str = id_higher_str,
                                                corr_higher_str_tf = 
                                                  corr_higher_str_tf,
                                                org_priors_initials_agrs,
                                                set_env,
                                                ...) {
      
      custom_order_prior_str <- c(paste0(set_nlpar_, "_prior_cor"))
      
        temp_gr_str_priors <- list()
        temp_gr_str_stanvars <- c()
        temp_gr_str_inits <- c()

        set_priors_initials_agrs_str <- org_priors_initials_agrs 
        
        gr_str_id <- id_higher_str
        counter_start_from_one_for_prior <- 0
        for (istrx in 2:length(gr_str_id)) {
          counter_start_from_one_for_prior <- counter_start_from_one_for_prior + 1
          get_corr_higher_str_tf <- corr_higher_str_tf[istrx]
          
          if(get_corr_higher_str_tf) {
            if(set_nlpar_ == 'sigma') {
              assign('sigma_arg_groupvar', gr_str_id[istrx], 
                     envir = set_env_what)
              set_priors_initials_agrs_str $ sigma_prior_cor <- 
                set_prior[counter_start_from_one_for_prior]
            } else {
              assign('group_arg_groupvar', gr_str_id[istrx], 
                     envir = set_env_what)
              set_priors_initials_agrs_str $ gr_prior_cor <- 
                set_prior[counter_start_from_one_for_prior]
            }
            prior_args_internal_str <- list()
            prior_args_internal_str <- mget(prior_args_internal_names, 
                                            envir = set_env_what)
            set_priors_initials_agrs_str $ prior_args_internal <- 
              prior_args_internal_str
            
            set_priors_initials_agrs_str $ custom_order_prior_str <-
              custom_order_prior_str
            
            bpriors_str <- do.call(set_priors_initials, 
                                   set_priors_initials_agrs_str, 
                                   envir = set_env_what)
            stanvars_str <- attr(bpriors_str, "stanvars")
            initials_str <- attr(bpriors_str, "initials")
            temp_gr_str_stanvars <- c(temp_gr_str_stanvars, stanvars_str)
            temp_gr_str_priors[[istrx]] <- bpriors_str
            
            # If prior is not evaluated, set it to NULL
            # This was encountered for cor for univariate_by and multivariate
            # See, set_priors_initials function, line 2086
            # currently brms does not allow setting separate ljk prior for
            # subset and multivariate
            # removing resp leads to duplicate priors, so need to set only once
            
            if(temp_gr_str_priors[[istrx]][,1] != "") {
              temp_gr_str_priors <- temp_gr_str_priors %>% do.call(rbind, .)
            } else if(temp_gr_str_priors[[istrx]][,1] == "") {
              temp_gr_str_priors <- temp_gr_str_stanvars <- NULL
              temp_gr_str_inits <- NULL
            }
            
          } # if(get_corr_higher_str_tf) {
          if(!get_corr_higher_str_tf) {
            temp_gr_str_priors <- temp_gr_str_stanvars <- NULL
              temp_gr_str_inits <- NULL
          }
          
        } # for (istrx in 2:length(gr_str_id)) {
       
        out <- list(temp_gr_str_priors = temp_gr_str_priors,
                    temp_gr_str_stanvars = temp_gr_str_stanvars,
                    temp_gr_str_inits = temp_gr_str_inits)
       
      out
    } # evaluate_higher_level_sd_priors
    

    temp_gr_str_priors_corr <- list()
    temp_gr_str_stanvars_corr <-  temp_gr_str_inits_corr <- c()
    for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
      set_nlpar_what <- set_randomsi_higher_levsli
      set_env_what   <- environment()

      id_higher_str  <- eval(parse(text = paste0(set_nlpar_what, 
                                                 "_str_unique_id")), 
                             envir = set_env_what)
      
      n_higher_str   <- length(id_higher_str)
      n_higher_str   <- n_higher_str - 1
      corr_higher_str_tf <- eval(parse(text = paste0(set_nlpar_what, 
                                                     "_str_corr_tf")),
                                 envir = set_env_what)

      # Evaluate hierarchy structure only if levels are 3 or more 
      # This because for the second level, the first argument of _str already 
      # evaluated with _sd priors
      
      if(n_higher_str > 0) {
        # assign _prior
        set_assign_prior_what <- '_prior'
        check_prior_ifp <- 
          extract_prior_str_lv(ept(paste0(set_nlpar_what, 
                                          paste0(set_assign_prior_what, 
                                                 "_", set_class_what, 
                                                 "_strsi"))))
        check_prior_ifp_true_false <- FALSE
        if(length(check_prior_ifp) == 1) {
          if(!is.null(ept(paste0(set_nlpar_what, 
                                 paste0(set_assign_prior_what, "_", 
                                        set_class_what, "_strsi") ))[[1]][1])) {
            check_prior_ifp_true_false <- TRUE 
          } else {
            check_prior_ifp_true_false <- FALSE
          }
          if(ept(paste0(set_nlpar_what, 
                        paste0(set_assign_prior_what, "_", 
                               set_class_what, "_strsi") )) != "NULL") {
            check_prior_ifp_true_false <- TRUE
          } else {
            check_prior_ifp_true_false <- FALSE
          }
        } else if(length(check_prior_ifp) > 1) {
          check_prior_ifp_true_false <- TRUE
        }
        if(check_prior_ifp_true_false) {
          assign(paste0(set_nlpar_what, paste0(set_assign_prior_what, "_",
                                               set_class_what, "si")),  
                 check_prior_ifp)
        }
        
        set_prior_cor_what <- ept(paste0(set_nlpar_what, 
                                     paste0(set_assign_prior_what, "_", 
                                            set_class_what, "si") ))
        
        
        if(length(set_prior_cor_what) > 1 & 
           length(set_prior_cor_what) != n_higher_str) {
          stop("Length of prior elements for random effect parameter ",
               "'", set_prior_cor_what, "'",
               " \n",
               "  specified by using the argument ", 
               "'", paste0(set_prior_cor_what, paste0(set_assign_prior_what, 
                                                      "_", 
                                                      set_class_what, "_str")), 
               "'",  " ",
               " \n",
               "  should be one or same as the levels of hierarchy minus one.",
               " \n",
               "  (minus one because prior for the second level of hierarchy",
               " \n",
               "  is taken from the ", 
               "'", paste0(set_prior_cor_what, 
                           paste0(set_assign_prior_what, "_", 
                                  set_class_what, "")), "'"
          )
        } else if(length(set_prior_cor_what) == 1) {
          set_prior_cor_what <- rep(set_prior_cor_what, n_higher_str)
        }
        # end assign _prior
        

        out2 <- evaluate_higher_level_corr_priors(set_nlpar_ = set_nlpar_what, 
                                          set_class  = set_class_what,
                                          set_prior = set_prior_cor_what,
                                          id_higher_str = id_higher_str,
                                          corr_higher_str_tf = corr_higher_str_tf,
                                          set_env = set_env_what,
                                          org_priors_initials_agrs = 
                                            set_org_priors_initials_agrs_what)
        
        temp_gr_str_priors_corr[[set_randomsi_higher_levsli]] <- 
          out2 $ temp_gr_str_priors
        temp_gr_str_stanvars_corr <- 
          c(temp_gr_str_stanvars_corr, out2 $ temp_gr_str_stanvars)
        temp_gr_str_inits_corr <- 
          c(temp_gr_str_inits_corr,    out2 $ temp_gr_str_inits)
      } # end if(n_higher_str > 0) {
    } # end for (set_randomsi_higher_levsli in set_randomsi_higher_levsl) {
    
    
    # Add temp_gr_str_priors_sd to bpriors above 
    
    higher_level_priors <- temp_gr_str_priors_corr %>% do.call(rbind, .)
    bpriors             <- rbind(bpriors, higher_level_priors)
    
    # Add temp_gr_str_stanvars_sd to stanvar_priors to above 
    
    if(length(temp_gr_str_stanvars_corr) > 0) {
      stanvar_priors_c <- temp_gr_str_stanvars_corr_c <- c()
      for (i in 1:length(stanvar_priors)) {
        stanvar_priors_c <- c(stanvar_priors_c, stanvar_priors[i])
      }
      for (i in 1:length(temp_gr_str_stanvars_corr)) {
        temp_gr_str_stanvars_corr_c <- c(temp_gr_str_stanvars_corr_c, 
                                         temp_gr_str_stanvars_corr[i])
      }
      stanvar_priors <- c(stanvar_priors_c, temp_gr_str_stanvars_corr_c)
    } # if(length(temp_gr_str_stanvars_sd_c) > 0) {
    
    # Add temp_gr_str_inits_sd to initials to above 
    # But this is not good because it sets initials only for the outset level
    # Instead, when 3 or more levels, set init = '0' or  init = 'random' or else 
    # vcov_init_0 = TRUE
    
    if(length(temp_gr_str_inits_corr) > 0) {
      initials_c <- temp_gr_str_inits_corr_c <- c()
      for (i in 1:length(initials)) {
        initials_c <- c(initials_c, initials[i])
      }
      for (i in 1:length(temp_gr_str_inits_corr)) {
        temp_gr_str_inits_corr_c <- c(temp_gr_str_inits_corr_c, 
                                      temp_gr_str_inits_corr[i])
      }
      initials <- c(initials_c, temp_gr_str_inits_corr) 
    } # if(length(temp_gr_str_inits_sd) > 0) {
    
    
    
    
    
    
    #######################################
    
   
    
    priorlist <- rbind(priorlist, bpriors)
    
    # Already added above
    # stanvar_priors <- attr(bpriors, "stanvars")
    
    stanvar_priors_names <- names(stanvar_priors)
    
    if(!"stanvars" %in% attr(stanvar_priors, 'class')) {
      attr(stanvar_priors, 'class') <- c("stanvars", attr(stanvar_priors, 'class'))
    }
    
    prior_stanvarlist[[ii]] <- stanvar_priors 
  
    # Already added above
    # initials <- attr(bpriors, "initials")
    

    
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
      # funlist_r_name <- 'funlist_r'
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
      # funlist_r_name <- paste0("funlist_r", "_", ysi)
    }
    
    # No need for response specific funlist_r_name because already are named
    # Therefore, ignoring above funlist_r_name and reassigning names 
    funlist_r_name <- 'funlist_r'
    funlist_rnamelist[[ii]] <- funlist_r_name
    funlist_rvaluelist[[ii]] <- funlist_r %>% unlist()
    
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
    sigma_groupvarvaluelist[[ii]] <- sigma_arg_groupvar
    
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
  
  # IMP - brms does not allow different lb conditions for sd parsm (e.e, all to be NA)
  # Error: Conflicting boundary information for coefficients of class 'sd'.
  # Because prior function automatically sets lb 0 for positive priors such as exponential
  # the following is need (again done at line 4753 )
  
  brmspriors <- brmspriors %>% 
    dplyr::mutate(lb = dplyr::if_else(class == 'sd', NA, lb))
  brmspriors <- brmspriors %>% 
    dplyr::mutate(ub = dplyr::if_else(class == 'sd', NA, ub))
  
  
  

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
    brms::stanvar(scode = paste(funlist, collapse = "\n"), block = "function")
  
  prior_stanvarlistlist <- c()
  for (i in 1:nys) {
    prior_stanvarlistlist[i] <- paste0("prior_stanvarlist[[", i, "]]")
  }
  
  bstanvars <-
    bstanvars + eval(parse(text = paste(prior_stanvarlistlist, 
                                        collapse = "+")))
  
  
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
      # to exclude student_nu _sd_nu parameters
      brmsinits_names <- brmsinits_names[!grepl('^_nu$|sd_nu', 
                                                brmsinits_names)]
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
      # to exclude student_nu _sd_nu parameters
      brmsinits_names <- brmsinits_names[!grepl('^_nu$|sd_nu', brmsinits_names)]
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
      # to exclude student_nu _sd_nu parameters
      brmsinits_names <- brmsinits_names[!grepl('^_nu$|sd_nu', brmsinits_names)]
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
        # to exclude student_nu _sd_nu parameters
        brmsinits_names <- brmsinits_names[!grepl('^_nu$|sd_nu', 
                                                  brmsinits_names)]
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
    
    
    # this is the right place to replace sd, z and L initials to 0
    # this is for univariate by model. for univariate and multivariate,
    # the already done in set_prior_initials worked. but all now here
    
    temp_stancode2 <- brms::make_stancode(formula = bformula,
                                    stanvars = bstanvars,
                                    prior = brmspriors,
                                    data = brmsdata)
    temp_standata2 <- brms::make_standata(formula = bformula,
                                    stanvars = bstanvars,
                                    prior = brmspriors,
                                    data = brmsdata)
    
    
    if(vcov_init_0e) {
      initialsx2 <- brmsinits
      for (initialsi in names(initialsx2)) {
        if(grepl("sd_", initialsi)) {
          # to exclude student_nu distribution parameter
          if(!grepl("sd_nu", initialsi, fixed = T)) {
            initialsx2[[initialsi]] <- NULL
            newinits <- set_init_gr_effects(temp_stancode2, 
                                            temp_standata2, what = 'sd')
            initialsx2 <- c(initialsx2, newinits)
          }
        }
        if(grepl("L_", initialsi)) {
          initialsx2[[initialsi]] <- NULL
          newinits <- set_init_gr_effects(temp_stancode2, 
                                          temp_standata2, what = 'L')
          initialsx2 <- c(initialsx2, newinits)
        }
        if(grepl("z_", initialsi)) {
          initialsx2[[initialsi]] <- NULL
          newinits <- set_init_gr_effects(temp_stancode2, 
                                          temp_standata2, what = 'z')
          initialsx2 <- c(initialsx2, newinits)
        }
      }
      uni_name <- unique(names(initialsx2))
      initialsx2 <- initialsx2[uni_name] 
      brmsinits <- initialsx2
    } # if(vcov_init_0e) {
    
    
    
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
          setarguments$threads <-  brms::threading(max.threads)
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
  
  # Remove these from the brms args passed via ... -> see also line 1388
  for (collect_dot_namesi in collect_dot_names) {
    if(!is.null(brmsdots_[[collect_dot_namesi]])) 
      brmsdots_[[collect_dot_namesi]] <- NULL
  }
 
  
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
      setdf_1 %>% 
      dplyr::mutate(index = interaction(class, coef, group, nlpar)) %>% 
      dplyr::mutate(order = row_number()) %>% 
      dplyr::arrange(index)
    setdf_2 <- 
      setdf_2 %>% 
      dplyr::mutate(index = interaction(class, coef, group, nlpar)) %>% 
      dplyr::mutate(order = row_number()) %>% 
      dplyr::arrange(index)
    vi_1 <- setdf_1 %>% 
      dplyr::mutate(valid = ifelse(!(class == 'sd' & coef == ""), 1, 0)) %>% 
      data.frame() %>% 
      dplyr::filter(valid == 1) %>% dplyr::select(index) %>% unlist() %>% 
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
      zz[[i]] <- setdf_2 %>% 
        dplyr::mutate(prior = dplyr::if_else(class == getx[['class']] &
                                               coef == getx[['coef']] &
                                               group == getx[['group']] &
                                               resp == getx[['resp']] &
                                               dpar == getx[['dpar']] &
                                               nlpar == getx[['nlpar']],
                                             getx$prior,
                                             setdf_2$prior)) %>% 
        dplyr::filter(class == getx[['class']] &
                        coef == getx[['coef']] &
                        group == getx[['group']] &
                        resp == getx[['resp']] &
                        dpar == getx[['dpar']] &
                        nlpar == getx[['nlpar']])
      
      cc[[i]] <- setdf_2 %>% 
        dplyr::mutate(prior = dplyr::if_else(class != getx[['class']] &
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
    
    brmsfit_sdcor <- do.call(brms::get_prior, brm_args) %>% 
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
     get_stancode & get_standata & get_formula & get_stanvars) {
    stop("Amongst 'get_priors' 'get_set_priors', 'validate_priors' ",
         "\n ",
         "'get_stancode', 'get_standata', 'get_formula', 'get_stanvars' ",
         "\n ",
         " arguments, only one can be set to TRUE at a time")
  }
  
  
  
  exe_model_fit <- TRUE
  if(get_stancode |
     get_standata |
     get_formula |
     get_stanvars |
     get_priors |
     get_set_priors |
     validate_priors |
     get_set_init) {
    exe_model_fit <- FALSE
  }
  
  
  
  # IMP - brms does not allow different lb conditions for sd parsm (all to NA)
  # Error: Conflicting boundary information for coefficients of class 'sd'.
  # Because prior function automatically sets lb 0 for positive priors such as 
  # exponential the following is need (again done at line 4002)
  
  lbbb_ <- ubbb_ <- NULL
  tempprior_hold <- brmspriors # brm_args$prior 
  setpriornamesorder <- colnames(tempprior_hold)
  tempprior_hold$lbbb_ <- tempprior_hold$lb
  tempprior_hold$ubbb_ <- tempprior_hold$ub
  tempprior_hold$lb <- tempprior_hold$ub <- NULL
  tempprior_hold <- tempprior_hold %>% 
    dplyr::mutate(lbbb_ = dplyr::if_else(class == 'sd', NA, lbbb_))
  tempprior_hold <- tempprior_hold %>% 
    dplyr::mutate(ubbb_ = dplyr::if_else(class == 'sd', NA, ubbb_))
  tempprior_hold$lb <- tempprior_hold$lbbb_
  tempprior_hold$ub <- tempprior_hold$ubbb_
  tempprior_hold$lbbb_ <- tempprior_hold$ubbb_ <- NULL
  tempprior_hold <- tempprior_hold %>% 
    dplyr::relocate(dplyr::all_of(setpriornamesorder))
  brmspriors <-   tempprior_hold
  
  
  
  if(!is.null(set_self_priors) & is.null(set_replace_priors)) {
    brmspriors <- set_self_priors
  }
  
  
  
  # if(!is.null(set_replace_priors) & is.null(set_self_priors)) {
  #   brmspriors <- brmspriors %>%
  #     dplyr::filter(source == 'user') %>%
  #     dplyr::bind_rows(., set_replace_priors)
  # }
  
  
  if(!is.null(set_replace_priors) & is.null(set_self_priors)) {
    if(!set_same_priors_hierarchy) {
      drop_old_sd_cor_groups <- unique(setup_new_priors_add[['group']])
      brmspriors <- brmspriors %>% filter(! group %in% drop_old_sd_cor_groups)
      brmspriors <- brmspriors %>%
        dplyr::filter(source == 'user') %>%
        dplyr::bind_rows(., set_replace_priors)
    } else if(set_same_priors_hierarchy) {
      brmspriors <- brmspriors %>%
        dplyr::filter(source == 'user')
    }
  }
  
  
  
  
  if(is.null(set_self_priors) & is.null(set_replace_priors)) {
    brmspriors <- brmspriors
  }
  
  
  brm_args$prior <- brmspriors
  
    
  if(!exe_model_fit) {
    if(get_priors) {
      options(mc.cores = mc.cores_restore)
      return(do.call(brms::get_prior, brm_args))
    } else if(get_standata) {
      options(mc.cores = mc.cores_restore)
      return(do.call(brms::make_standata, brm_args))
    } else if(get_stancode) {
      options(mc.cores = mc.cores_restore)
      return(do.call(brms::make_stancode, brm_args))
    } else if(get_set_priors) {
      options(mc.cores = mc.cores_restore)
      return(brm_args$prior)
    } else if(validate_priors) {
      options(mc.cores = mc.cores_restore)
      return(do.call(brms::validate_prior, brm_args))
    } else if(get_set_init) {
      return(brm_args$init)
    } else if(get_formula) {
      return(brm_args$formula)
    } else if(get_stanvars) {
      return(brm_args$stanvars)
    }
  } 
  
  
  
  
  # Fit model if get_set_priors get_priors get_standata get_stancode -> FALSE
  # this if(exe_model_fit) { is closed just before the end of the bsitar 
  
  if(exe_model_fit) {
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
      } else if(!is.list(brm_args$init) & length(brm_args$init) == 1) {
        if(brm_args$init == "0") {
          init_custom <- NULL
        } else if(brm_args$init == "random") {
          init_custom <- NULL
        } else if(brm_args$init == 0) {
          init_custom <- NULL
        } else {
          init_custom <- init_custom
        }
      } else if(!is.null(brm_args$init)) {
        init_custom <- init_custom
      }
    } # if(brm_args$backend == "cmdstanr") {
    
    
    
    
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
               " be a single named list (e.g., custom_init = list(x= 2,xx=5)) ",
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
    
    
    
    brmsfit <- do.call(brms::brm, brm_args)
  
    # Add model info for post-processing
    
    model_info <- list()
    
    for (i in 1:length(funlist_rnamelist)) {
      model_info[[funlist_rnamelist[[i]]]] <- funlist_rvaluelist[[i]]
    }
    
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
      model_info[[hierarchicalvarnamelist[[i]]]] <- 
        hierarchicalvarvaluelist[[i]]
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
    
    # model_info[[SplineFun_name]] <- SplineFun_name
    model_info[['StanFun_name']] <- SplineFun_name
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
    
    model_info[['brms_arguments_list']] <- brms_arguments_list
    
    model_info[['select_model']] <- select_model
    
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
      
      brmsfit <- expose_bsitar_functions(brmsfit, expose = TRUE)
      brmsfit$model_info[['expose_method']] <- 'S'
    } # if (expose_function) {
    
    
    
    if (!expose_function) {
      brmsfit <- expose_bsitar_functions(brmsfit, expose = FALSE, 
                                         select_model = select_model)
      brmsfit$model_info[['expose_method']] <- 'R'
    } # if (!expose_function) {
    
    
    
    
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


