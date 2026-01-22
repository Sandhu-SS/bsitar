


   # devtools::load_all()

##############################################################################
# set options
##############################################################################

# For current settings (model and draw_ids) min size should be at least 980MB
oopts <- options(future.globals.maxSize = 2.0 * 1e9)  ## 2.0 GB
on.exit(options(oopts), add = TRUE)

estimate_center   <- "mean"
estimate_interval <- "eti"
ec_ <- getOption("marginaleffects_posterior_center")
options("marginaleffects_posterior_center" = estimate_center)
on.exit(options("marginaleffects_posterior_center" = ec_), add = TRUE)

ei_ <- getOption("marginaleffects_posterior_interval")
options("marginaleffects_posterior_interval" = estimate_interval)
on.exit(options("marginaleffects_posterior_interval" = ei_), add = TRUE)

ec_agg <- getOption("marginaleffects_posterior_center")
ei_agg <- getOption("marginaleffects_posterior_interval")

if(is.null(ec_agg)) ec_agg <- "mean"
if(is.null(ei_agg)) ei_agg <- "eti"


##############################################################################
# set options
##############################################################################

# draw_ids <- 1:50

test_tolerance <- 0.01

if(test_univariate_fit_cov) {
  fit               = readRDS(testthat::test_path("models", 
                                                  "univariate_fit_cov.rds")) 
  resp              = uvar_resp
} else if(test_multivariate_fit_cov) {
  fit               = readRDS(testthat::test_path("models", 
                                                  "multivariate_fit_cov.rds"))  
  resp              = mvar_resp
} else {
  skip(message = 
         "Both test_univariate_fit_cov and test_multivariate_fit_cov FALSE")
}

# Need to re-assign functions to this local test environment
fit <- bsitar::expose_model_functions(fit, expose = F)

##############################################################################
# set options
##############################################################################

# brms 
draw_ids          = draw_ids
ndraws            = NULL
re_formula        = NA
dpar              = NULL
nlpar             = NULL
incl_autocor      = TRUE
allow_new_levels  = FALSE
sample_new_levels = "uncertainty"

##############################################################################
# set options
##############################################################################

# Key arguments -> by, cov, deriv, variables, condition, comparison, use_d1
xvar              = 'age'
cov               = NULL
deriv             = 1 
variables         = NULL
by                = FALSE
byfun             = NULL
comparison        = NULL # deriv = 0 = 'difference'; deriv > 0 = 'dydx'

# Keep it NULL / FALSE to match marginaleffects
use_d1            = NULL # TRUE use _d1, FALSE USES d0 - dydx, NULL == FALSE

##############################################################################
# set options
##############################################################################

model             <- fit
newdata           <- bsitar:::get.newdata(model)

# marginal future
set_future_method  <- 'future' 
set_future_session <- 'sequential' 
set_future_splits  <- TRUE
set_future_cores   <- 2

brms_args_names <- c('draw_ids', 'ndraws', 're_formula', 
                     'resp', 'dpar', 'nlpar', 'incl_autocor',
                     'allow_new_levels', 'sample_new_levels')


print_console     <- FALSE # set it TRUE outside testthat 


###############################################################################
###############################################################################
# marginaleffects::comparisons vs marginal_comparisons - average = FALSE
###############################################################################
###############################################################################

test_str_cat <- 
  "marginaleffects::comparisons vs marginal_comparisons -> average = FALSE"


htest_funcall                <- hypothesis_test.bgmfit
brms_funcall                 <- brms:::hypothesis.brmsfit
marginaleffects_funcall      <- marginaleffects::hypotheses
bayestestr_funcall           <- bayestestR:::equivalence_test.brmsfit
bayestestrparameters_funcall <- bayestestR:::bayesfactor_parameters.brmsfit



stop()







#############################################################################

# 
# devtools::load_all()
# 
# 
# equivalence_range <- get_test_range_null(
#   parameter = list(apgv = 1, pgv = 0),
#   sex       = list(Male = NULL, Female=1),
#   study     = list(S1 = 10, S2 = 20),
#   cohort    = list(C1 = 100, C2 = 200),
#   what = 'range',
#   data = NULL,
#   full_frame = T,
#   verbose = FALSE
# )
# 
# equivalence_range <- get_test_range_null(
#   parameter = list(apgv = 1, pgv = 0),
#   sex       = list(Male = NULL, Female=NULL),
#   full_frame = T,
#   verbose = FALSE
# )
# 
# 
# # 
# # 
# # res <- create_range_lists(
# #   parameter = list(apgv = 1, pgv = 0.5),
# #   sex       = list(Male = c(1.0, 0.5), Female=c(1.0, 0.5)),
# #   study     = list(S1 = c(0.2, 8, 3, 4), S2 = NULL),
# #   cohort    = list(C1 = NULL, C2 = NULL)
# #   )
# options(marginaleffects_safe=FALSE)
# aa <- 
#   marginal_growthparameters(fit,
#                             variables = c('sex'),
#                             by = c('sex'),
#                             comparison = "difference", 
#                             parameter = c('apgv', 'pgv'), 
#                             # method = 'custom', hypothesis = 'pairwise',
#                             hypothesis = difference ~ pairwise ,
#                             # equivalence = c(0,0),
#                             equivalence_test =
#                               list(
#                                 # range = list(c(0, 0), c(2, 3.0)),
#                                 range = equivalence_range,
#                                 # by = NULL,
#                                 inline = F, # if TRUE, evaluted inline, limited
#                                 format = NULL, # if TRUE, merge range hdi
#                                 get_form = NULL, # if TRUE, return range null format NA
#                                 get_value = NULL, # if TRUE, return range null format values
#                                 digits = 2,
#                                 as_percent = TRUE,
#                                 ci = 0.95),
#                             p_direction =
#                               list(
#                                 method = "direct",
#                                 inline = F, # if TRUE, evaluted inline, limited
#                                 null = 0,
#                                 as_p = FALSE,
#                                 as_percent = TRUE,
#                                 remove_na = TRUE,
#                                 rvar_col = NULL,
#                                 digits = NULL, 
#                                 get_form = NULL, # if TRUE, return range null format
#                                 get_value = NULL, # if TRUE, return range null format values
#                                 format = NULL # if TRUE, merge range hdi
#                               ),
#                             draw_ids = draw_ids)



##############################################################################

















if(exists('marginaleffects_args'))      rm('marginaleffects_args')
if(exists('marginal_args'))             rm('marginal_args')

if(exists('marginaleffects_out'))              rm('marginaleffects_out')
if(exists('marginal_out'))                     rm('marginal_out')
if(exists('marginal_out_custom_mdT'))          rm('marginal_out_custom_mdT')
if(exists('marginal_out_custom_mdF'))          rm('marginal_out_custom_mdF')
if(exists('marginal_out_custom_mdT_future'))   rm('marginal_out_custom_mdT_future')
if(exists('marginal_out_custom_mdF_future'))   rm('marginal_out_custom_mdF_future')
if(exists('marginal_out_custom_mdT_dofuture')) rm('marginal_out_custom_mdT_dofuture')
if(exists('marginal_out_custom_mdF_dofuture')) rm('marginal_out_custom_mdF_dofuture')


###############################################################################
# Set up arguments
###############################################################################

all_args_names <- NULL
all_args_names <- formalArgs(brms_funcall)
all_args_names <- c(all_args_names, brms_args_names)

marginaleffects_args <- list()
for (i in all_args_names) {
  if(i != "...") {
    if(!exists(i)) {
      stop("arg ", i, " not defined in the workspace")
    }
    marginaleffects_args[[i]] <- get(i)
  } # if(i != "...") {
} # for (i in all_args_names) {



###############################################################################
# hypothesis_test - brms:::hypothesis.brmsfit
###############################################################################





###############################################################################
# hypothesis_test - marginaleffects::hypotheses
###############################################################################



aa <- marginaleffects::comparisons(fit,
                                   variables = 'sex',
                                   by = 'sex',
                                   comparison = "difference",
                                   draw_ids = draw_ids)

mfx <- marginaleffects::get_draws(aa, shape = "DxP")

if(!is.null(attr(mfx, 'class'))) {
  attr(mfx, 'class') <- c(attr(mfx, 'class'), 'marginaleffects')
}


hypothesis_test(model = aa, 
                variables = "sex", by = "sex",
                hypothesis = ("a_Intercept - a_sexFemale > .2"),
                equivalence_test = 
                  list(
                    range = list(c(0, 0), c(2, 3.0)),
                    by = c("contrast"),
                    ci = 0.95),
                engine = 'data.frame',
                verbose = TRUE) # %>% str()


###############################################################################
# hypothesis_test - bayestestR:::equivalence_test.brmsfit
###############################################################################




###############################################################################
# hypothesis_test - bayestestR:::bayesfactor_parameters.brmsfit
###############################################################################




###############################################################################
# hypothesis_test - 
###############################################################################




###############################################################################
# hypothesis_test - 
###############################################################################







