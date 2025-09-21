# bsitar (development version)

### Breaking changes

A major refactoring of the internal code to streamline the integration of modelling the distributional parameter sigma, aligning it more closely with the location parameter mu. This update marks the first step in a series of improvements aimed at enabling robust modelling of location-scale models. Additionally, the core functions are significantly rewritten to enhance the speed and efficiency of model fitting. Note: This release should be considered experimental. Future versions may introduce further significant changes. Additionally, documentation has not yet been fully updated to reflect these modifications.

Due to the above changes, previously saved model objects will need to be refitted for post-processing to work correctly.
No other changes are required from the user’s perspective—only that the model should be re-run.

Major change part - 2

## New features/Additions

Added experimental support for estimating model-based individual growth parameters, such as age at peak growth velocity (APGV), as well as distance and velocity at APGV. For details, refer to the function ``modelbased_growthparameters()``. Note that the function name and its arguments may change in future releases.

Support has been added to allow the use of external functions, such as ``splines::ns()``, for modeling the distributional parameter ``sigma``. Users can now specify functions for both the fixed ``sigma_formula`` and random ``sigma_formula_gr`` effects. Further, different functions can be used in the fixed and random effects formulas. For example. For example. ``sigma_formula = ~ 1 + splines::ns(age, df = 3)``, and ``sigma_formula_gr = ~ 1 + stats::poly(age, degree = 2)``

Added support for ``tag`` feature implementing parameter specific prior sensitivity analysis in ``priorsense`` package via ``brms``.


### Minor changes/Enhancements

The efficiency of the post-processing function has been improved (average improvement in speed ~ 2x).

The ``marginal_comparison()`` function has been renamed as ``marginal_comparisons()`` to better reflect the wide range of comparisons offered. Hence, the name has been changed from **marginal_comparison** to **marginal_comparisons** (plural).

Also, the ``growthparameters_comparison()`` function has been renamed to ``marginal_growthparameters()`` making it clear that this function, like ``marginal_comparisons()`` and ``marginal_draws()``, is also based on the ``marginaleffects`` package.

For backward compatibility, the old functions ``marginal_comparison()`` and ``growthparameters_comparison()`` will be included as aliases for the new function names, i.e., ``marginal_comparisons()`` and ``marginal_growthparameters()``.

The default ``stype`` set to ``nsk`` instead of ``nsp``. 

 

### Bugfixes

The function block for the multivariate model did not render properly.

The random initial values for the residual correlation parameter of the multivariate model resulted in an empty list warning.


# bsitar 0.3.2


### New features/additions

  The ``bsitar`` now supports three different types of splines: ``'rcs'``, ``'nsp'``, and ``'nsk'``. 
  While ``'rcs'`` constructs the spline design matrix using the truncated power basis, both ``'nsp'`` 
  and ``'nsk'`` implement a B-spline-based natural cubic spline basis. The truncated power basis method,
  often referred to as Harrell's method, is implemented in the ``rcspline.eval()`` function of the 
  ``Hmisc`` package. The B-spline-based implementations of ``'nsp'`` and ``'nsk'`` are the same as 
  those described in the ``splines2`` package. Previously, only ``'rcs'`` was available. Now, the 
  default method is ``'nsp'``.

  The ``bsitar`` package allows for fitting the SITAR model with a four-parameter formulation ``(a + b + c + d)``, 
  where a fourth parameter, ``d``, is added. The parameter ``d``, the age slope, allows the adult part  of the growth 
  curve to vary in slope. The age slope ``d`` represents the regression coefficient of ``y`` on ``x``.
  Like the ``sitar`` package, the ``bsitar`` offers two parameterizations: one in which ``x`` is adjusted for the random 
  effects ``b`` and ``c`` i.e., ``(x-b)*exp(c)``, and the other in which the unadjusted ``x`` is used. This is controlled 
  by the argument ``d_adjusted`` provided in the ``bsitar::bsitar()`` function. When ``d_adjusted = TRUE``, the first 
  version is applied (``x`` adjusted for random effects ``b`` and ``c``), which means that individual developmental age, 
  rather than chronological age, is used in the slope regression. This makes parameter ``d`` more sensitive to the timing 
  of puberty in individuals. When ``d_adjusted = FALSE`` (the default), the unadjusted ``x`` is used. The default choice, 
  ``d_adjusted = FALSE``, is primarily to match the behavior of the ``sitar`` package, which sets ``d.adjusted = FALSE``.

  Added support for computing and comparing growth curves using the ``marginaleffects`` package as 
  the back-end (see ``marginal_draws()``, ``marginal_comparison()``, and ``growthparameters_comparison()``). 
  This enables the use of the computational flexibility offered by the ``marginaleffects`` package to 
  estimate various quantities of interest, such as adjusted growth curves (distance and velocity), and 
  growth parameters like age at peak growth velocity. All three functions support parallel computation via 
  the ``future`` and ``doFuture`` packages.

  The ``optimize_model()`` function now allows users to specify custom functions in ``optimize_x`` and ``optimize_y`` 
  when optimizing the Bayesian SITAR model. For example, it is now possible to use 
  ``optimize_x = list(function(x) log(x + 3/4))``. Thanks to Tim Cole for suggesting this feature. 
  This update greatly enhances the flexibility of ``optimize_model()`` and enables users to search for a 
  range of optimal ``x`` and ``y`` transformations.
   
    
  An experimental feature has been added to use ``$pathfinder()`` based initial values for the MCMC sampling
  ``$sample()`` (via the argument ``pathfinder_init = TRUE``, default is FALSE). The arguments for
  ``$pathfinder()`` can be specified as a named list using ``pathfinder_args``. Note that this feature is 
  only available when ``backend = 'cmdstanr'``.

   Added a new vignette comparing growth curves and growth parameters obtained from the frequentist 
   (``sitar`` package) and Bayesian (``bsitar`` package) versions of the SITAR model applied to the 
   height data. 

### Minor changes

 The prior distribution for each parameter has been changed from ``student_t()`` to ``normal()``.
 The prior distribution for all parameters, including regression coefficients as well as the standard
  deviation (sd) for the group-level random effects and the distributional parameter (sigma), has been 
  changed to ``normal()``. Previously, the distribution for regression coefficients and the sd for the 
  group-level random effects was ``student_t()``, while the distribution for the sd of the distributional
  parameter (sigma) was ``exponential()``. Note that the same location and scale parameters used earlier 
  for ``student_t()`` are now used for the ``normal()`` distribution. For example, prior specified
  earlier as ``student_t(3, 0, 15)`` (3 degree of freedom, location 0 and scale) has been revised as 
  ``normal(0, 15)``.

 The default setting for initial values is now ``random`` except for the population average parameters: 
 size (``a_init_beta``), timing (``b_init_beta``), intensity (``c_init_beta``), and spline coefficients 
 (``s_init_beta``). For size and spline coefficient parameters, the initial values are derived from the 
 linear regression fit and specified as ``a_init_beta = lm`` and ``s_init_beta = lm``. The initial 
 values for both timing and intensity parameters are set to '0', i.e., ``b_init_beta = 0`` and 
 ``c_init_beta = 0``.

 Improved documentation

### Bugfixes

 The ``sigma_cov_init_beta = random`` argument was setting incorrect initial values for the covariates 
 included in the ``sigma`` formula. The initial values for the Intercept (``sigma_init_beta``) were 
 incorrectly applied to the covariates as well.



### Miscellaneous
Users no longer need to set the environment to ``globalenv()``, i.e., ``envir = globalenv()``, for 
post-processing functions to work properly. The environment is now automatically set to match the 
environment of the exposed functions. It is important to note that manually setting the environment 
(via the ``envir`` argument) may result in errors. The ``envir`` argument is now primarily for internal
use, which is needed during tests.

Minor corrections and changes have been made to improve the efficiency of the R code.



# bsitar 0.2.1


### New feature

- Added ```add_model_criterion()``` function to compute fit criteria such as "loo", "waic", "kfold", 
 "loo_subsample", "bayes_R2" (Bayesian R-squared), "loo_R2" (LOO-adjusted R-squared), and "marglik" 
 (log marginal likelihood). The computed fit criteria are added to the model object for later use including 
 comparison of models.  The ```add_model_criterion()``` is a wrapper around the ```add_criterion()``` 
 function available from the  **brms** package


### Bugfixes

- bsitar(): The 'sigma_cov_init_beta = random' argument was setting wrong initial values for the covariates  
  included in the 'sigma' formula. The initials for Intercept ('sigma_init_beta') were used for covariates too.



### Miscellaneous
- Replaced Bayesian SITAR model fit shown as example from 'berkeley_mfit' applied to 20 randomly selected males   between 6 and 20 years of age ('berkeley_mdata') with 'berkeley_exfit' that is fit to 70 females between 8 and 20 years of age ('berkeley_exdata'). This is done to use the example model ('berkeley_exfit') in vignette that provide a detailed comparison between non Bayesian SITAR model fit (using the 'sitar' package) and Bayesian SITAR model fit (using the 'bsitar' package). The vignette included in the 'sitar' package analysed the exact same data (70 females between 8 and 20 years of age). 
- The ``bsitar::bsitar()```received options ```file```, ``file_refit```, and ``file_compress``` to save and retreive fitted objects. See  ```brms::brm``` help file for details. 
- Minor corrections/changes to make R code more efficient.


# bsitar 0.1.1


---

### New feature

- Added 'optimize_model' function to perform model optimization by fitting 
  model with varying degree of freedom and by transforming 'x' (predictor) 
  and 'y' (outcome) varibales. The allowed transformations for 'x' and 'y' 
  variables are 'log' (logarithmic) and 'sqrt' (square root) transformation. 
  The 'optimize_model' performs comparison of resulting model fits based on 
  user specified  criteria such as the Watanabe–Akaike information criterion
  ('waic'), leave-one-out cross-validation ('loo') and the Bayesian R square 
  ('bayes_R2'). Please see help file of 'optimize_model' to see documentation.

### Minor changes

- Updated documentation
- Added vignette about the Bayesian SITAR model fit to height data (univariate model)
- Added issue tracker url: https://github.com/Sandhu-SS/bsitar/issues


### Bugfixes

- plot_curves(): Fixed bug to remove warning "Duplicated aesthetics after name 
standardisation: group in when plotting together unadjusted and adjusted curves." 
- plot_curves(): Fixed bug to remove palettes error when tried plotting all 
four curves together (distance, velocity, adjusted and unadjusted)
- plot_curves(): and growthparameters() fixed issues relating to not using 
'dplyr::all_of()' when within the 'dplyr::select()'


### Miscellaneous
- Added more utility functions for internal use (for consistency and efficiency) 
- Minor corrections/changes to make R code more efficient and consistent across sub modules.



