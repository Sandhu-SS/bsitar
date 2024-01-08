## version 0.1.0.9000

---

### New feature

- Added feature to perform model optimization by fitting model with varying
  degree of freedom and by transforming 'x' (predictor such as age) and 'y' 
  (outcome such repeated height measurements on individuals) varibales. The    
  allowed transformations for 'x' and 'y' variables are the 'log' and 'sqrt'  
  transformation. The 'optimize_model' function implements the optimization.  
  Please see help file of 'optimize_model' to see documentation.

### Minor changes

- Updated documentation
- Added issue tracker url: https://github.com/Sandhu-SS/bsitar/issues


### Bugfixes

- plot_curves() Fixed bug that lead to Warning: Duplicated aesthetics after name 
standardisation: group in when plotting together unadjusted and adjusted curves. 
- plot_curves() Fixed bug that caused palettes error when tried plotting all 
four curves together (distance, velocity, adjusted and unadjusted)
- plot_curves() and growthparameters() fixed issues relating to not using 
'dplyr::all_of()' when within the 'dplyr::select()'


### Miscellaneous
- Minor corrections/changes to make R code more efficient and consistent across sub modules.



