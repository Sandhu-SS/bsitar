## version 0.1.0.9000

---


### Minor changes

- Updated documentation
- Typo correction in help files ('.rd' files) 


### Bugfixes

- plot_curves() Fixed bug that lead to Warning: Duplicated aesthetics after name 
standardisation: group in when plotting together unadjusted and adjusted curves. 
- plot_curves() Fixed bug that caused palettes error when tried plotting all 
four curves together (distance, velocity, adjusted and unadjusted)
- plot_curves() and growthparameters() fixed issues relating to not using 
'dplyr::all_of()' when within the 'dplyr::select()'



