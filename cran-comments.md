
This is a new submission, therefore 1 note across all test environments.

# Test environments (a total four)
 1. Windows 11 x64 (build 22621)
     	on devtools::check()
 2. Windows Server 2022 x64 (build 20348)
	on rhub::check_for_cran()
 3. Ubuntu 20.04.6 LTS, GCC
	on rhub::check_for_cran()
 4. Fedora Linux 36 (Container Image), clang, gfortran
	on rhub::check_for_cran()



# Test environments 1 - via devtools::check()
# Windows 11 x64 (build 22621) 
# platform: x86_64-w64-mingw32 (64-bit)

R CMD check results
0 errors | 0 warnings | 0 notes



# Test environments 2 - via rhub::check_for_cran()
# Windows Server 2022 x64 (build 20348)
# platform: x86_64-w64-mingw32

Status: success
Three notes:
 Note 1. New submission - Maintainer: 'Satpal Sandhu <satpal.sandhu@bristol.ac.uk>'
 Note 2. Found the following files/directories: ''NULL''
 Note 3. Found the following files/directories: 'lastMiKTeXException'

Note 2 is possibaly a bug/crash in MiKTeX as reported on R-hub [issue #503](https://github.com/r-hub/rhub/issues/503), and therefore can be ignored 

Note 3 can be ignored as it seems to be an Rhub issue [issue #560](https://github.com/r-hub/rhub/issues/560)



# Test environments 3 - via rhub::check_for_cran()
# Ubuntu 20.04.6 LTS 
# platform: x86_64-pc-linux-gnu (64-bit)

Status: success
Two notes:
 Note 1. New submission - Maintainer: 'Satpal Sandhu <satpal.sandhu@bristol.ac.uk>'
 Note 2. Skipping checking HTML validation: no command 'tidy' found

Note 2 seems to be a recurring issue on Rhub R-hub [issue #548](https://github.com/r-hub/rhub/issues/548), and therefore can be ignored 



# Test environments 4 - via rhub::check_for_cran()
# Fedora Linux 36 (Container Image)
# platform: x86_64-pc-linux-gnu

Two notes:
 Note 1. New submission - Maintainer: 'Satpal Sandhu <satpal.sandhu@bristol.ac.uk>'

 Note 2. Skipping checking HTML validation: no command 'tidy' found

Note 2 seems to be a recurring issue on Rhub [issue #548](https://github.com/r-hub/rhub/issues/548), and therefore can be ignored  




# Reviewer's comments
1. Comment: 
     "Please do not start the description with "This package", package name,
      title or similar."

   Response: 
     Edited, as suggested

2. Comment: 
      "You are setting options(warn=-1) in your function. This is not allowed.
	Please rather use suppressWarnings() if really needed. -> R/plot_curves.R

	Please make sure that you do not change the user's options, par or
	working directory. If you really have to do so within functions, please
	ensure with an *immediate* call of on.exit() that the settings are reset
	when the function is exited.
	e.g.:...
	old <- options() # code line i
	on.exit(options(old))# code line i+1
	...
	options(...)=# somewhere after

	e.g.: R/bgm.R; R/plot_curves.R
	If you're not familiar with the function, please check ?on.exit. This
	function makes it possible to restore options before exiting a function
	even if the function breaks. Therefore it needs to be called immediately
	after the option change within a function."

    Response: 
        Changed, as suggested. Now using on.exit(options(oldopts)) to 
        reset options for both R/bgm.R; R/plot_curves.R  


3. Comment: 
	"Please do not modifiy the .GlobalEnv. This is not allowed by the CRAN
	policies.
	->R/fitted_draws.R, R/growthparameters.R, R/loo_validation.R,
	R/plot_conditional_effects.R, R/plot_curves.R, R/plot_ppc.R,
	R/predict_draws.R, R/utils-helper-2.R, R/utils-helper-4.R"

   Response: 
        Corrected. 

4. Comment: 
	"Please do not install packages in your functions, examples or vignette.
	This can make the functions,examples and cran-check very slow. ->
	R/utils-helper-1.R"

   Response: 
	Corrected, the packages are not installed in function (R/utils-helper-1.R).

