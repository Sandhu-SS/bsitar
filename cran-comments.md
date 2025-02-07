
# The 'bsitar' package (version 0.3.2) was pre-test archived on 6 February 2025 because of the following 
Package size greater than 5 MB 
One example run took mopre than 5 sec

Both these issues have been fixed.


# The 'bsitar' package (version 0.1.1) was archived on 10 March 2024 because of the following 
error when running model2 <- update_model(model, df = 5, sample_prior = 'only') on
R Under development (unstable) 

Error in .make_numeric_version(x, strict, .standard_regexps()$valid_numeric_version) : 
  invalid non-character version specification 'x' (type: double)
Calls: update_model ... as.numeric_version -> numeric_version -> .make_numeric_version
Execution halted

The issue has beeb resolved and the package has beeb tested again (please see details below) 


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
 Windows 11 x64 (build 22621) 
 platform: x86_64-w64-mingw32 (64-bit)

R CMD check results
0 errors | 0 warnings | 0 notes



# Test environments 2 - via rhub::check_for_cran()
 Windows Server 2022 x64 (build 20348)
 platform: x86_64-w64-mingw32

Status: success
Notes:
 Found the following files/directories: ''NULL''
 Found the following files/directories: 'lastMiKTeXException'

Note about ''NULL'' is possibaly a bug/crash in MiKTeX as reported on R-hub [issue #503](https://github.com/r-hub/rhub/issues/503), and therefore can be ignored 

Note about 'lastMiKTeXException' can be ignored as it seems to be an Rhub issue [issue #560](https://github.com/r-hub/rhub/issues/560)



# Test environments 3 - via rhub::check_for_cran()
 Ubuntu 20.04.6 LTS 
 platform: x86_64-pc-linux-gnu (64-bit)

Status: success
Notes:
 Skipping checking HTML validation: no command 'tidy' found

Note about 'tidy' seems to be a recurring issue on Rhub R-hub [issue #548](https://github.com/r-hub/rhub/issues/548), and therefore can be ignored 



# Test environments 4 - via rhub::check_for_cran()
 Fedora Linux 36 (Container Image)
 platform: x86_64-pc-linux-gnu

Notes:
 Skipping checking HTML validation: no command 'tidy' found

Note about 'tidy' seems to be a recurring issue on Rhub R-hub [issue #548](https://github.com/r-hub/rhub/issues/548), and therefore can be ignored 
