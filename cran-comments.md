
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
