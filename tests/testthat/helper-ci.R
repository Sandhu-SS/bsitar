



# custom skip helper that runs tests only on CI and not on local R CMD Checks

skip_if_not_ci <- function() {
  ci <- Sys.getenv("CI")
  if (identical(ci, "") || identical(ci, "false") || identical(ci, "0")) {
    testthat::skip("Skipping: not on CI")
  }
}


# Usage 
# skip_if_not_ci()
# 
# test_that("big model behaviour", {
#   # heavy fits / large objects tested here
# })



