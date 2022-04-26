## ====================================================
## R packages needed not available via mini-conda
## ====================================================

## ____________________________________________________
## Setting-up the source repository
local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cloud.r-project.org"
  options(repos = r)
})
## ____________________________________________________

## ____________________________________________________
## Additional R packages needed by the user (CRAN)
install.packages("bookdown")
install.packages("mvtnorm")
install.packages("xtable")
install.packages("ranger")
install.packages("parallel")
## ____________________________________________________

## ____________________________________________________
## Additional R packages needed by the user ()
## remotes::install_github("user/package")
## ____________________________________________________
