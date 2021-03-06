# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}

install.packages(c("devtools", "tidyverse", "fs"))

?use_mit_license
usethis::use_mit_license("Yue Zhou")

usethis::use_package("mclust", type = "Imports")
usethis::use_package("dplyr", "Suggests")

usethis::use_roxygen_md()

usethis::use_readme_rmd()
devtools::build_readme()

devtools::load_all()
devtools::document()


