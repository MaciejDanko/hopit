devtools::document()
roxygen2::roxygenise()
#RcppEigen::RcppEigen.package.skeleton()

load("~/Documents/R-PRJ/gotm/data/healthsurvey.rda")
devtools::use_data(healthsurvey, overwrite = TRUE)