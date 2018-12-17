#devtools::document()

#RcppEigen::RcppEigen.package.skeleton()

load("~/Documents/R-PRJ/hopit/data/healthsurvey.rda")
devtools::use_data(healthsurvey, overwrite = TRUE)
roxygen2::roxygenise()
healthsurvey
