#

#RcppEigen::RcppEigen.package.skeleton()

load("~/Documents/R-PRJ/hopit/data/healthsurvey.rda")
devtools::use_data(healthsurvey, overwrite = TRUE)
roxygen2::roxygenise()
devtools::document()
devtools::source_package()
healthsurvey
browseVignettes("hopit")

R CMD Rd2dvi --pdf --title='Test of hopit' -o /tmp/hopit.pdf man/*.Rd
