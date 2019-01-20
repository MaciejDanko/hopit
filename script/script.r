#
library(tinytex)
library(roxygen2)
library(Rdpack)
#RcppEigen::RcppEigen.package.skeleton()
#tinytex::install_tinytex(force=TRUE)
#tinytex::tlmgr(c("info", "--list", "--only-installed", "--data", "name"))

{load("~/Documents/R-PRJ/hopit/data/healthsurvey.rda")
healthsurvey$contHM <- NULL
devtools::use_data(healthsurvey, overwrite = TRUE)
roxygen2::roxygenise()
}

#R CMD INSTALL --no-multiarch --with-keep.source hopit
devtools::clean_vignettes()
devtools::build_vignettes()

#CMD build '/home/maciej/Documents/R-PRJ/hopit' --no-resave-data
#tinytex::tlmgr()javascript:;
devtools::check()
devtools::document()
#devtools::source_package()
#healthsurvey
browseVignettes("hopit")

#R CMD Rd2dvi --pdf --title='Test of hopit' -o /tmp/hopit.pdf man/*.Rd

R CMD Rd2pdf --title="Package `hopit'" -o ~/Documents/R-PRJ/hopit/inst/doc/hopit.pdf man/*.Rd
