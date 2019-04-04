#' Calculation of the variance-covariance matrix for a specified survey design (experimental function)
#'
#' @param vcovMat a variance-covariance matrix.
#' @param estfun a gradient function of the log-likelihood function.
#' @param design a \code{survey.design} object.
#' @description
#' This function is an equivalent of \code{survey:::svy.varcoef}. In the original approach \code{estfun} is calculated from
#' glm's working residuals:\cr
#' \code{estfun <- model.matrix(glm.object) * resid(glm.object, "working") * glm.object$weights}\cr
#' In the hopit package, estfun is directly calculated as a gradient (vector of partial derivatives) of the log likelihood function.
#' Depending on detected design an appropriate \code{survey} function is called.
#' @seealso
#' \code{\link[survey]{svydesign}}
#' \code{\link{hopit}}
#' @importFrom survey svyrecvar twophasevar twophase2var svyCprod
svy.varcoef_hopit <- function (vcovMat, estfun, design) {
  x <- estfun %*% vcovMat
  if (inherits(design, "survey.design2"))
    varcoef <- survey::svyrecvar(x = x,
                                 clusters = design$cluster,
                                 stratas = design$strata,
                                 fpcs = design$fpc,
                                 postStrata = design$postStrata)
  else if (inherits(design, "twophase"))
    varcoef <- survey::twophasevar(x, design)
  else if (inherits(design, "twophase2"))
    varcoef <- survey::twophase2var(x, design)
  else if (inherits(design, "pps"))
    stop(hopit_msg(89),call.=NULL)
  else varcoef <-survey::svyCprod(x = x,
                                  strata = design$strata,
                                  psu = design$cluster[[1]],
                                  fpc = design$fpc,
                                  nPSU = design$nPSU,
                                  certainty = design$certainty,
                                  postStrata = design$postStrata)
  varcoef
}
