#' survey:::htvar.matrix clone
#'
#' @details clone taken from the survey v3.35 package.
#' @keywords internal
#' @param xcheck,Dcheck internal parameters.
#' @author Thomas Lumley
#' @importFrom Matrix crossprod
clone.of.htvar.matrix <- function (xcheck, Dcheck) {
  if (is.null(dim(xcheck)))
    xcheck <- as.matrix(xcheck)
  rval <- apply(xcheck, 2, function(xicheck) apply(xcheck, 2, function(xjcheck)
    as.matrix(Matrix::crossprod(xicheck, Dcheck %*% xjcheck))))
  if (is.null(dim(rval)))
    dim(rval) <- c(1, 1)
  rval
}

#' survey:::ygvar.matrix clone
#'
#' @details clone taken from the survey v3.35 package.
#' @param xcheck,Dcheck internal parameters.
#' @keywords internal
#' @author Thomas Lumley
clone.of.ygvar.matrix <-function (xcheck, Dcheck)
{
  ht <- clone.of.htvar.matrix(xcheck, Dcheck)
  if (is.null(dim(xcheck))) {
    corr <- sum(Dcheck %*% (xcheck * xcheck))
  }
  else {
    corr <- apply(xcheck, 2, function(xicheck)
      apply(xcheck, 2, function(xjcheck) sum(Dcheck %*% (xicheck * xjcheck))))
  }
  rval <- ht - corr
}

#' survey:::ppsvar clone
#'
#' @details clone taken from the survey v3.35 package.
#' @param x,design internal parameters.
#' @keywords internal
#' @author Thomas Lumley
clone.of.ppsvar<-function (x, design)
{
  postStrata <- design$postStrata
  est <- design$variance
  if (!is.null(postStrata)) {
    for (psvar in postStrata) {
      if (inherits(psvar, "greg_calibration")) {
        if (psvar$stage == 0) {
          y <- qr.resid(psvar$qr, y/psvar$w) * psvar$w
        }
        else {
          stop("calibration within clusters not yet available for PPS designs")
        }
      }
      else {
        psw <- attr(psvar, "weights")
        postStrata <- as.factor(psvar)
        psmeans <- rowsum(y/psw, psvar, reorder = TRUE)/
          as.vector(table(factor(psvar)))
        x <- y - psmeans[match(psvar, sort(unique(psvar))),
                         ] * psw
      }
    }
  }
  dcheck <- design$dcheck
  if (length(dcheck) != 1)
    stop("Multistage not implemented yet")
  rval <- switch(est,
                 HT = clone.of.htvar.matrix(
                   rowsum(x, dcheck[[1]]$id, reorder = FALSE),
                   dcheck[[1]]$dcheck),
                 YG = clone.of.ygvar.matrix(
                   rowsum(x, dcheck[[1]]$id, reorder = FALSE),
                   dcheck[[1]]$dcheck),
                 stop("can't happen"))
  rval
}

#' Calculation of variance-covariance matrix for specified survey design (experimental function)
#'
#' @param Ainv a variance-covariance matrix.
#' @param estfun a gradient function of the log-likelihood function.
#' @param design a \code{survey.design} object.
#' @description
#' This is a modification of \code{survey:::svy.varcoef}. In the original approach \code{estfun} is calcualted from
#' glm's working residuals:\cr
#' \code{estfun <- model.matrix(glm.object) * resid(glm.object, "working") * glm.object$weights}\cr
#' In the hopit package estfun is directly calculated as a gradient (vector of partial derivatives) of log likelihood function.
#' @seealso
#' \code{\link[survey]{svydesign}}
#' \code{\link{hopit}}
#' @details Based on the survey v3.35 package.
#' @importFrom survey svyrecvar twophasevar twophase2var svyCprod
#' @author Thomas Lumley, modified by Maciej J. Danko
svy.varcoef_hopit <- function (Ainv, estfun, design) {
  if (inherits(design, "survey.design2"))
    V <- survey::svyrecvar(estfun %*% Ainv, design$cluster, design$strata,
                           design$fpc, postStrata = design$postStrata)
  else if (inherits(design, "twophase"))
    V <- survey::twophasevar(estfun %*% Ainv, design)
  else if (inherits(design, "twophase2"))
    V <- survey::twophase2var(estfun %*% Ainv, design)
  else if (inherits(design, "pps"))
    V <- clone.of.ppsvar(estfun %*% Ainv, design)
  else V <-survey::svyCprod(estfun %*% Ainv,
                            design$strata,
                            design$cluster[[1]],
                            design$fpc,
                            design$nPSU,
                            design$certainty,
                            design$postStrata)
  V
}
