# INTERNAL : updating a model using new set ofvariable coefficients
#
# @details
# Gradient and Hessian of log-likelihood function as well as variance-covariance matrix are not re-calcualted in this function!
# @param model a fitted \code{hopit} model.
# @param newregcoef a new set of variable coefficients. A vector of length lower or equal to the length of \code{model$coef}.
# The coeficients are replaced starting from the first one.
# @param data a data used to fit original model.
# @author Maciej J. Danko
#' @noRd
#' @keywords internal
update.latent <- function(model, newregcoef, data){
  coefnames <- names(model$coef)
  thresh.names <- colnames(model$thresh.mm)
  model$coef[seq_along(newregcoef)] <- newregcoef
  class(model) <- 'hopit'
  names(model$coef) <- coefnames
  colnames(model$thresh.mm) <- thresh.names
  p <- hopit_ExtractParameters(model)
  model$alpha <- hopit_Threshold(p$thresh.lambda, p$thresh.gamma, model)
  model$y_latent_i <- hopit_Latent(p$latent.params, model)
  model$maxobservedlatentrange <-  range(model$y_latent_i)
  model$Ey_i <- factor(colSums(sapply(1L : model$N, function(k) model$alpha[k,]<model$y_latent_i[k])),levels=1L:model$J)
  levels(model$Ey_i) <- levels(model$y_i)
  model$coef.ls <- p
  model$deviance <- -2 * model$LL

  if (!length(model$design)) {
    k <- 2
    model$AIC <- model$deviance + k * (length(model$coef.ls$latent.params)+
                                         length(model$coef.ls$thresh.lambda)+
                                         length(model$coef.ls$thresh.gamma)+model$hasdisp)

  } else model$AIC <- NA

  return(model)
}

#' Bootstrapping hopit model
#'
#' \code{boot_hopit} performs the bootstrap of a function dependent on a fitted model.
#' In each of the bootstrap repetitions, a set of new model coefficients is drawn from the multivariate normal distribution,
#' assuming the originally estimated model coefficients (see \code{\link{coef.hopit}})
#' as a mean and using the model variance-covariance matrix (see \code{\link{vcov.hopit}}).
#' The drawn coefficients are then used to calculate the measure of interest using a function delivered by the \code{func} parameter.
#' @param model a fitted \code{hopit} model.
#' @param data data used to fit the model.
#' @param func a function to be bootstrapped of the form \code{func(model, ...)}.
#' @param nboot a number of bootstrap replicates.
#' @param unlist a logical indicating whether to unlist the boot object.
#' @param boot.only.latent a logical indicating whether to perform the bootstrap on latent variables only.
#' @param robust.vcov see \code{\link{vcov.hopit}}.
#' @param ... other parameters passed to the \code{func}.
#' @importFrom MASS mvrnorm
#' @author Maciej J. Danko
#' @return a list with bootstrapped elements.
#' @export
#' @seealso \code{\link{percentile_CI}}, \code{\link{getLevels}}, \code{\link{getCutPoints}}, \code{\link{latentIndex}}, \code{\link{standardiseCoef}}, \code{\link{hopit}}.
#' @examples
#' \donttest{
#' # DATA
#' data(healthsurvey)
#'
#' # the order of response levels decreases from the best health to
#' # the worst health; hence the hopit() parameter decreasing.levels
#' # is set to TRUE
#' levels(healthsurvey$health)
#'
#' # fit a model
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # Example 1 ---------------------
#' # bootstrapping cut-points
#'
#' # a function to be bootstrapped
#' cutpoints <-  function(model) getCutPoints(model)$cutpoints
#' B <- boot_hopit(model = model1, func = cutpoints, nboot = 100)
#'
#' # calculate lower and upper bounds using the percentile method
#' cutpoints.CI <- percentile_CI(B)
#'
#' # print estimated cutpoints and their confidence intervals
#' cutpoints(model1)
#' cutpoints.CI
#'
#' # Example 2 ---------------------
#' # bootstrapping differences in health levels
#'
#' # a function to be bootstrapped
#' diff_BadHealth <- function(model) {
#'   hl <- getLevels(model = model, formula=~ sex + ageclass, sep=' ')
#'   hl$original[,1] + hl$original[,2] - hl$adjusted[,1]- hl$adjusted[,2]
#' }
#'
#' # estimate the difference
#' est.org <- diff_BadHealth(model = model1)
#'
#' # perform the bootstrap
#' B <- boot_hopit(model = model1, func = diff_BadHealth, nboot = 100)
#'
#' # calculate lower and upper bounds using the percentile method
#' est.CI <- percentile_CI(B)
#'
#' # plot the difference and its (asymmetrical) confidence intervals
#' pmar <- par('mar'); par(mar = c(9.5,pmar[2:4]))
#' m <- max(abs(est.CI))
#' pos <- barplot(est.org, names.arg = names(est.org), las = 3,
#'                ylab = 'Original - Adjusted',
#'                ylim=c(-m, m), density = 20, angle = c(45, -45),
#'                col = c('blue', 'orange'))
#' for (k in seq_along(pos)) lines(c(pos[k,1],pos[k,1]),
#'                                 est.CI[,k], lwd = 2, col = 2)
#' abline(h = 0); box(); par(mar = pmar)
#' }
boot_hopit<-function(model, func, data=model$frame, nboot = 500, unlist = TRUE,
                     boot.only.latent = TRUE, robust.vcov, ...){
  if (missing(robust.vcov)) {
    if (length(model$design)) robust.vcov <- FALSE else robust.vcov <- TRUE
  }
  data <- model$na.action(data)
  if (model$control$transform.latent != 'none')
    data <- transform.data(model$latent.formula, data,
                           model$control$transform.latent)
  if (model$control$transform.thresh != 'none')
    data <- transform.data(model$thresh.formula, data,
                           model$control$transform.thresh)

  VCOV <- vcov.hopit(model, robust.vcov)
  if (boot.only.latent) N <- seq_len(model$parcount[1]) else N <- nrow(VCOV)
  if (length(VCOV) < 2) stop(call. = NULL, hopit_msg(23))
  bootsample <- MASS::mvrnorm(nboot, mu = model$coef[N], Sigma = VCOV[N,N])
  boots <- lapply(seq_len(nboot), function(k)
    func(model = update.latent(model,
                 bootsample[k,N],data = data), ...))
  if (unlist[1]) {
    boots <- sapply(boots,'[')
    class(boots) <- 'hopit.boot'
  } else class(boots) <- c('hopit.boot', 'list')
  boots
}

#' Calculating the confidence intervals of the bootstrapped function using the percentile method
#'
#' Calculate the confidence intervals of the bootstrapped function using the percentile method.
#' @param boot a matrix or a list of vectors with bootstrapped elements. If it is list, then each element of the list is one replication.
#' @param alpha a significance level.
#' @param bounds which bounds to return; one of \code{"both"}, \code{"lo"}, \code{"up"}.
#' @author Maciej J. Danko
#' @seealso \code{\link{boot_hopit}}, \code{\link{getLevels}}, \code{\link{getCutPoints}}, \code{\link{latentIndex}}, \code{\link{standardiseCoef}}, \code{\link{hopit}}.
#' @export
#' @examples
#' # see examples in boot_hopit() function.
percentile_CI <- function(boot, alpha = 0.05, bounds = c('both', 'lo', 'up')){
  bounds <- tolower(bounds[1])
  if (inherits(boot,'list')) boot <- sapply(boot,'[')
  probs <- switch(bounds,
                  up = 1-alpha/2,
                  lo = alpha/2,
                  both = c(alpha/2, 1-alpha/2))

  apply(boot, 1, stats::quantile, probs = probs)
}
