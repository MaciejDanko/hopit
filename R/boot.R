#' INTERNAL : updating model uses new set of latent variable coeficients
#'
#' @param model a fitted \code{hopit} model.
#' @param newregcoef a new set of latent variable coeficients. Vector of length of first element of model$parcount.
#' @param data used to fit original model.
#' @param hessian logical indicating if to calculate estfun, hessian and var-cov matrix.
#' @author Maciej J. Danko
#' @keywords internal
update.latent <- function(model, newregcoef, data, hessian=FALSE){
  coefnames <- names(model$coef)
  thresh.names <- colnames(model$thresh.mm)
  model$coef[seq_len(model$parcount[1])] <- newregcoef
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

  if(hessian) {
    hes <- my.grad(fn = hopit_derivLL, par = model$coef, model=model, eps = 1e-4, collapse = TRUE, negative=FALSE)
    if (model$hasdisp && remove.theta) {
      hes <- hes[-nrow(hes),-ncol(hes)] #remove theta from vcov
      model$coef <- model$coef[-length(model$coef)] #remove from coef
    }
    model$hessian <- hes

    model$vcov.basic <- try(base::solve(-hes), silent = FALSE)
    model$vcov.basic <- check_vcov(model$vcov.basic)

    if (model$hasdisp && remove.theta) COEF <- c(model$coef,model$coef.ls$logTheta) else COEF <- model$coef
    model$estfun <- hopit_derivLL(COEF, model, collapse = FALSE)
    if (remove.theta) model$estfun <- model$estfun[,-ncol(model$estfun)]

    if (length(model$design)) {
      model$vcov <- svy.varcoef.hopit(model$vcov.basic, model$estfun, design)
    } else {
      model$vcov <- model$vcov.basic
    }
  }
  if (!length(model$design)) {
    k <- 2
    model$AIC <- model$deviance + k * (length(model$coef.ls$latent.params)+
                                         length(model$coef.ls$thresh.lambda)+
                                         length(model$coef.ls$thresh.gamma)+model$hasdisp)

  } else model$AIC <- NA

  return(model)
}

#' Bootstraping hopit model
#'
#' \code{boot_hopit} performs bootstrap of a function dependent on fitted model.
#' In each of the bootstrap repetitions a set of new model coefficients is drawn from the multivariate normal distribution,
#' assuming originally estimated model coefficients (see \code{\link{coef.hopit}})
#' as a mean and using model variance-covariance matrix (see \code{\link{vcov.hopit}}).
#' The drawn coefficients are then used to calculate the measure of interest using a function delivered by \code{func} parameter.
#' @param model a fitted \code{Hopit} model.
#' @param data data used to fit the model.
#' @param func function to be bootstrapped of the form \code{func(model, data, ...)}.
#' @param nboot number of bootstrap replicates.
#' @param unlist logical indicting if to unlist boot object.
#' @param boot.only.latent logical indicating if to perform the bootstrap only on latent variables.
#' @param robust.vcov see \code{\link{vcov.hopit}}.
#' @param ... other parameters passed to the \code{func}.
#' @importFrom MASS mvrnorm
#' @author Maciej J. Danko
#' @return a list with bootstraped elements.
#' @export
#' @seealso \code{\link{boot_hopit_CI}}, \code{\link{getLevels}}, \code{\link{getCutPoints}}, \code{\link{latentIndex}}, \code{\link{standardiseCoef}}, \code{\link{hopit}}.
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # the order of response levels is decreasing (from the best health to the worst health)
#' levels(healthsurvey$health)
#'
#' # fitting a model
#' model1 <- hopit(latent.formula = health ~ hypertenssion + high_cholesterol +
#'                 heart_atack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # Example 1 ---------------------
#' # Bootstraping cutpoints
#'
#' # Function to be bootstraped
#' cutpoints <-  function(model, data) getCutPoints(model, plotf = FALSE)$cutpoints
#' B <- boot_hopit(model = model1, data = healthsurvey,
#'                 func = cutpoints, nboot = 100)
#'
#' # Calcualte lower and upper bounds using percentile method
#' cutpoints.CI <- boot_hopit_CI(B)
#'
#' # print estimated cutpoints and their confidence intervals
#' cutpoints(model1, healthsurvey)
#' cutpoints.CI
#'
#' # Example 2 ---------------------
#' # Bootstraping health levels differences
#'
#' # Function to be bootstraped
#' diff_BadHealth <- function(model, data) {
#'   hl <- getLevels(model = model, formula=~ sex + ageclass, data = data,
#'                   sep=' ', plotf=FALSE)
#'   hl$original[,1] + hl$original[,2] - hl$adjusted[,1]- hl$adjusted[,2]
#' }
#'
#' # Estimate of the difference
#' est.org <- diff_BadHealth(model = model1, data = healthsurvey)
#'
#' # Perform the bootstrap
#' B <- boot_hopit(model = model1, data = healthsurvey,
#'                 func = diff_BadHealth, nboot = 100)
#'
#' # Calcualte lower and upper bounds using percentile method
#' est.CI <- boot_hopit_CI(B)
#'
#' # Plotting the difference and its (assymetrical) confidence intervals
#' pmar <- par('mar'); par(mar = c(9.5,pmar[2:4]))
#' m <- max(abs(est.CI))
#' pos <- barplot(est.org, names.arg = names(est.org), las = 3, ylab = 'Orginal - Adjusted',
#'                ylim=c(-m, m), density = 20, angle = c(45, -45), col = c('blue', 'orange'))
#' for (k in seq_along(pos)) lines(c(pos[k,1],pos[k,1]), est.CI[,k], lwd = 2, col = 2)
#' abline(h = 0); box(); par(mar = pmar)
boot_hopit<-function(model, data, func, nboot=500, unlist = TRUE, boot.only.latent = TRUE, robust.vcov = TRUE, ...){
  VCOV <- vcov.hopit(model, robust.vcov)
  if (boot.only.latent) N <- seq_len(model$parcount[1]) else N <- nrow(VCOV)
  if (length(VCOV) < 2) stop(call. = NULL, hopit_msg(23))
  bootsample <- MASS::mvrnorm(nboot, mu = model$coef[N], Sigma = VCOV[N,N])
  boots <- lapply(seq_len(nboot), function(k) func(model = update.latent(model,
                                                   bootsample[k,N],data = data), data = data, ...))
  if (unlist) {
    boots <- sapply(boots,'[')
    class(boots) <- 'hopit.boot'
  } else class(boots) <- c('hopit.boot', 'list')
  boots
}

#' Calculating Confidence Intervals using percentile method
#'
#' @param boot boot object calculated by \code{\link{boot_hopit}} .
#' @param alpha significance level.
#' @param bounds one of \code{"both"}, \code{"lo"}, \code{"up"}.
#' @author Maciej J. Danko
#' @seealso \code{\link{boot_hopit}}, \code{\link{getLevels}}, \code{\link{getCutPoints}}, \code{\link{latentIndex}}, \code{\link{standardiseCoef}}, \code{\link{hopit}}.
#' @export
#' @examples
#' # see examples in boot_hopit() function.
boot_hopit_CI <- function(boot, alpha = 0.05, bounds=c('both','lo','up')){
  if (!inherits(boot,'hopit.boot')) stop(call.=NULL, hopit_msg(22))
  bounds <- tolower(bounds[1])
  if (inherits(boot,'list')) boot <- sapply(boot,'[')
  probs <- switch(bounds,
                         up = 1-alpha/2,
                         lo = alpha/2,
                         both = c(alpha/2, 1-alpha/2))

  apply(boot, 1, quantile, probs = probs)
}
