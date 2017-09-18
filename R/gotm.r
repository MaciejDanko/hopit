#' INTERNAL: Converts a vector of an categorical variable into a matrix with dummies in columns
#'
#' @param V a vector of categories.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
Vector2DummyMat<-function(V) sapply(levels(as.factor(V)), function(k) as.factor(V) == k)*1L

#' INTERNAL: Converts a matrix with dummies in columns into categorical vector
#'
#' @param D a matrix of dummies.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
DummyMat2Vector<-function(D) D %*% ((1L : dim(D)[2L]) -1L)

#' INTERNAL: Do cumsum() in each row of the matrix
#'
#' @param mat a matrix.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
cumsum_row<-function(mat) t(apply(as.matrix(mat), 1L, cumsum))

#' INTERNAL: Column path in a matrix
#'
#' For each row of a matrix \code{mat} extracts a value corresponding to a column stored in vector \code{y}.
#' @param mat a matrix
#' @param y a vector with column indexes corresponding to each row in the matrix \code{mat}.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
col_path<-function(mat, y) mapply(function (x, y) mat[x, y], seq_along(y), y)

#' INTERNAL: Calculation of cut-points (threshold)
#'
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_Threshold<-function(thresh.lambda, thresh.gamma, model = NULL, control = list(),
                         N = model$N, J = model$J, offset = model$thresh.offset,
                         mod.mat = model$thresh.mm, method = model$thresh.method,
                         use.cov = !model$thresh.no.cov,
                         extended.output = FALSE){

  control <- do.call("gotm.control", control)
  fn <- control$thresh.fun

  thresh.lambda <- t(as.matrix(thresh.lambda))
  thresh.offset <- matrix(offset, N, J - 1L)

  if (!use.cov){
    Lin.Tresh.mat <- t(matrix(thresh.lambda, J - 1L, N)) +
      thresh.offset
  } else {
    C <- dim(model$thresh.mm)[2L]
    thresh.gamma <- matrix(as.matrix(thresh.gamma), C, J - 1L)
    Lin.Tresh.mat <- t(matrix(thresh.lambda, J - 1L, N)) +
      mod.mat %*% thresh.gamma + thresh.offset
  }

  a <- matrix(NA, N, J + 1L)
  if (method == 'classic') {
    a[,1L] <- -Inf
    a[,2L] <- Lin.Tresh.mat[,1L]
  } else if (method == 'hopit') {
    a[,1L] <- 0L
    a[,2L] <- fn(Lin.Tresh.mat[,1L])
  } else stop('Unknown threshold metod.')
  a[,J + 1L] <- Inf
  b  <-  a
  if (J>=3L) {
    a[,3L : J]  <-  fn(Lin.Tresh.mat[,2L : (J - 1L)])
    b  <-  a
    a[,2L : J]  <-  cumsum_row(a[,2L : J])
  }
  if (extended.output){
    list(a = a, b = b, thresh.lambda = thresh.lambda, thresh.offset = thresh.offset, Lin.Tresh.mat = Lin.Tresh.mat)
  } else {
    a
  }
}

#' INTERNAL: Anticipate initial values of the model
#'
#' @param model \code{gotm} object.
#' @param data data.frame with data used to fit the model.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @importFrom ordinal clm
#' @keywords internal
gotm_GetInitial <- function(model, data = NULL){
  w <- as.numeric(as.vector(model$weights))
  data$w <- w
  fit0 <- clm(formula = model$reg.formula, data = data, weights = w, link = model$link)
  model$reg.start <- coef(fit0)[-(1L : (model$J - 1L))]
  alphas <- coef(fit0)[(1L : (model$J - 1L))]
  lambda1 <- c(alphas[1L], log(diff(alphas))) #classic

  sqf <- function(param){
    cpc <- cumsum(model$parcount[-1L])
    thresh.lambda <- param[1L : cpc[1L]]
    if (model$parcount[3L]) thresh.gamma <- param[(cpc[1L] + 1L) : cpc[2L]] else thresh.gamma <- 0L
    sum((t(gotm_Threshold(thresh.lambda, thresh.gamma, model)[,-1L][,1L : length(alphas)]) - alphas)^2L, na.rm = T)
  }

  sqf2 <- function(param){
    thresh.lambda <- rep(param[1L], model$parcount[2L])
    if (model$parcount[3L]) thresh.gamma <- rep(param[2L], model$parcount[3L]) else thresh.gamma <- 0L
    sum((t(gotm_Threshold(thresh.lambda, thresh.gamma, model)[,-1L][,1L : length(alphas)]) - alphas)^2L, na.rm = T)
  }

  param <- optim(par = rep(mean(lambda1), 2L), fn = sqf2)$par
  param <- c(lambda1, rep(param[2L], model$parcount[3L]))
  oldv <- 0L
  for (j in 1L : 10L){
    smodel <- optim(par = param, fn = sqf)
    if (abs(oldv - smodel$value) < 1e-12) break
    param <- smodel$par
    oldv <- smodel$value
  }

  cpc <- cumsum(model$parcount[-1L])
  thresh.lambda <- param[1L : cpc[1L]]
  if (model$parcount[3L]) thresh.gamma <- param[(cpc[1L] + 1L) : cpc[2L]] else thresh.gamma <- NULL

  model$lambda.start <- thresh.lambda
  model$gamma.start <- thresh.gamma
  model$start <- c(model$reg.start, model$lambda.start, model$gamma.start)
  model
}

#' INTERNAL: Calculate a latent variable
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_Latent <- function(reg.params, model = NULL, mod.mat = model$reg.mm, N = model$N, offset = model$reg.offset)
  mod.mat %*% (as.matrix(reg.params)) +
  matrix(offset, N, 1L)

#' INTERNAL: Extract model parameters in a form of list
#'
#' Extract model parameters in a form of a list
#' @param model \code{gotm} object
#' @param parameters model parameters (optional). If not delivered then taken from \code{model$coef}
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_ExtractParameters <- function(model, parameters){
  if (!length(model$parcount)) stop('Missing parcount in model object.')
  if (missing(parameters)) {
    parameters <- model$coef
    if (!length(parameters)) stop('Missing estimated parameters.')
  }
  if (length(parameters) != sum(model$parcount)) stop('Wrong number of parameters.')

  reg.params <- parameters[1L : model$parcount[1L]]
  cpc <- cumsum(model$parcount)

  if (model$parcount[2L]) {
    thresh.lambda <- parameters[(cpc[1L] + 1L) : cpc[2L]]
  } else {
    stop('Lambda must be given.')
  }

  if (model$parcount[3L]) {
    thresh.gamma <- parameters[(cpc[2L] + 1L) : cpc[3L]]
  } else {
    thresh.gamma <- NULL
  }
  list(reg.params = reg.params, thresh.lambda = thresh.lambda, thresh.gamma = thresh.gamma)
}


#' INTERNAL: The log likelihood function
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_negLL <- function(parameters, model, collapse = TRUE){
  p <- gotm_ExtractParameters(model, parameters)
  a <- gotm_Threshold(thresh.lambda = p$thresh.lambda, thresh.gamma = p$thresh.gamma,
                      model = model)
  b <- gotm_Latent(p$reg.params, model)
  y <- as.numeric(unclass(model$y_i))
  A2 <- pmax(col_path(a, y) - b, -20L)
  A1 <- pmin(col_path(a, y + 1) - b, 20L)
  P <- model$link.func(A1)-model$link.func(A2)
  cond <- all(P > 0L)
  if(collapse) {
    if (cond) -sum(model$weights * log(P)) else Inf
  } else -log(P)
}

#' INTERNAL: Fit \code{gotm}
#'
#' @param model \code{gotm} object
#' @param start starting parameters
#' @param control a list with control parameters. See \code{\link{gotm.control}}.
#' @importFrom maxLik maxNR
#' @importFrom numDeriv grad
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_fitter <- function(model, start = model$start, control = list()){

  control <- do.call("gotm.control", control)

  gotm_derivLL <- function(parameters = model$coef, model) grad(gotm_negLL, x = parameters, model = model) #numDeriv::

  refit <- function(fit, model, control){
    oldfit <- fit$value
    for (k in 1L : control$max.reiter) {
      try({fit <- optim(par = fit$par, fn = gotm_negLL, gr = gotm_derivLL, model = model, method = 'BFGS')}, silent = TRUE)
      fit <- optim(par = fit$par, fn = gotm_negLL, gr = gotm_derivLL, model = model)
      if (abs(fit$value - oldfit)<control$tol.reiter) break
      oldfit <- fit$value
      cat('.')
    }
    cat('\n')
    list(fit = fit, converged = (k<control$max.reiter))
  }

  z <- try({
    nmfit <- optim(par = start, fn = gotm_negLL, model = model)
    tmp <- refit(nmfit, model, control)
    nmFit <- tmp$fit
    fit <- nmFit
    if(!tmp$converged) {
      message('Convergence has not been reached yet, changing the method ...')
      control$fit.NR <- TRUE
    }
  }, silent = FALSE)
  if (class(z) == "try-error") stop('Impossible to find initial values.')
  if (control$fit.NR){
    nrFit <- maxNR(fn = function(p, model) -gotm_negLL(p, model), start = fit$par, model = model) #maxLik::
    nrFit$par <- nrFit$estimate
    nrFit$value <- -nrFit$maximum
    nrFit <- refit(nrFit, model, control)
    if (!nrFit$converged)
      warning(call. = FALSE,
              'The model probably did not converge. Try to increase "max.reiter" and/or set "forced.DEoptim" == TRUE.')
    nrFit <- nrFit$fit
    fit <- nrFit
  } else nrFit <- NULL

  params <- rbind(nmFit$par, nrFit$par)
  test <- c(nm = nmFit$value, nr = nrFit$value)
  ind <- which.min(test)
  Popt <- params[ind,]
  Vmax <- test[ind]
  model$coef <- Popt
  model$LL <- unname(-Vmax)
  model
}

#' Auxiliary for controlling \code{gotm} fitting
#'
#' @description
#' Auxiliary function for controlling \code{gotm} fitting. Use this function to set control
#' parameters of the \code{\link{gotm}} and other related functions.
#'
#' @param fit.NR logical indicating if to use Newton-Raphson algorithm after main optimization.
#' The Newton-Raphson algorithm is used anyway if the convergence for standard optimization was not found.
#' @param forced.DEoptim logical indicating if use DEoptim to find initial parameters.
#' @param max.reiter maximum number of repeats of the standard optimization procedure if optimum was not found.
#' @param tol.reiter the maximal tolerated difference between log-likelihoods of two
#' consequtive runs of standard optimization.
#' @param grad.eps epsilon for gradient function.
#' @param thresh.fun one parameter function used to calculate thresholds.
#' It can be either a user defined function or character string.
#' The default value is \code{'exp'}.
#' Other accepable character strings include \code{'identity' (or 'id')}, \code{'abs'}, and \code{'sqr' (or '^2')}.
#' @seealso \code{\link{gotm}}
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @export
gotm.control<-function(fit.NR = FALSE,
                       forced.DEoptim = FALSE,
                       max.reiter = 75L,
                       tol.reiter = 1e-8,
                       grad.eps = 1e-8,
                       thresh.fun = 'exp'){
  if (class(thresh.fun) == 'character') {
    if (tolower(thresh.fun) == 'exp') {
      thresh.fun <- exp
    } else if (tolower(thresh.fun) %in% c('identity', 'id')) {
      thresh.fun <- identity
    } else if (tolower(thresh.fun) == 'abs') {
      thresh.fun <- abs
    } else if (tolower(thresh.fun) %in% c('sqr', '^2')) {
      thresh.fun <- function(x) x^2L
    } else stop('Unknown "thresh.fun".')
  }
  list(fit.NR = fit.NR, thresh.fun = thresh.fun, forced.DEoptim = forced.DEoptim,
       max.reiter = max.reiter, tol.reiter = tol.reiter, grad.eps = grad.eps)
}

#frequency weight potraktowac jak w  clm czyli bez PSU
#' Auxiliary for setting a simple survey design for \code{gotm}
#'
#' @param PWeights Probability weights (the inverse of \code{FWeights})
#' @param FWeights Probability of an observation being selected into the sample.
#' Either \code{PWeight} or \code{FWeight} can be delivered (but not both simultaneously)
#' @param PSU Identificator of the PSU unit. Each P- and F- weight should correspond to exactly one PSU.
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
gotm.design<-function(PWeights = NULL, FWeights = NULL, PSU = NULL){
  tmp <- list(PWeights = PWeights, FWeights = FWeights, PSU = PSU)
  if (sum(sapply(tmp, length))){
    if(!length(PSU)) stop('"PSU" must be given.')
    if (!(length(PWeights) + length(FWeights))) stop('P- or F- weights must be given')
    if ((length(PWeights) && length(FWeights))) stop('Please deliver either "Pweights" or "FWeights".')
    if (length(PWeights) && (length(PSU) != length(PWeights))) stop('"PWeights" and "PSU" must have the same length.')
    if (length(FWeights) && (length(PSU) != length(FWeights))) stop('"FWeights" and "PSU" must have the same length.')
  }
  tmp
}

#' Fit Generelaized Ordered Choice Threshold Model
#'
#' @param reg.formula formula used to model latent process.
#' @param thresh.formula formula used to model threshold variable.
#' Any dependent variable (left side of "~") will be ignored.
#' @param data a data frame including all modeled variables.
#' @param reg.offset an offset for the latent variable. A vector with length equals the number of observations.
#' @param thresh.offset an offset for the threshold variable. A matrix of the size N x (J-1),
#' where N is the number of observations and J is the number of categorical responses.
#' @param survey an optional survey a survey design. Empty list indicates no survey design. See \code{\link{gotm.design}}.
#' @param link the link function. The possible values are \code{"probit"} (default) and \code{"logit"}.
#' @param lambda.est.method
#' \itemize{ an assumption about \code{lambda} parameter:
#' \item{\code{"multi"} - each of J-1 modeled thresholds depends on different \code{lambda}. Default option.}
#' \item{\code{"single"} - \code{lambda} is the same for each threshold.}
#' }
#' @param gamma.est.method
#' \itemize{ an assumption about \code{gamma} parameter:
#' \item{\code{"multi"} - each of J-1 modeled thresholds depends on different \code{gamma}. Default option.}
#' \item{\code{"single"} - \code{gamma} is the same for each threshold.}
#' }
#' @param thresh.method
#' \itemize{ Define how the threshold is calculated:
#' \item{\code{"classic"} - Assumes that: \code{alpha(0) = -Inf}, \code{alpha(1)} is a
#' linear function of threshold covariates + \code{lambda}, \code{alpha(2..J-1) = thresh.fun(lambda, gamma, threshold covariates)},
#' and \code{alpha(J) = Inf}. See also [1]. Default option.}
#' \item{\code{"hopit"} - Assumes that: \code{alpha(0) = 0},  \code{alpha(1..J-1) = thresh.fun(lambda, gamma, threshold covarites)},
#' and \code{alpha(J) = Inf}. See also [2].}
#' }
#' @param start starting values in the form \code{c(latent_parameters, threshold_lambdas, threshold_gammas)}
#' @param doFit logical indicating if to run the optimizationa algorithm. If \code{FALSE} starting values must be delivered.
#' @param control a list with control parameters. See \code{\link{gotm.control}}.
#' @importFrom DEoptim DEoptim
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
gotm<- function(reg.formula,
                thresh.formula = as.formula('~1'),
                data = NULL,
                reg.offset = NULL,
                thresh.offset = NULL,
                survey = list(),
                link = c('probit', 'logit'),
                lambda.est.method = c('multi', 'single'),
                gamma.est.method = c('multi', 'single'),
                thresh.method = c('classic', 'hopit'),
                start = NULL,
                doFit = TRUE,
                control = list()){

  control <- do.call("gotm.control", control)
  survey <- do.call("gotm.design", survey)

  model <- NULL
  model$control <- control

  lambda.est.method <- tolower(lambda.est.method[1L])
  gamma.est.method <- tolower(gamma.est.method[1L])
  thresh.method <- tolower(thresh.method[1L])

  link <- link[1L]
  model$link <- link
  if (tolower(link) == 'probit') {
    model$link.func <- pnorm
    model$distr.func <- dnorm
  } else if (tolower(link) == 'logit'){
    model$link.func <- function(x) exp(x)/(1L + exp(x))
    model$distr.func <- function (x) exp(-x)/((1L + exp(-x))^2L)
  } else stop('Unknown link function.')
  if (length(thresh.formula)>2L){
    warning(call. = FALSE, 'The treshold formula should be in the form "~ threshold variables ( +/- 1)".')
    thresh.formula[[2]] <- NULL
  }

  reg.formula <- update.formula(reg.formula, '~.+1')
  model$reg.formula <- reg.formula
  model$reg.mm <- as.matrix(model.matrix(reg.formula, data = data))
  reg.names <- colnames(model$reg.mm)
  grpi <- grepl('(Intercept)', colnames(model$reg.mm), fixed = TRUE)
  model$reg.mm <- as.matrix(model$reg.mm[,!grpi])
  reg.names <- reg.names[!grpi]
  model$reg.offset <- reg.offset

  thresh.formula <- update.formula(thresh.formula, '~.+1')
  model$thresh.formula <- thresh.formula
  model$thresh.mm <- as.matrix(model.matrix(thresh.formula, data = data))
  thresh.names <- colnames(model$thresh.mm)
  grpi <- grepl('(Intercept)', colnames(model$thresh.mm), fixed = TRUE)
  model$thresh.mm <- as.matrix(model$thresh.mm[,!grpi])
  thresh.names <- thresh.names[!grpi]
  if (any(dim(model$thresh.mm) == 0L)) {
    gamma.est.method <- 'single'
    model$thresh.no.cov <- TRUE
  } else {
    model$thresh.no.cov <- FALSE
  }
  model$thresh.method <- thresh.method
  model$thresh.offset <- thresh.offset
  model$lambda.est.method <- lambda.est.method
  model$gamma.est.method <- gamma.est.method

  if (any(thresh.names %in% reg.names) && (thresh.method == 'classic'))
    warning(call. = FALSE, 'The defined model is probably unidentifiable. Please use "hopit" method.')

  model$y_i <- model.frame(reg.formula, data = data)[,all.vars(reg.formula[[2]])]
  if (!is.factor(model$y_i)) stop('Response must be a factor with ordered levels.')
  model$y_latent_i <- NA# latent
  model$Ey_i <- NA# ordinal classified utput
  model$J <- length(levels(model$y_i))
  model$N <- length(model$y_i)
  if (model$J<3L) stop ('Response must have 3 or more levels.')

  if (!length(model$reg.offset)) model$reg.offset <- 0L
  if (!(length(model$reg.offset) %in% c(1L, model$N))) stop('Wrong "reg.offset" length.')
  model$reg.offset <- matrix(model$reg.offset, model$N, 1L)
  if (!length(model$thresh.offset)) {
    model$thresh.offset <- matrix(0L, model$N, model$J - 1L)
  } else if ((length(model$thresh.offset) == 1L) || (length(model$thresh.offset) == model$N)) {
    model$thresh.offset <- matrix(model$thresh.offset, model$N, model$J - 1L)
  } else if (length(model$thresh.offset) != (model$N*(model$J - 1L))) stop('Wrong "thresh.offset" length.')

  if (sum(sapply(survey, length))){
    model$design <- survey
    if (length(survey$PWeights)) model$weights <- 1L/survey$PWeights else model$weights <- survey$FWeights
    if (length(model$weights) != model$N) {
      stop('Vector of survey weights must be of the same length as data.')
    }
  } else {
    model$design <- NULL
    model$weights <- rep(1L, model$N)
  }

  Cr <- dim(model$reg.mm)[2L]
  Ct <- dim(model$thresh.mm)[2L]
  if (model$thresh.no.cov) Ct <- 0L
  metho <- paste(lambda.est.method, gamma.est.method, sep = '')
  model$parcount = switch (metho,
                           multimulti = c(Cr, model$J - 1L, Ct*(model$J - 1L)),
                           multisingle = c(Cr, model$J - 1L, Ct),
                           singlemulti = c(Cr, 1L, Ct*(model$J - 1L)),
                           singlesingle = c(Cr, 1L, Ct))
  model$parcount[3L] <- model$parcount[3L]*(model$thresh.no.cov == FALSE)
  interce <- paste(1L : (model$J - 1L), 2L : (model$J), sep = '|')
  if (model$thresh.no.cov){
    tmp <- NULL
    tmp2 <- NULL
  } else {
    tmp <- as.matrix(expand.grid(thresh.names, interce, KEEP.OUT.ATTRS = FALSE))
    tmp <- paste('(G)', apply(tmp, 1L, paste, sep = '', collapse = '.'), sep = '.')
    tmp2 <- paste('(G)', thresh.names, sep = '.')
  }

  coefnames <- switch (metho,
                       multimulti = c(reg.names, paste('(L)', interce, sep = '.'), tmp),
                       multisingle = c(reg.names, paste('(L)', interce, sep = '.'), tmp2),
                       singlemulti = c(reg.names, '(L)', tmp),
                       singlesingle = c(reg.names, '(L)', tmp2))

  model$weights <- as.vector(matrix(model$weights, 1L, model$N))

  if (!length(start)) {
    if (!doFit) stop('Starting values must be given for "doFit" == TRUE.')
    if (!control$forced.DEoptim) z <- try({gotm_GetInitial(model, data)}, silent = FALSE) else z <- NULL
    if (control$forced.DEoptim||(class(z) == "try-error")) {
      message('Initial values failed, using DEoptim to find them.')
      DEfit <- DEoptim(fn = gotm_negLL, #DEoptim::
                       lower = rep(-10L, sum(model$parcount)),
                       upper = rep(10L, sum(model$parcount)),
                       model = model, control = list(itermax = 500L, trace = TRUE))
      model$start <- DEfit$optim$bestmem
    } else model$start <- z$start
  } else model$start <- start

  if (doFit){
    model <- gotm_fitter(model, start = model$start, model$control)
  } else {
    model$coef <- model$start
    model$LL <- gotm_negLL(parameters = model$start, model)
  }
  names(model$coef) <- coefnames
  colnames(model$thresh.mm) <- thresh.names
  p <- gotm_ExtractParameters(model)
  model$alpha <- gotm_Threshold(p$thresh.lambda, p$thresh.gamma, model)
  model$y_latent_i <- gotm_Latent(p$reg.params, model)
  model$Ey_i <- levels(model$y_i)[colSums(sapply(1L : model$N, function(k) model$alpha[k,]<model$y_latent_i[k]))]
  class(model) <- 'gotm'
  return(model)
}

#' Calculate health index from the fitted \code{gotm} object
#'
#' @param model fitted \code{gotm} object.
HealthIndex.gotm <- function(model) {
  H <- model$y_latent_i
  (H - min(H)) / (max(H) - min(H))
}

#' Extracting coefficients of fitted \code{gotm} object
#'
#' @param object \code{gotm} object.
#' @param standardized logical indicating if to standardize the coefficients to get disability weights.
#' @param aslist logical indicating if model coefficients should be returned as a list of three vectors
#' related to latent variable, threshold lambdas, and threshold gammas.
#' @export
#' @keywords internal
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @aliases coefficients.gotm
coef.gotm <- function(object, standardized = FALSE, aslist = FALSE, ...){
  #stand only for probit
  params <- object$coef
  if (standardized) params <- params/HealthIndex.gotm(object)
  if (aslist) gotm_ExtractParameters(object, params) else  params
}

#' Printing basic information about fitted gotm
#'
#' @param x \code{gotm} object.
#' @export
#' @keywords internal
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
print.gotm<-function(x, ...){
  p <- gotm_ExtractParameters(x)
  cat("Formula (latent variables):", deparse(x$reg.formula), fill = TRUE)
  cat("Formula (threshold variables):", deparse(x$thresh.formula), fill = TRUE)
  cat('Link:', x$link, fill = TRUE)
  cat('Number of cases:', x$N, fill = TRUE)
  cat('Response levels:', toString(levels(x$y_i)), fill = TRUE)
  cat('Number of thresholds:', x$J - 1, fill = TRUE)
  cat('Method of threshold calculation: "', x$thresh.method, '"\n', fill = TRUE, sep = '')
  if (x$lambda.est.method == "multiple") cat('Assumption: Lambdas differ for different thresholds\n')
  if (x$lambda.est.method == "single") cat('Assumption: Lambdas are the same for each thresholds\n')
  if(length(p$thresh.gamma)){
    if (x$gamma.est.method == "multiple") cat('Assumption: Gammas differ between for different thresholds\n')
    if (x$gamma.est.method == "single") cat('Assumption: Gammas are the same for each thresholds\n')
  } else cat('Assumption: Thresholds are independent of covariates\n')
  cat('\nCoefficients (latent variables):\n')
  print(p$reg.params)
  cat('\nThreshold coefficents (Lambda):\n')
  print(p$thresh.lambda)
  if(length(p$thresh.gamma)){
    cat('\nThreshold coefficients (Gamma):\n')
    print(p$thresh.gamma)
  }
  invisible(NULL)
}

#' Extracting variance-covariance matrix from the fitted gotm
#'
#' @param object \code{gotm} object.
#' @param robust.vcov logical indicating if to use sandwich estimator to calculate variance-covariance matrix.
#' If survey deign is detected than this option is ignored.
#' @param control a list with control parameters. See \code{\link{gotm.control}}.
#' @importFrom numDeriv hessian
#' @importFrom survey svyrecvar
#' @export
#' @keywords internal
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
vcov.gotm<-function(object, robust.vcov, control = list(), ...){
  my.grad <- function(fn, par, eps, ...){
    sapply(1L : length(par), function(k){
      epsi <- rep(0L, length(par))
      epsi[k] <- eps
      (fn(par + epsi, ...) - fn(par - epsi, ...))/2/eps
    })
  }
  control <- do.call("gotm.control", control)
  hes <- hessian(gotm_negLL, object$coef, model = object) #numDeriv::
  z <- try(solve(hes), silent = T)
  if (class(z) == 'try-error') {
    z <- NA*hes
    warning(call. = FALSE, 'Model is probably unidentifiable, vcov cannot be computed. Please try to use a "hopit" model.')
  }
  if (length(object$design)){
    if (!missing(robust.vcov)) {
      warning('"robust.vcov" ignored, survey design was detected.')
      robust.vcov <- NA
    }
    gra <- my.grad(fn = gotm_negLL, par = object$coef, eps = control$grad.eps, model = object, collapse = FALSE)
    estfun <- gra # * object$weights : these weights are already included in gotm_negLL
    z <- svyrecvar(estfun %*% z, #survey::
                   data.frame(PSU = object$design$PSU),
                   data.frame(rep(1, object$N)),
                   list(popsize = NULL,
                        sampsize = matrix(length(unique(object$design$PSU)), object$N, 1L)))
  } else {
    if (missing(robust.vcov)) robust.vcov <- FALSE
    if(!robust.vcov){
    } else {
      gra <- my.grad(fn = gotm_negLL, par = object$coef, model = object, collapse = FALSE)
      z <- abs(z %*% crossprod(as.matrix(gra)) %*% z)
    }
  }
  attr(z, 'survey.design') <- (length(object$design) > 0L)
  attr(z, 'robust.vcov') <- robust.vcov
  class(z) <- 'vcov.gotm'
  z
}

#' Print object calculated by \code{\link{vcov.gotm}}
#'
#' @param x \code{gotm} object
#' @keywords internal
#' @export
print.vcov.gotm <- function(x, digits = 3L, ...){
  cat('Variance-covariance matrix:\n')
  print(round(x, digits))
  if (attr(x, 'survey.design')) cat('\nVariance-covariance matrix adjusted for survey design.\n')
  if (!is.na(attr(x, 'robust.vcov')) && attr(x, 'robust.vcov')) cat('\nVariance-covariance matrix based on sandwich estimmator.\n')
  invisible(NULL)
}

#' Calcualte model summary
#'
#' @param object \code{gotm} object.
#' @param robust.se logical indicating if to use robust standard errors based on the sandwich estimator.
#' If survey deign is detected than this option is ignored.
#' @param control a list with control parameters. See \code{\link{gotm.control}}.
#' @export
#' @keywords internal
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
summary.gotm <- function(object, robust.se, control = list(), ...){
  control <- do.call("gotm.control", control)
  varcov <- vcov.gotm(object, robust.se, control)
  SE <- suppressWarnings(sqrt(diag(varcov)))
  if (length(object$design)){
    cat('Survey weights detected. Standard errors was adjusted for survey design.\n')
  }
  if ((!robust.se) && (any(is.na(SE))) && !(length(object$design)))
    warning(call. = FALSE, 'Problem with some standard errors, please try option "robust.se" == TRUE, nd consider to use the "hopit" model.')
  testse <- abs(SE/object$coef)
  testse <- testse[!is.na(testse)]
  if (any(testse > 50L)) warning(call. = FALSE, 'Huge standard errors may suggest a problem with object identifiability. Please try tu use "hopit" model.')

  tstat <-  object$coef/SE
  pvalue <- pnorm(-abs(tstat))  * 2L
  table1 <- data.frame(Estimate = object$coef, 'Std. Error' = SE, 'z value' = tstat, 'Pr(>|z|)' = pvalue, check.names = FALSE)
  tmp <- list(table = table1, vcov = varcov, model = object, robust.se = robust.se)
  class(tmp) <- 'summary.gotm'
  tmp
}

#' Print object calculated by \code{\link{summary.gotm}}
#'
#' @param x object created with \code{\link{summary.gotm}}
#' @keywords internal
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
print.summary.gotm <- function(x, ...){
  model <- x$model
  p <- gotm_ExtractParameters(model)
  cat("Formula (latent variables):", deparse(model$reg.formula), fill = TRUE)
  cat("Formula (threshold variables):", deparse(model$thresh.formula), fill = TRUE)
  cat('\nLink:', model$link, fill = TRUE)
  cat('Number of cases:', model$N, fill = TRUE)
  cat('Response levels:', toString(levels(model$y_i)), fill = TRUE)
  cat('Number of thresholds:', model$J - 1, fill = TRUE)
  cat('Method of threshold calcualation: "', model$thresh.method, '"\n', fill = TRUE, sep = '')
  if (model$lambda.est.method == "multiple") cat('Assumption: (L)ambdas differ for different thresholds\n')
  if (model$lambda.est.method == "single") cat('Assumption: (L)ambdas are the same for each thresholds\n')
  if(length(p$thresh.gamma)){
    if (model$gamma.est.method == "multiple") cat('Assumption: (G)ammas differ for different thresholds\n')
    if (model$gamma.est.method == "single") cat('Assumption: (G)ammas are the same for each thresholds\n')
  } else cat('Assumption: Thresholds are independent of any covariate\n')
  if(x$robust.se) cat('\nRobust SE were used (sandwich estimator of varcov).\n')
  cat('\n')
  printCoefmat(x = x$table, P.values = TRUE, has.Pvalue = TRUE, digits = 4L, dig.tst = 2L)
  cat('\nLog-likelihood:', model$LL, fill = TRUE)
  cat('AIC:', AIC.gotm(model), fill = TRUE)
  invisible(NULL)
}

#' Extracts log likelihood of the fitted model
#'
#' @param object \code{gotm} object.
#' @keywords internal
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
logLik.gotm<-function(object, ...) object$LL

#' Extracts Akaike Information Criterion from the fitted model
#'
#' @param object \code{gotm} object.
#' @param k	numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @keywords internal
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
AIC.gotm<-function(object, ..., k = 2L) -2L*object$LL + length(object$coef)

#' Likelihood Ratio Tests
#'
#' @keywords internal
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
anova.gotm<-function(object, ..., method = c('sequential', 'with.first'), direction = c('decreasing', 'increasing')){
  method <- tolower(method[1L])
  direction <- tolower(direction[1L])
  if (length(list(object, ...)) > 1L) {
    objects <- list(object, ...)
    tmp <- deparse(substitute(list(object, ...)))
    ob.nam <- gsub(' ', '', strsplit(substring(tmp, 6L, nchar(tmp) - 1L), ',', fixed = TRUE)[[1L]])
  } else {
    if ((class(object) != 'list') || (length(object)<2L)) stop('At least two objects must be listed.')
    objects <- object
    ob.nam <- as.character(seq_along(objects))
  }
  if (length(objects) == 2L){
    if(length(objects[[1L]]$coef)>length(objects[[2L]]$coef)) {
      return(lrt.gotm(objects[[1L]], objects[[2L]]))
    } else {
      return(lrt.gotm(objects[[2L]], objects[[1L]]))
    }
  } else {
    out <- NULL
    rna <- NULL
    if (direction == 'increasing') objects <- objects[length(objects) : 1L] else if (direction != 'decreasing') stop('Unknown direction.')
    for (k in 1L : (length(objects) - 1L)) {
      if (tolower(method) == 'sequential'){
        tmp <- lrt.gotm(objects[[k]], objects[[k + 1L]]) # try models mut be of decreasing complexity, silent = F
        rna <- c(rna, paste(ob.nam[k], 'vs.', ob.nam[k + 1L], sep = ' '))
      } else if (tolower(method) == 'with.first') {
        tmp <- lrt.gotm(objects[[1L]], objects[[k + 1L]]) # the first model must be the most complex,  silent = F
        rna <- c(rna, paste(ob.nam[1L], 'vs', ob.nam[k + 1L], sep = ''))
      } else stop('Unknown method.')
      out <- rbind(out, c('Chi^2' = tmp$chisq, df = tmp$df, 'Pr(>Chi^2)' = tmp$pval))
    }
    rownames(out) <- rna
    if (direction == 'increasing') out <- out[dim(out)[1L] : 1L,]
  }
  out <- list(table = out, objets = objects, names = ob.nam, method = method)
  class(out) <- 'anova.gotm'
  out
}

#' Print object calcuated by \code{\link{anova.gotm}}
#'
#' @param x object generated by \code{\link{anova.gotm}}
#' @keywords internal
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
print.anova.gotm <- function(x, ...){
  cat('Anova (LRTs):\n')
  cat('Method: "', x$method, '"\n\n', sep = '')
  printCoefmat(x$table, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE, digits = 5L, dig.tst = 3L, tst.ind = 1L)
  invisible(NULL)
}

#' Likelihood ratio test for a pair of models
#'
#' @param full,nested {Models to be compared.}
#' @keywords internal
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
lrt.gotm <- function(full, nested){
  if (length(full$coef) <= length(nested$coef)) stop('The "full" model must have more parameters than the "nested" one.')
  if (full$LL - nested$LL < 0L) warning(call. = FALSE, 'The "nested" model has the higher likelihood than the "full" model. Try to improve the fit of the models.')
  if (ncol(full$reg.mm) < ncol(nested$reg.mm)) {
    cat('full model:\n')
    cat("-- Formula (latent variables):", deparse(full$reg.formula), fill = TRUE)
    cat('\nnested model:\n')
    cat("-- Formula (latent variables):", deparse(nested$reg.formula), fill = TRUE)
    stop('The latent formulas are not nested.')
  }
  if (ncol(full$thresh.mm) < ncol(nested$thresh.mm)) {
    cat('full model:\n')
    cat("-- Formula (threshold variables):", deparse(full$thresh.formula), fill = TRUE)
    cat('\nnested model:\n')
    cat("-- Formula (threshold variables):", deparse(nested$thresh.formula), fill = TRUE)
    stop('The threshold formulas are not nested.')
  }
  if ((ncol(full$reg.mm)) &&  (ncol(nested$reg.mm)))
    if (!(all(colnames(nested$reg.mm) %in% colnames(full$reg.mm)))) warning(call. = FALSE, 'Models use probably different (non-nested) data sets (latent variable formula).')
  if ((ncol(full$thresh.mm)) &&  (ncol(nested$thresh.mm)))
    if (!(all(colnames(nested$thresh.mm) %in% colnames(full$thresh.mm)))) warning(call. = FALSE, 'Models use probably different (non-nested) data sets (latent variable formula).')
  if ((full$lambda.est.method == 'single') && (nested$lambda.est.method == 'multiple')) stop('Threshold models are not nested.')
  if ((full$gamma.est.method == 'single') && (nested$gamma.est.method == 'multiple')) stop('Threshold models are not nested.')
  if (full$thresh.method != nested$thresh.method) stop('Methods of calcultion of thresholds are different.')

  stat <- 2L*(logLik.gotm(full) - logLik.gotm(nested))
  df.diff <- length(full$coef) - length(nested$coef)
  p <- 1L - pchisq(stat, df.diff)
  z <- list(chisq = stat, df = df.diff, pval = p, full = full, nested = nested)
  class(z) <- 'lrt.gotm'
  z
}

#' Print object calculated by \code{\link{lrt.gotm}}
#'
#' @param x object obtained from \code{\link{lrt.gotm}}
#' @keywords internal
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
print.lrt.gotm <- function(x, ...){
  cat('Likelihood ratio test:\n')
  cat('full model:\n')
  cat("-- Formula (latent variables):", deparse(x$full$reg.formula), fill = TRUE)
  cat("-- Formula (threshold variables):", deparse(x$full$thresh.formula), fill = TRUE)
  cat('\nnested model:\n')
  cat("-- Formula (latent variables):", deparse(x$nested$reg.formula), fill = TRUE)
  cat("-- Formula (threshold variables):", deparse(x$nested$thresh.formula), fill = TRUE)
  #uzyc signif
  cat('\nLikelihood ratio test:\n')
  out <- t(as.matrix(c('Chi^2' = x$chisq, df = x$df, 'Pr(>Chi^2)' = x$pval)))
  row.names(out) <- ''
  printCoefmat(out, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE, digits = 5L, dig.tst = 3L, tst.ind = 1L)
  invisible(NULL)
}

#' Model predictions
#'
#' @param object \code{gotm} object.
#' @param type the type of prediction required. The default \code{"link"}
#' is on the scale of linear predictors (latent variable). The alternative \code{"response"}
#' is on the scale of categorical response variable. The \code{"threshold"}
#' gives the latent threshold variables for each observation.
#' @param standardized logical indicating if to use a standardization to calculate disability weight. See [1].
#' @param strata stratification variable used during standardization.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
predict.gotm <- function(object, type = c('link', 'response', 'threshold'),
                         standardized = FALSE, strata = NULL, ...){
  type <- tolower(type[1L])
  if (!(type %in% c('link', 'response', 'threshold'))) stop('Unknown type.')
  H <- switch(type,
              link = object$y_latent_i,
              response = object$Ey_i,
              threshold = object$alpha)

  standardized <- standardized[1L]
  if (standardized && object$link != 'probit') {
    standardized <- FALSE
    warning(call. = FALSE, 'Standardization omitted. It makes sense only for probit models.')
  } else if (standardized && type != 'link') {
    standardized <- FALSE
    warning(call. = FALSE, 'Standardization omitted. It makes sense only for latent variable scale.')
  }

  if (!standardized) {
    if (length(strata) != 0) warning(call. = FALSE, 'The "strata" ignored. It makes only sense for "standardized" == TRUE.')
    return(H)
  } else {
    HealthIndex = function(H) (H - min(H)) / (max(H) - min(H))
    if (!length(strata)) {
      return(HealthIndex(H))
    } else {
      cH <- H*NA
      strata <- as.character(strata)
      U <- unique(strata)
      for(k in seq_along(U)){
        cH[strata == U[k]] <- HealthIndex(H[strata == U[k]])
      }
      return(cH)
    }
  }
}

#' #' Simulation model output
#' #'
#' #' Given a data and model parameters simulate the categorical response.
#' #'
#' #' @export
#' #' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' gotm_simulate <- function(reg.formula,
#'                         thresh.formula = as.formula('~1'),
#'                         J,
#'                         reg.params,
#'                         th.lambda,
#'                         th.gamma,
#'                         offset = 0L,
#'                         fn = exp,
#'                         method = 'classic',
#'                         add.epsilon = TRUE,
#'                         data){
#'   response <- as.character(reg.formula[[2L]])
#'   reg.mat <- model.matrix(reg.formula[-2L], data = data)[,-1L]
#'   thr.mat <- model.matrix(thresh.formula, data = data)[,-1L]
#'   data$y <- gotm_Latent(reg.params = reg.params, mod.mat = reg.mat, N = dim(data)[1L], offset = offset)
#'   if (add.epsilon) data$y <- data$y + rnorm(dim(data)[1L])
#'   thresh.alpha <- gotm_Threshold(thresh.lambda = th.lambda, thresh.gamma = th.gamma, model = NULL,
#'                               N = dim(data)[1L], J = J, offset = offset,
#'                               mod.mat = model.matrix(thresh.formula, data = data)[,-1L], method = method,
#'                               use.cov = length(thr.mat) > 0L, list(thresh.fun = fn))
#'   data <- cbind(data, as.factor(colSums(sapply(1L : dim(data)[1L], function(k) (data$y[k]>thresh.alpha[k,])))))
#'   colnames(data)[dim(data)[2L]] <- response
#'   params <- c(reg.params, th.lambda, th.gamma)
#'   list(data = data, a = thresh.alpha, par = params)
#' }


