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
#' @param y a vector with column indices corresponding to each row in the matrix \code{mat}.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @export
col_path<-function(mat, y) mapply(function (x, y) mat[x, y], seq_along(y), y)

#' @export
unravel <-function(mat, freq)  {
  mat <- cbind(mat, freq)
  FreqInd <- NCOL(mat)
  ls <- apply(mat,1,function(k) ((rep_row(as.matrix(k[-FreqInd]),k[FreqInd]))))
  do.call('rbind', ls)
}  

#' Convert individual data to frequency table of unique combination of dependent and independent variables
#' 
#' @param formula formula indicating, which variables will be used to construct new database.
#' @param data data.frame including all variables listed in formula.
#' @param FreqNam name of the column with frequencies.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @export
data2freq<-function(formula, data, FreqNam='Freq'){
  what <- c(deparse(formula[[2]]),attr(terms(formula),"term.labels"))
  tmp <- data[,which(names(data)%in%what)]
  V <- as.numeric(as.factor(apply(tmp,1,paste,collapse='',sep='')))
  tmp <- cbind(tmp,V)
  tmp <- tmp[order(V),]
  V <- V[order(V)]
  V2 <- sapply(unique(V),function(k) sum(V==k))
  newd <- tmp[match(unique(V), V),]
  newd[,NCOL(newd)] <- V2
  colnames(newd)[NCOL(newd)] <- FreqNam
  newd
}

#' INTERNAL: Calculation of cut-points (threshold)
#'
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_reg_thresh<-function(thresh.lambda,thresh.gamma, model){
  tmp <- t(matrix(thresh.lambda, model$J - 1L, model$N)) 
  for (k in seq_len(NCOL(model$thresh.mm))) tmp <- tmp + 
      model$thresh.mm[, k, drop=FALSE] %*% thresh.gamma[(2:model$J)-1+ (k-1)*(model$J-1)] 
  tmp
}

#' INTERNAL: Calculation of cut-points (threshold)
#'
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_Threshold<-function(thresh.lambda, thresh.gamma, model = NULL, 
                         extended.output = FALSE){

  control <- model$control
  fn <- control$thresh.fun
  J <- model$J
  N <- model$N
  thresh.lambda <- t(as.matrix(thresh.lambda))

  if (model$thresh.no.cov){
    Lin.Tresh.mat <- t(matrix(thresh.lambda, J - 1L, N)) 
  } else {
    # C <- dim(model$thresh.mm)[2L]
    # thresh.gamma <- t(matrix(as.matrix(thresh.gamma),  J - 1L, C))
    # Lin.Tresh.mat <- t(matrix(thresh.lambda, J - 1L, N)) +
    #   model$thresh.mm %*% thresh.gamma 
    Lin.Tresh.mat <- gotm_reg_thresh(thresh.lambda, thresh.gamma, model)
  }

  a <- matrix(NA, N, J + 1L)
  if (model$thresh.method == 'classic') {
    a[,1L] <- if (control$alpha_0 == 'auto') -Inf else control$alpha_0
    a[,2L] <- Lin.Tresh.mat[,1L]
  } else if (model$thresh.method == 'hopit') {
    a[,1L] <- if (control$alpha_0 == 'auto') 0L else control$alpha_0
    a[,2L] <- fn(Lin.Tresh.mat[,1L]) + a[,1L]
  } else stop('Unknown threshold metod.')
  a[,J + 1L] <- Inf
  b  <-  a
  if (J>=3L) {
    a[,3L : J]  <-  fn(Lin.Tresh.mat[,2L : (J - 1L)])
    b  <-  a
    a[,2L : J]  <-  cumsum_row(a[,2L : J])
  }
  if (extended.output){
    list(a = a, b = b, thresh.lambda = thresh.lambda, Lin.Tresh.mat = Lin.Tresh.mat)
  } else {
    a
  }
}

#' INTERNAL: Fit vglm to get parameters of the model
#'
#' @param model \code{gotm} object.
#' @param data data.frame with data used to fit the model.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @importFrom VGAM vglm
#' @keywords internal
get.vglm.start<-function(model, data, start = NULL){
  subs <- function (a, k, z) {a[k] <- z; a}
  reg.formula <- model$reg.formula
  thresh.formula <- model$thresh.formula
  thresh.method <- model$thresh.method
  if (length(thresh.formula)>2) thresh.formula[[2]] <- NULL
  thrf <- deparse(thresh.formula[[2]])
  big.formula <- update(reg.formula, paste('~ ', thrf,' + . + 1'))
  Y <<- Vector2DummyMat(data[,paste(reg.formula[[2]])])
  big.formula[[2]] <- as.name('Y')
  small.formula <- formula(paste('FALSE ~', thrf))
  cat('Running vglm...')
  mv2<-switch(model$link,
              probit = vglm(big.formula, weights = model$weights, data = data, coefstart = start,
                            family = cumulative(parallel = small.formula, link = 'probit')), #direct substitution of link doesn't work
              logit = vglm(big.formula, weights = model$weights, data = data, coefstart = start,
                           family = cumulative(parallel = small.formula, link = 'logit')))
  cat(' done\nRecalculaing parameters...')
  rm(Y, envir = .GlobalEnv)
  mv2.par <- c(coef(mv2)[-(1:sum(model$parcount[2:3]))], coef(mv2)[(1:sum(model$parcount[2:3]))])
  par.ls <- gotm_ExtractParameters(model, mv2.par)
  if (length(par.ls$thresh.gamma))
    alpha.thresh <- gotm_reg_thresh(thresh.lambda = par.ls$thresh.lambda, 
                                    thresh.gamma = par.ls$thresh.gamma, 
                                    model = model)
    #alpha.thresh<-t(matrix(par.ls$thresh.lambda, model$J - 1L, model$N)) + model$thresh.mm %*% par.ls$thresh.gamma
  par.ls$thresh.lambda <- c(par.ls$thresh.lambda[1], diff(par.ls$thresh.lambda))
  ini.gamma<-par.ls$thresh.gamma
  if (length(par.ls$thresh.gamma)) par.ls$thresh.gamma <- c(par.ls$thresh.gamma[1], diff(par.ls$thresh.gamma))  
  if (!identical(deparse(model$control$thresh.fun), deparse(identity))) {
    if (thresh.method == 'classic'){
      par.ls$thresh.lambda <- c(par.ls$thresh.lambda[1],log(par.ls$thresh.lambda[2:length(par.ls$thresh.lambda)]))
      if (length(par.ls$thresh.gamma)) {
        for (k in 2:(model$J-1)){  # probably analitical solution also exists
          fn <- function(x) sum((alpha.thresh[,k]-gotm_Threshold(subs(par.ls$thresh.lambda,k,x[1]),
                                                                 subs(par.ls$thresh.gamma,k,x[2]), model)[,k+1])^2)
          oo<-optim(par=c(par.ls$thresh.lambda[k],par.ls$thresh.gamma[k]),fn=fn)
          oo<-optim(par=c(oo$par[1],oo$par[2]),fn=fn) # one more time
          par.ls$thresh.lambda[k]<-oo$par[1]
          par.ls$thresh.gamma[k]<-oo$par[2]
        }
      }                      
    } else if (thresh.method == 'hopit'){
      par.ls$thresh.lambda <- log(par.ls$thresh.lambda[1:length(par.ls$thresh.lambda)])
      if (length(par.ls$thresh.gamma)) {
        for (k in 1:(model$J-1)){
          fn <- function(x) sum((alpha.thresh[,k]-gotm_Threshold(subs(par.ls$thresh.lambda,k,x[1]),
                                                                 subs(par.ls$thresh.gamma,k,x[2]), model)[,k+1])^2)
          oo<-optim(par=c(par.ls$thresh.lambda[k],par.ls$thresh.gamma[k]),fn=fn)
          par.ls$thresh.lambda[k]<-oo$par[1]
          par.ls$thresh.gamma[k]<-oo$par[2]
        }
      }
    }
  }
  cat(' done\n')
  model$vglm.start <- mv2.par
  model$vglm.summary <- summary(mv2)
  model$reg.start <- par.ls$reg.params
  model$lambda.start <- par.ls$thresh.lambda
  model$gamma.start <-par.ls$thresh.gamma
  model$start <- c(model$reg.start, model$lambda.start, model$gamma.start)
  model
}

#as.matrix(coef(model.simple.exp))
#as.matrix(par.ls$thresh.gamma)
#' INTERNAL: Calculate a latent variable
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_Latent <- function(reg.params, model = NULL) model$reg.mm %*% (as.matrix(reg.params)) 

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
gotm_negLL <- function(parameters, model, collapse = TRUE, include.weights = TRUE){
  p <- gotm_ExtractParameters(model, parameters)
  a <- gotm_Threshold(thresh.lambda = p$thresh.lambda, thresh.gamma = p$thresh.gamma,
                      model = model)
  b <- gotm_Latent(p$reg.params, model)
  y <- as.numeric(unclass(model$y_i))
  A2 <- pmax(col_path(a, y) - b, -20L)
  A1 <- pmin(col_path(a, y + 1L) - b, 20L)
  P <- model$link.func(A1)-model$link.func(A2)
  cond <- all(P > 0L)
  if (include.weights) w <- model$weights else w <- 1L
  if(collapse) {
    if (cond) -sum(w * log(P)) else Inf
  } else -log(P) * w
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
gotm_fitter <- function(model, start = model$start){

  control <- model$control

  gotm_derivLL <- function(parameters = model$coef, model) grad(gotm_negLL, x = parameters, model = model) #numDeriv::

  refit <- function(fit, model){
    oldfit <- fit$value
    for (k in 1L : control$max.reiter) {
      if (k>3) try({fit <- optim(par = fit$par, fn = gotm_negLL, gr = gotm_derivLL, model = model, method = 'BFGS')}, silent = TRUE)
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
    tmp <- refit(nmfit, model)
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
    nrFit <- refit(nrFit, model)
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
#' @param max.reiter maximum number of repeats of the standard optimization procedure if optimum was not found.
#' @param tol.reiter the maximal tolerated difference between log-likelihoods of two
#' consequtive runs of standard optimization.
#' @param grad.eps epsilon for gradient function.
#' @param thresh.fun one parameter function used to calculate thresholds.
#' It can be either a user defined function or character string.
#' The default value is \code{'exp'} for exponential function.
#' Other accepable character strings include \code{'identity'} or \code{'id'} for identity function.
#' @param alpha_0 value for "zero" threshold. If 'auto' then alpha_0 is set specificaly 
#' to the 'classic' or 'hopit' model. See \code{\link{gotm}}.
#' @seealso \code{\link{gotm}}
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @export
gotm.control<-function(fit.NR = FALSE,
                       max.reiter = 75L,
                       tol.reiter = 1e-6,
                       grad.eps = 1e-8,
                       thresh.fun = c('exp','identity','id'),
                       alpha_0 = 'auto'){
  thresh.fun <- match.arg(thresh.fun)
  if (tolower(alpha_0) != 'auto'){
    if (!is.numeric(alpha_0)) stop('"alpha_0" must be a numeric or equal "auto".')
  } else alpha_0 <- tolower(alpha_0)
  if (tolower(thresh.fun) == 'exp') {
    thresh.fun <- exp
  } else if (tolower(thresh.fun) %in% c('identity', 'id')) {
    thresh.fun <- identity
  } 
  list(fit.NR = fit.NR, thresh.fun = thresh.fun, 
       max.reiter = max.reiter, tol.reiter = tol.reiter, 
       grad.eps = grad.eps, alpha_0 = alpha_0)
}

#frequency weight potraktowac jak w  clm czyli bez PSU
#' Auxiliary for setting a simple survey design for \code{gotm}
#'
#' @param PWeights Probability weights (the inverse of Probability of an observation being selected into the sample)
#' @param FWeights Frequency weights. 
# Either \code{PWeight} or \code{FWeight} can be delivered (but not both simultaneously)
#' @param PSU Identificator of the PSU unit. Each P- and F- weight should correspond to exactly one PSU.
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
# if ((length(PWeights) && length(FWeights))) stop('Please deliver either "PWeights" or "FWeights".')
gotm.design<-function(PWeights = NULL, FWeights = NULL, PSU = NULL){
  tmp <- list(PWeights = PWeights, FWeights = FWeights, PSU = PSU)
  if (sum(sapply(tmp, length))){
    if (!length(PSU) && length(PWeights)) stop('"PSU" must be given.')
    if (length(PSU) && length(unique(PSU)) == 1) stop('There is only one "PSU". "PWeights" cannot be used.')
    if (!(length(PWeights) + length(FWeights))) stop('P- or F- weights must be given')
    if (length(PWeights) && (length(PSU) != length(PWeights))) stop('"PWeights" and "PSU" must have the same length.')
  }
  tmp
}

#' @export
'%!in%' <-function(x, table) match(x, table, nomatch = 0L) == 0L
#' @export
'%c%' <-function(x, table) all(match(x, table, nomatch = 0L))
#' @export
'%!c%' <- function(x, table) !all(match(x, table, nomatch = 0L))
#' @export
'%<=>%' <- function(x, y) identical(x,y)

#' @export
rep_row <- function(mat, times) t(matrix(t(mat), NCOL(mat), NROW(mat) * times))

#' Get starting parameters from less or more complicated hierarchical models
#' 
#' @export
get.start.gotm <- function(object, reg.formula, thresh.formula, data, asList = FALSE){
  old.rf <- object$reg.formula
  old.tf <- object$thresh.formula
  if (deparse(object$reg.formula[[2]])!=deparse(reg.formula[[2]])) stop('Models have different dependent variables')
  if (length(thresh.formula)>2L){
    warning(call. = FALSE, 'The treshold formula should be given without dependent variable.')
    thresh.formula[[2]] <- NULL
  }
  old.rt<-attr(terms(old.rf),"term.labels")
  old.tt<-attr(terms(old.tf),"term.labels")
  new.rt<-attr(terms(as.formula(reg.formula)),"term.labels")
  new.tt<-attr(terms(as.formula(thresh.formula)),"term.labels")
  pr <- gotm_ExtractParameters(object) 
  pr.new <-pr
  if ((old.rt %c% new.rt) && (new.rt %!c% old.rt)) {
    reg.mm <- model.matrix(reg.formula,data)[,-1]
    pr.new$reg.params <- rep(0, NCOL(reg.mm))
    old.ind <- which(colnames(reg.mm)%in%colnames(object$reg.mm))
    pr.new$reg.params[old.ind] <- pr$reg.params
    names(pr.new$reg.params)[old.ind] <- names(pr$reg.params) 
    new.ind <- which(colnames(reg.mm)%!in%colnames(object$reg.mm))
    nnam <- colnames(reg.mm)[new.ind]
    names(pr.new$reg.params)[new.ind] <- nnam
  } else if ((new.rt %c% old.rt) && (old.rt %!c% new.rt)) {
    reg.mm <- model.matrix(reg.formula,data)[,-1]
    rm.ind <- which(colnames(object$reg.mm)%!in%colnames(reg.mm))
    tmp <- pr.new$reg.params[rm.ind]
    pr.new$reg.params <- pr.new$reg.params[-rm.ind]
    pr.new$thresh.lambda[1] <- pr.new$thresh.lambda[1] - 
      mean(object$reg.mm[,rm.ind] * rep_row(tmp, NROW(reg.mm))) * length(tmp)
  }
  if ((old.tt %c% new.tt) && (new.tt %!c% old.tt)){
    thresh.mm <- model.matrix(thresh.formula,data)[,-1]
    pr.new$thresh.params <- rep(0, NCOL(thresh.mm))
    old.ind <- which(colnames(thresh.mm)%in%colnames(object$thresh.mm))
    pr.new$thresh.gamma[old.ind] <- pr$thresh.gamma
    names(pr.new$thresh.gamma)[old.ind] <- names(pr$thresh.gamma) 
    new.ind <- which(colnames(thresh.mm)%!in%colnames(object$thresh.mm))
    nnam <- colnames(thresh.mm)[new.ind]
    names(pr.new$thresh.gamma)[new.ind] <- nnam
  } else if ((new.tt %c% old.tt) && (old.tt %!c% new.tt)) {
    thresh.mm <- model.matrix(thresh.formula,data)[,-1]
    rm.ind <- which(colnames(object$thresh.mm)%!in%colnames(thresh.mm))
    tmp <- pr.new$thresh.gamma[rm.ind]
    pr.new$thresh.gamma <- pr.new$thresh.gamma[-rm.ind]
    pr.new$thresh.lambda[1] <- pr.new$thresh.lambda[1] - 
      mean(object$thresh.mm[,rm.ind] * rep_row(tmp, NROW(thresh.mm))) * length(tmp)
  }
  if (asList) return(pr.new) else 
    return(c(pr.new$reg.params, pr.new$thresh.lambda, pr.new$thresh.gamma))
}

#' Fit Generelaized Ordered Choice Threshold Model
#'
#' @param reg.formula formula used to model latent process.
#' @param thresh.formula formula used to model threshold variable.
#' Any dependent variable (left side of "~") will be ignored.
#' @param data a data frame including all modeled variables.
#' @param survey an optional survey a survey design. Empty list indicates no survey design. See \code{\link{gotm.design}}.
#' @param link the link function. The possible values are \code{"probit"} (default) and \code{"logit"}.
#' @param thresh.method
#' \itemize{ Define how the threshold is calculated:
#' \item{\code{"classic"} - Assumes that: \code{alpha(0) = -Inf}, \code{alpha(1)} is a
#' linear function of threshold covariates + \code{lambda}, \code{alpha(2..J-1) = thresh.fun(lambda, gamma, threshold covariates)},
#' and \code{alpha(J) = Inf}. See also [1]. Default option.}
#' \item{\code{"hopit"} - Assumes that: \code{alpha(0) = 0},  \code{alpha(1..J-1) = thresh.fun(lambda, gamma, threshold covarites)},
#' and \code{alpha(J) = Inf}. See also [2].}
#' }
#' @param start starting values in the form \code{c(latent_parameters, threshold_lambdas, threshold_gammas)}
#' @param doFit character string, \code{'full'} perform the full fit, \code{'vglm'} use starting values obteined from vglm, 
#' but will not improve the fit, and \code{'no'} use external starting values, which must be delivered.
#' @param control a list with control parameters. See \code{\link{gotm.control}}.
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
gotm<- function(reg.formula,
                thresh.formula = as.formula('~ 1'),
                data,
                survey = list(),
                link = c('probit', 'logit'),
                thresh.method = c('classic', 'hopit'),
                start = NULL,
                doFit = c('full','vglm','no'),
                control = list()){
  
  my.grad <- function(fn, par, eps, ...){
    sapply(1L : length(par), function(k){
      epsi <- rep(0L, length(par))
      epsi[k] <- eps
      (fn(par + epsi, ...) - fn(par - epsi, ...))/2/eps
    })
  }
  
  if (missing(data)) data <- environment(reg.formula)
  
  doFit <- match.arg(doFit)
  thresh.method <- match.arg(thresh.method)
  link <- match.arg(link)
  control <- do.call("gotm.control", control)
  survey <- do.call("gotm.design", survey)

  if (length(start) && class(start) == 'gotm'){
    if ((thresh.method != start$thresh.method) || 
        (link != start$link)) {
      warning ('Model in "start" is not compatible and will be not used.')
      start <- NULL
    } else {
      tmp <- deparse(substitute(start))
      start <- get.start.gotm(object = start, reg.formula = reg.formula, 
                         thresh.formula = thresh.formula, 
                         data = data, asList = FALSE)
      cat('Model "',tmp,'" was used to get starting values.\n',sep='')
    } 
  } else if (length(start) && !is.double(start)) stop('Wrong format of "start".')
  
  model <- NULL
  model$control <- control
  model$link <- link
  
  if (tolower(link) == 'probit') {
    model$link.func <- pnorm
    model$distr.func <- dnorm
  } else if (tolower(link) == 'logit'){
    model$link.func <- function(x) exp(x)/(1L + exp(x))
    model$distr.func <- function (x) exp(-x)/((1L + exp(-x))^2L)
  } 
  
  if (length(thresh.formula)>2L){
    warning(call. = FALSE, 'The treshold formula should be given without dependent variable.')
    thresh.formula[[2]] <- NULL
  }

  reg.formula <- update.formula(reg.formula, '~.+1')
  if (any(grepl('offset(',as.character(reg.formula[[3]]),fixed=TRUE))) stop('Offset not supported.')
  model$reg.formula <- reg.formula
  model$reg.mm <- as.matrix(model.matrix(reg.formula, data = data))
  model$reg.lev<-lapply(model.frame(model$reg.formula, data = data), function(k) if (is.factor(k)) as.matrix(table(k)) else 'Not a facor')
  reg.names <- colnames(model$reg.mm)
  grpi <- grepl('(Intercept)', colnames(model$reg.mm), fixed = TRUE)
  model$reg.mm <- as.matrix(model$reg.mm[,!grpi])
  reg.names <- reg.names[!grpi]

  thresh.formula <- update.formula(thresh.formula, '~.+1')
  if (any(grepl('offset(',as.character(thresh.formula[[2]]),fixed=TRUE))) stop('Offset not supprted.')
  model$thresh.formula <- thresh.formula
  model$thresh.mm <- as.matrix(model.matrix(thresh.formula, data = data))
  model$thresh.lev <- lapply(model.frame(model$thresh.formula, data = data), function(k) if (is.factor(k)) as.matrix(table(k)) else 'Not a facor')
  thresh.names <- colnames(model$thresh.mm)
  grpi <- grepl('(Intercept)', colnames(model$thresh.mm), fixed = TRUE)
  model$thresh.mm <- as.matrix(model$thresh.mm[,!grpi])
  thresh.names <- thresh.names[!grpi]
  if (any(dim(model$thresh.mm) == 0L)) {
    model$thresh.no.cov <- TRUE
  } else {
    model$thresh.no.cov <- FALSE
  }
  model$thresh.method <- thresh.method
  model$y_i <- model.frame(reg.formula, data = data)[,all.vars(reg.formula[[2]])]
  if (!is.factor(model$y_i)) stop('Response must be a factor with ordered levels.')
  model$y_latent_i <- NA# latent
  model$Ey_i <- NA# ordinal classified utput
  model$J <- length(levels(model$y_i))
  model$N <- length(model$y_i)
  if (model$J<3L) stop ('Response must have 3 or more levels.')

  if (sum(sapply(survey, length))){
    model$design <- survey
    if (length(survey$PWeights) && (length(survey$FWeights))) {
      model$weights <- survey$FWeights/survey$PWeights 
    } else if (length(survey$PWeights)) {
      model$weights <- 1L/survey$PWeights 
    } else if (length(survey$FWeights)){
      model$weights <- survey$FWeights
    } else stop('This should not happen.')
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
  model$parcount <- c(Cr, model$J - 1L, Ct*(model$J - 1L))
  model$parcount[3L] <- model$parcount[3L]*(model$thresh.no.cov == FALSE)
  interce <- paste(1L : (model$J - 1L), 2L : (model$J), sep = '|')
  if (model$thresh.no.cov){
    tmp <- NULL
  } else {
    tmp <- as.matrix(expand.grid(interce, thresh.names, KEEP.OUT.ATTRS = FALSE))
    tmp <- tmp[,c(2,1)]
    tmp <- paste('(G)', apply(tmp, 1L, paste, sep = '', collapse = '.'), sep = '.')
  }

  coefnames <-  c(reg.names, paste('(L)', interce, sep = '.'), tmp)

  model$weights <- as.vector(matrix(model$weights, 1L, model$N))

  if (!length(start) || (doFit == 'vglm')) {
    if (doFit == 'no') stop('Starting values must be given.')
    z <- try({get.vglm.start(model, data)}, silent = FALSE)
    if (class(z) == "try-error") stop('Initial values failed.')
    model$start <- z$start
  } else model$start <- start

  if (doFit == 'full'){
    cat('Improving fit...')
    model <- gotm_fitter(model, start = model$start)
    cat(' done\n')
  } else {
    model$coef <- model$start
    model$LL <- gotm_negLL(parameters = model$start, model)
  }
  names(model$coef) <- coefnames
  colnames(model$thresh.mm) <- thresh.names
  p <- gotm_ExtractParameters(model)
  model$alpha <- gotm_Threshold(p$thresh.lambda, p$thresh.gamma, model)
  model$y_latent_i <- gotm_Latent(p$reg.params, model)
  model$Ey_i <- factor(colSums(sapply(1L : model$N, function(k) model$alpha[k,]<model$y_latent_i[k])),levels=1L:model$J)
  levels(model$Ey_i) <- levels(model$y_i)
  cat('Calcualting hessaian...')
  hes <- numDeriv::hessian(gotm_negLL, model$coef, model = model) #numDeriv::
  model$vcov <- try(solve(hes), silent = T)
  if (class(z) == 'try-error') {
    z <- NA*hes
    warning(call. = FALSE, 'Model is probably unidentifiable, vcov cannot be computed. Please try to use a "hopit" model.')
  }
  cat(' done\nCalcualting estfun...')
  
  model$estfun <- my.grad(fn = gotm_negLL, par = model$coef,
                          eps = control$grad.eps, model = model, collapse = FALSE)
  cat(' done\n')
  class(model) <- 'gotm'
  return(model)
}

# #' Calculate health index from the fitted \code{gotm} object
# #'
# #' @param model fitted \code{gotm} object.
# HealthIndex.gotm <- function(model) {
#   H <- model$y_latent_i
#   (H - min(H)) / (max(H) - min(H))
# }

#' Extracting coefficients of fitted \code{gotm} object
#'
#' @param object \code{gotm} object.
#' @param standardized logical indicating if to standardize the coefficients to get disability weights.
#' @param aslist logical indicating if model coefficients should be returned as a list of three vectors
#' related to latent variable, threshold lambdas, and threshold gammas.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @keywords internal
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @aliases coefficients.gotm
coef.gotm <- function(object, standardized = FALSE, aslist = FALSE, ...){
  #stand only for probit
  params <- object$coef
  if (standardized) params <- -params/(max(object$y_latent_i) - min(object$y_latent_i))
  if (aslist) gotm_ExtractParameters(object, params) else  params
}

#' Printing basic information about fitted gotm
#'
#' @param x \code{gotm} object.
#' @param ...	further arguments passed to or from other methods.
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

#' #' Extracting residuals from the fitted gotm
#' #'
#' #' @param object \code{gotm} object.
#' #' @param type type of residuals
#' #' @param ...	further arguments passed to or from other methods.
#' #' @aliases residuals
#' #' @keywords internal
#' #' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' #' @export
#' resid.gotm<-function(object, type = c("working", "response"), ...){
#'   warning("resid.gotm() is an experimentaly function, proceed with caution.")
#'   Y <- unclass(object$y_i)
#'   Yhat <- unclass(object$Ey_i)
#'   res <- Y - Yhat
#'   type <- tolower(type[1])
#'   if (!(type %in%  c("working", "response"))) stop('Unknown type.')
#'   if (type == 'working') {
#'     stop('Don't know how to calculate)
#'     #cate <- predict(object, type = 'threshold_link')
#'     Ygrad <- dnorm(object$y_latent_i)
#'     res <- res / Ygrad
#'   }
#'   res
#' }

#' Extracting variance-covariance matrix from the fitted gotm
#'
#' @param object \code{gotm} object.
#' @param robust.vcov logical indicating if to use sandwich estimator to calculate variance-covariance matrix.
#' If survey deign is detected than this option is ignored.
#' @param control a list with control parameters. See \code{\link{gotm.control}}.
# @param robust.method method of calculation of log-likelihood gradient. 
# Set \code{"grad"} (default) for numerical gradient or \code{"working"} 
# for the method based on working residuals. 
# The latter one is an experimental method so warning will apear.
# robust.method = c("grad","working")
#' @param ...	further arguments passed to or from other methods.
#' @importFrom numDeriv hessian
#' @importFrom survey svyrecvar
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
vcov.gotm<-function(object, robust.vcov, control = list(), ...){
  #robust.method <- tolower(robust.method[1])
  #if (!(robust.method %in% c("grad","working"))) stop('Unknown method.')
  control <- do.call("gotm.control", control)
  z <- object$vcov 
  if (length(object$design$PSU)){
    if (!missing(robust.vcov) && (robust.vcov)) {
      warning(call. = FALSE, '"robust.vcov" ignored, survey design was detected.')
      robust.vcov <- NA
    }
    z <- svyrecvar(object$estfun %*% z, #survey::
                   data.frame(PSU = object$design$PSU),
                   data.frame(rep(1, object$N)),
                   list(popsize = NULL,
                        sampsize = matrix(length(unique(object$design$PSU)), object$N, 1L)))
  } else {
    if (missing(robust.vcov)) robust.vcov <- FALSE
    if (robust.vcov) z <- (z %*% t(object$estfun) %*% (object$estfun/object$design$FWeights) %*% (z)) 
  }
  attr(z, 'survey.design') <- (length(object$design) > 0L)
  attr(z, 'robust.vcov') <- robust.vcov
  class(z) <- 'vcov.gotm'
  z
}

#' Print object calculated by \code{\link{vcov.gotm}}
#'
#' @param x \code{gotm} object
#' @param ...	further arguments passed to or from other methods.
#' @keywords internal
#' @export
print.vcov.gotm <- function(x, digits = 3L, ...){
  cat('Variance-covariance matrix:\n')
  print(round(x, digits))
  if (attr(x, 'survey.design')) cat('\nVariance-covariance matrix adjusted for survey design.\n')
  if (!is.na(attr(x, 'robust.vcov')) && attr(x, 'robust.vcov')) cat('\nVariance-covariance matrix based on sandwich estimmator.\n')
  invisible(NULL)
}

#' Calculate model summary
#'
#' @param object \code{gotm} object.
#' @param robust.se logical indicating if to use robust standard errors based on the sandwich estimator.
#' If survey deign is detected than this option is ignored.
#' @param control a list with control parameters. See \code{\link{gotm.control}}.
# @param robust.method method of calculation of log-likelihood gradient. 
# Set \code{"grad"} (default) for numerical gradient or \code{"working"} 
# for the method based on working residuals. 
# The latter one is an experimental method so warning will apear.
# robust.method = c("grad","working")
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
summary.gotm <- function(object, robust.se = FALSE, control = list(), ...){
  control <- do.call("gotm.control", control)
  varcov <- vcov(object, robust.se, control, ...)
  SE <- suppressWarnings(sqrt(diag(varcov)))
  if (length(object$design)){
    cat('Survey weights detected. Standard errors was adjusted for survey design.\n')
  }
  if ((!robust.se) && (any(is.na(SE))) && !(length(object$design)))
    warning(call. = FALSE, 'Problem with some standard errors, please try option "robust.se" == TRUE, nd consider to use the "hopit" model.')
  testse <- abs(SE/object$coef)
  testse <- testse[!is.na(testse)]
  if (any((testse > 50L)&(SE > 20L))) warning(call. = FALSE, 'Huge standard errors may suggest a problem with object identifiability. Please try tu use "hopit" model.')

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
#' @param ...	further arguments passed to or from other methods.
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
#' @param ...	additional objects of the same type.
#' @keywords internal
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
logLik.gotm<-function(object, ...) {
  objects <- list(object, ...)
  tmp <- deparse(substitute(list(object, ...)))
  ob.nam <- gsub(' ', '', strsplit(substring(tmp, 6L, nchar(tmp) - 1L), ',', fixed = TRUE)[[1L]])
  res <- sapply(objects,function(object) object$LL)
  names(res) <- ob.nam
  res
}
#' Extracts Akaike Information Criterion from the fitted model
#'
#' @param object \code{gotm} object.
#' @param k	numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @param ...	additional objects of the same type.
#' @keywords internal
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
AIC.gotm<-function(object, ..., k = 2L) {
  objects <- list(object, ...)
  tmp <- deparse(substitute(list(object, ...)))
  ob.nam <- gsub(' ', '', strsplit(substring(tmp, 6L, nchar(tmp) - 1L), ',', fixed = TRUE)[[1L]])
  res <- sapply(objects,function(object) -2L*object$LL + k * length(object$coef))
  names(res) <- ob.nam
  res
}

#' LRT Tables
#'
#' Compute likelihood rato test for two or more \code{gotm} objecs.
#' @param object an object containing the results returned by a \code{gotm}.
#' @param ...	additional objects of the same type.
#' @param method the method of model comparison. Choose \code{"sequential"} for 1-2, 2-3, 3-4, ... comparisons or
#' \code{"with.first"} for 1-2, 1-3, 1-4, ... comparisons.
#' @param direction determine if complexity of listed models is \code{"increasing"} or \code{"decreasing"} (default).
#' @keywords internal
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
anova.gotm<-function(object, ..., method = c('sequential', 'with.first'), direction = c('decreasing', 'increasing')){
  method <- match.arg(method)
  direction <- match.arg(direction)
  if (length(list(object, ...)) > 1L) {
    objects <- list(object, ...)
    tmp <- deparse(substitute(list(object, ...)))
    ob.nam <- gsub(' ', '', strsplit(substring(tmp, 6L, nchar(tmp) - 1L), ',', fixed = TRUE)[[1L]])
  } else  stop('At least two objects must be listed.')
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
      } else 
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
#' @param ...	further arguments passed to or from other methods.
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
  if (full$LL - nested$LL < -.Machine$double.eps) warning(call. = FALSE, 'The "nested" model has the higher likelihood than the "full" model. Try to improve the fit of the models.')
  if (ncol(full$reg.mm) < ncol(nested$reg.mm)) {
    cat('Full model:\n')
    cat("-- Formula (latent variables):", deparse(full$reg.formula), fill = TRUE)
    cat('\nNested model:\n')
    cat("-- Formula (latent variables):", deparse(nested$reg.formula), fill = TRUE)
    stop('The latent formulas are not nested.')
  }
  if (ncol(full$thresh.mm) < ncol(nested$thresh.mm)) {
    cat('Full model:\n')
    cat("-- Formula (threshold variables):", deparse(full$thresh.formula), fill = TRUE)
    cat('\nNested model:\n')
    cat("-- Formula (threshold variables):", deparse(nested$thresh.formula), fill = TRUE)
    stop('The threshold formulas are not nested.')
  }
  if ((ncol(full$reg.mm)) &&  (ncol(nested$reg.mm)))
    if (!(all(colnames(nested$reg.mm) %in% colnames(full$reg.mm)))) warning(call. = FALSE, 'Models use probably different (non-nested) data sets (latent variable formula).')
  if ((ncol(full$thresh.mm)) &&  (ncol(nested$thresh.mm)))
    if (!(all(colnames(nested$thresh.mm) %in% colnames(full$thresh.mm)))) warning(call. = FALSE, 'Models use probably different (non-nested) data sets (latent variable formula).')
  if (full$thresh.method != nested$thresh.method) stop('Methods of calcultion of thresholds are different.')

  stat <- 2L*(logLik.gotm(full) - logLik.gotm(nested))
  df.diff <- length(full$coef) - length(nested$coef)
  p <- 1L - pchisq(stat, df.diff)
  z <- list(chisq = stat, df = df.diff, pval = p, full = full, nested = nested)
  class(z) <- 'lrt.gotm'
  z
}

'%lrt%' <- lrt.gotm

#' Print object calculated by \code{\link{lrt.gotm}}
#'
#' @param x object obtained from \code{\link{lrt.gotm}}
#' @param ...	further arguments passed to or from other methods.
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
#' gives the thresholds for each observation, whereas the \code{"threshold_link"} gives meaningful thresholds
#' together with latent variable for each observation (a data.frame with fields \code{$left.boundary},
#' \code{$latent.variable}, and \code{$right.boundary}).
#' @param standardized logical indicating if to use a standardization to calculate disability weight. See [1].
#' @param strata stratification variable used during standardization.
#' @param unravelFreq logical indicating if to represent results on individual scale if FWeights were used.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
predict.gotm <- function(object, type = c('link', 'response', 'threshold', 'threshold_link'),
                         standardized = FALSE, strata = NULL, unravelFreq = TRUE, ...){
  
  if (length(object$design$FWeights) && unravelFreq) conv<-function(x) unravel(x,freq=object$design$FWeights) else conv<-identity
  type <- match.arg(type)
  if (type == 'latent') type <- 'link'
  H <- switch(type,
              link = conv(object$y_latent_i),
              response = conv(object$Ey_i),
              threshold = conv(object$alpha),
              threshold_link = data.frame(left.boundary=conv(col_path(object$alpha, unclass(object$Ey_i)+1)),
                                          latent.variable=conv(object$y_latent_i),
                                          right.boundary=conv(col_path(object$alpha, unclass(object$Ey_i)+2))))

  standardized <- standardized[1L]
  if (standardized && (object$link != 'probit')) {
    standardized <- FALSE
    warning(call. = FALSE, 'Standardization omitted. It makes sense only for probit models.')
  } else if (standardized && (type == 'response')) {
    standardized <- FALSE
    warning(call. = FALSE, 'Standardization omitted. It makes sense only for latent variable scale.')
  }

  if (!standardized) {
    if (length(strata) != 0) warning(call. = FALSE, 'The "strata" ignored. It makes only sense for "standardized" == TRUE.')
    return(H)
  } else {
    L <- conv(object$y_latent_i)
    HealthIndex <- function(H, L) (H - min(L)) / (max(L) - min(L))
    if (!length(strata)) {
      if (NCOL(H)>1) cH <- apply(H,2,function(k) HealthIndex(k,L)) else cH <- HealthIndex(H,L)
    } else {
      cH <- H*NA
      strata <- as.character(conv(strata))
      U <- unique(strata)
      for(k in seq_along(U)){
        ind <- strata == U[k]
        if (NCOL(H)>1) {
          cH[ind,] <- apply(H[ind,],2,function(k) HealthIndex(k,L[ind]))
        } else {
          cH[ind] <- HealthIndex(H[ind],L[ind])
        }
      }
    }
    cH[cH < 0] <- 0
    cH[cH > 1] <- 1
    cH
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

factor.mat.gotm<-function(object, by.formula = object$thresh.formula){
  if (length(by.formula)==3) by.formula[[2]] <- NULL
  by.what<-ex.lb(by.formula)
  rm.interaction <- grep(':',by.what,fixed=TRUE)
  if (length(rm.interaction)) by.what <- by.what[-rm.interaction]
  if (any(by.what %!in% c(ex.lb(object$reg.formula),ex.lb(object$thresh.formula)))) stop('Some elemnts of "by.formula" not present in the data.')
  R. <- (object$reg.mm[,unlist(sapply(by.what, grep, x = colnames(object$reg.mm), fixed =TRUE))])
  if(!length(R.)) R. <- NULL
  T. <- object$thresh.mm[,unlist(sapply(by.what, grep, x = colnames(object$thresh.mm), fixed =TRUE))]
  if (!length(T.)) T. <- NULL        
  if (any(colnames(T.) %in% colnames(R.))) T. <- T.[,-which(colnames(T.) %in% colnames(R.))]
  vari <- cbind(T., R.)
  get.reg.ref <- sapply(object$reg.lev,function(k) rownames(k)[1])[-1]
  get.thresh.ref <- sapply(object$thresh.lev,function(k) rownames(k)[1])
  get.ref <- c(get.reg.ref, get.thresh.ref)
  miss.ref <- get.ref[names(get.ref)%in%by.what]
  miss.col <- sapply(names(miss.ref), function(k) 1-rowSums(as.matrix(vari[,grep(k,colnames(vari),fixed=TRUE)])))
  miss.nam <- paste(names(miss.ref),miss.ref,sep='')
  colnames(miss.col) <- miss.nam
  fullmat<-cbind(vari,miss.col)
  fullmat<-fullmat[,order(colnames(fullmat))]
  fullmat
}