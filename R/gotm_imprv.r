#' INTERNAL: Calculation of cut-points (threshold)
#'
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib gotm
#' @importFrom Rcpp evalCpp
gotm_Threshold<-function(thresh.lambda, thresh.gamma, model){
  if (model$method == 0) {
    getThresholds(model$thresh.mm, thresh.lambda, thresh.gamma, model$thresh.no.cov,
                  thresh_start = model$control$thresh.start, thresh_1_exp = model$control$thresh.1.exp) #RcppEigen
  } else {
    VgetThresholds(model$thresh.mm, thresh.lambda, thresh.gamma, model$thresh.no.cov,
                   thresh_start = model$control$thresh.start) #RcppEigen
  }
}


#' INTERNAL: Calculate latent variable
#' @author Maciej J. Danko
#' @keywords internal
gotm_Latent <- function(reg.params, model = NULL) model$reg.mm %*% (as.matrix(reg.params))


#' INTERNAL: Calculate maximum possible latent range
#' @param model a fitted \code{gotm} model.
#' @param data a data used to fit the model
#' @keywords internal
gotm_latentrange <- function (model, data) {
  #Function must be checked for variables that are not factors

  cfm <- model$coef[seq_len(model$parcount[1])]
  ttr <- terms.formula(model$reg.formula)
  ttr <- delete.response(ttr)
  tt <- attr(ttr,'variables')
  ttn <- attr(ttr,'term.labels')
  li <- lapply(eval(tt, data), function(k) if (class(k) == 'factor') levels(k) else range(k, na.rm=TRUE))
  names(li) <- ttn
  L=sapply(li, length)-1
  cfm_neg <- cfm * (cfm<0)
  cfm_pos <- cfm * (cfm>0)
  pos=c(0,cumsum(L))+1
  cfm_neg_ls <- lapply(1:(length(pos)-1),function(k) cfm_neg[pos[k]:(pos[k+1]-1)])
  cfm_pos_ls <- lapply(1:(length(pos)-1),function(k) cfm_pos[pos[k]:(pos[k+1]-1)])
  c(sum(sapply(cfm_neg_ls, min, na.rm=TRUE)),
    sum(sapply(cfm_pos_ls, max, na.rm=TRUE)))
}


#' INTERNAL: Extract model parameters in a form of list
#'
#' Extract model parameters in a form of a list
#' @param model \code{gotm} object
#' @param parameters model parameters (optional). If not delivered then taken from \code{model$coef}
#' @author Maciej J. Danko
#' @keywords internal
gotm_ExtractParameters <- function(model, parameters, parcount = model$parcount){
  if (!length(parcount)) stop('Missing parcount in model object.')
  if (missing(parameters)) {
    parameters <- model$coef
    if (!length(parameters)) stop('Missing (estimated) parameters.')
  }
  if (length(parameters) != sum(parcount)) stop('Wrong number of parameters.')

  reg.params <- parameters[1L : parcount[1L]]
  cpc <- cumsum(parcount)

  if (parcount[2L]) {
    thresh.lambda <- parameters[(cpc[1L] + 1L) : cpc[2L]]
  } else {
    stop('Lambda must be given.')
  }

  if (parcount[3L]) {
    thresh.gamma <- parameters[(cpc[2L] + 1L) : cpc[3L]]
  } else {
    thresh.gamma <- NULL
  }
  list(reg.params = reg.params, thresh.lambda = thresh.lambda, thresh.gamma = thresh.gamma)
}


#' INTERNAL: The log likelihood function
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib gotm
#' @importFrom Rcpp evalCpp
gotm_negLL <- function(parameters=model$coef, model, collapse = TRUE, use_weights = TRUE, negative = TRUE){
  link = gotm_c_link(model)
  if (collapse) {
    LLFunc(parameters, yi=as.numeric(unclass(model$y_i)),reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
           link=link,thresh_no_cov=model$thresh.no.cov, negative=negative, thresh_1_exp = model$control$thresh.1.exp,
           weights=model$weights,use_weights = use_weights, thresh_start=model$control$thresh.start, out_val = model$control$LL_out_val,
           method = model$method)
  } else {
    LLFuncIndv(parameters, yi=as.numeric(unclass(model$y_i)),reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
               link=link,thresh_no_cov=model$thresh.no.cov, negative=negative, thresh_1_exp = model$control$thresh.1.exp,
               weights=model$weights, thresh_start = model$control$thresh.start, use_weights = use_weights,
               method = model$method)
  }
}


#' INTERNAL: The gradient of the log likelihood function
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib gotm
#' @importFrom Rcpp evalCpp
gotm_derivLL <- function(parameters=model$coef, model,
                         collapse = TRUE, use_weights = TRUE, negative = FALSE){
  link = gotm_c_link(model)
  #be compatible with older versions
  if ((!length(model$YYY1)) || (!length(model$YYY1))){
    model <- calcYYY(model)
  }
  if (collapse) {
    LLGradFunc(parameters, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2,
               YYY3=model$YYY3[,-model$J],YYY4=model$YYY3[,-1],
               reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
               link=link,thresh_no_cov=model$thresh.no.cov,negative=negative, thresh_1_exp = model$control$thresh.1.exp,
               weights=model$weights, thresh_start = model$control$thresh.start, use_weights = use_weights,
               method = model$method)
  } else {
    LLGradFuncIndv(parameters, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2, YYY3=model$YYY3[,-model$J], YYY4=model$YYY3[,-1],
                   reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
                   link=link,thresh_no_cov=model$thresh.no.cov, negative=negative, thresh_1_exp = model$control$thresh.1.exp,
                   weights=model$weights,thresh_start = model$control$thresh.start, use_weights = use_weights,
                   method = model$method)
  }
}


#' INTERNAL: Fit \code{gotm}
#'
#' Fit the model.
#' @param model \code{gotm} object
#' @param start starting parameters
#' @param control a list with control parameters. See \code{\link{gotm.control}}.
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib gotm
#' @importFrom Rcpp evalCpp
gotm_fitter <- function(model, start = model$start, use_weights = TRUE){

  #be compatible with older versions
  if ((!length(model$YYY1)) || (!length(model$YYY1))){
    model <- calcYYY(model)
  }
  control <- model$control

  link <- gotm_c_link(model)

  LLgr <- function(par, neg = TRUE) LLGradFunc(par, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2, YYY3=model$YYY3[,-model$J],
                                               YYY4=model$YYY3[,-1],
                                               reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
                                               link=link,thresh_no_cov=model$thresh.no.cov, negative=neg, thresh_1_exp = model$control$thresh.1.exp,
                                               weights=model$weights, thresh_start=model$control$thresh.start, use_weights = use_weights,
                                               method = model$method)
  LLfn <- function(par, neg = TRUE) LLFunc(par, yi=as.numeric(unclass(model$y_i)),reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
                                           link=link, thresh_no_cov=model$thresh.no.cov, negative=neg, thresh_1_exp = model$control$thresh.1.exp,
                                           weights=model$weights,use_weights = use_weights, thresh_start=model$control$thresh.start,
                                           out_val = model$control$LL_out_val, method = model$method)
  if (model$method==0) {

    #Two methods one after each othe
    fastgradfit <- function(fit){
      #BFGS and CG method
      try({
        fit <- optim(par = fit$par, fn = LLfn, gr = LLgr,
                     method = 'BFGS', hessian = FALSE, control=list(maxit=model$control$bgfs.maxit, reltol=model$control$bgfs.reltol))
        fit <- optim(par = fit$par, fn = LLfn, gr = LLgr,
                     method = 'CG', hessian = FALSE, control=list(maxit=model$control$cg.maxit, reltol=model$control$cg.reltol))
      }, silent = FALSE)

      return(fit)
    }

    z <- try({
      fit <- list(par = start)
      fit <- fastgradfit(fit)
      if (! model$control$quick.fit) {
        fit <- suppressWarnings(nlm(f = LLfn, p=fit$par, gradtol = model$control$nlm.gradtol, steptol = model$control$nlm.steptol, hessian = FALSE, iterlim=model$control$nlm.maxit))
        fit <- list(par=fit$estimate, value=fit$minimum)
      }
    }, silent = FALSE)
    if (class(z) == "try-error") stop('Optimization cannot continue.')

    model$coef <- fit$par
    model$LL <- unname(-fit$value)
  } else {
    model$coef <- start
    model$LL <- LLfn(start)
  }
  model
}


#' Auxiliary for controlling \code{gotm} fitting
#'
#' @description
#' Auxiliary function for controlling \code{gotm} fitting. Use this function to set control
#' parameters of the \code{\link{gotm}} and other related functions.
#' @param grad.eps epsilon for hessian function.
#' @param quick.fit logical, if TRUE extensive optimization methods are ignored and only BFGS and CG methods are run.
#' @param trace logical, if to trace model fitting
#' @param thresh.start internal parameter (threshold 0), under developement, do not change.
#' @param thresh.1.exp internal parameter (logical if to exponentiate threshold function at threshold 1), under developement, do not change.
#' @param LL_out_val internal parameter (LL value returned if LL cannot be computed), under developement, do not change.
#' @seealso \code{\link{gotm}}
#' @author Maciej J. Danko
#' @export
gotm.control<-function(grad.eps = 3e-5,
                       bgfs.maxit = 1e4,
                       cg.maxit = 1e4,
                       nlm.maxit = 150,
                       bgfs.reltol = 5e-10,
                       cg.reltol = 5e-10,
                       nlm.gradtol = 1e-7,
                       nlm.steptol = 1e-7,
                       quick.fit = TRUE,
                       trace = TRUE,
                       thresh.start = -Inf,
                       thresh.1.exp = FALSE,
                       LL_out_val = -Inf){

  list(grad.eps = grad.eps,
       bgfs.maxit = bgfs.maxit,
       cg.maxit = cg.maxit,
       nlm.maxit = nlm.maxit,
       bgfs.reltol = bgfs.reltol,
       cg.reltol = cg.reltol,
       nlm.gradtol = nlm.gradtol,
       nlm.steptol = nlm.steptol,
       quick.fit = quick.fit,
       trace = trace,
       thresh.start = thresh.start,
       thresh.1.exp = thresh.1.exp,
       LL_out_val = LL_out_val)
}


#frequency weight potraktowac jak w  clm czyli bez PSU
#' Auxiliary for setting a simple survey design for \code{gotm}
#'
#' @param PWeights Probability weights (the inverse of Probability of an observation being selected into the sample)
#' @param FWeights Frequency weights.
# Either \code{PWeight} or \code{FWeight} can be delivered (but not both simultaneously)
#' @param PSU Identificator of the PSU unit. Each P- and F- weight should correspond to exactly one PSU.
#' @param CountryID PSU are typically prescribed to only one country.
#' @export
#' @author Maciej J. Danko
# if ((length(PWeights) && length(FWeights))) stop('Please deliver either "PWeights" or "FWeights".')
gotm.design<-function(PWeights = NULL, FWeights = NULL, PSU = NULL, CountryID = NULL){
  if (length(PSU) && length(CountryID)) {
    if (length(PSU)!= length(CountryID)) stop('"PSU" and "CountryID" not of equal lengths.',call.=NULL)
    PSU <- paste(CountryID,PSU,sep='_')
  }
  tmp <- list(PWeights = PWeights, FWeights = FWeights, PSU = PSU)
  if (sum(sapply(tmp, length))){
    if (!length(PSU) && length(PWeights)) stop('"PSU" must be given.',call.=NULL)
    if (length(PSU) && length(unique(PSU)) == 1) stop('There is only one "PSU". "PWeights" cannot be used.',call.=NULL) #Correct this!!!!
    if (!(length(PWeights) + length(FWeights))) stop('P- or F- weights must be given',call.=NULL)
    if (length(PWeights) && (length(PSU) != length(PWeights))) stop('"PWeights" and "PSU" must have the same length.',call.=NULL)
  }
  tmp
}


#' Fit Generelaized Ordered Choice Threshold Model
#'
#'-- use_weights is unused....
#' @param reg.formula formula used to model latent process.
#' @param thresh.formula formula used to model threshold variable.
#' Any dependent variable (left side of "~") will be ignored.
#' @param data a data frame including all modeled variables.
#' @param survey an optional survey a survey design. Empty list indicates no survey design. See \code{\link{gotm.design}}.
#' @param link the link function. The possible values are \code{"probit"} (default) and \code{"logit"}.
#' @param start starting values in the form \code{c(latent_parameters, threshold_lambdas, threshold_gammas)}
#' @param hessain logical indicating if to calculate hessian matrix
#' but will not improve the fit, and \code{'no'} use external starting values, which must be delivered.
#' @param control a list with control parameters. See \code{\link{gotm.control}}.
#' @export
#' @author Maciej J. Danko
gotm<- function(reg.formula,
                thresh.formula = as.formula('~ 1'),
                data,
                method = c('linear','exp'),
                #overdispersion = FALSE,
                survey = list(),
                link = c('probit', 'logit'),
                start = NULL,
                hessian = TRUE,
                control = list()){

  my.grad <- function(fn, par, eps, ...){
    sapply(1L : length(par), function(k){
      epsi <- rep(0L, length(par))
      epsi[k] <- eps
      (fn(par + epsi, ...) - fn(par - epsi, ...))/2/eps
    })
  }

  if (missing(data)) data <- environment(reg.formula)

  method <- tolower(method[1])
  if (method=='linear') method=1 else if (method=='exp') method=0 else stop('Unknown method')
  link <- match.arg(link)
  control <- do.call("gotm.control", control)
  survey <- do.call("gotm.design", survey)

  if (length(start) && class(start) == 'gotm'){
    if (link != start$link) {
      stop ('Model in "start" is not compatible and will not be used.', call.=NULL)
    } else {
      tmp <- deparse(substitute(start))
      start <- get.start.gotm(object = start, reg.formula = reg.formula,
                              thresh.formula = thresh.formula,
                              data = data, asList = FALSE)
      cat('Model "',tmp,'" was used to get starting values.\n',sep='')
    }
  } else if (length(start) && !is.double(start)) stop('Wrong format of "start".', call.=NULL)

  model <- NULL
  model$control <- control
  model$link <- link
  model$method <- method

  if (length(thresh.formula)>2L){
    warning(call. = FALSE, 'The treshold formula should be given without dependent variable.')
    thresh.formula[[2]] <- NULL
  }

  reg.formula <- update.formula(reg.formula, '~.+1')
  if (any(grepl('offset(',as.character(reg.formula[[3]]),fixed=TRUE))) stop('Offset not supported.', call.=NULL)
  model$reg.formula <- reg.formula
  model$reg.mm <- as.matrix(model.matrix(reg.formula, data = data))
  model$reg.lev<-lapply(model.frame(model$reg.formula, data = data), function(k) if (is.factor(k)) as.matrix(table(k)) else 'Not a facor')
  reg.names <- colnames(model$reg.mm)
  grpi <- grepl('(Intercept)', colnames(model$reg.mm), fixed = TRUE)
  model$reg.mm <- as.matrix(model$reg.mm[,!grpi])
  reg.names <- reg.names[!grpi]

  thresh.formula <- update.formula(thresh.formula, '~.+1')
  if (any(grepl('offset(',as.character(thresh.formula[[2]]),fixed=TRUE))) stop('Offset not supprted.', call.=NULL)
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

  model$y_i <- model.frame(reg.formula, data = data)[,all.vars(reg.formula[[2]])]
  if (!is.factor(model$y_i)) stop('Response must be a factor with ordered levels.', call.=NULL)
  model$y_latent_i <- NA# latent
  model$Ey_i <- NA# ordinal classified utput
  model$J <- length(levels(model$y_i))
  model$N <- length(model$y_i)
  if (model$J<3L) stop ('Response must have 3 or more levels.', call.=NULL)

  model$thresh.extd <- matrix(rep_row(model$thresh.mm, model$J-1),model$N, NCOL(model$thresh.mm)*(model$J-1))

  if (sum(sapply(survey, length))){
    model$design <- survey
    if (length(survey$PWeights) && (length(survey$FWeights))) {
      model$weights <- survey$FWeights/survey$PWeights
    } else if (length(survey$PWeights)) {
      model$weights <- 1L/survey$PWeights
    } else if (length(survey$FWeights)){
      model$weights <- survey$FWeights
    } else stop('This should not happen.', call.=NULL)
    if (length(model$weights) != model$N) {
      stop('Vector of survey weights must be of the same length as data.', call.=NULL)
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

  #calculate special matrices for gradient calaculation
  model <- calcYYY(model)

  if (!length(start)) {
    if (model$control$trace) cat('Calculating starting parameters...')
    model <- suppressWarnings(get.vglm.start(model, data))
  } else {
    model$start <- start
  }

  if (model$control$trace && !model$method) cat(' done\nFitting the model...')

  model <- gotm_fitter(model, start = model$start)
  # if (model$method) {
  #   model$coef <- model$vglm.start
  #   model$coef.ls <- model$vglm.start.ls
  # }
  class(model) <- 'gotm'
  colnames(model$thresh.mm) <- thresh.names
  names(model$coef) <- coefnames

  if (model$control$trace) cat(' done\nCalculating maximum latent range...')

  model$Ey_i <- classify.ind(model)
  model$y_latent_i <- gotm_Latent(model$coef[seq_len(model$parcount[1])], model)
  p <- gotm_ExtractParameters(model)
  model$coef.ls <- p
  model$maxlatentrange <- sort(gotm_latentrange(model=model, data=data))
  model$maxobservedlatentrange <-  sort(range(gotm_Latent(p$reg.params,model)))
  if (model$control$trace) cat(' done\n')

  if (hessian) {
    if (model$control$trace) cat('Calculating hessian...')
    hes <- my.grad(fn = gotm_derivLL, par = model$coef, model=model, eps = model$control$grad.eps, collapse = TRUE, negative=FALSE)
    model$hessian <- hes
    model$vcov <- try(solve(-hes), silent = T)
    if (class(model$vcov) == 'try-error')
      warning(call. = FALSE, 'Model is probably unidentifiable, $vcov (variance-covariance matrix) cannot be computed.')
    if (model$control$trace) cat(' done\nCalculating estfun...')
    model$estfun <- gotm_derivLL(model$coef, model, collapse = FALSE)
    if (model$control$trace) cat(' done\n')
  }

  return(model)
}


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
  if ((old.rt %c% new.rt) && (new.rt %notc% old.rt)) {
    reg.mm <- model.matrix(reg.formula,data)[,-1]
    pr.new$reg.params <- rep(0, NCOL(reg.mm))
    old.ind <- which(colnames(reg.mm)%in%colnames(object$reg.mm))
    pr.new$reg.params[old.ind] <- pr$reg.params
    names(pr.new$reg.params)[old.ind] <- names(pr$reg.params)
    new.ind <- which(colnames(reg.mm)%notin%colnames(object$reg.mm))
    nnam <- colnames(reg.mm)[new.ind]
    names(pr.new$reg.params)[new.ind] <- nnam
  } else if ((new.rt %c% old.rt) && (old.rt %notc% new.rt)) {
    reg.mm <- model.matrix(reg.formula,data)[,-1]
    rm.ind <- which(colnames(object$reg.mm)%notin%colnames(reg.mm))
    tmp <- pr.new$reg.params[rm.ind]
    pr.new$reg.params <- pr.new$reg.params[-rm.ind]
    pr.new$thresh.lambda[1] <- pr.new$thresh.lambda[1] -
      mean(object$reg.mm[,rm.ind] * rep_row(tmp, NROW(reg.mm))) * length(tmp)
  }
  if ((old.tt %c% new.tt) && (new.tt %notc% old.tt)){
    thresh.mm <- model.matrix(thresh.formula,data)[,-1]
    pr.new$thresh.params <- rep(0, NCOL(thresh.mm))
    old.ind <- which(colnames(thresh.mm)%in%colnames(object$thresh.mm))
    pr.new$thresh.gamma[old.ind] <- pr$thresh.gamma
    names(pr.new$thresh.gamma)[old.ind] <- names(pr$thresh.gamma)
    new.ind <- which(colnames(thresh.mm)%notin%colnames(object$thresh.mm))
    nnam <- colnames(thresh.mm)[new.ind]
    names(pr.new$thresh.gamma)[new.ind] <- nnam
  } else if ((new.tt %c% old.tt) && (old.tt %notc% new.tt)) {
    thresh.mm <- model.matrix(thresh.formula,data)[,-1]
    rm.ind <- which(colnames(object$thresh.mm)%notin%colnames(thresh.mm))
    tmp <- pr.new$thresh.gamma[rm.ind]
    pr.new$thresh.gamma <- pr.new$thresh.gamma[-rm.ind]
    pr.new$thresh.lambda[1] <- pr.new$thresh.lambda[1] -
      mean(object$thresh.mm[,rm.ind] * rep_row(tmp, NROW(thresh.mm))) * length(tmp)
  }
  if (asList) return(pr.new) else
    return(c(pr.new$reg.params, pr.new$thresh.lambda, pr.new$thresh.gamma))
}
