#' INTERNAL: Calculation of cut-points (threshold)
#'
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_Threshold<-function(thresh.lambda, thresh.gamma, model){
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
hopit_Latent <- function(reg.params, model = NULL) model$reg.mm %*% (as.matrix(reg.params))


# #' INTERNAL: Calculate maximum possible latent range
# #' @param model a fitted \code{hopit} model.
# #' @param data a data used to fit the model
# #' @keywords internal
# hopit_latentrange <- function (model, data) {
#
#   cfm <- model$coef[seq_len(model$parcount[1])]
#
#   L <- apply(model$reg.mm,2,function(x)length(unique(x))) -1
#   FACTORS <- as.logical(L==1)
#   if (any(FACTORS)) {
#     d <- cfm[FACTORS]
#     r <- c(sum(d[d<0]), sum(d[d>0]))
#   }
#   if (any(!FACTORS)) {
#     #message('Continous covariate detected, model may not give realistic results.',appendLF = FALSE)
#     d <- matrix(cfm[!FACTORS],length(cfm[!FACTORS]),model$N) * t(model$reg.mm[,!FACTORS])
#     r <- c(min(d[d<0],r), max(r,d[d>0]))
#   }
#   r
# }


#' INTERNAL: Extract model parameters in a form of list
#'
#' Extract model parameters in a form of a list
#' @param model \code{hopit} object
#' @param parameters model parameters (optional). If not delivered then taken from \code{model$coef}
#' @author Maciej J. Danko
#' @keywords internal
hopit_ExtractParameters <- function(model, parameters, parcount = model$parcount){
  logTheta <- 0
  if (!length(parcount)) stop('Missing parcount in model object.')
  if (missing(parameters)) {
    parameters <- model$coef
    if (!length(parameters)) stop('Missing (estimated) parameters.')
  }
  if (length(parameters) != sum(parcount) + model$hasdisp) stop('Wrong number of parameters.')

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

  if (model$hasdisp) {
    list(reg.params = reg.params, thresh.lambda = thresh.lambda, thresh.gamma = thresh.gamma, logTheta = parameters[length(parameters)])
  } else {
    list(reg.params = reg.params, thresh.lambda = thresh.lambda, thresh.gamma = thresh.gamma, logTheta = logTheta)
  }
}


#' INTERNAL: The log likelihood function
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_negLL <- function(parameters=model$coef, model, collapse = TRUE, use_weights = model$use.weights, negative = TRUE){
  link = hopit_c_link(model)
  if (collapse) {
    LL <- LLFunc(parameters, yi=as.numeric(unclass(model$y_i)),reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
                 hasdisp = model$hasdisp,
                 link=link,thresh_no_cov=model$thresh.no.cov, negative=negative, thresh_1_exp = model$control$thresh.1.exp,
                 weights=model$weights,use_weights = use_weights, thresh_start=model$control$thresh.start, out_val = model$control$LL_out_val,
                 method = model$method)
  } else {
    LL <- LLFuncIndv(parameters, yi=as.numeric(unclass(model$y_i)),reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
                     hasdisp = model$hasdisp,
                     link=link,thresh_no_cov=model$thresh.no.cov, negative=negative, thresh_1_exp = model$control$thresh.1.exp,
                     weights=model$weights, thresh_start = model$control$thresh.start, use_weights = use_weights,
                     method = model$method)
  }
  if (use_weights) {
    LL <- LL / sum(model$weights) * model$N #scale likelihood
  }
  LL
}


#' INTERNAL: The gradient of the log likelihood function
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_derivLL <- function(parameters=model$coef, model,
                          collapse = TRUE, use_weights = model$use.weights, negative = FALSE){
  link = hopit_c_link(model)
  #be compatible with older versions
  if ((!length(model$YYY1)) || (!length(model$YYY1))){
    model <- calcYYY(model)
  }
  if (collapse) {
    LLgr <- LLGradFunc(parameters, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2,
                       YYY3=model$YYY3[,-model$J],YYY4=model$YYY3[,-1], hasdisp = model$hasdisp,
                       reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
                       link=link,thresh_no_cov=model$thresh.no.cov,negative=negative, thresh_1_exp = model$control$thresh.1.exp,
                       weights=model$weights, thresh_start = model$control$thresh.start, use_weights = use_weights,
                       method = model$method)
  } else {
    LLgr <- LLGradFuncIndv(parameters, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2, YYY3=model$YYY3[,-model$J], YYY4=model$YYY3[,-1],
                           reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
                           hasdisp = model$hasdisp,
                           link=link,thresh_no_cov=model$thresh.no.cov, negative=negative, thresh_1_exp = model$control$thresh.1.exp,
                           weights=model$weights,thresh_start = model$control$thresh.start, use_weights = use_weights,
                           method = model$method)
  }

  if (use_weights) {
    LLgr <- LLgr / sum(model$weights) * model$N #scale likelihood
  }
  LLgr
}


#' INTERNAL: Fit \code{hopit}
#'
#' Fit the model.
#' @param model \code{hopit} object
#' @param start starting parameters
#' @param control a list with control parameters. See \code{\link{hopit.control}}.
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_fitter <- function(model, start = model$start, use_weights = model$use.weights){

  #be compatible with older versions
  if ((!length(model$YYY1)) || (!length(model$YYY1))){
    model <- calcYYY(model)
  }
  control <- model$control

  link <- hopit_c_link(model)

  LLgr <- function(par, neg = TRUE) LLGradFunc(par, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2, YYY3=model$YYY3[,-model$J],
                                               YYY4=model$YYY3[,-1],hasdisp = model$hasdisp,
                                               reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
                                               link=link,thresh_no_cov=model$thresh.no.cov, negative=neg, thresh_1_exp = model$control$thresh.1.exp,
                                               weights=model$weights, thresh_start=model$control$thresh.start, use_weights = use_weights,
                                               method = model$method)
  LLfn <- function(par, neg = TRUE) LLFunc(par, yi=as.numeric(unclass(model$y_i)),reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
                                           link=link, thresh_no_cov=model$thresh.no.cov, negative=neg, thresh_1_exp = model$control$thresh.1.exp,
                                           weights=model$weights,use_weights = use_weights, thresh_start=model$control$thresh.start,
                                           out_val = model$control$LL_out_val, hasdisp = model$hasdisp, method = model$method)
  if (model$method==0) {

    #Two methods one after each othe
    fastgradfit <- function(fit, meto){
      #BFGS and CG method
      try({
        if (meto == 'BFGS') {
          fit <- optim(par = fit$par, fn = LLfn, gr = LLgr,
                       method = 'BFGS', hessian = FALSE, control=list(maxit=model$control$bgfs.maxit, reltol=model$control$bgfs.reltol))
        } else if (meto == 'CG'){
          fit <- optim(par = fit$par, fn = LLfn, gr = LLgr,
                       method = 'CG', hessian = FALSE, control=list(maxit=model$control$cg.maxit, reltol=model$control$cg.reltol))
        }
      }, silent = FALSE)

      return(fit)
    }

    z <- try({
      fit <- list(par = start)
      # sstart=start
      # sstart[length(sstart)]=1
      # LLFunc(sstart, yi=as.numeric(unclass(model$y_i)),reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
      #        link=link, thresh_no_cov=model$thresh.no.cov, negative=TRUE, thresh_1_exp = model$control$thresh.1.exp,
      #        weights=model$weights,use_weights = use_weights, thresh_start=model$control$thresh.start,
      #        out_val = model$control$LL_out_val, hasdisp = TRUE, method = model$method)
      # start=model$start
      # LLGradFunc(start, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2, YYY3=model$YYY3[,-model$J],
      #            YYY4=model$YYY3[,-1],hasdisp = TRUE,
      #            reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
      #            link=link,thresh_no_cov=model$thresh.no.cov, negative=TRUE, thresh_1_exp = model$control$thresh.1.exp,
      #            weights=model$weights, thresh_start=model$control$thresh.start, use_weights = use_weights,
      #            method = model$method)
      # my.grad(LLfn,start,1e-6)

      if ('CG' %in% control$fit.methods) fit <- fastgradfit(fit, meto = 'CG')
      if ('BFGS' %in% control$fit.methods) fit <- fastgradfit(fit, meto = 'BFGS')
      if (!model$control$quick.fit || !length(fit$par)) {
        if (!length(fit$par)) fit$par <- start
        if (model$control$trace) cat(' done\nImproving fit with nlm...')
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
  # cat('VGLM LL:',model$vglm.LL,'\n')
  # cat('START LL:',LLfn(start, neg=FALSE),'\n')
  # cat('FITTED LL:',LLfn(model$coef, neg=FALSE),'\n')
  # cat('FITTED LL:',model$LL,'\n')
  if (use_weights) {
    model$LL <- model$LL / sum(model$weights) * model$N #scale likelihood
  }
  model
}


#' survey:::htvar.matrix clone
#' @keywords internal
#' @author Thomas Lumley
#' @importFrom Matrix crossprod
clone.of.htvar.matrix <- function (xcheck, Dcheck) {
  if (is.null(dim(xcheck)))
    xcheck <- as.matrix(xcheck)
  rval <- apply(xcheck, 2, function(xicheck) apply(xcheck, 2, function(xjcheck) as.matrix(Matrix::crossprod(xicheck,
                                                                                                            Dcheck %*% xjcheck))))
  if (is.null(dim(rval)))
    dim(rval) <- c(1, 1)
  rval
}

#' survey:::ygvar.matrix clone
#' @keywords internal
#' @author Thomas Lumley
clone.of.ygvar.matrix <-function (xcheck, Dcheck)
{
  ht <- clone.of.htvar.matrix(xcheck, Dcheck)
  if (is.null(dim(xcheck))) {
    corr <- sum(Dcheck %*% (xcheck * xcheck))
  }
  else {
    corr <- apply(xcheck, 2, function(xicheck) apply(xcheck,
                                                     2, function(xjcheck) sum(Dcheck %*% (xicheck * xjcheck))))
  }
  rval <- ht - corr
}

#' survey:::ppsvar clone
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
        psmeans <- rowsum(y/psw, psvar, reorder = TRUE)/as.vector(table(factor(psvar)))
        x <- y - psmeans[match(psvar, sort(unique(psvar))),
                         ] * psw
      }
    }
  }
  dcheck <- design$dcheck
  if (length(dcheck) != 1)
    stop("Multistage not implemented yet")
  rval <- switch(est, HT = clone.of.htvar.matrix(rowsum(x, dcheck[[1]]$id,
                                                        reorder = FALSE), dcheck[[1]]$dcheck),
                 YG = clone.of.ygvar.matrix(rowsum(x, dcheck[[1]]$id, reorder = FALSE), dcheck[[1]]$dcheck),
                 stop("can't happen"))
  rval
}

#' Calculation of variance-covariance matrix for specified survey design
#' @keywords internal
#' @param object a hopit object
#' @param design a survey.design object
#' @author Thomas Lumley, modified by Maciej J. Dańko
svy.varcoef.hopit <- function (Ainv, estfun, design) {
  # Ainv <- object$vcov #summary(glm.object)$cov.unscaled
  # estfun <- object$estfun #model.matrix(glm.object) * resid(glm.object, "working") *
  if (inherits(design, "survey.design2"))
    V <- survey::svyrecvar(estfun %*% Ainv, design$cluster, design$strata,
                           design$fpc, postStrata = design$postStrata)
  else if (inherits(design, "twophase"))
    V <- survey::twophasevar(estfun %*% Ainv, design)
  else if (inherits(design, "twophase2"))
    V <- survey::twophase2var(estfun %*% Ainv, design)
  else if (inherits(design, "pps"))
    V <- clone.of.ppsvar(estfun %*% Ainv, design)
  else V <-survey::svyCprod(estfun %*% Ainv, design$strata, design$cluster[[1]],
                            design$fpc, design$nPSU, design$certainty, design$postStrata)
  if (inherits(design, "svyrep.design")) {
    Vtest <- vcov(survey::svymean(estfun %*% Ainv*design$prob, design)) * length(design$prob)^2
    if (max(abs(Vtest-V)) > 1e-4) warning('Something wrong with survey package.')
  }
}


#' Auxiliary for controlling \code{hopit} fitting
#'
#' @description
#' Auxiliary function for controlling \code{hopit} fitting. Use this function to set control
#' parameters of the \code{\link{hopit}} and other related functions.
#' @param grad.eps epsilon for hessian function.
#' @param quick.fit logical, if TRUE extensive optimization methods are ignored and only BFGS and CG methods are run.
#' @param trace logical, if to trace model fitting
#' @param thresh.start internal parameter (threshold 0), under developement, do not change.
#' @param thresh.1.exp internal parameter (logical if to exponentiate threshold function at threshold 1), under developement, do not change.
#' @param LL_out_val internal parameter (LL value returned if LL cannot be computed), under developement, do not change.
#' @seealso \code{\link{hopit}}
#' @author Maciej J. Danko
#' @export
hopit.control<-function(grad.eps = 3e-5,
                        bgfs.maxit = 1e4,
                        cg.maxit = 1e4,
                        nlm.maxit = 150,
                        bgfs.reltol = 5e-10,
                        cg.reltol = 5e-10,
                        nlm.gradtol = 1e-7,
                        nlm.steptol = 1e-7,
                        fit.methods = c('CG','BFGS'),
                        quick.fit = TRUE,
                        trace = TRUE,
                        thresh.start = -Inf,
                        thresh.1.exp = FALSE,
                        LL_out_val = -Inf){

  if (!length(fit.methods)) stop('Please specify fit method',call. = NULL) else fit.methods <- toupper(fit.methods)
  if (any(fit.methods %notin% c('CG','BFGS'))) stop ('Unknown fit method.',call.=NULL)

  list(grad.eps = grad.eps,
       bgfs.maxit = bgfs.maxit,
       cg.maxit = cg.maxit,
       nlm.maxit = nlm.maxit,
       bgfs.reltol = bgfs.reltol,
       cg.reltol = cg.reltol,
       nlm.gradtol = nlm.gradtol,
       nlm.steptol = nlm.steptol,
       fit.methods = fit.methods,
       quick.fit = quick.fit,
       trace = trace,
       thresh.start = thresh.start,
       thresh.1.exp = thresh.1.exp,
       LL_out_val = LL_out_val)
}


#' Extract Theta parameter from the model
#'
#' @param moodel fitted \code{hopit} model.
#' @export
#' @author Maciej J. Danko
gettheta <- function(model) unname(exp(model$coef.ls$logTheta))


#' Fit Generelaized Ordered Choice Threshold Model
#'
#' @param reg.formula formula used to model latent process.
#' @param thresh.formula formula used to model threshold variable.
#' Any dependent variable (left side of "~") will be ignored.
#' @param data a data frame including all modeled variables.
#' @param design an optional survey design. Use \code{\link[survey]{svydesign}} function to specify the design.
#' @param weights an optional weights. Use design to construct survey weights.
#' @param link the link function. The possible values are \code{"probit"} (default) and \code{"logit"}.
## @param start starting values in the form \code{c(latent_parameters, threshold_lambdas, threshold_gammas)}
#' @param control a list with control parameters. See \code{\link{hopit.control}}.
#' @export
#' @author Maciej J. Danko
hopit<- function(reg.formula,
                 thresh.formula = as.formula('~ 1'), # ~1 not tested!!!
                 data,
                 method = c('hopit','vglm'),
                 overdispersion = FALSE,
                 design = list(),
                 weights = NULL,
                 link = c('probit', 'logit'),
                 start.method = c('glm','vglm'),
#                 start = NULL,
                 control = list()){

  if (!overdispersion) remove.theta = FALSE else remove.theta = TRUE
  if (missing(data)) data <- environment(reg.formula)
  start.method <- tolower(start.method[1])
  method <- tolower(method[1])
  if (method=='vglm') method <- 1 else if (method=='hopit') method <- 0 else stop('Unknown method')
  link <- match.arg(link)
  control <- do.call("hopit.control", control)

  # if (length(start) && class(start) == 'hopit'){
  #   if (link != start$link) {
  #     stop ('Model in "start" is not compatible and will not be used.', call.=NULL)
  #   } else {
  #     tmp <- deparse(substitute(start))
  #     start <- get.start.hopit(object = start, reg.formula = reg.formula,
  #                              thresh.formula = thresh.formula,
  #                              data = data, asList = FALSE)
  #     cat('Model "',tmp,'" was used to get starting values.\n',sep='')
  #   }
  # } else if (length(start) && !is.double(start)) stop('Wrong format of "start".', call.=NULL)

  model <- NULL
  model$control <- control
  model$link <- link
  model$method <- as.logical(method)
  model$start.method <- start.method
  model$hasdisp <- overdispersion
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
  if (any(grepl('offset(',tolower(as.character(thresh.formula[[2]])),fixed=TRUE))) stop('Offset not supprted.', call.=NULL)
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

  model$reg.terms <- attr(terms(as.formula(reg.formula)),"term.labels")
  model$thresh.terms <-attr(terms(as.formula(thresh.formula)),"term.labels")

  model$y_i <- model.frame(reg.formula, data = data)[,all.vars(reg.formula[[2]])]
  if (!is.factor(model$y_i)) stop('Response must be a factor with ordered levels.', call.=NULL)
  model$y_latent_i <- NA# latent
  model$Ey_i <- NA# ordinal classified utput
  model$J <- length(levels(model$y_i))
  model$N <- length(model$y_i)
  if (model$J<3L) stop ('Response must have 3 or more levels.', call.=NULL)

  model$thresh.extd <- matrix(rep_row(model$thresh.mm, model$J-1),model$N, NCOL(model$thresh.mm)*(model$J-1))

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
  if (model$hasdisp) coefnames <- c(coefnames, 'logTheta')

  model$weights <- NULL
  if (length(weights) && length(design)) stop('Multiple weights specification detected. Please use either design or weights parameter.', call.=NULL)
  if (length(design)) model$weights <- design$prob else if (length(weights)) model$weights <- weights
  if (!length(model$weights)) {
    model$weights <- rep(1, model$N)
    model$use.weights <- FALSE
  } else model$use.weights <- TRUE

  if (length(model$weights) != model$N) {
    print(length(model$weights))
    print(model$N)
    stop('Vector of survey weights must be of the same length as data.', call.=NULL)
  }
  model$weights <- as.vector(matrix(model$weights, 1L, model$N))
  #scaling weights
  model$weights <- model$N * model$weights / sum(model$weights)

  #calculate special matrices for gradient calaculation
  model <- calcYYY(model)

  # if (!length(start)) {
  if (model$control$trace) cat('Calculating starting parameters...')
  model <- suppressWarnings(get.vglm.start(model, data))
  # } else {
  #   model$start <- start
  # }

  if (model$control$trace && !model$method) cat(' done\nFitting the model...')

  model <- hopit_fitter(model, start = model$start)
  model$vglm <- NULL
  model$vglm.LL <- NULL

  class(model) <- 'hopit'
  colnames(model$thresh.mm) <- thresh.names
  names(model$coef) <- coefnames

  if (model$control$trace) cat(' done\nCalculating maximum latent range...')
  p <- hopit_ExtractParameters(model)
  model$alpha <- hopit_Threshold(p$thresh.lambda, p$thresh.gamma, model)
  model$Ey_i <- classify.ind(model)
  model$y_latent_i <- hopit_Latent(model$coef[seq_len(model$parcount[1])], model)
  model$coef.ls <- p
  #model$maxlatentrange <- sort(hopit_latentrange(model=model, data=data)) #waht is the difference in the context of continuous variables
  model$maxobservedlatentrange <-  range(model$y_latent_i)
  if (model$control$trace) cat(' done\n')
  model$deviance <- -2 * model$LL
  k <- 2

  if (model$control$trace) cat('Calculating hessian...')

  hes <- my.grad(fn = hopit_derivLL, par = model$coef, model=model, eps = 1e-4, collapse = TRUE, negative=FALSE)
  if (model$hasdisp && remove.theta) {
    hes <- hes[-nrow(hes),-ncol(hes)] #remove theta from vcov
    model$coef <- model$coef[-length(model$coef)] #remove from coef
  }
  model$hessian <- hes

  model$vcov.basic <- try(base::solve(-hes), silent = FALSE)
  if (class(model$vcov) == 'try-error') {
    warning(call. = FALSE, 'Model is probably unidentifiable, $vcov (variance-covariance matrix) cannot be computed.')
    model$vcov.basic <- NA
  }
  if (model$control$trace) cat(' done\nCalculating estfun...')
  if (model$hasdisp && remove.theta) COEF <- c(model$coef,model$coef.ls$logTheta) else COEF <- model$coef
  model$estfun <- hopit_derivLL(COEF, model, collapse = FALSE)
  if (remove.theta) model$estfun <- model$estfun[,-ncol(model$estfun)]
  if (model$control$trace) cat(' done\n')

  if (length(model$design)) {
    if (model$control$trace) cat('Including survey design...')
    model$vcov <- svy.varcoef.hopit(model$vcov.basic, model$estfun, design)

    # model$misspec <- try(base::eigen(base::solve(model$vcov.basic) %*% model$vcov, only.values = TRUE)$values, silent = TRUE)
    # if (class(model$misspec) == 'try-error' || is.na(model$vcov)) {
    #   warning(call. = FALSE, 'Cannot estimate AIC using survey design.')
    #   model$misspec <- model$deltabar <- model$eff.p <- model$AIC <- NA
    # } else {
    #   # model$deltabar <- mean(model$misspec)
    #   # model$eff.p <- sum(model$misspec)
    #   # model$AIC <- model$deviance + k * sum(model$misspec)
    # }
    #model$df.residual <- survey::degf(design) + 1 - length(model$coef[!is.na(model$coef)])
    model$misspec <- model$deltabar <- model$eff.p <- model$AIC <- NA
    if (model$control$trace) cat(' done\n')
  } else {
    model$vcov <- model$vcov.basic
    model$AIC <- model$deviance + k * (length(model$coef.ls$reg.params)+
                                       length(model$coef.ls$thresh.lambda)+
                                       length(model$coef.ls$thresh.gamma)+model$hasdisp)
    model$misspec <- model$deltabar <- model$eff.p <- NA
  }
  return(model)
}


#' #' Get starting parameters from less or more complicated hierarchical models
#' #'
#' #' @export
#' get.start.hopit <- function(object, reg.formula, thresh.formula, data, asList = FALSE){
#'   old.rf <- object$reg.formula
#'   old.tf <- object$thresh.formula
#'   if (deparse(object$reg.formula[[2]])!=deparse(reg.formula[[2]])) stop('Models have different dependent variables')
#'   if (length(thresh.formula)>2L){
#'     warning(call. = FALSE, 'The treshold formula should be given without dependent variable.')
#'     thresh.formula[[2]] <- NULL
#'   }
#'   old.rt<-attr(terms(old.rf),"term.labels")
#'   old.tt<-attr(terms(old.tf),"term.labels")
#'   new.rt<-attr(terms(as.formula(reg.formula)),"term.labels")
#'   new.tt<-attr(terms(as.formula(thresh.formula)),"term.labels")
#'   pr <- hopit_ExtractParameters(object)
#'   pr.new <-pr
#'   if ((old.rt %c% new.rt) && (new.rt %notc% old.rt)) {
#'     reg.mm <- model.matrix(reg.formula,data)[,-1]
#'     pr.new$reg.params <- rep(0, NCOL(reg.mm))
#'     old.ind <- which(colnames(reg.mm)%in%colnames(object$reg.mm))
#'     pr.new$reg.params[old.ind] <- pr$reg.params
#'     names(pr.new$reg.params)[old.ind] <- names(pr$reg.params)
#'     new.ind <- which(colnames(reg.mm)%notin%colnames(object$reg.mm))
#'     nnam <- colnames(reg.mm)[new.ind]
#'     names(pr.new$reg.params)[new.ind] <- nnam
#'   } else if ((new.rt %c% old.rt) && (old.rt %notc% new.rt)) {
#'     reg.mm <- model.matrix(reg.formula,data)[,-1]
#'     rm.ind <- which(colnames(object$reg.mm)%notin%colnames(reg.mm))
#'     tmp <- pr.new$reg.params[rm.ind]
#'     pr.new$reg.params <- pr.new$reg.params[-rm.ind]
#'     pr.new$thresh.lambda[1] <- pr.new$thresh.lambda[1] -
#'       mean(object$reg.mm[,rm.ind] * rep_row(tmp, NROW(reg.mm))) * length(tmp)
#'   }
#'   if ((old.tt %c% new.tt) && (new.tt %notc% old.tt)){
#'     thresh.mm <- model.matrix(thresh.formula,data)[,-1]
#'     pr.new$thresh.params <- rep(0, NCOL(thresh.mm))
#'     old.ind <- which(colnames(thresh.mm)%in%colnames(object$thresh.mm))
#'     pr.new$thresh.gamma[old.ind] <- pr$thresh.gamma
#'     names(pr.new$thresh.gamma)[old.ind] <- names(pr$thresh.gamma)
#'     new.ind <- which(colnames(thresh.mm)%notin%colnames(object$thresh.mm))
#'     nnam <- colnames(thresh.mm)[new.ind]
#'     names(pr.new$thresh.gamma)[new.ind] <- nnam
#'   } else if ((new.tt %c% old.tt) && (old.tt %notc% new.tt)) {
#'     thresh.mm <- model.matrix(thresh.formula,data)[,-1]
#'     rm.ind <- which(colnames(object$thresh.mm)%notin%colnames(thresh.mm))
#'     tmp <- pr.new$thresh.gamma[rm.ind]
#'     pr.new$thresh.gamma <- pr.new$thresh.gamma[-rm.ind]
#'     pr.new$thresh.lambda[1] <- pr.new$thresh.lambda[1] -
#'       mean(object$thresh.mm[,rm.ind] * rep_row(tmp, NROW(thresh.mm))) * length(tmp)
#'   }
#'   if (asList) return(pr.new) else
#'     return(c(pr.new$reg.params, pr.new$thresh.lambda, pr.new$thresh.gamma))
#' }

