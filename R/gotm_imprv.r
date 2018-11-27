#' INTERNAL: Converts a vector of an categorical variable into a matrix with dummies in columns
#'
#' @param V a vector of categories.
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @keywords internal
Vector2DummyMat<-function(V) sapply(levels(as.factor(V)), function(k) as.factor(V) == k)*1L

#' INTERNAL: Converts a matrix with dummies in columns into categorical vector
#'
#' @param D a matrix of dummies.
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @keywords internal
DummyMat2Vector<-function(D) D %*% ((1L : dim(D)[2L]) -1L)

#' INTERNAL: Do cumsum() in each row of a matrix
#'
#' @param mat a matrix.
#' @keywords internal
cumsum_row<-function(mat) t(apply(as.matrix(mat), 1L, cumsum))

#' INTERNAL: Column path in a matrix
#'
#' For each row of a matrix \code{mat} extracts a value corresponding to a column stored in a vector \code{y}.
#' @param mat a matrix
#' @param y a vector with column indices corresponding to each row in the matrix \code{mat}.
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @keywords internal
col_path<-function(mat, y, offset = 0) colpath(mat, v, offset) # RcppEigen

#' Convert individual data to frequency table of unique combination of dependent and independent variables
#'
#' @param formula formula indicating, which variables will be used to construct new database.
#' @param data data.frame including all variables listed in formula.
#' @param FreqNam name of the column with frequencies.
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
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
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_reg_thresh<-function(thresh.lambda,thresh.gamma, model)  ind_reg_thresh(model$thresh.mm, thresh.lambda, thresh.gamma) #RcppEigen

#' @keywords internal
gotm_c_init<-function(model){
  model$link <- tolower(model$link)
  model$thresh.method <- tolower(model$thresh.method)
  #model$control$thresh.fun <- tolower(model$control$thresh.fun)

  #CHANGE (identical(deparse(model$control$thresh.fun), deparse(identity))) thresh func only as "exp" "ident"

  # if (model$control$thresh.fun %in% c('exp','identity')){
  #   if (model$control$thresh.fun=='exp') thresh.fun=0 else thresh.fun=1
  # } else stop(paste('Unknown threshold function:', model$control$thresh.fun), call. = NULL)

  if (model$link %in% c('probit','logit')){
    if (model$link=='probit') link=0 else link=1
  } else stop(paste('Unknown link function:',model$link),call. = NULL)

  if (model$thresh.method %in% c('jurges','hopit','vglm')){
    if (model$thresh.method=='vglm') stop('Function works only for "Jurges" and "Hopit models.',call. = NULL)
    if (model$thresh.method=='jurges') thresh_method <- 0 else thresh_method <-1
  } else stop(paste('Unknown threshold method:',model$thresh_method),call. = NULL)

  if (model$control$alpha_0 == 'auto') {
    use_alpha <- 0
    alpha_0 <- 0
    } else {
    if (is.numeric(model$control$alpha_0)) {
      use_alpha <- 1
      alpha_0 <- model$control$alpha_0
    } else stop(paste('Cannot interpret control$alpha_0 =', model$control$alpha_0),call. = NULL)
  }
  list(use_alpha=use_alpha, thresh_method=thresh_method, #thresh_fun=thresh.fun,
       link=link, alpha_0=alpha_0)
}

#' INTERNAL: Calculation of cut-points (threshold)
#'
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_Threshold<-function(thresh.lambda, thresh.gamma, model = NULL){
  z <- gotm_c_init(model)
  getThresholds(model$thresh.mm, thresh.lambda, thresh.gamma, model$thresh.no.cov,
                z$thresh_method, #z$thresh_fun,
                z$use_alpha, z$alpha_0) #RcppEigen

}

#' INTERNAL: Fit vglm to get parameters of the model
#'
#' @importFrom VGAM vglm
#' @keywords internal
fit.vglm <-function(model, data, start=NULL){
  reg.formula <- model$reg.formula
  thresh.formula <- model$thresh.formula
  thresh.method <- model$thresh.method
  #thresh.fun <- model$control$thresh.fun
  if (length(thresh.formula)>2) thresh.formula[[2]] <- NULL
  thrf <- deparse(thresh.formula[[2]])
  big.formula <- update(reg.formula, paste('~ ', thrf,' + . + 1'))
  Y <<- Vector2DummyMat(data[,paste(reg.formula[[2]])])
  w <- model$weights
  data$w <- w
  big.formula[[2]] <- as.name('Y')
  small.formula <- formula(paste('FALSE ~', thrf))
  cat('Running vglm...')
  mv2<-switch(model$link,
              probit = vglm(big.formula, weights = w, data = data, coefstart = start,
                            family = cumulative(parallel = small.formula, link = 'probit')), #direct substitution of link doesn't work
              logit = vglm(big.formula, weights = w, data = data, coefstart = start,
                           family = cumulative(parallel = small.formula, link = 'logit')))
  rm(Y, envir = .GlobalEnv)
  mv2.par <- c(coef(mv2)[-(1:sum(model$parcount[2:3]))], coef(mv2)[(1:sum(model$parcount[2:3]))])
  model$vglm <- mv2
  model$vglm.LL<-logLik(mv2)
  model$vglm.start.ls <- gotm_ExtractParameters(model, mv2.par)
  model$vglm.start <- mv2.par
  model
}


#' INTERNAL: Fit vglm to get parameters of the model
#'
#' @param model \code{gotm} object.
#' @param data data.frame with data used to fit the model.
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @importFrom VGAM vglm
#' @keywords internal
get.vglm.start<-function(model, data, start = NULL){
  subs <- function (a, k, z) {a[k] <- z; a}
  gmat <- function (x) matrix(x, model$J-1L, as.integer(round(length(x) / (model$J-1L))))
  subsgmat <- function (a, k, z) {
    an <- names(a)
    a <- gmat(a)
    a[k,] <- z
    a <- as.vector(a)
    names(a) <- an
    a
  }
  extrgmat<-function (a, k){
    an <- names(a)
    a <- gmat(a)
    an <- gmat(an)
    z <- as.vector(a[k,])
    names(z) <- as.vector(an[k,])
    z
  }

  model <- fit.vglm(model, data, start)
  mv2.par <- model$vglm.start
  par.ls <- model$vglm.start.ls
  cat(' done\n')
  if (thresh.method == 'vglm'){
    model$vglm.start <- mv2.par
    model$reg.start <- par.ls$reg.params
    model$lambda.start <- par.ls$thresh.lambda
    model$gamma.start <-par.ls$thresh.gamma
    model$start <- c(model$reg.start, model$lambda.start, model$gamma.start)
    model$start.LL <- model$vglm.LL
    model$LL <- model$vglm.LL
    model$coef <- model$start
  } else {
    cat('Recalculating parameters...')

    if (thresh.method == 'jurges'){
      z <- vglm2gotm_jurges_exp(par.ls$reg.params, par.ls$thresh.lambda, par.ls$thresh.gamma)#, (thresh.method != 'jurges')*1)
      if (model$control$alpha_0 == 0) {
        Z=min(z$thresh_lambda[1],z$reg_params)
        if (Z<0) {
          z$thresh_lambda[1] <- z$thresh_lambda[1] - Z
          z$reg_params <- z$reg_params - Z
        }
      }
      model$vglm.start <- mv2.par
      model$reg.start <- z$reg_params
      model$lambda.start <- z$thresh_lambda
      model$gamma.start <-z$thresh_gamma
      model$start <- z$coef
      model$start.LL <- gotm_negLL(model$start, model, negative = FALSE)
    } else stop(paste('The thresh.method =',thresh.method,'not implemented yet.'),call. = NULL)
    # par.ls$reg.params <- -par.ls$reg.params
    # #par.ls$thresh.lambda <- c(par.ls$thresh.lambda[1], diff(par.ls$thresh.lambda))
    #
    # if (length(par.ls$thresh.gamma))
    #   alpha.thresh <- gotm_reg_thresh(thresh.lambda = par.ls$thresh.lambda,
    #                                   thresh.gamma = par.ls$thresh.gamma,
    #                                   model = model)
    #
    # #alpha.thresh<-t(matrix(par.ls$thresh.lambda, model$J - 1L, model$N)) + model$thresh.mm %*% par.ls$thresh.gamma
    # #get lambda: get intercepts
    # par.ls$thresh.lambda <- c(par.ls$thresh.lambda[1], diff(par.ls$thresh.lambda)) # moved up does not work
    #
    # #get gamma : threshold params
    # #  get ini 1#
    # ini.gamma<-par.ls$thresh.gamma
    # #  get ini 2#
    # if (length(par.ls$thresh.gamma)) {
    #   tmpg <- gmat(par.ls$thresh.gamma)
    #   tmpg2 <- apply(tmpg,2,diff)
    #   tmpg3 <- as.vector(rbind(tmpg[1,],tmpg2))
    #   names(tmpg3)<-names(par.ls$thresh.gamma)
    #   par.ls$thresh.gamma <- tmpg3
    # }
    #
    # if (thresh.fun!='identity') {
    #   if (thresh.method == 'jurges'){
    #     # update lambda according to Jurges method
    #     par.ls$thresh.lambda <- c(par.ls$thresh.lambda[1],log(par.ls$thresh.lambda[2:length(par.ls$thresh.lambda)]))
    #
    #     # calculate gamma
    #
    #     if (length(par.ls$thresh.gamma)) {
    #       fn <- function(x) sum((alpha.thresh-gotm_Threshold(x[seq_len(model$J-1)],x[-seq_len(model$J-1)],model)[,-c(1,model$J+1)])^2)
    #       oo<-optim(par=c(par.ls$thresh.lambda,par.ls$thresh.gamma),fn=fn)
    #       for (i in 1:50) {
    #         oldv=oo$value
    #         oo<-optim(par=c(oo$par),fn=fn)
    #         cat(oo$value,',')
    #         if (abs(oo$value-oldv)<1e3) break
    #       }
    #       par.ls$thresh.lambda<-oo$par[seq_len(model$J-1)]
    #       par.ls$thresh.gamma<-oo$par[-seq_len(model$J-1)]
    #       cat('\n')
    #     }
    #   } else if (thresh.method == 'hopit'){
    #     par.ls$thresh.lambda <- log(par.ls$thresh.lambda[1:length(par.ls$thresh.lambda)])
    #     if (length(par.ls$thresh.gamma)) {
    #       for (k in 1:(model$J-1)){
    #         fn <- function(x) sum(abs(alpha.thresh[,k] - gotm_Threshold(subs(par.ls$thresh.lambda,k,x[1]),
    #                                                                   subs(par.ls$thresh.gamma,k,x[2]), model)[,k+1]))
    #         oo<-optim(par=c(par.ls$thresh.lambda[k],par.ls$thresh.gamma[k]),fn=fn)
    #         par.ls$thresh.lambda[k]<-oo$par[1]
    #         par.ls$thresh.gamma[k]<-oo$par[2]
    #       }
    #     }
    #   }
    # }
    # cat(' done\n')
    # model$vglm.start <- mv2.par
    # model$reg.start <- par.ls$reg.params
    # model$lambda.start <- par.ls$thresh.lambda
    # model$gamma.start <-par.ls$thresh.gamma
    # model$start <- c(model$reg.start, model$lambda.start, model$gamma.start)
    # #gotm_Threshold(par.ls$thresh.lambda, par.ls$thresh.gamma, model)
    # model$start.LL <- -gotm_negLL(model$start, model)
  }
  model
}

#as.matrix(coef(model.simple.exp))
#as.matrix(par.ls$thresh.gamma)
#' INTERNAL: Calculate a latent variable
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_Latent <- function(reg.params, model = NULL) model$reg.mm %*% (as.matrix(reg.params))

#' INTERNAL: Calculate maximum possible latent range
#' @param model a fitted \code{gotm} model.
#' @param data a data used to fit the model
#' @keywords internal
gotm_latentrange <- function (model, data) {
  cfm <- model$coef[seq_len(model$parcount[1])]
  ttr <- terms.formula(model$reg.formula)
  ttr <- delete.response(ttr)
  tt <- attr(ttr,'variables')
  ttn <- attr(ttr,'term.labels')
  li <- lapply(eval(tt, data), function(k) if (class(k) == 'factor') levels(k) else range(k, na.rm=TRUE))
  names(li) <- ttn
  # new.data <- expand.grid(li)
  # mm <- stats::model.matrix(ttr, data = new.data)[,-1]
  # V <- mm %*% as.matrix(cfm)
  # range(V)
  L=sapply(li, length)-1
  cfm_neg <- cfm * (cfm<0)
  cfm_pos <- cfm * (cfm>0)
  pos=c(0,cumsum(L))+1
  cfm_neg_ls <- lapply(1:(length(pos)-1),function(k) cfm_neg[pos[k]:(pos[k+1]-1)])
  cfm_pos_ls <- lapply(1:(length(pos)-1),function(k) cfm_pos[pos[k]:(pos[k+1]-1)])
  c(sum(sapply(cfm_neg_ls, min, na.rm=TRUE)),
    sum(sapply(cfm_pos_ls, max, na.rm=TRUE)))
  #must be checked for variables that are not factors
}

#' INTERNAL: Extract model parameters in a form of list
#'
#' Extract model parameters in a form of a list
#' @param model \code{gotm} object
#' @param parameters model parameters (optional). If not delivered then taken from \code{model$coef}
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
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
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_negLL <- function(parameters=model$coef, model, collapse = TRUE, use_weights = TRUE, negative = TRUE){
  z = gotm_c_init(model)
  if (collapse) {
    LLFunc(parameters, yi=as.numeric(unclass(model$y_i)),reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
           link=z$link,thresh_no_cov=model$thresh.no.cov*1,thresh_method=z$thresh_method,
           #thresh_func=z$thresh_fun,
           use_alpha = z$use_alpha, alpha_0=z$alpha_0, negative=1*negative,
           weights=model$weights,use_weights = 1*use_weights, out_val = model$control$LL_out_val)
  } else {
    LLFuncIndv(parameters, yi=as.numeric(unclass(model$y_i)),reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
               link=z$link,thresh_no_cov=model$thresh.no.cov*1,thresh_method=z$thresh_method,
               #thresh_func=z$thresh_fun,
               use_alpha = z$use_alpha, alpha_0=z$alpha_0, negative=1*negative,
               weights=model$weights,use_weights = 1*use_weights)
  }
}

#' Calculate special matrices for gradient calaculation
#' @keywords internal
calcYYY<-function(model){
  y <- as.numeric(unclass(model$y_i))
  dY <- Vector2DummyMat(y)
  YY2 <- (1-cumsum_row(dY))
  YY1 <- YY2 + dY
  YY2 <- YY2 + dY * (y == NCOL(YY1))
  YY1 <- YY1[, -NCOL(YY1)]
  YY2 <- YY2[, -NCOL(YY2)]
  YY1bound <- matrix(1 * (y!=max(y)), model$N, model$parcount[2])
  YY2bound <- matrix(1 * (y!=1), model$N, model$parcount[2])
  model$YYY1 <- YY1 * YY1bound
  model$YYY2 <- YY2 * YY2bound
  model
}

#' INTERNAL: The gradient of the log likelihood function
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_derivLL <- function(parameters=model$coef, model, collapse = TRUE, use_weights = TRUE, negative = FALSE){
  z = gotm_c_init(model)
  #be compatible with older versions
  if ((!length(model$YYY1)) || (!length(model$YYY1))){
    model <- calcYYY(model)
  }
  if (collapse) {
    LLGradFunc(parameters, yi=as.numeric(unclass(model$y_i)), YYY1=as.matrix(unname(model$YYY1)), YYY2=as.matrix(unname(model$YYY2)),
               reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
               link=z$link,thresh_no_cov=model$thresh.no.cov*1,thresh_method=z$thresh_method,
               #thresh_func=z$thresh_fun,
               use_alpha = z$use_alpha, alpha_0=z$alpha_0, negative=1*negative,
               weights=model$weights, use_weights = 1*use_weights)
  } else {
    LLGradFuncIndv(parameters, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2,
                   reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
                   link=z$link,thresh_no_cov=model$thresh.no.cov*1,thresh_method=z$thresh_method,
                   #thresh_func=z$thresh_fun,
                   use_alpha = z$use_alpha, alpha_0=z$alpha_0, negative=1*negative,
                   weights=model$weights, use_weights = 1*use_weights)
  }
}

#' INTERNAL: Fit \code{gotm}
#'
#' Fitter for Jurges and Hopit models
#' @param model \code{gotm} object
#' @param start starting parameters
#' @param control a list with control parameters. See \code{\link{gotm.control}}.
#' @importFrom maxLik maxNR
#' @importFrom numDeriv grad
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @keywords internal
gotm_fitter <- function(model, start = model$start, use_weights = TRUE){

  #be compatible with older versions
  if ((!length(model$YYY1)) || (!length(model$YYY1))){
    model <- calcYYY(model)
  }
  control <- model$control

  z <- gotm_c_init(model)

  LLgr <- function(par, neg=1) LLGradFunc(par, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2,
                                   reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
                                   link=z$link,thresh_no_cov=model$thresh.no.cov*1,thresh_method=z$thresh_method,
                                   #thresh_func=z$thresh_fun,
                                   use_alpha = z$use_alpha, alpha_0=z$alpha_0, negative=neg,
                                   weights=model$weights, use_weights = use_weights*1)
  LLfn <- function(par, neg=1) LLFunc(par, yi=as.numeric(unclass(model$y_i)),reg_mm=model$reg.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
                               link=z$link,thresh_no_cov=model$thresh.no.cov*1,thresh_method=z$thresh_method,
                               #thresh_func=z$thresh_fun,
                               use_alpha = z$use_alpha, alpha_0=z$alpha_0, negative=neg,
                               weights=model$weights,use_weights = use_weights*1, out_val = model$control$LL_out_val)

  refit <- function(fit, model){
    oldfit <- fit$value
    for (k in 1L : control$max.reiter) {
      if (k>0) try({fit <- optim(par = fit$par, fn = LLFunc, gr = LLgr,
                                 method = 'BFGS', hessian = FALSE, control=list(maxit=1000))}, silent = TRUE)

      fit <- optim(par = fit$par, fn = LLfn, gr = LLgr, hessian = FALSE, control=list(maxit=1000))

      lldiff <- abs(fit$value - oldfit)
      if (lldiff < control$tol.reiter) break
      oldfit <- fit$value
      if (k > 1L) cat(', ')
      cat(-fit$value,'(',lldiff,')',sep='')
    }
    list(fit = fit, converged = (k<control$max.reiter))
  }

  z <- try({
    nmfit <- optim(par = start, fn = LLfn)
    tmp <- refit(nmfit, model)
    nmFit <- tmp$fit
    fit <- nmFit
    if(!tmp$converged) {
      message('\nConvergence has not been reached yet (try to increase control$max.reiter), changing the method ...')
      control$fit.NR <- TRUE
    }
  }, silent = FALSE)
  if (class(z) == "try-error") stop('Impossible to find initial values.')
  if (control$fit.NR){

    nrFit <- maxNR(fn = LLfn, grad = LLgr, neg=0, start = fit$par) #maxLik::
    nrFit$par <- nrFit$estimate
    nrFit$value <- -nrFit$maximum
    nrFit <- refit(nrFit, model)
    if (!nrFit$converged)
      warning(call. = FALSE,
              '\nThe model probably did not converge. Try to increase "max.reiter" and/or set "forced.DEoptim" == TRUE.')
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
# @param thresh.fun function used to calculate thresholds. Possible values are \code{'exp'}(default) and \code{'identity'}.
#' @param alpha_0 value for "zero" threshold. If 'auto' then alpha_0 is set specificaly
#' to the 'jurges' or 'hopit' model. See \code{\link{gotm}}.
#' @seealso \code{\link{gotm}}
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @export
gotm.control<-function(fit.NR = FALSE,
                       max.reiter = 100L, #increase if faster optimizaton is developed
                       tol.reiter = 5e-5,
                       grad.eps = 1e-7,
                       LL_out_val = -Inf,
                       #thresh.fun = c('exp','identity','id','lin','linear'),
                       alpha_0 = 'auto'){

  # if (class(thresh.fun)=='character') {
  #   thresh.fun <- thresh.fun[1]
  #   if (tolower(thresh.fun) == 'exp') {
  #     thresh.fun <- exp
  #   } else if (tolower(thresh.fun) %in% c('identity', 'id')) {
  #     thresh.fun <- identity
  #   } else stop('Unknown threshold function.')
  # } else {
  #   if (!(identical(deparse(thresh.fun), deparse(exp)) || identical(deparse(thresh.fun), deparse(identity))))
  #     stop('Unknown threshold function.')
  # }
  # thresh.fun <- tolower(thresh.fun)[1]
  # if (thresh.fun %in% c('exp','identity','id','lin','linear')){
  #   if (thresh.fun %in% c('id','lin','linear','identity')) {
  #     thresh.fun == 'identity'
  #     stop('Currently only thresh.fun = "exp" is possible',call. = NULL)
  #   }
  # } else stop('Unknown threshold function.', call.=NULL)

  alpha_0 <- alpha_0[1]
  if (tolower(alpha_0) != 'auto'){
    #stop('Currently only alpha_0 = "auto" is possible.', call.=NULL)
    if (!is.numeric(alpha_0)) stop('"alpha_0" must be a numeric (-Inf or 0) or equal "auto".', call.=NULL) else if (alpha_0 %!in% c(-Inf,0)) stop('"alpha_0" can take only -Inf or 0 values.', call.=NULL)
  } else alpha_0 <- tolower(alpha_0)

  list(fit.NR = fit.NR, # thresh.fun = thresh.fun,
       max.reiter = max.reiter, tol.reiter = tol.reiter,
       grad.eps = grad.eps, alpha_0 = alpha_0, LL_out_val = LL_out_val)
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
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
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

#' Not %in% function
#'
#' @export
'%!in%' <-function(x, table) match(x, table, nomatch = 0L) == 0L

#' Check if one set is a subset of an another subset
#'
#' @export
'%c%' <-function(x, table) all(match(x, table, nomatch = 0L))

#' Not %c% function
#'
#' @export
'%!c%' <- function(x, table) !all(match(x, table, nomatch = 0L))

#' REMINDER: PUT HERE c++ FUNCTION!!
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
#' \item{\code{"jurges"} - Assumes that: \code{alpha(0) = -Inf}, \code{alpha(1)} is a
#' linear function of threshold covariates + \code{lambda}, \code{alpha(2..J-1) = thresh.fun(lambda, gamma, threshold covariates)},
#' and \code{alpha(J) = Inf}. See also [1]. Default option.}
#' \item{\code{"hopit"} - Assumes that: \code{alpha(0) = 0},  \code{alpha(1..J-1) = thresh.fun(lambda, gamma, threshold covarites)},
#' and \code{alpha(J) = Inf}. See also [2].}
#' }
#' @param start starting values in the form \code{c(latent_parameters, threshold_lambdas, threshold_gammas)}
#' @param doFit character string, \code{'full'} perform the full fit, \code{'vglm'} use starting values obteined from vglm,
#' @param hessain logical indicating if to calculate hessian matrix
#' but will not improve the fit, and \code{'no'} use external starting values, which must be delivered.
#' @param control a list with control parameters. See \code{\link{gotm.control}}.
#' @export
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
gotm<- function(reg.formula,
                thresh.formula = as.formula('~ 1'),
                data,
                survey = list(),
                link = c('probit', 'logit'),
                thresh.method = c('jurges', 'hopit','vglm'),
                start = NULL,
                doFit = c('full','vglm','no'),
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

  doFit <- match.arg(doFit)
  thresh.method <- match.arg(thresh.method)
  link <- match.arg(link)
  control <- do.call("gotm.control", control)
  survey <- do.call("gotm.design", survey)

  if (length(start) && class(start) == 'gotm'){
    if ((thresh.method != start$thresh.method) ||
        (link != start$link)) {
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

  # if (link == 'probit') {  #change to c++ or remove
  #   model$link.func <- pnorm
  #   model$distr.func <- dnorm
  # } else if (link == 'logit'){
  #   model$link.func <- function(x) exp(x)/(1L + exp(x))
  #   model$distr.func <- function (x) exp(-x)/((1L + exp(-x))^2L)
  # }

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
  model$thresh.method <- thresh.method[1]
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

  if (!length(start) || (doFit == 'vglm') || (model$thresh.method == 'vglm')) {  # was !=
    if (doFit == 'no') stop('Starting values must be given.', call.=NULL)
    z <- suppressWarnings(get.vglm.start(model, data))
    #z$start
    cat('VGLM logLik:',z$vglm.LL,'\n')
    if (model$thresh.method != 'vglm') cat('Recalculated VGLM (gotm start) logLik:',
                                           gotm_negLL(parameters = z$start, model=z, negative = FALSE),'\n')
    # z <- try({get.vglm.start(model, data)}, silent = FALSE)
    # if (class(z) == "try-error") stop('Initial values failed.')
    model <- z
  } else model$start <- start

  if ((doFit == 'full') && (model$thresh.method != 'vglm')){
    cat('Improving fit...')
    model <- gotm_fitter(model, start = model$start)
    cat(' done\nGotm logLik:', gotm_negLL(parameters = model$coef, model),'\n')
  } else {
    model$coef <- model$start
    if (model$thresh.method != 'vglm') model$LL <- gotm_negLL(parameters = model$start, model)
  }
  #add class vglm if thr.mth = vglm
  class(model) <- 'gotm'
  names(model$coef) <- coefnames
  colnames(model$thresh.mm) <- thresh.names
  p <- gotm_ExtractParameters(model)

  # try to calculate this for vglm anyway
  if (model$thresh.method != 'vglm') model$alpha <- gotm_Threshold(p$thresh.lambda, p$thresh.gamma, model) else model$alpha <- 'Not implemented'
  model$y_latent_i <- gotm_Latent(p$reg.params, model)
  if (model$thresh.method != 'vglm') {
    model$Ey_i <- factor(colSums(sapply(1L : model$N, function(k) model$alpha[k,]<model$y_latent_i[k])),levels=1L:model$J)
    levels(model$Ey_i) <- levels(model$y_i)
  } else {
    model$Ey_i <- 'Not implemented'
  }

  cat('Calculating maximum latent range...')
  if (model$thresh.method != 'vglm') k <- 1 else k <- -1
  model$maxlatentrange <- sort(k*gotm_latentrange(model=model, data=data))
  model$maxobservedlatentrange <-  sort(k*range(gotm_Latent(p$reg.params,model)))
  #if (class(model$maxlatentrange) == 'try-error') model$maxlatentrange=NA
  cat(' done\n')

  if (hessian && (model$thresh.method != 'vglm')) {
    cat('Calculating hessian...')
    # system.time(hes <- numDeriv::hessian(gotm_negLL, model$coef, model = model)) #numDeriv::)
    # system.time(hes2 <- my.grad(fn = gotm_derivLL, par = model$coef, model=model, eps = model$control$grad.eps, collapse = TRUE))
    # system.time(hes3 <- numDeriv::jacobian(func = gotm_derivLL, x = model$coef, model=model))
    # system.time(hes4 <- numDeriv::jacobian(func = gotm_derivLL, x = model$coef,
    #                                        model=model, method='simple',
    #                                        method.args=list(eps=model$control$grad.eps/2)))
    # dim(hes)
    # dim(hes2)
    # dim(hes3)
    # sum(abs(hes+hes4))
    # sum(abs(hes+hes2))
    # sum(hes+hes3)
    # sum(abs(hes2-hes4))
    hes <- my.grad(fn = gotm_derivLL, par = model$coef, model=model, eps = model$control$grad.eps, collapse = TRUE, negative=FALSE)
    model$vcov <- try(solve(-hes), silent = T)
    if (class(z) == 'try-error') {
      z <- NA*hes
      warning(call. = FALSE, 'Model is probably unidentifiable, vcov cannot be computed. Please try to use a "hopit" model.')
    }
    cat(' done\n')
    cat('Calculating estfun...')

    # model$estfun <- -my.grad(fn = gotm_negLL, par = model$coef,
    #                         eps = control$grad.eps, model = model, collapse = FALSE)

    model$estfun <- gotm_derivLL(model$coef, model, collapse = FALSE)

    #sum(model$estfun2-model$estfun)
    cat(' done\n')
  }
  return(model)
}

#' Extracting coefficients of fitted \code{gotm} object
#'
#' @param object \code{gotm} object.
#' @param aslist logical indicating if model coefficients should be returned as a list of three vectors
#' related to latent variable, threshold lambdas, and threshold gammas.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @keywords internal
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
#' @aliases coefficients.gotm
coef.gotm <- function(object, aslist = FALSE, ...){
  params <- object$coef
  if (aslist) gotm_ExtractParameters(object, params) else params
}

#' Printing basic information about fitted gotm
#'
#' @param x \code{gotm} object.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @keywords internal
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
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
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
vcov.gotm<-function(object, robust.vcov, control = list(), ...){
  #robust.method <- tolower(robust.method[1])
  #if (!(robust.method %in% c("grad","working"))) stop('Unknown method.')
  control <- do.call("gotm.control", control)
  z <- object$vcov
  if (!length(z)) stop('Hessian was not calculated.')
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
    if (length(object$design$FWeights)) divw <- object$design$FWeights else divw <- 1
    if (robust.vcov) z <- (z %*% t(object$estfun) %*% (object$estfun/divw) %*% (z))
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
  print.default(x)
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
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
summary.gotm <- function(object, robust.se = FALSE, control = object$control, ...){
  control <- do.call("gotm.control", control)
  if (object$thresh.method != 'vglm') {
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
  } else {
    if (length(object$design$PSU)) warning(call. = FALSE, 'The vglm method currently does not include survey design (except weights), SE may be biased.')
    if (robust.se) {
      warning(call. = FALSE, 'Calcualtion of robust SE is currently not supported for vglm method.')
      robust.se <- FALSE
    }
    varcov <- vcov(object$vglm)
    SE <- suppressWarnings(sqrt(diag(varcov)))
  }
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
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
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
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
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
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
AIC.gotm<-function(object, ..., k = 2L) {
  if (length(object$design$PSU)) warning(call. = FALSE, 'The AIC is currently not supported for survey design, the value may be biased.')
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
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
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
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
print.anova.gotm <- function(x, ...){
  cat('Anova (LRTs):\n')
  cat('Method: "', x$method, '"\n\n', sep = '')
  printCoefmat(x$table, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE, digits = 5L, dig.tst = 3L, tst.ind = 1L)
  invisible(NULL)
}

#' Likelihood ratio test for a pair of models
#'
#' @param full, nested {Models to be compared.}
#' @keywords internal
#' @export
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
lrt.gotm <- function(full, nested){
  if (length(full$design$PSU)) warning(call. = FALSE, 'The LRT is currently not supported for survey design, the value may be biased.')
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


#' Print object calculated by \code{\link{lrt.gotm}}
#'
#' @param x object obtained from \code{\link{lrt.gotm}}
#' @param ...	further arguments passed to or from other methods.
#' @keywords internal
#' @export
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
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
  out <- t(as.matrix(c('Chi^2' = unname(x$chisq), df = unname(x$df), 'Pr(>Chi^2)' = unname(x$pval))))
  row.names(out) <- ''
  printCoefmat(out, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE, digits = 5L, dig.tst = 3L, tst.ind = 1L)
  invisible(NULL)
}


#' @keywords internal
unravel <-function(mat, freq)  {
  mat <- cbind(mat, freq)
  FreqInd <- NCOL(mat)
  ls <- apply(mat,1,function(k) ((rep_row(as.matrix(k[-FreqInd]),k[FreqInd]))))
  do.call('rbind', ls)
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
#' @param unravelFreq logical indicating if to represent results on individual scale if FWeights were used.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @author Maciej J. Danko  <\email{maciej.danko@gmail.com}>
predict.gotm <- function(object, newdata=NULL, type = c('link', 'response', 'threshold', 'threshold_link'),
                         unravelFreq = TRUE, ...){
  if (length(newdata)) stop('"new data" not implemented.')
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
  return(H)
}


# factor.mat.gotm<-function(object, by.formula = object$thresh.formula){
#   if (length(by.formula)==3) by.formula[[2]] <- NULL
#   by.what<-ex.lb(by.formula)
#   rm.interaction <- grep(':',by.what,fixed=TRUE)
#   if (length(rm.interaction)) by.what <- by.what[-rm.interaction]
#   if (any(by.what %!in% c(ex.lb(object$reg.formula),ex.lb(object$thresh.formula)))) stop('Some elemnts of "by.formula" not present in the data.')
#   R. <- (object$reg.mm[,unlist(sapply(by.what, grep, x = colnames(object$reg.mm), fixed =TRUE))])
#   if(!length(R.)) R. <- NULL
#   T. <- object$thresh.mm[,unlist(sapply(by.what, grep, x = colnames(object$thresh.mm), fixed =TRUE))]
#   if (!length(T.)) T. <- NULL
#   if (any(colnames(T.) %in% colnames(R.))) T. <- T.[,-which(colnames(T.) %in% colnames(R.))]
#   vari <- cbind(T., R.)
#   get.reg.ref <- sapply(object$reg.lev,function(k) rownames(k)[1])[-1]
#   get.thresh.ref <- sapply(object$thresh.lev,function(k) rownames(k)[1])
#   get.ref <- c(get.reg.ref, get.thresh.ref)
#   miss.ref <- get.ref[names(get.ref)%in%by.what]
#   miss.col <- sapply(names(miss.ref), function(k) 1-rowSums(as.matrix(vari[,grep(k,colnames(vari),fixed=TRUE)])))
#   miss.nam <- paste(names(miss.ref),miss.ref,sep='')
#   colnames(miss.col) <- miss.nam
#   fullmat<-cbind(vari,miss.col)
#   fullmat<-fullmat[,order(colnames(fullmat))]
#   fullmat
# }

#' @export
plot.gotm<-function(x,...){
  tmp1<-predict.gotm(x, type='link',unravelFreq = FALSE)
  tmp2<-unclass(predict.gotm(x, type='response',unravelFreq = FALSE))
  tmp<-predict.gotm(x, type='threshold',unravelFreq = FALSE)
  K<-range(tmp,na.rm=TRUE,finite=TRUE)
  oo<-order(apply(cbind(tmp[,-c(1,x$J+1)])-K[1]+1,1,sum))
  #oo<-order(tmp2*tmp1)
  tmp<-tmp[oo,]
  tmp1<-tmp1[oo,]
  tmp2<-tmp2[oo]
  plot(NA,NA,xlim=c(0,NROW(tmp)),ylim=K,xlab='observation',ylab='latent variable')
  lines(tmp1,type='p',cex=0.2,col=adjustcolor(tmp2,alpha.f=0.5))
  for (k in 1:NCOL(tmp)) {
    lines(tmp[,k],col=k,lwd=2)
  }

}
