#' Extracting coefficients of fitted \code{hopit} object
#'
#' @param object \code{hopit} object.
#' @param aslist logical indicating if model coefficients should be returned as a list of three vectors
#' related to latent variable, threshold lambdas, and threshold gammas.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @usage \method{coef}{hopit}(object, aslist = FALSE, ...)
#' @keywords internal
#' @author Maciej J. Danko
#' @aliases coefficients.hopit
coef.hopit <- function(object, aslist = FALSE, ...)  if (aslist) object$coef.ls else object$coef


#' Printing basic information about fitted hopit
#'
#' @param x \code{hopit} object.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @usage \method{print}{hopit}(x, ...)
#' @keywords internal
#' @author Maciej J. Danko
print.hopit<-function(x, ...){
  cat("Formula (latent variables):", deparse(x$reg.formula), fill = TRUE)
  cat("Formula (threshold variables):", deparse(x$thresh.formula), fill = TRUE)
  cat('Link:', x$link, fill = TRUE)
  cat('Number of cases:', x$N, fill = TRUE)
  cat('Response levels:', toString(levels(x$y_i)), fill = TRUE)
  if (x$hasdisp) cat('Dispersion parameter (Theta):', x$coef[length(x$coef)], fill = TRUE)
  cat('\nCoefficients of the latent variable:\n')
  print(x$coef.ls$reg.params)
  cat('\nThreshold coefficents (Lambda):\n')
  print(x$coef.ls$thresh.lambda)
  if(length(x$coef.ls$thresh.gamma)){
    cat('\nThreshold coefficients (Gamma):\n')
    print(x$coef.ls$thresh.gamma)
  }
  if(length(x$coef.ls$theta)){
    if (x$hasdisp) cat('\nTheta:\n') else cat('\nFixed Theta:\n')
    cat(x$coef.ls$theta,'\n')
  }
  invisible(NULL)
}


#' Extracting variance-covariance matrix from the fitted hopit
#'
#' @param object \code{hopit} object.
#' @param robust.vcov logical indicating if to use sandwich estimator to calculate variance-covariance matrix.
#' If survey deign is detected than this option is ignored.
#' @param control a list with control parameters. See \code{\link{hopit.control}}.
# @param robust.method method of calculation of log-likelihood gradient.
# Set \code{"grad"} (default) for numerical gradient or \code{"working"}
# for the method based on working residuals.
# The latter one is an experimental method so warning will apear.
# robust.method = c("grad","working")
#' @param ...	further arguments passed to or from other methods.
#' @importFrom survey svyrecvar
#' @export
#' @usage \method{vcov}{hopit}(object, robust.vcov, ...)
#' @keywords internal
#' @author Maciej J. Danko
vcov.hopit<-function(object, robust.vcov, ...){
  #robust.method <- tolower(robust.method[1])
  #if (!(robust.method %in% c("grad","working"))) stop('Unknown method.')
  z <- object$vcov
  if (class(z) == "try-error") stop(paste('Cannot compute variance-covariance matrix:\n',attr(z,"condition"),sep=''),call.=NULL)
  if (!length(z)) stop('Hessian was not calculated.')
  if (length(object$design)){
    if (!missing(robust.vcov) && (robust.vcov)) {
      warning(call. = FALSE, '"robust.vcov" ignored, survey design was detected.')
      robust.vcov <- NA
    }
  } else {
    if (missing(robust.vcov)) robust.vcov <- TRUE
    if (length(object$weights)) divw <- object$weights else divw <- 1
    if (robust.vcov) z <- (z %*% t(object$estfun) %*% (object$estfun/divw) %*% (z)) #check how weights work here, they must be standardized.
  }
  attr(z, 'survey.design') <- (length(object$design) > 0L)
  attr(z, 'robust.vcov') <- robust.vcov
  class(z) <- 'vcov.hopit'
  z
}


#' Print object calculated by \code{\link{vcov.hopit}}
#'
#' @param x \code{hopit} object
#' @param ...	further arguments passed to or from other methods.
#' @usage \method{print}{vcov.hopit}(x, digits = 3L, ...)
#' @keywords internal
#' @export
print.vcov.hopit <- function(x, digits = 3L, ...){
  cat('Variance-covariance matrix:\n')
  print.default(x)
  if (attr(x, 'survey.design')) cat('\nVariance-covariance matrix adjusted for survey design.\n')
  if (!is.na(attr(x, 'robust.vcov')) && attr(x, 'robust.vcov')) cat('\nVariance-covariance matrix based on sandwich estimmator.\n')
  invisible(NULL)
}


#' Calculate model summary
#'
#' @param object \code{hopit} object.
#' @param robust.se logical indicating if to use robust standard errors based on the sandwich estimator.
#' If survey deign is detected than this option is ignored.
#' @param control a list with control parameters. See \code{\link{hopit.control}}.
# @param robust.method method of calculation of log-likelihood gradient.
# Set \code{"grad"} (default) for numerical gradient or \code{"working"}
# for the method based on working residuals.
# The latter one is an experimental method so warning will apear.
# robust.method = c("grad","working")
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @author Maciej J. Danko
#' @useDynLib hopit
#' @usage \method{summary}{hopit}(object, robust.se = TRUE, ...)
#' @keywords internal
#' @importFrom Rcpp evalCpp
summary.hopit <- function(object, robust.se = TRUE, ...){

  varcov <- vcov(object, robust.se, ...)
  SE <- suppressWarnings(sqrt(diag(varcov)))
  if (length(object$design)){
    cat('Survey weights detected. Standard errors was adjusted for survey design.\n')
  }
  if ((!robust.se) && (any(is.na(SE))) && !(length(object$design)))
    warning(call. = FALSE, 'Problem with some standard errors, please try option "robust.se" == TRUE.')
  if (length(object$coef) != length(SE)) stop('Something wrong.',call.=NULL)
  tstat <-  object$coef/SE
  pvalue <- pstdnorm(-abs(tstat))  * 2L
  table1 <- data.frame(Estimate = object$coef, 'Std. Error' = SE, 'z value' = tstat, 'Pr(>|z|)' = pvalue, check.names = FALSE)
  tmp <- list(coef = table1, vcov = varcov, model = object, robust.se = robust.se)
  class(tmp) <- 'summary.hopit'
  tmp
}


#' Print object calculated by \code{\link{summary.hopit}}
#'
#' @param x object created with \code{\link{summary.hopit}}
#' @param ...	further arguments passed to or from other methods.
#' @keywords internal
#' @export
#' @usage \method{print}{summary.hopit}(x, ...)
#' @author Maciej J. Danko
print.summary.hopit <- function(x, ...){
  model <- x$model
  cat("Formula (latent variable):", deparse(model$reg.formula), fill = TRUE)
  cat("Formula (threshold):", deparse(model$thresh.formula), fill = TRUE)
  cat('\nLink:', model$link, fill = TRUE)
  cat('Number of cases:', model$N, fill = TRUE)
  cat('Response levels:', toString(levels(model$y_i)), fill = TRUE)
  if (model$hasdisp) cat('Dispersion parameter (Theta):', model$coef[length(model$coef)], fill = TRUE)
  if(x$robust.se) cat('\nRobust SE were used (sandwich estimator of varcov).\n')
  cat('\n')
  printCoefmat(x = x$coef, P.values = TRUE, has.Pvalue = TRUE, digits = 4L, dig.tst = 2L)
  cat('\nTheta:', model$coef.ls$theta, fill = TRUE)
  cat('\nLog-likelihood:', model$LL, fill = TRUE)
  cat('\nDeviance:', model$deviance, fill = TRUE)
  if (!length(model$design)) cat('AIC:', AIC.hopit(model), fill = TRUE)
  cat('\n')
  invisible(NULL)
}


#' Extracts log likelihood of the fitted model
#'
#' @param object \code{hopit} object.
#' @param ...	additional objects of the same type.
#' @keywords internal
#' @export
#' @usage \method{logLik}{hopit}(object, ...)
#' @author Maciej J. Danko
logLik.hopit<-function(object, ...) {
  #if (length(object$design$PSU)) warning(call. = FALSE, 'The LogLik function is currently not supported for survey design, the value may be biased.')
  objects <- list(object, ...)
  tmp <- deparse(substitute(list(object, ...)))
  ob.nam <- gsub(' ', '', strsplit(substring(tmp, 6L, nchar(tmp) - 1L), ',', fixed = TRUE)[[1L]])
  res <- sapply(objects,function(object) object$LL)
  names(res) <- ob.nam
  res
}

#' Extracts Akaike Information Criterion from the fitted model
#'
#' @param object \code{hopit} object.
#' @param k	numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @param ...	additional objects of the same type.
#' @keywords internal
#' @export
#' @usage \method{AIC}{hopit}(object, ..., k = 2L)
#' @author Maciej J. Danko
AIC.hopit<-function(object, ..., k = 2L) {
  objects <- list(object, ...)
  tmp <- deparse(substitute(list(object, ...)))
  ob.nam <- gsub(' ', '', strsplit(substring(tmp, 6L, nchar(tmp) - 1L), ',', fixed = TRUE)[[1L]])
  res <- sapply(objects,function(object) if (!length(object$design)) object$AIC else
    stop('AIC for models with survey design not implemented yet.', call=NULL))
  names(res) <- ob.nam
  res
}

#' LRT Tables
#'
#' Compute likelihood rato test for two or more \code{hopit} objecs.
#' @param object an object containing the results returned by a \code{hopit}.
#' @param ...	additional objects of the same type.
#' @param method the method of model comparison. Choose \code{"sequential"} for 1-2, 2-3, 3-4, ... comparisons or
#' \code{"with.first"} for 1-2, 1-3, 1-4, ... comparisons.
#' @param direction determine if complexity of listed models is \code{"increasing"} or \code{"decreasing"} (default).
# @keywords internal
#' @usage \method{anova}{hopit}(object, ..., method = c("sequential", "with.first"),
#' direction = c("decreasing", "increasing"))
#' @export
#' @author Maciej J. Danko
anova.hopit<-function(object, ..., method = c('sequential', 'with.first'), direction = c('decreasing', 'increasing')){

  method <- match.arg(method)
  direction <- match.arg(direction)
  if (length(list(object, ...)) > 1L) {
    objects <- list(object, ...)
    tmp <- deparse(substitute(list(object, ...)))
    ob.nam <- gsub(' ', '', strsplit(substring(tmp, 6L, nchar(tmp) - 1L), ',', fixed = TRUE)[[1L]])
  } else  stop('At least two objects must be listed.')
  if (length(objects) == 2L){
    if(length(objects[[1L]]$coef)>length(objects[[2L]]$coef)) {
      return(lrt.hopit(objects[[1L]], objects[[2L]]))
    } else {
      return(lrt.hopit(objects[[2L]], objects[[1L]]))
    }
  } else {
    out <- NULL
    rna <- NULL
    if (direction == 'increasing') objects <- objects[length(objects) : 1L] else if (direction != 'decreasing') stop('Unknown direction.')
    for (k in 1L : (length(objects) - 1L)) {
      if (tolower(method) == 'sequential'){
        tmp <- lrt.hopit(objects[[k]], objects[[k + 1L]]) # try models mut be of decreasing complexity, silent = F
        rna <- c(rna, paste(ob.nam[k], 'vs.', ob.nam[k + 1L], sep = ' '))
      } else if (tolower(method) == 'with.first') {
        tmp <- lrt.hopit(objects[[1L]], objects[[k + 1L]]) # the first model must be the most complex,  silent = F
        rna <- c(rna, paste(ob.nam[1L], 'vs', ob.nam[k + 1L], sep = ''))
      } else
        out <- rbind(out, c('Chi^2' = tmp$chisq, df = tmp$df, 'Pr(>Chi^2)' = tmp$pval))
    }
    rownames(out) <- rna
    if (direction == 'increasing') out <- out[dim(out)[1L] : 1L,]
  }
  out <- list(table = out, objets = objects, names = ob.nam, method = method)
  class(out) <- 'anova.hopit'
  out
}


#' Print object calcuated by \code{\link{anova.hopit}}
#'
#' @param x object generated by \code{\link{anova.hopit}}
#' @param ...	further arguments passed to or from other methods.
#' @keywords internal
#' @export
#' @author Maciej J. Danko
#' @usage \method{print}{anova.hopit}(x, ...)
print.anova.hopit <- function(x, ...){
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
#' @author Maciej J. Danko
lrt.hopit <- function(full, nested){

  if (!identical(full$design, nested$design)) stop('Models have different survey designs.',call. = NULL)
  if (length(full$coef) <= length(nested$coef)) stop('The "full" model must have more parameters than the "nested" one.',call. = NULL)
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
  if ((full$hasdisp) < (nested$hasdisp)) stop('Theta params are not nested.')

  if ((ncol(full$reg.mm)) &&  (ncol(nested$reg.mm)))
    if (!(all(colnames(nested$reg.mm) %in% colnames(full$reg.mm)))) warning(call. = FALSE, 'Models use probably different (non-nested) data sets (latent variable formula).')
  if ((ncol(full$thresh.mm)) &&  (ncol(nested$thresh.mm)))
    if (!(all(colnames(nested$thresh.mm) %in% colnames(full$thresh.mm)))) warning(call. = FALSE, 'Models use probably different (non-nested) data sets (latent variable formula).')

  stat <- 2L*( logLik.hopit(full) - logLik.hopit(nested))
  #df.diff <- length(full$coef) - length(nested$coef) + length(full$coef.ls$theta) - length(nested$coef.ls$theta)

  if (!length(full$design)) {
    df.diff <- length(full$coef.ls$reg.params) - length(nested$coef.ls$reg.params) +
      length(full$coef.ls$thresh.lambda) - length(nested$coef.ls$thresh.lambda) +
      length(full$coef.ls$thresh.gamma) - length(nested$coef.ls$thresh.gamma) +
      (full$hasdisp) - (nested$hasdisp)
    p <- 1L - pchisq(stat, df.diff)
    scalef <- NULL
  } else {
    stop('LRT for models with survey design not implemented yet.', call=NULL)
    # ....
    # misspec <- eigen(solve(V0) %*% V, only.values = TRUE)$values
    #
    # p <- survey::pchisqsum(stat, rep(1, length(misspec)), misspec,
    #                        method = "sad", lower.tail = FALSE)
    # print(c(p1,p2))
    # df.diff <- NULL
    # scalef <- full$misspec/mean(full$misspec)
  }

  z <- list(chisq = stat, df = df.diff, pval = p, scalef<-scalef, full = full, nested = nested)
  class(z) <- 'lrt.hopit'
  z
}


#' Print object calculated by \code{\link{lrt.hopit}}
#'
#' @param x object obtained from \code{\link{lrt.hopit}}
#' @param ...	further arguments passed to or from other methods.
#' @keywords internal
#' @export
#' @usage \method{print}{lrt.hopit}(x, ...)
#' @author Maciej J. Danko
print.lrt.hopit <- function(x, ...){
  cat('Likelihood ratio test:\n')
  cat('full model:\n')
  cat("-- Formula (latent variables):", deparse(x$full$reg.formula), fill = TRUE)
  cat("-- Formula (threshold variables):", deparse(x$full$thresh.formula), fill = TRUE)
  cat("-- Estimated theta:",x$full$hasdisp, fill=TRUE)
  cat('\nnested model:\n')
  cat("-- Formula (latent variables):", deparse(x$nested$reg.formula), fill = TRUE)
  cat("-- Formula (threshold variables):", deparse(x$nested$thresh.formula), fill = TRUE)
  cat("-- Estimated theta:",x$nested$hasdisp, fill=TRUE)
  #uzyc signif
  cat('\nLikelihood ratio test:\n')
  if (length(x$df)) {
    out <- t(as.matrix(c('Chi^2' = unname(x$chisq), df = unname(x$df), 'Pr(>Chi^2)' = unname(x$pval))))
    out2 <- NULL
  } else {
    out <- t(as.matrix(c('Chi^2' = unname(x$chisq), 'Pr(>Chi^2)' = unname(x$pval))))
    out2 <- x$scalef
  }
  row.names(out) <- ''
  printCoefmat(out, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE, digits = 5L, dig.tst = 3L, tst.ind = 1L)
  if (length(out2)) print(paste('Scale factors',out2))
  invisible(NULL)
}


#' Model predictions
#'
#' @param object \code{hopit} object.
#' @param type the type of prediction required. The default \code{"link"}
#' is on the scale of linear predictors (latent variable). The alternative \code{"response"}
#' is on the scale of categorical response variable. The \code{"threshold"}
#' gives the thresholds for each observation, whereas the \code{"threshold_link"} gives meaningful thresholds
#' together with latent variable for each observation (a data.frame with fields \code{$left.boundary},
#' \code{$latent.variable}, and \code{$right.boundary}).
#' @param unravelFreq logical indicating if to represent results on individual scale if FWeights were used.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @usage \method{predict}{hopit}(object, newdata=NULL,
#' type = c('link', 'response', 'threshold', 'threshold_link'),
#' unravelFreq = TRUE, ...)
#' @author Maciej J. Danko
predict.hopit <- function(object, newdata=NULL, type = c('link', 'response', 'threshold', 'threshold_link'),
                         unravelFreq = TRUE, ...){
  if (length(newdata)) stop('"new data" not implemented yet.')
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

#' Calculate log likelihood profile for fitted hopit model
#'
#' @param fitted \code{hopit} object. Fitted model.
#' @param scope value (fraction) defining plotting range for a coeficient. The range is \code{c(coef \* (1-scope), coef \* (1+scope))}.
#' @param steps at how many equaly spaced points calcualte the log likelihood function for each coeficient.
#' @param ... unused now.
#' @export
#' @keywords internal
#' @author Maciej J. Danko
#' @usage \method{profile}{hopit}(fitted, ..., scope = 0.15, steps = 101)
profile.hopit<-function(fitted, ..., scope=0.15, steps=101){
  steps <- floor(steps/2)*2+1
  if (fitted$hasdisp) COEF <- c(fitted$coef, fitted$coef.ls$theta) else COEF <- fitted$coef
  sub <- function(x,y) if (x==1) c(y,COEF[2:length(COEF)]) else if (x==length(COEF)) c(COEF[-length(COEF)],y) else
    c(COEF[1:(x-1)],y,COEF[(x+1):length(COEF)])
  lo <- COEF*(1-scope)
  hi <- COEF*(1+scope)
  GG <- function(x) sapply(seq(lo[x],hi[x],length.out=steps),function(y) hopit_negLL(parameters=sub(x,y),fitted,negative = FALSE))
  val <- sapply(seq_along(COEF), function(x) GG(x))
  attr(val,'scope') <- scope
  attr(val,'steps') <- steps
  attr(val,'lo') <- lo
  attr(val,'hi') <- hi
  colnames(val) <- names(COEF)
  class(val) <- c("profile.hopit", "profile")
  val
}


#' Plot log likelihood profile for profile.hopit object
#'
#' Plot method for profile.hopit object.
#' @param x \code{profile.hopit} object.
#' @param leg.cex character expansion factor relative to current \code{par("cex")} (see \code{\link{legend}}).
#' @param leg.col the color used for the legend text.
#' @param ... arguments to be passed to \code{plot}() function (see \code{\link{par}}).
#' @export
#' @keywords internal
#' @usage \method{plot}{profile.hopit}(x, ..., leg.cex = 0.85, leg.col = 'blue4')
#' @author Maciej J. Danko
plot.profile.hopit<-function(x, ..., leg.cex = 0.85, leg.col = 'blue4'){
  z <- sqrt(ncol(x))
  zy <- round(z)
  zx <- ceiling(z)
  spar <- par(c('mfrow','mar'))
  par(mfrow=c(zx,zy),mar=c(0,0,0,0))
  for (j in seq_len(ncol(x))) {
    plot(x[,j],type='l',axes='F', ...)
    abline(v=floor(nrow(x)/2)+1,col=2,lty=2)
    legend('bottom',colnames(x)[j], bty='n',cex=leg.cex, text.col=leg.col)
    box()
  }
  suppressWarnings(par(spar))
}


#' Print method for profile.hopit object
#'
#' @param x \code{profile.hopit} object.
#' @param plotf \code{hopit} object.
#' @param ... arguments to be passed to \code{plot}() function (see \code{\link{plot.profile.hopit}}).
#' @export
#' @keywords internal
#' @usage \method{print}{profile.hopit}(x, ..., plotf = TRUE)
#' @author Maciej J. Danko
print.profile.hopit<-function(x, ..., plotf = TRUE){
  test <- apply(x,2,which.max)==floor(nrow(x)/2)+1
  if(plotf) plot.profile.hopit(x, ...)
  if (any(!test)) {
    message('Log likelihood maximum not reached.')
    message(paste('Problem in:',paste(names(test)[!test],sep='',collapse = ',  ')))
  } else {
    cat('All parameters seem to be at arg.max (optimum).\n')
  }
}

