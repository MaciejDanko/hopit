#' Extracting coefficients of fitted \code{gotm} object
#'
#' @param object \code{gotm} object.
#' @param aslist logical indicating if model coefficients should be returned as a list of three vectors
#' related to latent variable, threshold lambdas, and threshold gammas.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @keywords internal
#' @author Maciej J. Danko
#' @aliases coefficients.gotm
coef.gotm <- function(object, aslist = FALSE, ...)  if (aslist) object$coef.ls else object$coef


#' Printing basic information about fitted gotm
#'
#' @param x \code{gotm} object.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @keywords internal
#' @author Maciej J. Danko
print.gotm<-function(x, ...){
  cat("Formula (latent variables):", deparse(x$reg.formula), fill = TRUE)
  cat("Formula (threshold variables):", deparse(x$thresh.formula), fill = TRUE)
  cat('Link:', x$link, fill = TRUE)
  cat('Number of cases:', x$N, fill = TRUE)
  cat('Response levels:', toString(levels(x$y_i)), fill = TRUE)
  cat('Number of thresholds:', x$J - 1, fill = TRUE)
  cat('\nCoefficients (latent variables):\n')
  print(model$coef.ls$reg.params)
  cat('\nThreshold coefficents (Lambda):\n')
  print(model$coef.ls$thresh.lambda)
  if(length(model$coef.ls$thresh.gamma)){
    cat('\nThreshold coefficients (Gamma):\n')
    print(model$coef.ls$thresh.gamma)
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
#' @importFrom survey svyrecvar
#' @export
#' @author Maciej J. Danko
vcov.gotm<-function(object, robust.vcov, ...){
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
    if (missing(robust.vcov)) robust.vcov <- FALSE
    if (length(object$weights)) divw <- object$weights else divw <- 1
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
#' @author Maciej J. Danko
summary.gotm <- function(object, robust.se = FALSE, ...){

  varcov <- vcov(object, robust.se, ...)
  SE <- suppressWarnings(sqrt(diag(varcov)))
  if (length(object$design)){
    cat('Survey weights detected. Standard errors was adjusted for survey design.\n')
  }
  if ((!robust.se) && (any(is.na(SE))) && !(length(object$design)))
    warning(call. = FALSE, 'Problem with some standard errors, please try option "robust.se" == TRUE, nd consider to use the "hopit" model.')
  testse <- abs(SE/object$coef)
  testse <- testse[!is.na(testse)]
  #if (any((testse > 50L)&(SE > 20L))) warning(call. = FALSE, 'Huge standard errors may suggest a problem with model identifiability.')

  #varcov <- vcov(object$vglm)
  #SE <- suppressWarnings(sqrt(diag(varcov)))

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
#' @author Maciej J. Danko
print.summary.gotm <- function(x, ...){
  model <- x$model
  cat("Formula (latent variables):", deparse(model$reg.formula), fill = TRUE)
  cat("Formula (threshold variables):", deparse(model$thresh.formula), fill = TRUE)
  cat('\nLink:', model$link, fill = TRUE)
  cat('Number of cases:', model$N, fill = TRUE)
  cat('Response levels:', toString(levels(model$y_i)), fill = TRUE)
  cat('Number of thresholds:', model$J - 1, fill = TRUE)
  if(x$robust.se) cat('\nRobust SE were used (sandwich estimator of varcov).\n')
  cat('\n')
  printCoefmat(x = x$table, P.values = TRUE, has.Pvalue = TRUE, digits = 4L, dig.tst = 2L)
  cat('\nLog-likelihood:', model$LL, fill = TRUE)
  cat('\nDeviance:', model$deviance, fill = TRUE)
  cat('AIC:', AIC.gotm(model), fill = TRUE)
  cat('\n')
  invisible(NULL)
}


#' Extracts log likelihood of the fitted model
#'
#' @param object \code{gotm} object.
#' @param ...	additional objects of the same type.
#' @keywords internal
#' @export
#' @author Maciej J. Danko
logLik.gotm<-function(object, ...) {
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
#' @param object \code{gotm} object.
#' @param k	numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @param ...	additional objects of the same type.
#' @keywords internal
#' @export
#' @author Maciej J. Danko
AIC.gotm<-function(object, ..., k = 2L) {
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
#' Compute likelihood rato test for two or more \code{gotm} objecs.
#' @param object an object containing the results returned by a \code{gotm}.
#' @param ...	additional objects of the same type.
#' @param method the method of model comparison. Choose \code{"sequential"} for 1-2, 2-3, 3-4, ... comparisons or
#' \code{"with.first"} for 1-2, 1-3, 1-4, ... comparisons.
#' @param direction determine if complexity of listed models is \code{"increasing"} or \code{"decreasing"} (default).
# @keywords internal
#' @export
#' @author Maciej J. Danko
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
#' @author Maciej J. Danko
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
#' @author Maciej J. Danko
lrt.gotm <- function(full, nested){

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
  if ((ncol(full$reg.mm)) &&  (ncol(nested$reg.mm)))
    if (!(all(colnames(nested$reg.mm) %in% colnames(full$reg.mm)))) warning(call. = FALSE, 'Models use probably different (non-nested) data sets (latent variable formula).')
  if ((ncol(full$thresh.mm)) &&  (ncol(nested$thresh.mm)))
    if (!(all(colnames(nested$thresh.mm) %in% colnames(full$thresh.mm)))) warning(call. = FALSE, 'Models use probably different (non-nested) data sets (latent variable formula).')

  stat <- 2L*(logLik.gotm(full) - logLik.gotm(nested))
  df.diff <- length(full$coef) - length(nested$coef)

  if (!length(full$design)) {
    df.diff <- length(full$coef) - length(nested$coef)
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
  class(z) <- 'lrt.gotm'
  z
}


#' Print object calculated by \code{\link{lrt.gotm}}
#'
#' @param x object obtained from \code{\link{lrt.gotm}}
#' @param ...	further arguments passed to or from other methods.
#' @keywords internal
#' @export
#' @author Maciej J. Danko
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
  if (length(x$df.diff)) {
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
#' @author Maciej J. Danko
predict.gotm <- function(object, newdata=NULL, type = c('link', 'response', 'threshold', 'threshold_link'),
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
