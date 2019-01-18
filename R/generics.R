#' Extracting coefficients of the fitted \code{hopit} model
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
  cat(hopit_msg(65), deparse(x$latent.formula), fill = TRUE)
  cat(hopit_msg(66), deparse(x$thresh.formula), fill = TRUE)
  cat(hopit_msg(72), x$link, fill = TRUE)
  cat(hopit_msg(73), x$N, fill = TRUE)
  cat(hopit_msg(74), toString(levels(x$y_i)), fill = TRUE)
  if (x$hasdisp) cat(hopi_msg(77), x$coef[length(x$coef)], fill = TRUE)
  cat(hopit_msg(78))
  print(x$coef.ls$latent.params)
  cat(hopit_msg(79))
  print(x$coef.ls$thresh.lambda)
  if(length(x$coef.ls$thresh.gamma)){
    cat(hopit_msg(80))
    print(x$coef.ls$thresh.gamma)
  }
  #if(length(x$coef.ls$logTheta)){
    if (x$hasdisp) cat(hopit_msg(82)) else cat(hopit_msg(81))
    cat(exp(x$coef.ls$logTheta),'\n')
  #}
  invisible(NULL)
}


#' Extracting variance-covariance matrix from the fitted hopit
#'
#' @param object \code{hopit} object.
#' @param robust.vcov logical indicating if to use sandwich estimator to calculate variance-covariance matrix.
#' If survey deign is detected than this option is ignored.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @usage \method{vcov}{hopit}(object, robust.vcov, ...)
#' @keywords internal
#' @author Maciej J. Danko
vcov.hopit<-function(object, robust.vcov, ...){
  z <- object$vcov
  if (class(z) == "try-error") stop(paste(hopit_msg(37),attr(z,"condition"),sep=''),call.=NULL)
  if (!length(z)) stop(hopit_msg(38),call.=NULL)
  if (length(object$design)){
    if (!missing(robust.vcov) && (robust.vcov)) {
      warning(call. = FALSE, hopit_msg(39))
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
#' @param digits see \code{\link{print.default}}
#' @param ... further arguments passed to or from other methods.
#' @usage \method{print}{vcov.hopit}(x, digits = 3L, ...)
#' @keywords internal
#' @export
print.vcov.hopit <- function(x, digits = 3L, ...){
  cat(hopit_msg(40))
  print.default(x, digits = digits, ...)
  if (attr(x, 'survey.design')) cat(hopit_msg(41))
  if (!is.na(attr(x, 'robust.vcov')) && attr(x, 'robust.vcov')) cat(hopit_msg(42))
  invisible(NULL)
}


#' Calculate model summary
#'
#' @param object \code{hopit} object.
#' @param robust.se logical indicating if to use robust standard errors based on the sandwich estimator.
#' If survey deign is detected than this option is ignored.
#' @param control a list with control parameters. See \code{\link{hopit.control}}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @author Maciej J. Danko
#' @useDynLib hopit
#' @usage \method{summary}{hopit}(object, robust.se = TRUE, ...)
#' @keywords internal
#' @importFrom Rcpp evalCpp
summary.hopit <- function(object, robust.se = TRUE, ...){
  varcov <- vcov(object, robust.se, ...)
  dvcov <- diag(varcov)
  if (any(dvcov<0)) warning(hopit_msg(43),call.=NA)
  SE <- suppressWarnings(sqrt(abs(dvcov)))
  if (length(object$design)){
    cat(hopit_msg(44))
  }
  if ((!robust.se) && (any(is.na(SE))) && !(length(object$design)))
    warning(call. = FALSE, hopit_msg(45))
  if (length(object$coef) != length(SE)) stop(hopit_msg(46),call.=NULL)
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
  cat(hopit_msg(65), deparse(model$latent.formula), fill = TRUE)
  cat(hopit_msg(66), deparse(model$thresh.formula), fill = TRUE)
  cat(hopit_msg(72), model$link, fill = TRUE)
  cat(hopit_msg(73), model$N, fill = TRUE)
  cat(hopit_msg(74), toString(levels(model$y_i)), fill = TRUE)
  if(x$robust.se) cat(hopit_msg(71))
  cat('\n')
  printCoefmat(x = x$coef, P.values = TRUE, has.Pvalue = TRUE, digits = 4L, dig.tst = 2L)
  cat(hopit_msg(70), exp(model$coef.ls$logTheta), fill = TRUE)
  cat(hopit_msg(75), model$LL, fill = TRUE)
  cat(hopit_msg(76), model$deviance, fill = TRUE)
  if (!length(model$design)) cat('AIC:', AIC.hopit(model), fill = TRUE)
  cat('\n')
  invisible(NULL)
}


#' Extracts log likelihood of the fitted model
#'
#' @param object \code{hopit} object.
#' @param ... additional objects of the class \code{hopit}.
#' @keywords internal
#' @export
#' @usage \method{logLik}{hopit}(object, ...)
#' @author Maciej J. Danko
logLik.hopit<-function(object, ...) {
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
#' @param k a penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @param ... additional objects of the class \code{hopit}.
#' @keywords internal
#' @export
#' @usage \method{AIC}{hopit}(object, ..., k = 2L)
#' @author Maciej J. Danko
AIC.hopit<-function(object, ..., k = 2L) {
  objects <- list(object, ...)
  tmp <- deparse(substitute(list(object, ...)))
  ob.nam <- gsub(' ', '', strsplit(substring(tmp, 6L, nchar(tmp) - 1L), ',', fixed = TRUE)[[1L]])
  res <- sapply(objects,function(object) if (!length(object$design)) object$AIC else
    stop(hopit_msg(47), call=NULL))
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
  } else  stop(hopit_msg(48))
  if (length(objects) == 2L){
    if(length(objects[[1L]]$coef)+objects[[1L]]$hasdisp>length(objects[[2L]]$coef)+objects[[2L]]$hasdisp) {
      return(lrt.hopit(objects[[1L]], objects[[2L]]))
    } else {
      return(lrt.hopit(objects[[2L]], objects[[1L]]))
    }
  } else {
    out <- NULL
    rna <- NULL
    if (direction == 'increasing') objects <- objects[length(objects) : 1L] else if (direction != 'decreasing')
      stop(call.=NULL, hopit_msg(50))
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
  cat(hopit_msg(49), x$method, '"\n\n', sep = '')
  printCoefmat(x$table, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE, digits = 5L, dig.tst = 3L, tst.ind = 1L)
  invisible(NULL)
}

#' Likelihood ratio test for a pair of models
#'
#' @param full,nested Models to be compared.
#' @keywords internal
#' @export
#' @author Maciej J. Danko
lrt.hopit <- function(full, nested){
  if (!identical(full$design, nested$design)) stop(hopit_msg(51),call. = NULL)
  if (length(full$coef) + full$hasdisp <= length(nested$coef)+ nested$hasdisp) stop(hopit_msg(52),call. = NULL)
  if (abs(full$LL - nested$LL) < .Machine$double.eps^0.6) message(hopit_msg(53)) else
    if (full$LL - nested$LL < -.Machine$double.eps^0.6) warning(call. = FALSE, hopit_msg(54))
  if (ncol(full$latent.mm) < ncol(nested$latent.mm)) {
    cat(hopit_msg(64))
    cat("--",hopit_msg(65), deparse(full$latent.formula), fill = TRUE)
    cat(hopit_msg(67))
    cat("--",hopit_msg(65), deparse(nested$latent.formula), fill = TRUE)
    stop(hopit_msg(68))
  }
  if (ncol(full$thresh.mm) < ncol(nested$thresh.mm)) {
    cat(hopit_msg(64))
    cat("--",hopit_msg(66), deparse(full$thresh.formula), fill = TRUE)
    cat(hopit_msg(67))
    cat("--",hopit_msg(66), deparse(nested$thresh.formula), fill = TRUE)
    stop(hopit_msg(69))
  }
  if ((full$hasdisp) < (nested$hasdisp)) stop(hopit_msg(55))

  if ((ncol(full$latent.mm)) &&  (ncol(nested$latent.mm)))
    if (!(all(colnames(nested$latent.mm) %in% colnames(full$latent.mm)))) warning(call. = FALSE, hopit_msg(56))
  if ((ncol(full$thresh.mm)) &&  (ncol(nested$thresh.mm)))
    if (!(all(colnames(nested$thresh.mm) %in% colnames(full$thresh.mm)))) warning(call. = FALSE, hopit_msg(57))

  stat <- 2L*( logLik.hopit(full) - logLik.hopit(nested))

  if (!length(full$design)) {
    df.diff <- length(full$coef.ls$latent.params) - length(nested$coef.ls$latent.params) +
      length(full$coef.ls$thresh.lambda) - length(nested$coef.ls$thresh.lambda) +
      length(full$coef.ls$thresh.gamma) - length(nested$coef.ls$thresh.gamma) +
      (full$hasdisp) - (nested$hasdisp)
    p <- 1L - pchisq(stat, df.diff)
  } else {
    stop(hopit_msg(58), call=NULL)
  }

  z <- list(chisq = stat, df = df.diff, pval = p, full = full, nested = nested)
  class(z) <- 'lrt.hopit'
  z
}


#' Print object calculated by \code{\link{lrt.hopit}}
#'
#' @param x object obtained from \code{\link{lrt.hopit}}
#' @param short logical, shortened description
#' @param ...	further arguments passed to or from other methods.
#' @keywords internal
#' @export
#' @usage \method{print}{lrt.hopit}(x, short = FALSE, ...)
#' @author Maciej J. Danko
print.lrt.hopit <- function(x, short = FALSE, ...){
  if (!short) {
    cat(hopit_msg(64))
    cat("--", hopit_msg(65), deparse(x$full$latent.formula), fill = TRUE)
    cat("--",hopit_msg(66), deparse(x$full$thresh.formula), fill = TRUE)
    cat("--",hopit_msg(70),x$full$hasdisp, fill=TRUE)
    cat(hopit_msg(67))
    cat("--",hopit_msg(65), deparse(x$nested$latent.formula), fill = TRUE)
    cat("--",hopit_msg(66), deparse(x$nested$thresh.formula), fill = TRUE)
    cat("--",hopit_msg(70),x$nested$hasdisp, fill=TRUE)
  }
  #uzyc signif
  cat(hopit_msg(59))
  if (length(x$df)) {
    out <- t(as.matrix(c('Chi^2' = unname(x$chisq), df = unname(x$df), 'Pr(>Chi^2)' = unname(x$pval))))
    out2 <- NULL
  } else {
    out <- t(as.matrix(c('Chi^2' = unname(x$chisq), 'Pr(>Chi^2)' = unname(x$pval))))
    out2 <- x$scalef
  }
  row.names(out) <- ''
  printCoefmat(out, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE, digits = 5L, dig.tst = 3L, tst.ind = 1L)
  if (length(out2)) print(paste(hopit_msg(63),out2))
  invisible(NULL)
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
  if (fitted$hasdisp) COEF <- c(fitted$coef, fitted$coef.ls$logTheta) else COEF <- fitted$coef
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
#' @param ylim see \code{\link{plot}}
#' @param relative logical indicating if \code{ylim} on each panel should be the same (\code{TRUE}) or not (\code{FALSE})
#' @param ... arguments to be passed to \code{\link{plot}}() function (see \code{\link{par}}).
#' @export
#' @keywords internal
#' @usage \method{plot}{profile.hopit}(x, ..., leg.cex = 0.85, leg.col = 'blue4')
#' @author Maciej J. Danko
plot.profile.hopit<-function(x, ..., ylim = NULL, relative = FALSE, leg.cex = 0.85, leg.col = 'blue4'){
  z <- sqrt(ncol(x))
  zy <- round(z)
  zx <- ceiling(z)
  spar <- par(c('mfrow','mar'))
  par(mfrow=c(zx,zy),mar=c(0,0,0,0))
  if (relative) ylim <- range(x)
  for (j in seq_len(ncol(x))) {
    plot(x[,j],type='l',axes='F', ylim = ylim, ...)
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
    message(hopit_msg(60))
    message(paste(hopit_msg(61),paste(names(test)[!test],sep='',collapse = ',  ')))
  } else {
    cat(hopit_msg(62))
  }
}

