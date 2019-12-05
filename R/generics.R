#' Extracting the model coefficients
#'
#' Extract the model coefficients from the \code{hopit} model.
#' @param object a \code{hopit} object.
#' @param aslist a logical indicating whether the model coefficients should be returned
#' as a list of three vectors
#' related to latent variables, threshold lambdas, and threshold gammas.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @usage \method{coef}{hopit}(object, aslist = FALSE, ...)
#' @keywords internal
#' @author Maciej J. Danko
#' @aliases coefficients.hopit
#' @method coef hopit
coef.hopit <- function(object, aslist = FALSE, ...)
  if (aslist[1]) object$coef.ls else object$coef


#' Printing basic information about fitted \code{hopit} model
#'
#' Print a \code{hopit} model.
#' @param x a \code{hopit} object.
#' @param ...	further arguments passed to or from other methods.
#' @export
#' @usage \method{print}{hopit}(x, ...)
#' @method print hopit
#' @keywords internal
#' @author Maciej J. Danko
print.hopit<-function(x, ...){
  cat(hopit_msg(65), deparse(x$latent.formula), fill = TRUE)
  cat(hopit_msg(66), deparse(x$thresh.formula), fill = TRUE)
  cat(hopit_msg(72), x$link, fill = TRUE)
  cat(hopit_msg(73), x$N, fill = TRUE)
  cat(hopit_msg(74), toString(levels(x$y_i)), fill = TRUE)
  if (x$hasdisp) cat(hopit_msg(77), x$coef[length(x$coef)], fill = TRUE)
  cat(hopit_msg(78))
  print(x$coef.ls$latent.params)
  cat(hopit_msg(79))
  print(x$coef.ls$thresh.lambda)
  if(length(x$coef.ls$thresh.gamma)){
    cat(hopit_msg(80))
    print(x$coef.ls$thresh.gamma)
  }
  #if(length(x$coef.ls$logSigma)){
  if (x$hasdisp) cat(hopit_msg(82)) else cat(hopit_msg(81))
  cat(exp(x$coef.ls$logSigma),'\n')
  #}
  invisible(NULL)
}


#' Variance-covariance matrix from the fitted model
#'
#' Returns the variance-covariance matrix of the main parameters of a fitted \code{hopit} model object.
#' @param object a \code{hopit} object.
#' @param robust.vcov a logical indicating whether to use the sandwich estimator to
#' calculate the variance-covariance matrix.
#' If a survey design is detected, then this option is ignored.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @usage \method{vcov}{hopit}(object, robust.vcov, ...)
#' @keywords internal
#' @author Maciej J. Danko
#' @method vcov hopit
vcov.hopit<-function(object, robust.vcov, ...){
  z <- object$vcov
  if ("try-error" %in% class(z)) stop(paste(hopit_msg(37),
                                        attr(z,"condition"),sep=''),call.=NULL)
  if (!length(z)) stop(hopit_msg(38),call.=NULL)
  if (length(object$design)){
    if (!missing(robust.vcov) && (robust.vcov)) {
      warning(call. = FALSE, hopit_msg(39))
      robust.vcov <- FALSE
    }
  } else {
    if (missing(robust.vcov)) robust.vcov <- TRUE
    if (length(object$weights)) divw <- object$weights else divw <- 1
    #check how weights work here, they must be standardized.
    if (robust.vcov) z <- (z %*% t(object$estfun) %*%
                             (object$estfun/divw) %*% (z))
  }
  attr(z, 'survey.design') <- (length(object$design) > 0L)
  attr(z, 'robust.vcov') <- robust.vcov
  class(z) <- 'vcov.hopit'
  z
}


# Printing the variance-covariance matrix
#
# Print the variance-covariance matrix calculated by the \code{\link{vcov.hopit}}.
# @param x a \code{vcov.hopit} object
# @param digits see \code{\link{print.default}}.
# @param ... further arguments passed to or from other methods.
# @usage \method{print}{vcov.hopit}(x, digits = 3L, ...)
#' @noRd
#' @method print vcov.hopit
#' @keywords internal
#' @export
print.vcov.hopit <- function(x, digits = 3L, ...){
  cat(hopit_msg(40))
  print.default(x, digits = digits, ...)
  if (attr(x, 'survey.design')) cat(hopit_msg(41))
  if (!is.na(attr(x, 'robust.vcov')) && attr(x, 'robust.vcov'))
    cat(hopit_msg(42))
  invisible(NULL)
}


#' Calculate the model summary
#'
#' Summarize a \code{hopit}  model.
#' @param object a \code{hopit} object.
#' @param robust.se a logical indicating whether to use robust standard errors based
#' on the sandwich estimator.
#' If a survey design is detected, then this option is ignored.
#' @param control a list with control parameters.
#' See \code{\link{hopit.control}}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @author Maciej J. Danko
#' @useDynLib hopit
#' @usage \method{summary}{hopit}(object, robust.se, ...)
#' @method summary hopit
#' @keywords internal
#' @importFrom Rcpp evalCpp
summary.hopit <- function(object, robust.se, ...){
  if (missing(robust.se)) {
    if (length(object$design)) robust.se <- FALSE else robust.se <- TRUE
  }
  varcov <- vcov.hopit(object, robust.se, ...)
  robust.se <- attr(varcov, 'robust.vcov')
  dvcov <- diag(varcov)
  if (any(dvcov<0)) warning(hopit_msg(43),call.=NA)
  SE <- suppressWarnings(sqrt(abs(dvcov)))
  if (length(object$design)){
    message(hopit_msg(44))
  }
  if ((!robust.se) && (any(is.na(SE))) && !(length(object$design)))
    warning(call. = FALSE, hopit_msg(45))
  if (length(object$coef) != length(SE)) stop(hopit_msg(46),call.=NULL)
  tstat <-  object$coef/SE
  pvalue <- pstdnorm(-abs(tstat))  * 2L
  table1 <- data.frame(Estimate = object$coef, 'Std. Error' = SE,
                  'z value' = tstat, 'Pr(>|z|)' = pvalue, check.names = FALSE)
  tmp <- list(coef = table1, vcov = varcov, model = object,
              robust.se = robust.se)
  class(tmp) <- 'summary.hopit'
  tmp
}


# Print an object calculated by \code{\link{summary.hopit}}
#
# @param x an object created with \code{\link{summary.hopit}}.
# @param ...	further arguments passed to or from other methods.
#' @keywords internal
#' @export
# @usage \method{print}{summary.hopit}(x, ...)
#' @method print summary.hopit
# @author Maciej J. Danko
#' @noRd
print.summary.hopit <- function(x, ...){
  model <- x$model
  cat(hopit_msg(65), deparse(model$latent.formula), fill = TRUE)
  cat(hopit_msg(66), deparse(model$thresh.formula), fill = TRUE)
  cat(hopit_msg(72), model$link, fill = TRUE)
  cat(hopit_msg(73), model$N, fill = TRUE)
  cat(hopit_msg(74), toString(levels(model$y_i)), fill = TRUE)
  if(x$robust.se) cat(hopit_msg(71))
  cat('\n')
  stats::printCoefmat(x = x$coef, P.values = TRUE,
                      has.Pvalue = TRUE, digits = 4L, dig.tst = 2L)
  cat(hopit_msg(70), exp(model$coef.ls$logSigma), fill = TRUE)
  cat(hopit_msg(75), model$LL, fill = TRUE)
  cat(hopit_msg(76), model$deviance, fill = TRUE)
  if (!length(model$design)) cat('AIC:', AIC.hopit(model), fill = TRUE)
  cat('\n')
  invisible(NULL)
}

#' @keywords internal
#' @export
#' @method print hopitDW
#' @noRd
print.hopitDW <- function(x, digits=4, quote=FALSE, ...){
  x1 <- as.matrix(format(round(as.vector(x),digits),nsmall=digits, digits=digits,
         scientific = FALSE))
  dimnames(x1) <- dimnames(x)
  print.default(x1, quote=quote, digits=digits, ...)
}

#' @keywords internal
#' @export
#' @method summary hopitDW
#' @noRd
summary.hopitDW<-function(object, ...){
  x <- data.frame(cbind('Std. coef'=as.vector(object), attr(object,'Ptab')),
                  check.names = FALSE)
  class(x) <- c('summary.hopitDW','data.frame')
  attr(x,'legend') <- attr(object,'legend')
  return(x)
}

#' @keywords internal
#' @export
#' @method print summary.hopitDW
#' @noRd
print.summary.hopitDW <- function(x, digits=4, show.coef.names=TRUE, ...){
  if (show.coef.names[1]) {
    g <- data.frame('Coefficient name' = x[,2],
                  'Standardized coeficient'=
                    format(round(x[,1],digits),nsmall=digits, digits=digits,
                           scientific = FALSE),
                  'Pr(>|z|)' = format(round(x[,3],digits),
                         nsmall=digits, digits=digits, scientific = FALSE),
                  x[,4],
                  fix.empty.names = FALSE, check.names = FALSE)
  } else {
    g <- data.frame('Standardized coeficient'=
                      format(round(x[,1],digits), nsmall=digits, digits=digits,
                             scientific = FALSE),
                    'Pr(>|z|)' = format(round(x[,3],digits),
                          nsmall=digits, digits=digits, scientific = FALSE),
                    x[,4],
                    fix.empty.names = FALSE, check.names = FALSE)
  }
  rownames(g)<-rownames(x)
  print(g)
  cat("---\nSignif. codes: ",attr(x,'legend'),'\n')
}

#' Plotting standardized coefficients
#' @keywords internal
#' @export
#' @method plot hopitDW
#' @param x a object generated by \code{\link{standardizeCoef}} function.
#' @param ordered a logical indicating whether to sort the disability weights.
#' @param show.signif show significance codes.
#' @param mar,oma graphic parameters, see \code{\link{par}}.
#' @param ylab,xlab,density,angle,col,las,... arguments passed to \code{\link{barplot}}.
#' @importFrom graphics text
#' @author Maciej J. Danko
plot.hopitDW <- function(x,
                         ordered = TRUE,
                         show.signif = TRUE,
                         mar = c(10, 4, 1.5, 1),
                         oma = c(0, 0, 0, 0),
                         ylab = "Disability weight",
                         xlab = "",
                         density = 20,
                         angle = 45,
                         col = 'orange',
                         las = 3, ...){
  x <- summary.hopitDW(x)
  if (ordered[1]) {
    oz <- order(x[,1], decreasing = TRUE)
    x <- x[oz,]
  }
  if (length(mar) || length(oma)) {
    opar <- graphics::par(c("mar", "oma"))
    if (length(mar)) graphics::par(mar = mar)
    if (length(oma)) graphics::par(oma = oma)
  }
  rr <- graphics::barplot(x[,1], las = las, names.arg = rownames(x), ylab = ylab,
                          angle=angle, col=col, density=density, ...)
  if (show.signif[1]) graphics::text(rr,x[,1],paste('',x[,4]),srt=90, xpd=NA, adj=0, col='darkred')
}

#' Plotting Latent Index
#' @keywords internal
#' @export
#' @method plot hopitHI
#' @param x a object generated by \code{\link{latentIndex}} function.
#' @param response X-axis plotting option; choose \code{'data'} for the raw responses and \code{'fitted'} for the responses reclassified by the model.
#' @param xlab a label of the x-axis.
#' @param ylab a label of the y-axis.
#' @param ... further parameters passed to the \code{\link{plot}} function.
#' @author Maciej J. Danko
plot.hopitHI<-function(x,
                       response = c('data','fitted'),
                       xlab = '', ylab = 'Latent index', ...){
  mi <- attr(x, 'model.info')
  response <- tolower(match.arg(response))
  if (response=='data') YY <- mi$y_i else
    if (response=='fitted') YY <- mi$Ey_i else
      stop(hopit_msg(83),call.=NULL)
  graphics::plot(YY, x ,las=3, ylab=ylab, xlab=xlab,...)
}

#' @keywords internal
#' @export
#' @method print hopitHI
#' @noRd
print.hopitHI <- function(x, digits = 2, quote=FALSE, ...){
  print.default(format(round(as.vector(x), digits),
                       scientific=FALSE, digits=digits,
                       nsmall=digits), quote=quote, ...)
}

#' @keywords internal
#' @export
#' @method summary hopitHI
#' @noRd
summary.hopitHI <- function(object, breaks = seq(0,1,length.out = 21), ...){
  z<-cut(object, include.lowest = FALSE, breaks=breaks)
  z<-table(z)
  dimnames(z)<-unname(dimnames(z))
  z
}

#' Plotting Cut-Points
#' @keywords internal
#' @export
#' @method plot hopitCP
#' @param decreasing.levels a logical indicating whether self-reported health classes are ordered in decreasing order.
#' @param plotf a logical indicating whether to plot the results.
#' @param XLab,XLab.cex a label of the x axis and it's size.
#' @param YLab,YLab.cex a label of the y axis and it's size.
#' @param border.lwd,border.lty,border.col graphic parameters for vertical lines used to plot cut-points.
#' @param mar,oma graphic parameters, see \code{\link{par}}.
#' @param group.labels.type a position of the legend. One of \code{middle}, \code{border}, or \code{none}.
#' @param ... further plotting arguments.
#' @author Maciej J. Danko
plot.hopitCP<-function(x,
                       decreasing.levels=x$decreasing.levels,
                       mar=c(4,4,1,1),
                       oma=c(0,0,0,0),
                       XLab='Health index',
                       XLab.cex=1.1,
                       YLab='Counts',
                       YLab.cex=1.1,
                       border.col=2,
                       border.lty=2,
                       border.lwd=1.5,
                       group.labels.type=c('middle','border','none'),
                       ...){
  if (decreasing.levels[1]) dorev <- rev else dorev <- identity
  lv <- dorev(as.character(levels(attr(x,'y_i'))))
  Nm <- paste(lv[-length(lv)],lv[-1],sep=' | ')
  Nm <- sapply(Nm, function(k) bquote(bold(.(k))))
  group.labels.type<-tolower(group.labels.type[1])
  if (group.labels.type %notin%  c('middle','border','none'))
    stop(hopit_msg(84),call.=NULL)
  opar <- graphics::par(c('mar','oma'))
  graphics::par(mar=mar, oma=oma)

  z<-graphics::hist(x$h.index, 100,xlab='',ylab='' ,
                    main='', yaxs='i', col=grDevices::grey(0.4, alpha = 0.5),
                    border=grDevices::grey(0.4, alpha = 0.5))
  if (group.labels.type == 'border') {
    for (j in seq_along(Nm)) graphics::text(x=R1[j],y=(1.1*max(z$counts))/2,
                                            labels=Nm[[j]],
                                            srt=90,pos=2,offset=0.67,col=2)
  } else if (group.labels.type == 'middle'){
    R1 <- x$cutpoints
    R11=-diff(c(0,R1,1))/2+c(R1,1)+graphics::strheight('S',units='figure')/2
    for (j in seq_along(lv)) graphics::text(x=R11[j],
                                            y=(3*1.1*max(z$counts))/4,
                                            labels=lv[j],
                                            srt=90,
                                            pos=3,
                                            offset=0.67,
                                            col=2)
  }
  graphics::box()
  graphics::abline(v=R1,lwd=border.lwd,col=border.col,lty=border.lty)
  graphics::mtext(XLab, 1, cex=XLab.cex, line = 2.5)
  graphics::mtext(YLab, 2, cex=YLab.cex, line = 2.5)
  suppressWarnings(graphics::par(opar))
}

#' Plotting getLevels object
#' @keywords internal
#' @export
#' @method plot hopitLV
#' @param x aobject generated by \code{\link{getLevels}}.
#' @param mar,oma graphic parameters, see \code{\link{par}}.
#' @param YLab,YLab.cex a label for the y-axis and it's size.
#' @param legbg a legend background color. See \code{bg} parameter in \code{\link{legend}}.
#' @param legbty a legend box type. See \code{bty} parameter in \code{\link{legend}}.
#' @param las the style of axis labels, see \code{\link{par}}.
#' @param ... further plotting arguments.
#' @author Maciej J. Danko
plot.hopitLV<-function(x,
                      mar = c(7,2,1.5,0.5), oma = c(0,3,0,0),
                      YLab = 'Fraction [%]',
                      YLab.cex = 1.1,
                      legbg = grDevices::adjustcolor('white',alpha.f=0.4),
                      legbty = 'o',
                      las = 3,
                      ...){

    opar <- graphics::par(c('mar','oma','mfrow'))
    graphics::par(mfrow=c(1,2))
    graphics::par(mar=mar,oma=oma)
    graphics::barplot(t(x$original),las=las,main='Original')
    graphics::barplot(t(x$adjusted),las=las,main='Adjusted', legend.text=TRUE,
                      args.legend = list(x='center', box.col=NA,
                                         bg=legbg, bty=legbty))
    graphics::par(mfrow=c(1,1))
    graphics::par(mar=mar,oma=rep(0,4))
    graphics::mtext(YLab,2,cex=YLab.cex)
    suppressWarnings(graphics::par(opar))
}

#' Extracting a log likelihood of the fitted model
#'
#' Extract the log likelihood of a \code{hopit} model.
#' @param object a \code{hopit} object.
#' @param ... additional objects of the class \code{hopit}.
#' @keywords internal
#' @export
#' @usage \method{logLik}{hopit}(object, ...)
#' @author Maciej J. Danko
#' @method logLik hopit
logLik.hopit<-function(object, ...) {
  objects <- list(object, ...)
  tmp <- deparse(substitute(list(object, ...)))
  ob.nam <- gsub(' ', '', strsplit(substring(tmp, 6L,
                                  nchar(tmp) - 1L), ',', fixed = TRUE)[[1L]])
  res <- sapply(objects,function(object) object$LL)
  names(res) <- ob.nam
  res
}

#' Extracting the Akaike Information Criterion from the fitted model
#'
#' Extract the Akaike Information Criterion (AIC) from a fitted \code{hopit} model.
#' @param object a \code{hopit} object.
#' @param k a penalty per parameter to be used; the default k = 2 is the
#' classical AIC.
#' @param ... additional objects of the class \code{hopit}.
#' @keywords internal
#' @export
#' @usage \method{AIC}{hopit}(object, ..., k = 2L)
#' @method AIC hopit
#' @author Maciej J. Danko
AIC.hopit<-function(object, ..., k = 2L) {
  objects <- list(object, ...)
  tmp <- deparse(substitute(list(object, ...)))
  ob.nam <- gsub(' ', '', strsplit(substring(tmp, 6L, nchar(tmp) - 1L), ',',
                                   fixed = TRUE)[[1L]])
  res <- sapply(objects,function(object)
    if (!length(object$design)) object$AIC else
      stop(hopit_msg(47), call=NULL))
  names(res) <- ob.nam
  res
}

#' Likelihood Ratio Test Tables
#'
#' Perform the likelihood ratio test(s) for two or more \code{hopit} objects.
#' @param object an object containing the results returned by a \code{hopit}.
#' @param ...	an additional object(s) of the same type.
#' @param method the method of ordered model comparisons. Choose \code{"sequential"}
#' for 1-2, 2-3, 3-4, ... comparisons or
#' \code{"with.most.complex"} for 1-2, 1-3, 1-4, ... comparisons,
#' where 1 is the most complex model (the least complex for \code{"with.least.complex"}).
#' @param direction determine if the complexity of listed models is
#' \code{"increasing"} or \code{"decreasing"} (default).
# @keywords internal
#' @usage \method{anova}{hopit}(object, ..., method = c("sequential",
#' "with.most.complex", 'with.least.complex'),
#' direction = c("decreasing", "increasing"))
#' @method anova hopit
#' @export
#' @return a vector or a matrix with the results of the test(s).
#' @author Maciej J. Danko
#' @seealso \code{\link{print.lrt.hopit}},
#' \code{\link{lrt.hopit}}, \code{\link{hopit}}.
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # the order of response levels decreases from the best health to
#' # the worst health; hence the hopit() parameter decreasing.levels
#' # is set to TRUE
#' levels(healthsurvey$health)
#'
#' # Example 1 ---------------------
#' \donttest{
#' # fitting two nested models
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # a model with an interaction between hypertension and high_cholesterol
#' model2 <- hopit(latent.formula = health ~ hypertension * high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # a likelihood ratio test
#' lrt1 <- anova(model1, model2)
#' lrt1
#'
#' # print results in a shorter form
#' print(lrt1, short = TRUE)
#'
#' # or equivalently
#' lrt.hopit(model2, model1)
#' }
#' # Example 2 ---------------------
#' \donttest{
#' # fitting additional nested models
#' model3 <- hopit(latent.formula = health ~ hypertension * high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese * diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' model4 <- hopit(latent.formula = health ~ hypertension * high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese * diabetes + other_diseases,
#'               thresh.formula = ~ sex * ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # sequential likelihood ratio tests
#' # model complexity increases so direction = "increasing"
#' anova(model1, model2, model3, model4,
#'       direction = "increasing", method = "sequential")
#'
#' # likelihood ratio tests of the most complex model with the rest of the models
#' anova(model1, model2, model3, model4,
#'       direction = "increasing", method = "with.most.complex")
#'
#' # likelihood ratio tests of the least complex model with the rest of the models
#' anova(model1, model2, model3, model4,
#'       direction = "increasing", method = "with.least.complex")
#' }
anova.hopit<-function(object, ..., method = c('sequential', 'with.most.complex',
                                              'with.least.complex'),
                      direction = c('decreasing', 'increasing')){

  method <- tolower(match.arg(method))
  direction <- tolower(match.arg(direction))
  if (length(list(object, ...)) > 1L) {
    objects <- list(object, ...)
    tmp <- deparse(substitute(list(object, ...)))
    ob.nam <- gsub(' ', '', strsplit(substring(tmp, 6L, nchar(tmp) - 1L), ',',
                                     fixed = TRUE)[[1L]])
  } else  stop(hopit_msg(48))
  if (length(objects) == 2L){
    if(length(objects[[1L]]$coef)+objects[[1L]]$hasdisp >
       length(objects[[2L]]$coef)+objects[[2L]]$hasdisp) {
      return(lrt.hopit(objects[[1L]], objects[[2L]]))
    } else {
      return(lrt.hopit(objects[[2L]], objects[[1L]]))
    }
  } else {
    out <- NULL
    rna <- NULL
    if (direction == 'increasing') {
      objects <- rev(objects)
      ob.nam  <- rev(ob.nam)
    } else if (direction != 'decreasing')
      stop(call.=NULL, hopit_msg(50))
    for (k in 1L : (length(objects) - 1L)) {
      if (method == 'sequential'){
        tmp <- lrt.hopit(objects[[k]], objects[[k + 1L]])
        rna <- c(rna, paste(ob.nam[k], 'vs.', ob.nam[k + 1L], sep = ' '))
      } else if (method == 'with.most.complex') {
        tmp <- lrt.hopit(objects[[1L]], objects[[k + 1L]])
        rna <- c(rna, paste(ob.nam[1L], 'vs.', ob.nam[k + 1L], sep = ' '))
      } else if (method == 'with.least.complex') {
        tmp <- lrt.hopit(objects[[k]], objects[[length(objects)]])
        rna <- c(rna, paste( ob.nam[k], 'vs.',
                             ob.nam[length(objects)], sep = ' '))
      }
      out <- rbind(out, c('Chi^2' = tmp$chisq, df = tmp$df,
                          'Pr(>Chi^2)' = tmp$pval))
    }
    rownames(out) <- rna
    if (direction == 'increasing') out <- out[dim(out)[1L] : 1L,]
  }
  out <- list(table = out, objets = objects, names = ob.nam, method = method)
  class(out) <- 'anova.hopit'
  out
}


# Print an object calculated by \code{\link{anova.hopit}}
#
# @param x an object generated by \code{\link{anova.hopit}}.
# @param ...	further arguments passed to or from other methods.
#' @keywords internal
#' @export
# @author Maciej J. Danko
# @usage \method{print}{anova.hopit}(x, ...)
# @seealso \code{\link{anova.hopit}}, \code{\link{hopit}}.
#' @method print anova.hopit
#' @noRd
print.anova.hopit <- function(x, ...){
  cat(hopit_msg(49), x$method, '"\n\n', sep = '')
  stats::printCoefmat(x$table, signif.stars = TRUE, P.values = TRUE,
                    has.Pvalue = TRUE, digits = 5L, dig.tst = 3L, tst.ind = 1L)
  invisible(NULL)
}


#' Likelihood ratio test for a pair of models
#'
#' @param full,nested models to be compared.
#' @keywords internal
#' @return a vector with the results of the test.
#' @export
#' @author Maciej J. Danko
#' @seealso \code{\link{print.lrt.hopit}}, \code{\link{anova.hopit}},
#' \code{\link{hopit}}.
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
#' # Example 1 ---------------------
#'
#' # fitting two nested models
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # model with an interaction between hypertension and high_cholesterol
#' model2 <- hopit(latent.formula = health ~ hypertension * high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # Likelihood ratio test
#' lrt1 <- lrt.hopit(full = model2, nested = model1)
#' lrt1
#'
#' # print the results in a shorter form
#' print(lrt1, short = TRUE)
#'
#' # equivalently
#' print(anova(model2, model1), short = TRUE)
#' }
lrt.hopit <- function(full, nested){
  if (!identical(full$design, nested$design)) stop(hopit_msg(51),call. = NULL)
  if (length(full$coef) + full$hasdisp <= length(nested$coef)+ nested$hasdisp)
    stop(hopit_msg(52),call. = NULL)
  if (abs(full$LL - nested$LL) < .Machine$double.eps^0.45) {
    message(hopit_msg(53))
  } else if (full$LL - nested$LL < -.Machine$double.eps^0.45){
    warning(call. = FALSE, hopit_msg(54))
  }
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
    if (!(all(colnames(nested$latent.mm) %in% colnames(full$latent.mm))))
      warning(call. = FALSE, hopit_msg(56))
  if ((ncol(full$thresh.mm)) &&  (ncol(nested$thresh.mm)))
    if (!(all(colnames(nested$thresh.mm) %in% colnames(full$thresh.mm))))
      warning(call. = FALSE, hopit_msg(57))

  stat <- 2L*( logLik.hopit(full) - logLik.hopit(nested))

  if (!length(full$design)) {
    df.diff <- length(full$coef.ls$latent.params) -
      length(nested$coef.ls$latent.params) +
      length(full$coef.ls$thresh.lambda) -
      length(nested$coef.ls$thresh.lambda) +
      length(full$coef.ls$thresh.gamma) -
      length(nested$coef.ls$thresh.gamma) +
      (full$hasdisp) - (nested$hasdisp)
    p <- 1L - stats::pchisq(stat, df.diff)
  } else {
    stop(hopit_msg(58), call=NULL)
  }

  z <- list(chisq = stat, df = df.diff, pval = p, full = full, nested = nested)
  class(z) <- 'lrt.hopit'
  z
}


#' Printing an object calculated by \code{\link{lrt.hopit}}
#'
#' @param x an object obtained from \code{\link{lrt.hopit}}.
#' @param short a logical indicating whether to show a shortened description.
#' @param ...	further arguments passed to or from other methods.
#' @keywords internal
#' @export
#' @usage \method{print}{lrt.hopit}(x, short = FALSE, ...)
#' @author Maciej J. Danko
#' @method print lrt.hopit
#' @seealso \code{\link{lrt.hopit}}, \code{\link{anova.hopit}},
#' \code{\link{hopit}}.
print.lrt.hopit <- function(x, short = FALSE, ...){
  if (!short[1]) {
    cat(hopit_msg(64))
    cat("--", hopit_msg(65), deparse(x$full$latent.formula), fill = TRUE)
    cat("--",hopit_msg(66), deparse(x$full$thresh.formula), fill = TRUE)
    cat("--",hopit_msg(70),x$full$hasdisp, fill=TRUE)
    cat(hopit_msg(67))
    cat("--",hopit_msg(65), deparse(x$nested$latent.formula), fill = TRUE)
    cat("--",hopit_msg(66), deparse(x$nested$thresh.formula), fill = TRUE)
    cat("--",hopit_msg(70),x$nested$hasdisp, fill=TRUE)
    cat('\n')
  }
  cat(hopit_msg(59))
  # if (length(x$df)) {
  out <- t(as.matrix(c('Chi^2' = unname(x$chisq), df = unname(x$df),
                       'Pr(>Chi^2)' = unname(x$pval))))
  #  out2 <- NULL
  # else {
  #   out <- t(as.matrix(c('Chi^2' = unname(x$chisq),
  #                        'Pr(>Chi^2)' = unname(x$pval))))
  #   out2 <- x$scalef
  # }
  row.names(out) <- ''
  stats::printCoefmat(out, signif.stars = TRUE, P.values = TRUE,
                  has.Pvalue = TRUE, digits = 5L, dig.tst = 3L, tst.ind = 1L)
  # if (length(out2)) print(paste(hopit_msg(63),out2))
  invisible(NULL)
}


#' Calculate the log likelihood profile for the fitted \code{hopit} model
#'
#' @param fitted a \code{hopit} object (a fitted model).
#' @param scope a value (fraction) defining the plotting range for a coefficient.
#' The range is \code{c(coef \* (1-scope), coef \* (1+scope))}.
#' @param steps at how many equally spaced points the log likelihood
#' function is calculated for each coefficient.
#' @param ... unused now.
#' @export
#' @keywords internal
#' @author Maciej J. Danko
#' @usage \method{profile}{hopit}(fitted, ..., scope = 0.15, steps = 101)
#' @method profile hopit
#' @seealso \code{\link{plot.profile.hopit}}, \code{\link{print.profile.hopit}},
#'  \code{\link{hopit}}
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
#' # Example 1 ---------------------
#'
#' # fitting the model:
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # check the fit using the profile function (at 51 points)
#' pr <- profile(model1, steps = 51)
#' print(pr, plotf = FALSE)
#'
#' # plot profile
#' plot(pr, relative = FALSE)
#'
#' # alternative plot
#' plot(pr, relative = TRUE)
#' }
profile.hopit<-function(fitted, ..., scope=0.15, steps=101){
  steps <- floor(steps/2)*2+1
  if (fitted$hasdisp) COEF <- c(fitted$coef, fitted$coef.ls$logSigma) else
    COEF <- fitted$coef
  sub <- function(x,y) if (x==1) c(y,COEF[2:length(COEF)]) else
    if (x==length(COEF)) c(COEF[-length(COEF)],y) else
       c(COEF[1:(x-1)],y,COEF[(x+1):length(COEF)])
  lo <- COEF*(1-scope)
  hi <- COEF*(1+scope)
  GG <- function(x) sapply(seq(lo[x],hi[x],length.out=steps),function(y)
    hopit_negLL(parameters=sub(x,y),fitted,negative = FALSE))
  val <- sapply(seq_along(COEF), function(x) GG(x))
  attr(val,'scope') <- scope
  attr(val,'steps') <- steps
  attr(val,'lo') <- lo
  attr(val,'hi') <- hi
  colnames(val) <- names(COEF)
  class(val) <- c("profile.hopit", "profile")
  val
}


#' Plot the log likelihood profile for a profile.hopit object
#'
#' Plot the method for a profile.hopit object.
#' @param x a \code{profile.hopit} object.
#' @param leg.cex a character expansion factor relative to the current
#' \code{par("cex")} (see \code{\link{legend}}).
#' @param leg.col a color used for the legend text.
#' @param ylim see \code{\link{plot}}.
#' @param relative a logical indicating whether \code{ylim} on each panel should be the
#'  same (\code{TRUE}) or not (\code{FALSE}).
#' @param ... arguments to be passed to the \code{\link{plot}}() function (see
#' \code{\link{par}}).
#' @export
#' @keywords internal
#' @usage \method{plot}{profile.hopit}(x, ..., ylim = NULL, relative = FALSE,
#' leg.cex = 0.85, leg.col = 'blue4')
#' @author Maciej J. Danko
#' @method plot profile.hopit
#' @seealso \code{\link{profile.hopit}}, \code{\link{print.profile.hopit}},
#' \code{\link{hopit}}.
plot.profile.hopit<-function(x, ..., ylim = NULL, relative = FALSE,
                             leg.cex = 0.85, leg.col = 'blue4'){
  z <- sqrt(ncol(x))
  zy <- round(z)
  zx <- ceiling(z)
  spar <- graphics::par(c('mfrow','mar'))
  graphics::par(mfrow=c(zx,zy),mar=c(0,0,0,0))
  if (relative[1]) ylim <- range(x)
  for (j in seq_len(ncol(x))) {
    graphics::plot(x[,j],type='l',axes='F', ylim = ylim, ...)
    graphics::abline(v=floor(nrow(x)/2)+1,col=2,lty=2)
    graphics::legend('bottom',colnames(x)[j], bty='n',cex=leg.cex,
                     text.col=leg.col)
    graphics::box()
  }
  suppressWarnings(graphics::par(spar))
}


#' Print method for a profile.hopit object
#'
#' @param x a \code{profile.hopit} object.
#' @param plotf a logical indicating whether to plot the profile.
#' @param ... arguments to be passed to the \code{plot}() function (see
#' \code{\link{plot.profile.hopit}}).
#' @export
#' @keywords internal
#' @usage \method{print}{profile.hopit}(x, ..., plotf = TRUE)
#' @author Maciej J. Danko
#' @method print profile.hopit
#' @seealso \code{\link{profile.hopit}}, \code{\link{plot.profile.hopit}},
#' \code{\link{hopit}}
print.profile.hopit<-function(x, ..., plotf = TRUE){
  test <- apply(x,2,which.max)==floor(nrow(x)/2)+1
  if(plotf[1]) plot.profile.hopit(x, ...)
  if (any(!test)) {
    message(hopit_msg(60))
    message(paste(hopit_msg(61),paste(names(test)[!test],sep='',
                                      collapse = ',  ')))
  } else {
    cat(hopit_msg(62))
  }
}


#' Extract the \code{Sigma} parameter from a \code{hopit} model
#'
#' Extract the \code{Sigma} parameter from a \code{hopit} model
#' @param model a fitted \code{hopit} model.
#' @usage \method{sigma}{hopit}(object, ...)
#' @keywords internal
#' @export
#' @author Maciej J. Danko
#' @method sigma hopit
sigma.hopit <- function(object, ...)
  unname(exp(object$coef.ls$logSigma))
