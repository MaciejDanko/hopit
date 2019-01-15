#' Calculate latent index
#' @description
#' Calcualte latent index from the latent variable. It takes values from 0 to 1, where
#' zero is prescribed to the worse possible state (maximum posible value fo the latent variable;
#' all conditions/diseases with negative effects on health are present) and 1 is prescribed
#' to the best possible health (calculated analogically).
# WRONG DESCRIPTION
#' @param model a fitted \code{hopit} model.
#' @param subset an optional vector specifying a subset of observations.
#' @param plotf logical indicating if to plot summary figure.
#' @param crude logical indicating if to calculate crude health measure based directly on self reported health levels.
#' @param healthlevelsorder order of self-reported healh levels. Possible values are \code{'increasing'} and \code{'decreasing'}
#' @param response X axis plotting option, choose \code{'data'} for raw responses and \code{'fitted'} for model reclassified responses
# #' @param method a method of calcualtion of of latent range and standardization of. For \code{'observed'} (default) health index is
# #' standardzized using the minimum and maximum observed latent variable, whereas for \code{'theoretical} the minimum of latent range is
# #' deffined as latent variable for a hypothetical individual heaving no "disabilities", while the maximum of latent range is defined as a
# #' latent variable for an hypothetical individual having all possible disabilities.
#' @param YLab a label of y axis.
#' @export
latentIndex <- function(model, subset=NULL, plotf = FALSE, crude = FALSE, #scaled = TRUE,
                        healthlevelsorder = c('decreasing','increasing'), #put this also in hopit and change to decreased if necessary,
                        #here it will be removed and read from model$responselevelorder, or even not
                        response = c('data','fitted'),
                        #method = c('observed','theoretical'), #leave only observed option in the future
                        YLab = 'Latent index') {
  #0 is the worse possible health, 1 is the best possible health
  #method <- tolower(method[1])
  response <- tolower(response[1])
  healthlevelsorder <- tolower(healthlevelsorder[1])
  if (response=='data') YY <- model$y_i else if (response=='fitted') YY <- model$Ey_i else stop('Unknown response')
  p <- hopit_ExtractParameters(model)
  #if (!length(model$maxlatentrange)) model$maxlatentrange <- sort(hopit_latentrange(model=model, data=data))
  if (!length(model$maxobservedlatentrange)) model$maxobservedlatentrange <- sort(range(hopit_Latent(p$reg.params,model)))
  if (length(subset) == 0) subset=seq_along(YY)
  if (crude) {
    hi <- as.numeric(unclass(YY))
    hi <- (hi - min(hi))/(diff(range(hi)))
    if (healthlevelsorder == 'decreasing') hi <- 1-hi else
      if (healthlevelsorder != 'increasing') stop('Unknown value for healthlevelsorder.')
    hi <- hi[subset]
    if (plotf) plot(YY[subset], hi,las=3, ylab=YLab)
  } else {
    #if (method=='theoretical') r <- model$maxlatentrange else
      #if (method=='observed')
    r <- model$maxobservedlatentrange #else
    #    stop ('Unknown method.',call.=NULL)
    hi <- (1 - ((model$y_latent_i - r[1]) / diff(r)))[subset]
    if (plotf) plot(YY[subset], hi,las=3, ylab=YLab)
  }
  if (plotf) invisible(hi) else return(hi)
}

#' @rdname latentIndex
healthIndex <- latentIndex

#' Standardization of coefficients
#' @description
#' Calculate standardized coeficients - disability weights
#' computed as the regression parameters from the generalised ordered probit model divided by the maximum possible range
#' of its linear prediction. The range is calcualted as difference between maximum and minimum possible value of the latent variable
#' given estimated parameters.
#' @seealso \code{\link{healthIndex}}
#' @param model a fitted \code{hopit} model.
# #' @param latent.method see \code{\link{latentindex}}.
#' @param plotf logical indicating if to plot results.
#' @param plotpval logical indicating if to plot p-values.
#' @param mar see \code{\link{par}}.
#' @param oma see \code{\link{par}}.
#' @param YLab,YLab.cex labale and size of the label for y axis.
#' @param namesf one argument function that modifies names of coeficients.
#' @param ... arguments passed to \code{\link{boxplot}}.
#' @name standardizeCoef
#' @export
standardizeCoef <- function (model, # latent.method = c('observed','theoretical'),
                               plotpval = FALSE,
                               plotf = TRUE, mar = c(15, 4, 1, 1), oma = c(0, 0, 0, 0),
                               YLab = "Disability weight",
                               YLab.cex = 1.25,
                               namesf = identity, ...) {
  #latent.method <- latent.method[1]
  p <- hopit_ExtractParameters(model)
  #if (!length(model$maxlatentrange)) stop('maxlatentenrange field not present in the model.') #model$maxlatentrange <- sort(k*hopit_latentrange(model=model, data=data))
  if (!length(model$maxobservedlatentrange)) stop('maxobservedlatentenrange field not present in the model.') #model$maxobservedlatentrange <- sort(k*range(hopit_Latent(p$reg.params,model)))
  z <- (model$coef)[seq_len(model$parcount[1])]

  if (class(namesf)=='function') names(z) <- namesf(names(z)) else names(z) <- namesf
  oz <- order(z, decreasing = TRUE)
  cfm <- z[oz]
  r <- model$maxobservedlatentrange

  #if (latent.method=='theoretical') r <- model$maxlatentrange else if (latent.method=='observed') r <- model$maxobservedlatentrange

  #if (!force.positive) {
    #res <- as.matrix((cfm - r[1])/diff(r))
    res <- as.matrix((cfm)/diff(r))
  #} else #if (force.positive) {
  #  res <- as.matrix((cfm - max(r[1],0))/min(diff(r),max(r)))
  #} else stop('Unknown method.', call. = FALSE)

  if (plotf) {

    opar <- par(c("mar", "oma"))
    par(mar = mar, oma = oma)
    rr <- barplot(t(res), las = 3, ...)
    if(plotpval) {
      y <- summary(model)
      pval <- format(round(y$table$`Pr(>|z|)`[seq_len(model$parcount[1])],4),digits=4,scientific=FALSE)[oz]
      yr <- res/2
      ind <- yr < max(res)*0.1
      yr[ind] <- (res+max(res)*0.1)[ind]
      text(rr,yr,paste('P =',pval),srt=90,col=c('white','black')[1+ind])
    }
    mtext(YLab, 2, cex = YLab.cex, line = 2.5)
    suppressWarnings(par(opar))
  }
  if (plotf)
    invisible(res)
  else return(res)
}

#'@rdname standardizeCoef
standardiseCoef<-standardizeCoef

#'@rdname standardizeCoef
disabilityWeights<-standardizeCoef

# #' @keywords internal
# rowplot <- function (M,xlim,ylim,xlab='',ylab='',
#                      labels=rownames(M),col=seq_len(NCOL(M)),
#                      pch=19,cex=1.5,...){
#   if (missing(xlim)) xlim <- c(floor(min(M,na.rm=TRUE)*10)/10,ceiling(max(M,na.rm=TRUE)*10)/10)
#   if (missing(ylim)) ylim <- c(1,NROW(M))
#   plot(NA,NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE,...)
#   abline(h=seq_len(NROW(M)),lty=3)
#   for (k in seq_len(NROW(M))) lines(M[k,],rep(k,NCOL(M)),type='p',col=col,pch=pch,cex=cex,...)
#   axis(1)
#   axis(2,at=seq_len(NROW(M)),labels=labels,las=1)
#   box()
# }


#' Calcualte threshold cut-points using Jurges' method
#' @description
#' Calcualte threshold cut-points using Jurges' method.
#' @param model a fitted \code{hopit} model.
#' @param subset an optional vector specifying a subset of observations.
#' @param plotf logical indicating if to plot the results.
#' @param revf logical indicating if self-reported health classes are ordered in increasing order. ?????
#' @param healthlevelsorder order of self-reported healh levels. Possible values are \code{'increasing'} and \code{'decreasing'}
#' @param mar see \code{\link{par}}.
#' @param oma see \code{\link{par}}.
#' @export
getCutPoints <- function(model, subset=NULL, plotf = TRUE, mar=c(4,4,1,1),oma=c(0,0,0,0),
                         XLab='Health index', XLab.cex=1.25, YLab='Counts', YLab.cex=1.25,
                           revf=NULL, group.labels.type=c('middle','border','none')){
  if (length(subset) == 0) subset=seq_along(model$y_i)
  Y <- model$y_i[subset]
  if (length(revf) == 0) {
    message(paste('Are the levels',toString(levels(Y)),'sorted in incresing order?'))
    message('if not please set revf to TRUE.')
    stop('"revf" must be given',call.=FALSE)
  }
  if (revf) dorev <- rev else dorev <- identity
  h.index <- healthIndex(model, subset, plotf = FALSE)
  tY <- table(Y)
  tmp <- dorev(as.vector(tY))
  invcs <- (cumsum(tmp)/sum(tmp))[-length(tmp)]
  R1 <- quantile(h.index, invcs)

  lv <- dorev(as.character(levels(model$y_i)))
  Nm <- paste(lv[-length(lv)],lv[-1],sep=' | ')
  Nm <- sapply(Nm, function(k) bquote(bold(.(k))))
  if (plotf) {
    group.labels.type<-tolower(group.labels.type[1])
    if(group.labels.type %notin%  c('middle','border','none')) stop ('Unknown group.labels.type.',call.=NULL)
    opar <- par(c('mar','oma'))
    par(mar=mar, oma=oma)
    z<-hist(h.index, 200,xlab='',ylab='' ,
            main='', yaxs='i', col=grey(0.8, alpha = 0.5),border=grey(0.4, alpha = 0.5))
    if (group.labels.type == 'border') {
    for (j in seq_along(Nm)) text(x=R1[j],y=(1.1*max(z$counts))/2,labels=Nm[[j]],
                                  srt=90,pos=2,offset=0.67,col=2)
    } else if (group.labels.type == 'middle'){
      R11=-diff(c(0,R1,1))/2+c(R1,1)+strheight('S',units='figure')/2
      for (j in seq_along(lv)) text(x=R11[j],y=(3*1.1*max(z$counts))/4,labels=lv[j],
                                    srt=90,pos=3,offset=0.67,col=2)
    }
    box()
    abline(v=R1,lwd=2,col=2)
    mtext(XLab, 1, cex=XLab.cex, line = 2.5)
    mtext(YLab, 2, cex=YLab.cex, line = 2.5)
    suppressWarnings(par(opar))
  }
  if (length(h.index)){
    CIN <- c(0,R1,1)
    if (anyNA(CIN)){
      ind=seq_along(CIN)[is.na(CIN)]
      for (k in ind) CIN[ind]=CIN[ind-1]
    }
    if (anyDuplicated(CIN)) {
      warning('Some cut-points are NA or duplicated.')
      CIN <- CIN + cumsum(duplicated(CIN)*CIN/1e7)
      CIN <- CIN / max(CIN)
    }
    adjused.health.levels<- cut(h.index, CIN,labels= dorev(levels(Y)))
  } else adjused.health.levels <- NA
  res <- list(cutpoints=R1, adjused.health.levels=adjused.health.levels)
  if (plotf) invisible(res) else return(res)
}

#' Calcualte adjusted health levels.
#' @description
#' Calcualte adjusted health levels according to th Jurges' method.
#' @param model a fitted \code{hopit} model.
#' @param formula a formula containing the variables. It is by default set to threshold formula.
#' @param data data used to fit the model.
#' @param plotf logical indicating if to plot the results.
#' @param sep separator for levels names.
#' @param revf logical indicating if self-reported health classes are ordered in increasing order.
#' @param mar see \code{\link{par}}.
#' @param oma see \code{\link{par}}.
#' @export
getLevels<-function(model, formula=model$thresh.formula,
                          data=environment(model$thresh.formula), revf = NULL,
                          sort.flag=FALSE,
                          plotf = TRUE, sep='_',mar=c(7,2,1.5,0.5),oma=c(0,3,0,0)){
  if (class(formula)=='formula') inte_ <- formula2classes(formula, data, sep=sep, return.matrix = TRUE) else stop('Not implemented.')
  inte <- inte_$x
  namind <- inte_$class.mat
  nam <- levels(inte)
  cpall<-getCutPoints(model, plotf = FALSE, revf = revf)
  TAB1 <- round(table(original=model$y_i, adjusted=cpall$adjused.health.levels)*100/length(model$y_i),2)
  tmp <- untable(t(table(factor(model$y_i,levels=levels(cpall$adjused.health.levels)), inte)))
  N1 <- tmp
  tmp <-tmp/rowSums(tmp)
  if (sort.flag) {
    oD1 <- order(tmp[,NCOL(tmp)]+tmp[,NCOL(tmp)-1])
    tmp <- tmp[oD1,]
    orignalind <- namind[oD1,]
  } else orignalind <- namind

  tmp2 <- untable(t(table(cpall$adjused.health.levels, inte)))
  N2 <- tmp2
  tmp2 <- tmp2/rowSums(tmp2)
  if (sort.flag) {
    oD2 <- order(tmp2[,NCOL(tmp2)]+tmp2[,NCOL(tmp2)-1])
    tmp2 <- tmp2[oD2,]
    adjustedind <- namind[oD2,]
  } else adjustedind <- namind

  if (plotf) {
    opar <- par(c('mar','oma'))
    par(mfrow=c(1,2))
    par(mar=mar,oma=oma)
    barplot(t(tmp),las=3,main='Original')
    barplot(t(tmp2),las=3,main='Adjusted', legend.text=TRUE,
            args.legend = list(x='center', box.col=NA,
                               bg=adjustcolor('white',alpha.f=0.4)))
    par(mfrow=c(1,1))
    par(mar=mar,oma=rep(0,4))
    mtext('Fraction [%]',2,cex=1.5)
    suppressWarnings(par(opar))
  }
  res <- list(original= tmp,
              adjusted= tmp2,
              tab= TAB1,
              N.original= N1,
              N.adjusted= N2,
              I.orignal= orignalind,
              I.adjusted= adjustedind,
              mat=cbind(inte_$mat,
                      original= model$y_i,
                      adjusted= cpall$adjused.health.levels))
  class(res) <- 'healthlevels'
  if (plotf) invisible(res) else return(res)
}


