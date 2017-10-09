#' Calculate health index
#' @description
#' Calcualte halth index from the latent variable. It takes values from 0 (the worse possible health) to
#' 1 (the best possible health).
#' @param model a fitted \code{gotm} model.
#' @param subset an optional vector specifying a subset of observations.
#' @param plotf logical indicating if to plot summary figure.
#' @export
healthindex <- function(model, subset=NULL, plotf = FALSE) {
  #0 is the worse possible health, 1 is the best possible health
  if (length(subset)==0) subset=seq_along(model$y_i)
  cfm <- coef(model)[seq_len(model$parcount[1])]
  minc <- min(cfm, 0)
  Zs <- sum((cfm * (cfm>0))) - minc
  hi <- (1 - ((model$y_latent_i - minc) / Zs))[subset]
  if (plotf) plot(model$y_i[subset], hi,las=3, ylab='Health index')
  if (plotf) invisible(hi) else return(hi)
}

#' Calculate disability weights
#' @description
#' Calculate disability weights.
#' @param model a fitted \code{gotm} model.
#' @param plotf logical indicating if to plot results.
#' @param mar see \code{\link{par}}.
#' @param oma see \code{\link{par}}.
#' @export
disabilityweights <-function(model, plotf = TRUE, mar=c(15,4,1,1),oma=c(0,0,0,0)){
  cfm <- sort(coef(model)[seq_len(model$parcount[1])], decreasing = TRUE)
  minc <- min(cfm,0)
  Zs <- sum((cfm*(cfm>0))) - minc
  res <- as.matrix((cfm - minc) / Zs)
  if (plotf){
    opar <- par(c('mar','oma'))
    par(mar=mar,oma=oma)
    barplot(t(res),las=3)
    mtext('Disability weight',2,cex=1.5,line=2.5)
    suppressWarnings(par(opar))
  }
  if (plotf) invisible(res) else return(res)
}

#' @keywords internal
rowplot <- function (M,xlim,ylim,xlab='',ylab='',
                     labels=rownames(M),col=seq_len(NCOL(M)),
                     pch=19,cex=1.5,...){
  if (missing(xlim)) xlim <- c(floor(min(M,na.rm=TRUE)*10)/10,ceiling(max(M,na.rm=TRUE)*10)/10)
  if (missing(ylim)) ylim <- c(1,NROW(M))
  plot(NA,NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE,...)
  abline(h=seq_len(NROW(M)),lty=3)
  for (k in seq_len(NROW(M))) lines(M[k,],rep(k,NCOL(M)),type='p',col=col,pch=pch,cex=cex,...)
  axis(1)
  axis(2,at=seq_len(NROW(M)),labels=labels,las=1)
  box()
}

#' @keywords internal
untable <- function(x) {
  names(attr(x, "dimnames")) <- c('','')
  as.matrix(x)
}

#' @keywords internal
formula2classes <- function(formula, data, sep='_'){
  tmp <- model.frame(formula, data)
  colnames(tmp)
  lv <- lapply(seq_len(NCOL(tmp)),function (k) levels(tmp[,k]))
  names(lv) <-colnames(tmp)
  tmp2 <- expand.grid(lv)
  tmp2 <- sapply(seq_len(NCOL(tmp2)), function (k) paste(colnames(tmp2)[k],'[',tmp2[,k],']',sep=''))
  nlv <- levels(interaction(as.data.frame(tmp2),sep=sep))
  tmp <- sapply(seq_len(NCOL(tmp)), function (k) paste(colnames(tmp)[k],'[',tmp[,k],']',sep=''))
  tmp <- interaction(as.data.frame(tmp),sep=sep)
  tmp <- factor(tmp, levels=nlv)
  tmp
}

#' Get health index quantiles with respect to specified vaiables
#' @description
#' Get health index quantiles with respect to specified vaiables.
#' @param model a fitted \code{gotm} model.
#' @param formula a formula containing the variables. It is by default set to threshold formula.
#' @param data used to fit the model.
#' @param plotf logical indicating if to plot the results.
#' @param sep separator for levls names.
#' @param mar see \code{\link{par}}.
#' @param oma see \code{\link{par}}.
#' @export
gethealthindexquantiles<-function(model, formula=model$thresh.formula, data=environment(model$thresh.formula),
                                  plotf = TRUE, sep='_',
                                  mar=c(4,8,1.5,0.5),oma=c(0,0,0,0)){
  if (class(formula)=='formula') tmp <- formula2classes(formula, data, sep=sep) else stop('Not implemented.')
  D <- t(sapply(levels(tmp),function(k) quantile(healthindex(model, tmp==k))))
  D0 <- quantile(healthindex(model))
  if (plotf){
    opar <- par(c('mar','oma'))
    par(mar=mar,oma=oma)
    rowplot(D[NROW(D):1,-c(1,model$J)], pch=c(21,19,21), col=1)
    abline(v=D0[-c(1,model$J)])
    mtext('Health index',1,cex=1.5,line = 2.5)
    suppressWarnings(par(opar))
  }
  res <- list(q.all=D0, q=D)
  if (plotf) invisible(res) else return(res)
}

#' Calcualte threshold cut-points using Jurges' method
#' @description
#' Calcualte threshold cut-points using Jurges' method.
#' @param model a fitted \code{gotm} model.
#' @param subset an optional vector specifying a subset of observations.
#' @param plotf logical indicating if to plot the results.
#' @param revf logical indicating if self-reported health classes are ordered in increasing order.
#' @param mar see \code{\link{par}}.
#' @param oma see \code{\link{par}}.
#' @keywords internal
basiccutpoints <- function(model, subset=NULL, plotf = TRUE, mar=c(4,4,1,1),oma=c(0,0,0,0), revf=NULL){

  if (length(subset) == 0) subset=seq_along(model$y_i)
  Y <- model$y_i[subset]
  if (length(revf) == 0) {
    message(paste('Are the levels',toString(levels(Y)),'sorted in incresing order?'))
    message('if not please set revf to TRUE.')
    stop('"revf" must be given',call.=FALSE)
  }
  if (revf) dorev <- rev else dorev <- identity
  h.index <- healthindex(model, subset, plotf = FALSE)
  tY <- table(Y)
  tmp <- dorev(as.vector(tY))
  invcs <- (cumsum(tmp)/sum(tmp))[-length(tmp)]
  R1 <- quantile(h.index, invcs)
  lv <- dorev(as.character(levels(model$y_i)))
  Nm <- paste(lv[-length(lv)],lv[-1],sep=' | ')
  Nm <- sapply(Nm, function(k) bquote(bold(.(k))))
  if (plotf) {
    opar <- par(c('mar','oma'))
    par(mar=mar, oma=oma)
    z<-hist(h.index, 200,xlab='',ylab='' ,
            main='', yaxs='i', col=grey(0.8, alpha = 0.5),border=grey(0.4, alpha = 0.5))
    for (j in seq_along(Nm)) text(x=R1[j],y=(1.1*max(z$counts))/2,labels=Nm[[j]],
                                  srt=90,pos=2,offset=0.67,col=2)
    box()
    abline(v=R1,lwd=2,col=2)
    mtext('Health index',1,cex=1.5,line = 2.5)
    mtext('Counts',2,cex=1.5, line=2.5)
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
    adjused.health.levels<- cut(h.index,CIN,labels= dorev(levels(Y)))
  } else adjused.health.levels <- NA
  res <- list(cutpoints=R1, adjused.health.levels=adjused.health.levels)
  if (plotf) invisible(res) else return(res)
}

#' Calcualte threshold cut-points using Jurges' method
#' @description
#' Calcualte threshold cut-points using Jurges' method.
#' @param model a fitted \code{gotm} model.
#' @param formula a formula containing the variables.
#' It is by default set to threshold formula.
#' If set to \code{NULL} then threshold cut-point are calcualted for the whole population.
#' @param data data used to fit the model.
#' @param plotf logical indicating if to plot the results.
#' @param sep separator for levels names.
#' @param revf logical indicating if self-reported health classes are ordered in increasing order.
#' @param mar see \code{\link{par}}.
#' @param oma see \code{\link{par}}.
#' @export
getcutpoints<-function(model, formula=model$thresh.formula,
                       data=environment(model$thresh.formula),
                       plotf = TRUE, sep='_',
                       revf=NULL, mar, oma){
  if (!length(formula)) {
    if (missing(oma)) oma=c(0,0,0,0)
    if (missing(mar)) mar=c(4,4,1,1)
    res <- basiccutpoints(model, subset=NULL, plotf = plotf, mar = mar, oma = oma, revf = revf)
  } else {
    if (missing(oma)) oma=c(0,0,0,0)
    if (missing(mar)) mar=c(4,8,1.5,0.5)
    Y <- model$y_i
    if (length(revf) == 0) {
      message(paste('Are the levels',toString(levels(Y)),'sorted in incresing order?'))
      message('if not please set revf to TRUE.')
      stop('"revf" must be given',call.=FALSE)
    }
    if (revf) dorev <- rev else dorev <- identity
    if (class(formula)=='formula') tmp <- formula2classes(formula, data, sep=sep) else stop('Not implemented.')
    D <- t(sapply(levels(tmp),function(k) basiccutpoints(model=model, tmp==k, plotf = FALSE, revf = revf)$cutpoints))
    D0 <- basiccutpoints(model, plotf = FALSE, revf = revf)$cutpoints
    lv <- dorev(as.character(levels(model$y_i)))
    Nm <- paste(lv[-length(lv)],lv[-1],sep=' | ')
    colnames(D) <- Nm
    Nm <- sapply(Nm, function(k) bquote(bold(.(k))))
    if (plotf){
      opar <- par(c('mar','oma'))
      par(mar=mar, oma=oma)
      rowplot(D[NROW(D):1,], pch=19, col=1)
      abline(v=D0, col=2,lwd=2)
      mtext('Health index cut points',1,cex=1.5,line=2.5)
      for (j in seq_along(Nm)) text(x=D0[j],y=NROW(D)/4+NROW(D)/2,labels=Nm[[j]],
                                    srt=90,pos=2,offset=0.9,col=2)
      par(new = TRUE)
      rowplot(D[NROW(D):1,], pch=19, col=1:model$J)
      par(new = TRUE)
      rowplot(D[NROW(D):1,], pch=21, col=1)
      suppressWarnings(par(opar))
    }
    res <- list(q.all=D0, q=D)
  }
  if (plotf) invisible(res) else return(res)
}

#' Calcualte adjusted health levels.
#' @description
#' Calcualte adjusted health levels according to th Jurges' method.
#' @param model a fitted \code{gotm} model.
#' @param formula a formula containing the variables. It is by default set to threshold formula.
#' @param data data used to fit the model.
#' @param plotf logical indicating if to plot the results.
#' @param sep separator for levels names.
#' @param revf logical indicating if self-reported health classes are ordered in increasing order.
#' @param mar see \code{\link{par}}.
#' @param oma see \code{\link{par}}.
#' @export
gethealthlevels<-function(model, formula=model$thresh.formula,
                          data=environment(model$thresh.formula), revf = NULL,
                          plotf = TRUE, sep='_',mar=c(7,2,1.5,0.5),oma=c(0,3,0,0)){
  if (class(formula)=='formula') inte <- formula2classes(formula, data, sep=sep) else stop('Not implemented.')
  cpall<-basiccutpoints(model, plotf = FALSE, revf = revf)
  TAB1 <- round(table(original=model$y_i, adjusted=cpall$adjused.health.levels)*100/length(model$y_i),2)
  tmp <- untable(t(table(factor(data$r.health,levels=levels(cpall$adjused.health.levels)), inte)))
  tmp <-tmp/rowSums(tmp)
  tmp2 <- untable(t(table(cpall$adjused.health.levels, inte)))
  tmp2 <- tmp2/rowSums(tmp2)
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
  res <- list(original= tmp, adjusted= tmp2, tab= TAB1)
  class(res) <- 'healthlevels'
  if (plotf) invisible(res) else return(res)
}

#' @keywords internal
getsq <- function(x, xf=1) c(ceiling(x/ceiling(xf*x/sqrt(x))),ceiling(xf*x/sqrt(x)))

#' Plot a comparison of original and adjusted health levels
#' @param object an object generated by \code{\link{gethealthlevels}}
#' @param pch,xlab,ylab,mar,oma common graphical parameters.
#' @param ratio an aspect ratio for panels composition.
#' @export
comparehealthlevels<-function(object, pch=19,xlab='Original frequency [%]',ylab='Adjusted frequency [%]',
                            mar=c(2.5, 1, 1.5, 1), oma=c(4, 4, .1, .1),
                            ratio = 1){
  if (class(object) != 'healthlevels') stop('The object must be of class: "healthlevels".')
  sq <- getsq(NROW(object$tab), ratio)
  opar <- par(c('mar','oma'))
  par(mfrow=sq,oma=oma,mar=mar)
  for(k in seq_len(NROW(object$tab))){
    plot(object$original[,k]*100,object$adjusted[,k]*100,pch=pch,xlab='',ylab='',main=colnames(object$adjusted)[k],
         xaxs='r',yaxs='r')
    lines(-20:120,-20:120)
    posi <- approx(x=range(object$original[,k]),y=c(4,2),xout=object$original[,k],method="constant")$y
    text(object$original[,k]*100,object$adjusted[,k]*100,labels=names(object$adjusted[,k]),cex=0.6,pos=posi,offset=0.5)
    legend('topleft','underrated     ',bg=gray(0.8,alpha = 0.5),inset=0.05,xjust=0.5,box.col=NA,y.intersp=0.5)
    legend('bottomright','overrated     ',bg=gray(0.8,alpha = 0.5),inset=0.05,xjust=0.5,box.col=NA,y.intersp=0.5)
  }
  par(mfrow=c(1,1))
  par(oma=c(0,0,0,0),mar=mar)
  mtext(xlab,1,cex=1.5)
  mtext(ylab,2,cex=1.5, line=-0.5)
  suppressWarnings(par(opar))
  invisible(object)
}

#Here is the place for OAXACA-Blinder decomposition of self reported helth beween countries, or genders