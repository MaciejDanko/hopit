#' INTERNAL: Update model according to new parameters
#'
#' @param model a fitted \code{gotm} model.
#' @param newregcoef new coeficients for health variable deffined in model$reg.formula
#' @param data data used to fit the model
#' @param crude logical indicating if to re-calculate hessian, variance-covariance matrix, and estfun.
#' @keywords internal
update.latent <-function(model, newregcoef, data, update.hessian=FALSE){
  coefnames <- names(model$coef)
  thresh.names <- colnames(model$thresh.mm)
  model$coef[seq_len(model$parcount[1])]=newregcoef
  class(model) <- 'gotm'
  names(model$coef) <- coefnames
  colnames(model$thresh.mm) <- thresh.names
  p <- gotm_ExtractParameters(model)
  model$maxlatentrange <- sort(gotm_latentrange(model=model, data=data))
  model$maxobservedlatentrange <-  sort(range(gotm_Latent(p$reg.params,model)))
  model$y_latent_i <- gotm_Latent(p$reg.params, model)
  model$Ey_i <- factor(colSums(sapply(1L : model$N, function(k) model$alpha[k,]<model$y_latent_i[k])),levels=1L:model$J)
  levels(model$Ey_i) <- levels(model$y_i)

  if (update.hessian) {
    my.grad <- function(fn, par, eps, ...){
      sapply(1L : length(par), function(k){
        epsi <- rep(0L, length(par))
        epsi[k] <- eps
        (fn(par + epsi, ...) - fn(par - epsi, ...))/2/eps
      })
    }
    hes <- my.grad(fn = gotm_derivLL, par = model$coef, model=model, eps = model$control$grad.eps, collapse = TRUE, negative=FALSE)
    model$hessian <- hes
    model$vcov <- try(solve(-hes), silent = T)
    if (class(model$vcov) == 'try-error')
      warning(call. = FALSE, 'Model is probably unidentifiable, $vcov (variance-covariance matrix) cannot be computed.')
    model$estfun <- gotm_derivLL(model$coef, model, collapse = FALSE)
  }
  model
}


#' Calculate health index
#' @description
#' Calcualte halth index from the latent variable. It takes values from 0 to 1, where
#' zero is prescribed to the worse possible health (maximum posible value fo the latent variable;
#' all conditions/diseases with negative effects on health are present) and 1 is prescribed
#' to the best possible health (calculated analogically).
#' @param model a fitted \code{gotm} model.
#' @param subset an optional vector specifying a subset of observations.
#' @param plotf logical indicating if to plot summary figure.
#' @param crude logical indicating if to calculate crude health measure based directly on self reported health levels.
#' @param healthlevelsorder order of self-reported healh levels. Possible values are \code{'increasing'} and \code{'decreasing'}
#' @param method the method of calcualtion of of latent range .........
#' @export
healthindex <- function(model, subset=NULL, plotf = FALSE, crude = FALSE, #scaled = TRUE,
                        healthlevelsorder = 'decreasing', method = c('observed','theoretical')) {
  #0 is the worse possible health, 1 is the best possible health
  method = method[1]
  p <- gotm_ExtractParameters(model)
  if (!length(model$maxlatentrange)) model$maxlatentrange <- sort(gotm_latentrange(model=model, data=data))
  if (!length(model$maxobservedlatentrange)) model$maxobservedlatentrange <- sort(range(gotm_Latent(p$reg.params,model)))
  if (length(subset) == 0) subset=seq_along(model$y_i)
  if (crude) {
    hi <- as.numeric(unclass(model$y_i))
    hi <- (hi - min(hi))/(diff(range(hi)))
    if (healthlevelsorder == 'decreasing') hi <- 1-hi else if (healthlevelsorder != 'increasing') stop('Unknown value for healthlevelsorder.')
    hi <- hi[subset]
    if (plotf) plot(model$y_i[subset], hi,las=3, ylab='Health index')
  } else {
    if (method=='theoretical') r <- model$maxlatentrange else if (method=='observed') r <- model$maxobservedlatentrange
    hi <- (1 - ((model$y_latent_i - r[1]) / diff(r)))[subset]
    if (plotf) plot(model$y_i[subset], hi,las=3, ylab='Health index')
  }
  if (plotf) invisible(hi) else return(hi)
}


contingencytables <- function(model, formula, data, names.reg=identity){
  if (class(names.reg)=='function') NN <- names.reg(colnames(model$reg.mm)) else NN <- names.reg
  if (class(formula)=='formula')
    tmp <- formula2classes(formula, data, sep=' ', return.matrix = TRUE) else stop('The formula parameter must be of class "formula".')
  colnames(model$reg.mm) <- NN
  M <- tmp$x
  cTAB <- sapply(levels(M), function(k) colSums(model$reg.mm[k==M,]))
  cTAB <- cbind(cTAB, ALL=rowSums(cTAB))
  fTAB.1 <- 100*cTAB / length(M)
  cTAB <- rbind(cTAB, 'Number of subjects'=c(table(tmp$x),length(M)))
  #fTAB.2 <- (100 * cTAB / t(matrix(cTAB[dim(cTAB)[1],],dim(cTAB)[2],dim(cTAB)[1])))[-dim(cTAB)[1],]
  fTAB.2 <- (100 * cTAB / matrix(cTAB[,dim(cTAB)[2]],dim(cTAB)[1],dim(cTAB)[2]))[,-dim(cTAB)[2]]
  list(countsTAB=cTAB, freqTAB1=fTAB.1, freqTAB2=fTAB.2)
}

#' Calculate disability weights
#' @description
#' Calculate disability weights
#' computed as the regression parameters from the generalised ordered probit model divided by the maximum possible range
#' of its linear prediction. The range is calcualted as difference between maximum and minimum possible value of the latent variable
#' given estimated parameters.
#' @seealso \code{\link{healthindex}}
#' @param model a fitted \code{gotm} model.
#' @param plotf logical indicating if to plot results.
#' @param mar see \code{\link{par}}.
#' @param oma see \code{\link{par}}.
#' @export
disabilityweights <- function (model, latent.method = c('observed','theoretical'),
                               plotpval = FALSE,
                               plotf = TRUE, mar = c(15, 4, 1, 1), oma = c(0, 0, 0, 0),
                               namesf = identity) {
  latent.method <- latent.method[1]
  #if (model$thresh.method != 'vglm') k <- 1 else k <- -1
  k=1
  p <- gotm_ExtractParameters(model)
  if (!length(model$maxlatentrange)) model$maxlatentrange <- sort(k*gotm_latentrange(model=model, data=data))
  if (!length(model$maxobservedlatentrange)) model$maxobservedlatentrange <- sort(k*range(gotm_Latent(p$reg.params,model)))
  z <- (k * model$coef)[seq_len(model$parcount[1])]

  if (class(namesf)=='function') names(z) <- namesf(names(z)) else names(z) <- namesf
  oz <- order(z, decreasing = TRUE)
  cfm <- z[oz]
  if (latent.method=='theoretical') r <- model$maxlatentrange else if (latent.method=='observed') r <- model$maxobservedlatentrange
  #if (!force.positive) {
    #res <- as.matrix((cfm - r[1])/diff(r))
    res <- as.matrix((cfm)/diff(r))
  #} else #if (force.positive) {
  #  res <- as.matrix((cfm - max(r[1],0))/min(diff(r),max(r)))
  #} else stop('Unknown method.', call. = FALSE)

  if (plotf) {

    opar <- par(c("mar", "oma"))
    par(mar = mar, oma = oma)
    rr <- barplot(t(res), las = 3)
    if(plotpval) {
      y <- summary(model)
      pval <- format(round(y$table$`Pr(>|z|)`[seq_len(model$parcount[1])],4),digits=4,scientific=FALSE)[oz]
      yr <- res/2
      ind <- yr < max(res)*0.1
      yr[ind] <- (res+max(res)*0.1)[ind]
      text(rr,yr,paste('P =',pval),srt=90,col=c('white','black')[1+ind])
    }
    mtext("Disability weight", 2, cex = 1.5, line = 2.5)
    suppressWarnings(par(opar))
  }
  if (plotf)
    invisible(res)
  else return(res)
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
formula2classes <- function(formula, data, sep='_', add.var.names = FALSE, return.matrix = FALSE){
  tmp <- model.frame(formula, data)
  mod.mat <- tmp
  lv <- lapply(seq_len(NCOL(tmp)),function (k) levels(tmp[,k]))
  names(lv) <-colnames(tmp)
  tmp2 <- expand.grid(lv)
  if (add.var.names) tmp2 <- sapply(seq_len(NCOL(tmp2)), function (k) paste(colnames(tmp2)[k],'[',tmp2[,k],']',sep=''))
  nlv <- levels(interaction(as.data.frame(tmp2),sep=sep))
  if (add.var.names) tmp <- sapply(seq_len(NCOL(tmp)), function (k) paste(colnames(tmp)[k],'[',tmp[,k],']',sep=''))
  tmp <- interaction(as.data.frame(tmp),sep=sep)
  tmp <- factor(tmp, levels=nlv)
  if (return.matrix) list(x = tmp, mat = mod.mat, class.mat = tmp2) else tmp
}

#' Get health index quantiles with respect to specified vaiables
#' @description
#' Get health index quantiles with respect to specified vaiables.
#' @param model a fitted \code{gotm} model.
#' @param formula a formula containing the variables. It is by default set to threshold formula.
#' @param data used to fit the model.
#' @param plotf logical indicating if to plot the results.
#' @param healthlevelsorder order of self-reported healh levels. Possible values are \code{'increasing'} and \code{'decreasing'}
#' @param sep separator for levls names.
#' @param mar see \code{\link{par}}.
#' @param oma see \code{\link{par}}.
#' @export
gethealthindexquantiles<-function(model, formula=model$thresh.formula, data=environment(model$thresh.formula),
                                  plotf = TRUE, sep='\n',sort.flag=FALSE,
                                  mar=c(4,8,1.5,0.5),oma=c(0,0,0,0), healthlevelsorder = 'decreasing'){
  if (class(formula)=='formula') tmp <- formula2classes(formula, data, sep=sep) else stop('Not implemented.')
  D <- t(sapply(levels(tmp),function(k) quantile(healthindex(model, tmp==k))))
  Jh <- floor(model$J/2)
  if (healthlevelsorder != 'decreasing') Jh <- model$J - Jh

  M.crude <- t(sapply(levels(tmp),function(k) sum(model$weights[tmp==k]*healthindex(model, tmp==k, crude = TRUE,
                                                             healthlevelsorder = healthlevelsorder))/sum(model$weights[tmp==k])))
  M.crude2 <- t(sapply(levels(tmp),function(k) {
    Hi <- as.numeric(unclass(model$y_i[tmp==k]))
    Mi <- Hi <= Jh
    W=model$weights[tmp==k]
    sum(W*Mi)/sum(model$weights[tmp==k])}
    ))

  M <- t(sapply(levels(tmp),
                function(k) sum(model$weights[tmp==k]*
                                  healthindex(model, tmp==k, crude = FALSE,
                                              healthlevelsorder = healthlevelsorder))/sum(model$weights[tmp==k])))
  if (sort.flag) {
    oD <- order(D[,3], decreasing = TRUE)
    D <- D[oD, ]
    M.crude <- M.crude[oD]
    M.crude2 <- M.crude2[oD]
    M <- M[oD]
  }

  D0 <- quantile(healthindex(model))
  IQR <- D[,4] - D[,2]
  if (plotf){
    opar <- par(c('mar','oma'))
    par(mar=mar,oma=oma)
    rowplot(D[NROW(D):1,-c(1,model$J)], pch=c(NA,19,NA), col=1)
    for (y in seq_len(NROW(D))) lines(c(D[y,2],D[y,4]),NROW(D)-c(y,y)+1,lwd=2)
    lines(M,rev(seq_along(M)),type='p',pch=4,col='red3',cex=1.5)
    #lines(M.crude2,rev(seq_along(M.crude2)),type='p',pch=15,col='blue3',cex=1.5)
    abline(v=D0[-c(1,model$J)])
    #text(x=1,y=rev(seq_along(M)),format(M.crude2,digits=3),col='blue3')
    mtext('Health index',1,cex=1.5,line = 2.5)
    suppressWarnings(par(opar))
  }
  res <- list(q.all=D0, q=D, IQR=IQR, mean = M, mean.crude = M.crude, frac.crude = M.crude2)
  if (plotf) invisible(res) else return(res)
}

#' Calcualte threshold cut-points using Jurges' method
#' @description
#' Calcualte threshold cut-points using Jurges' method.
#' @param model a fitted \code{gotm} model.
#' @param subset an optional vector specifying a subset of observations.
#' @param plotf logical indicating if to plot the results.
#' @param revf logical indicating if self-reported health classes are ordered in increasing order.
#' @param healthlevelsorder order of self-reported healh levels. Possible values are \code{'increasing'} and \code{'decreasing'}
#' @param mar see \code{\link{par}}.
#' @param oma see \code{\link{par}}.
#' @keywords internal
basiccutpoints <- function(model, subset=NULL, plotf = TRUE, mar=c(4,4,1,1),oma=c(0,0,0,0), revf=NULL, simple.plot=TRUE){

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
    if (!simple.plot) {
    for (j in seq_along(Nm)) text(x=R1[j],y=(1.1*max(z$counts))/2,labels=Nm[[j]],
                                  srt=90,pos=2,offset=0.67,col=2)
    } else {
      R11=-diff(c(0,R1,1))/2+c(R1,1)+strheight('S',units='figure')/2
      for (j in seq_along(lv)) text(x=R11[j],y=(3*1.1*max(z$counts))/4,labels=lv[j],
                                    srt=90,pos=3,offset=0.67,col=2)

    }
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
    adjused.health.levels<- cut(h.index, CIN,labels= dorev(levels(Y)))
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
                       plotf = TRUE, sep='_',sort.flag=FALSE,
                       revf=NULL, mar, oma, simple.names = TRUE){
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

    if (sort.flag) {
      oD <- order(D[,2], decreasing = FALSE)
      D <- D[oD, ]
    }

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
                          sort.flag=FALSE,
                          plotf = TRUE, sep='_',mar=c(7,2,1.5,0.5),oma=c(0,3,0,0)){
  if (class(formula)=='formula') inte_ <- formula2classes(formula, data, sep=sep, return.matrix = TRUE) else stop('Not implemented.')
  inte <- inte_$x
  namind <- inte_$class.mat
  nam <- levels(inte)
  cpall<-basiccutpoints(model, plotf = FALSE, revf = revf)
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

helatlevel.boot.plot<-function(object, model, namefunc = identity,
                               col=c('red4','blue4'),pch=c(1,0),
                               symlim = FALSE, pos.method = 2,
                               legtxt=c('Women','Men'),out1=NULL,out2=NULL,
                               myposi1 = NULL,
                               myposi2 = NULL) {
  par(mfrow=c(1,2),oma=c(3.5,2,0,0),mar=c(1,2.2,1,0.4))
  incl<-function(x,f) if (f) return (x) else return(NULL)
  plot(NA,NA,xlim=range(out1,object$org.lo[,1],object$org.hi[,1],incl(c(object$adj.lo[,1],object$adj.hi[,1]),symlim)),
       ylim=range(object$adj.lo[,1],object$adj.hi[,1],incl(c(out1,object$org.lo[,1],object$org.hi[,1]),symlim)),
       xlab='',ylab='',
       main=paste(levels(model$y_i)[model$J-c(0:1)],collapse='&'),xaxs='r',yaxs='r')
  lines(c(-100,100),c(-100,100),col='gray')
  lines(object$org.est[,1],object$adj.est[,1],type='p',pch=rep(pch, dim(object$org.est)[1]/2), col=rep(col, dim(object$org.est)[1]/2))
  for (j in seq_len(dim(object$org.est)[1])) lines(c(object$org.est[j,1],object$org.est[j,1]),c(object$adj.lo[j,1],object$adj.hi[j,1]),lwd=2, col=col[1+(j+1)%%2])

  if (pos.method == 1) {
    posi <- approx(x=range(object$org.est[,1]),y=c(4,2),xout=object$org.est[,1],method="linear")$y
    posi <- round(posi / 2) * 2
  } else if (pos.method == 3) {
    posi <- rep(0,length(object$org.est[,1]))
    posi[order(object$org.est[,1])] <- rep(c(2,4),length(object$org.est[,1])/2)
  } else if (pos.method == 2) {
    posi <- rep(0,length(object$org.est[,1]))
    posi[order(object$org.est[,1]*object$adj.est[,1])] <- rep(c(2,4),length(object$org.est[,1])/2)
  } else if (pos.method == 4) {
    if (!length(myposi1)) stop ('myposi should be given.')
    posi <-  myposi1
  }
  pos1 = posi
  text(object$org.est[,1],object$adj.est[,1],labels=namefunc(gsub('_',' ',names(object$adj.est[,1]))),cex=0.65,pos=posi,offset=0.3)
  if (length(legtxt)) legend('topleft',legtxt, pch=pch, col=col, bty='n')

  plot(NA,NA,xlim=range(out2,object$org.lo[,2],object$org.hi[,2],incl(c(object$adj.lo[,2],object$adj.hi[,2]),symlim)),
       ylim=range(object$adj.lo[,2],object$adj.hi[,2],incl(c(out2,object$org.lo[,2],object$org.hi[,2]),symlim)),xlab='',ylab='',
       main=paste(levels(model$y_i)[1:2],collapse='&'),xaxs='r',yaxs='r')
  lines(c(-100,100),c(-100,100),col='gray')
  lines(object$org.est[,2],object$adj.est[,2],type='p',pch=rep(pch, dim(object$org.est)[1]/2), col=rep(col, dim(object$org.est)[1]/2))
  for (j in seq_len(dim(object$org.est)[1])) lines(c(object$org.est[j,2],object$org.est[j,2]),c(object$adj.lo[j,2],object$adj.hi[j,2]),lwd=2, col=col[1+(j+1)%%2])

  if (pos.method == 1) {
    posi <- approx(x=range(object$org.est[,2]),y=c(4,2),xout=object$org.est[,2],method="linear")$y
    posi <- round(posi / 2) * 2
  } else if (pos.method == 3) {
    posi <- rep(0,length(object$org.est[,2]))
    posi[order(object$org.est[,2])] <- rep(c(2,4),length(object$org.est[,2])/2)
  } else if (pos.method == 2) {
    posi <- rep(0,length(object$org.est[,2]))
    posi[order(object$org.est[,2]*object$adj.est[,2])] <- rep(c(2,4),length(object$org.est[,2])/2)
  } else if (pos.method == 4) {
    if (!length(myposi2)) stop ('myposi should be given.')
    posi  <- myposi2
  }
  pos2 = posi

  text(object$org.est[,2],object$adj.est[,2],labels=namefunc(gsub('_',' ',names(object$adj.est[,2]))),cex=0.65,pos=posi,offset=0.3)
  par(mfrow=c(1,1),oma=c(0,0,0,0))
  mtext('Adjusted SRH [%]',2,line=1,cex=1.5)
  mtext('Original SRH [%]',1,line=-0.25,cex=1.5)
  list('p1'=pos1, 'p2'=pos2)
}

#' Bootstraping health levels
#'
#' @export
#' @author Maciej J. Danko
#' @importFrom MASS mvrnorm
gethealthlevels_boot<-function(model, formula=model$thresh.formula,
                          data=environment(model$thresh.formula), revf = NULL,
                          nboot=1000, alpha=0.05, plotF= FALSE, namefunc = identity,
                          col=c('red4','blue4'),pch=c(1,0)){

  N=seq_len(model$parcount[1])
  bootsample=MASS::mvrnorm(nboot,mu=model$coef[N],Sigma=summary(model)$vcov[N,N])
  GG0 <- gethealthlevels(model=model,formula=formula,data=data,revf=revf,plotf=FALSE)
  boots=sapply(seq_len(nboot), function(k) {
    bmodel <- update.latent(model,bootsample[k,N],data=data)
    GG <- gethealthlevels(model=bmodel,formula=formula,data=data,revf=revf,plotf=FALSE)
    tmp <- cbind(L=GG$adjusted[,1]+GG$adjusted[,2],H=GG$adjusted[,model$J-1]+GG$adjusted[,model$J])
    as.vector(tmp)
  })
  adj.est=100*cbind(L=GG0$adjusted[,1]+GG0$adjusted[,2],H=GG0$adjusted[,model$J-1]+GG0$adjusted[,model$J])
  org.est=100*cbind(L=GG0$original[,1]+GG0$original[,2],H=GG0$original[,model$J-1]+GG0$original[,model$J])
  orgN=cbind(L=GG0$N.original[,1]+GG0$N.original[,2],H=GG0$N.original[,model$J-1]+GG0$N.original[,model$J])
  adj.mean=100*matrix(apply(boots,1,mean),dim(org.est)[1],dim(org.est)[2]); dimnames(adj.mean)=dimnames(org.est)
  adj.lo=100*matrix(apply(boots,1, function(k) quantile(k,c(alpha/2))),dim(org.est)[1],dim(org.est)[2]); dimnames(adj.lo)=dimnames(org.est)
  adj.hi=100*matrix(apply(boots,1, function(k) quantile(k,1-c(alpha/2))),dim(org.est)[1],dim(org.est)[2]); dimnames(adj.hi)=dimnames(org.est)
  org.tmp=sapply(seq_len(length(org.est)), function(k) {
    tmp=replicate(nboot,mean(rbinom(n=as.vector(orgN)[k],size=1, prob=as.vector(org.est/100)[k])))
    quantile(tmp,c(alpha/2,1-c(alpha/2)))
    })
  org.lo=matrix(org.tmp[1,]*100,dim(org.est)[1],dim(org.est)[2]); dimnames(org.lo)=dimnames(org.est)
  org.hi=matrix(org.tmp[2,]*100,dim(org.est)[1],dim(org.est)[2]); dimnames(org.hi)=dimnames(org.est)
  res=list(org.est=org.est,org.lo=org.lo,org.hi=org.hi,adj.est=adj.est,adj.lo=adj.lo,adj.hi=adj.hi)

  if (plotF) {
    helatlevel.boot.plot(res, model, namefunc, col, pch)
  }
  return(res)
}


#' @keywords internal
getsq <- function(x, xf=1) c(ceiling(x/ceiling(xf*x/sqrt(x))),ceiling(xf*x/sqrt(x)))

#' Plot a comparison of original and adjusted health levels
#' @param object an object generated by \code{\link{gethealthlevels}}
#' @param simplified logical
#' @param pch,xlab,ylab,mar,oma common graphical parameters.
#' @param ratio an aspect ratio for panels composition.
#' @export
comparehealthlevels<-function(object,
                              simplified = TRUE,
                              pch=19,
                              xlab='Original frequency [%]',ylab='Adjusted frequency [%]',
                              mar=c(2.5, 1, 1.5, 1), oma=c(4, 4, .1, .1),
                              ratio = 1){
  if (class(object) != 'healthlevels') stop('The object must be of class: "healthlevels".')
  if (simplified) {
    if (NCOL(object$adjusted) %% 2) {
      sq=c(1,3)
      z <- floor(NCOL(object$adjusted)/2)
      f1 <- seq_len(z)
      f2 <- z + 1
      f3 <- sort(NCOL(object$adjusted) - f1) + 1
      h <- colnames(object$adjusted)
      h <- c(paste(h[f1],collapse='&'),h[f2],paste(h[f3],collapse='&'))
      object$adjusted <- cbind(rowSums(object$adjusted[,f1]),
                               (object$adjusted[,f2]),
                               rowSums(object$adjusted[,f3]))
      object$original <- cbind(rowSums(object$original[,f1]),
                               (object$original[,f2]),
                               rowSums(object$original[,f3]))
      colnames(object$original) <- colnames(object$adjusted) <- h
    } else {
      sq=c(1,2)
      z <- round(NCOL(object$adjusted)/2)
      f1 <- seq_len(z)
      f2 <- sort(NCOL(object$adjusted) - f1) + 1
      h <- colnames(object$adjusted)
      h <- c(paste(h[f1],collapse='&'),paste(h[f2],collapse='&'))
      object$adjusted <- cbind(rowSums(object$adjusted[,f1]),
                               rowSums(object$adjusted[,f2]))
      object$original <- cbind(rowSums(object$original[,f1]),
                               rowSums(object$original[,f2]))
      colnames(object$original) <- colnames(object$adjusted) <- h
    }
  } else sq <- getsq(NROW(object$tab), ratio)
  opar <- par(c('mar','oma'))
  par(mfrow=sq,oma=oma,mar=mar)
  U_ <- 'underrated     '
  O_ <- 'overrated     '
  z <- floor(NCOL(object$adjusted)/2)
  for(k in seq_len(NCOL(object$adjusted))){
    rr <- range(c(object$original[,k]*100,object$adjusted[,k]*100))
    plot(object$original[,k]*100,object$adjusted[,k]*100, pch=pch,
         xlab='', ylab='',main=colnames(object$adjusted)[k],
         xaxs='r', yaxs='r', xlim=rr, ylim=rr)
    lines(-20:120,-20:120)
    posi <- approx(x=range(object$original[,k]),y=c(4,2),xout=object$original[,k],method="linear")$y
    posi <- round(posi / 2) * 2
    text(object$original[,k]*100,object$adjusted[,k]*100,labels=names(object$adjusted[,k]),cex=0.6,pos=posi,offset=0.5)
    if (k <= z) {
      U <- O_
      O <- U_
      legend('topleft',U,bg=gray(0.8,alpha = 0.5),inset=0.05,xjust=0.5,box.col=NA,y.intersp=0.5)
      legend('bottomright',O,bg=gray(0.8,alpha = 0.5),inset=0.05,xjust=0.5,box.col=NA,y.intersp=0.5)
    } else if (k >= NCOL(object$adjusted)-z +1){
      U <- U_
      O <- O_
      legend('topleft',U,bg=gray(0.8,alpha = 0.5),inset=0.05,xjust=0.5,box.col=NA,y.intersp=0.5)
      legend('bottomright',O,bg=gray(0.8,alpha = 0.5),inset=0.05,xjust=0.5,box.col=NA,y.intersp=0.5)
    }

  }
  par(mfrow=c(1,1))
  par(oma=c(0,0,0,0),mar=mar)
  mtext(xlab,1,cex=1.5)
  mtext(ylab,2,cex=1.5, line=-0.5)
  suppressWarnings(par(opar))
  invisible(object)
}

