#' @export
healthindex<-function(model, subset=NULL, plotf = FALSE) {
  #0 is the worse possible health, 1 is the best possible health
  if (length(subset)==0) subset=seq_along(model$y_i)
  hi <- (1 - (model$y_latent_i/sum(coef(model)[1:model$parcount[1]])))[subset]
  if (plotf) plot(model$y_i[subset], hi,las=3, ylab='Health index')
  hi
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
untable<-function(x) {
  names(attr(x, "dimnames")) <- c('','')
  as.matrix(x)
}

#' @keywords internal
formula2classes<-function(formula, data, sep='_'){
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

#' @export
gethealthindexquantiles<-function(model, x=model$thresh.formula, data, plotf = TRUE, sep='_',
                                  mar=c(4,8,1.5,0.5),oma=c(0,0,0,0)){
  if (class(x)=='formula') tmp <- formula2classes(x, data, sep=sep) else stop('Not implemented.')
  D <- t(sapply(levels(tmp),function(k) quantile(healthindex(model, tmp==k))))
  D0 <- quantile(healthindex(model))
  if (plotf){
    par(mar=mar,oma=oma)
    rowplot(D[NROW(D):1,-c(1,model$J)], pch=c(21,19,21), col=1)
    abline(v=D0[-c(1,model$J)])
    mtext('Health index',1,cex=1.5,line = 2.5)
  }
  list(q.all=D0, q=D)
}

#' @export
cutpoints<-function(model, subset=NULL, plotf = TRUE, mar=c(4,4,1,1),oma=c(0,0,0,0), revf=NULL){
  
  if (length(subset)==0) subset=seq_along(model$y_i)
  Y=model$y_i[subset]
  if (length(revf)==0) { 
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
    par(mar=mar, oma=oma)
    z<-hist(h.index,200,xlab='',ylab='' ,
            main='', yaxs='i', col=grey(0.8, alpha = 0.5),border=grey(0.4, alpha = 0.5))
    for (j in seq_along(Nm)) text(x=R1[j],y=(1.1*max(z$counts))/2,labels=Nm[[j]],
                                  srt=90,pos=2,offset=0.67,col=2)
    box()
    abline(v=R1,lwd=2,col=2)
    mtext('Health index',1,cex=1.5,line = 2.5)
    mtext('Counts',2,cex=1.5, line=2.5)
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
  list(cutpoints=R1, adjused.health.levels=adjused.health.levels)
}

#' @export
getcutpoints<-function(model, x=model$thresh.formula, data, plotf = TRUE, sep='_',
                       revf=NULL, mar=c(4,8,1.5,0.5),oma=c(0,0,0,0)){
  Y=model$y_i
  if (length(revf)==0) { 
    message(paste('Are the levels',toString(levels(Y)),'sorted in incresing order?'))
    message('if not please set revf to TRUE.')
    stop('"revf" must be given',call.=FALSE)
  }
  if (revf) dorev <- rev else dorev <- identity
  if (class(x)=='formula') tmp <- formula2classes(x, data, sep=sep) else stop('Not implemented.')
  D <- t(sapply(levels(tmp),function(k) cutpoints(model=model, tmp==k, plotf = FALSE, revf = revf)$cutpoints))
  D0 <- cutpoints(model, plotf = FALSE, revf = revf)$cutpoints
  lv <- dorev(as.character(levels(model$y_i)))
  Nm <- paste(lv[-length(lv)],lv[-1],sep=' | ')
  colnames(D) <- Nm
  Nm <- sapply(Nm, function(k) bquote(bold(.(k))))
  
  if (plotf){
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
  }
  list(q.all=D0, q=D)
}

#' @export
gethealthlevels<-function(model, x=model$thresh.formula, 
                          data=environment(model$thresh.formula), 
                          plotf = TRUE, sep='_',mar=c(7,2,1.5,0.5),oma=c(0,3,0,0),
                          revf = NULL){
  if (class(x)=='formula') inte <- formula2classes(x, data, sep=sep) else stop('Not implemented.')
  cpall<-cutpoints(model, plotf = FALSE, revf = revf)
  TAB1 <- round(table(original=model$y_i, adjusted=cpall$adjused.health.levels)*100/length(model$y_i),2)
  tmp <- untable(t(table(factor(dta.ch$r.health,levels=levels(cpall$adjused.health.levels)), inte)))
  tmp <-tmp/rowSums(tmp)
  tmp2 <- untable(t(table(cpall$adjused.health.levels, inte)))
  tmp2 <- tmp2/rowSums(tmp2)
  if (plotf) {
    par(mfrow=c(1,2))
    par(mar=mar,oma=oma)
    barplot(t(tmp),las=3,main='Original')
    barplot(t(tmp2),las=3,main='Adjusted', legend.text=TRUE)
    par(mfrow=c(1,1))
    par(mar=mar,oma=rep(0,4))
    mtext('Fraction [%]',2,cex=1.5)
  }
  list(original= tmp, adjusted= tmp2, tab= TAB1)
}  

#' @keywords internal
getsq <- function(x,xf=1) c(ceiling(x/ceiling(xf*x/sqrt(x))),ceiling(xf*x/sqrt(x)))

#' @export
comparehealthlevels<-function(model, data, pch=19,xlab='Original',ylab='Adjusted',
                              revf=NULL,
                              mar=c(2.5, 1, 1.5, 1), oma=c(4, 4, .1, .1)){
  sq <- getsq(model$J)
  d <- gethealthlevels(model, data=data, plotf = FALSE, sep='\n', revf = revf)
  par(mfrow=sq,oma=oma,mar=mar)
  for(k in seq_len(model$J)){
    plot(d$original[,k]*100,d$adjusted[,k]*100,pch=pch,xlab='',ylab='',main=colnames(d$adjusted)[k],
         xaxs='r',yaxs='r')
    lines(-20:120,-20:120)
    posi <- approx(x=range(d$original[,k]),y=c(4,2),xout=d$original[,k],method="constant")$y
    text(d$original[,k]*100,d$adjusted[,k]*100,labels=names(d$adjusted[,k]),cex=0.6,pos=posi,offset=0.5)
    legend('topleft','underrated     ',bg=gray(0.8,alpha = 0.5),inset=0.05,xjust=0.5,box.col=NA,y.intersp=0.5)
    legend('bottomright','overrated     ',bg=gray(0.8,alpha = 0.5),inset=0.05,xjust=0.5,box.col=NA,y.intersp=0.5)
  }
  par(mfrow=c(1,1))
  par(oma=c(0,0,0,0),mar=mar)
  mtext(xlab,1,cex=1.5)
  mtext(ylab,2,cex=1.5)
  invisible(d)
}

