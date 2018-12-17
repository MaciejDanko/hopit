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
