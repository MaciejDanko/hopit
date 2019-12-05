#' Not \%in\% function
#'
#' @param x,y numeric vectors
#' @usage x \%notin\% y
#' @author Maciej J. Danko
#' @export
'%notin%' <-function(x, y) match(x, y, nomatch = 0L) == 0L


#' Check whether one set contains all elements of another set
#'
#' Check whether the y set contains all elements of set x
#' @param x,y numeric vectors
#' @usage x \%c\% y
#' @author Maciej J. Danko
#' @export
'%c%' <-function(x, y) all(match(x, y, nomatch = 0L))


#' Not \%c\% function
#'
#' Check whether the y set contains none of elements of the x set
#' @param x,y numeric vectors
#' @usage x \%notc\% y
#' @author Maciej J. Danko
#' @export
'%notc%' <- function(x, y) !all(match(x, y, nomatch = 0L))


#' @noRd
#' @keywords internal
rep_row <- function(mat, times) t(matrix(t(mat), NCOL(mat), NROW(mat) * times))


# INTERNAL: Calculate special matrices for gradient calculation
# @param model fitted model.
# @return updated model.
#' @keywords internal
# @author Maciej J. Danko
#' @noRd
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
  model$YYY3 <- dY
  model
}


# INTERNAL: numerical gradient
#' @keywords internal
# @param fn function.
# @param par parameters
# @param eps epsilon.
# @param ... other parameters passed to fn
# @author Maciej J. Danko
#' @noRd
my.grad <- function(fn, par, eps, ...){
  sapply(1L : length(par), function(k){
    epsi <- rep(0L, length(par))
    epsi[k] <- eps
    (fn(par + epsi, ...) - fn(par - epsi, ...))/2/eps
  })
}


#' @noRd
#' @keywords internal
hopit_c_link<-function(model){
  model$link <- tolower(model$link)[1]
  if (model$link %in% c('probit','logit')){
    if (model$link=='probit') link=0 else link=1
  } else stop(paste(hopit_msg(17),model$link),call. = NULL)
  link
}


# INTERNAL: Converts a vector of an categorical variable into a matrix with dummies in columns
#' @noRd
#' @keywords internal
Vector2DummyMat<-function(V) sapply(levels(as.factor(V)), function(k)
  as.factor(V) == k)*1L


#' @keywords internal
#' @noRd
decomposeinteractions<-function(allterms) lapply(allterms,
                                function(k) strsplit(k,split=':', fixed=TRUE)[[1]])

#' @keywords internal
#' @noRd
untable <- function(x) {
  names(attr(x, "dimnames")) <- c('','')
  as.matrix(x)
}


#' @keywords internal
#' @noRd
findintercept<-function(varnames) grepl('(Intercept)', varnames, fixed = TRUE)


#' @keywords internal
#' @noRd
formula2classes <- function(formula, data, sep='_', add.var.names = FALSE,
                            return.matrix = FALSE){
  tmp <- stats::model.frame(formula, data)
  mod.mat <- tmp
  lv <- lapply(seq_len(NCOL(tmp)),function (k) levels(as.factor(tmp[,k])))
  names(lv) <-colnames(tmp)
  tmp2 <- expand.grid(lv)
  if (add.var.names) tmp2 <- sapply(seq_len(NCOL(tmp2)), function (k)
    paste(colnames(tmp2)[k],'[',tmp2[,k],']',sep=''))
  nlv <- levels(interaction(as.data.frame(tmp2),sep=sep))
  if (add.var.names) tmp <- sapply(seq_len(NCOL(tmp)), function (k)
    paste(colnames(tmp)[k],'[',tmp[,k],']',sep=''))
  tmp <- interaction(as.data.frame(tmp),sep=sep)
  tmp <- factor(tmp, levels=nlv)
  if (return.matrix) list(x = tmp, mat = mod.mat, class.mat = tmp2) else tmp
}


#' @noRd
#' @keywords internal
cumsum_row<-function(mat) t(apply(as.matrix(mat), 1L, cumsum))


# INTERNAL: Clasify individuals according to the latent.params and calculated thresholds
#
#' @keywords internal
# @param model \code{hopit} object.
# @author Maciej J. Danko
#' @noRd
classify.ind<-function(model){
  p <- hopit_ExtractParameters(model, model$coef)
  a <- hopit_Threshold(thresh.lambda = p$thresh.lambda,
                       thresh.gamma = p$thresh.gamma,
                       model = model)
  b <- hopit_Latent(p$latent.params, model)
  a_0=a[,-1]
  a_J=a[,-ncol(a)]
  Ey_i <- sapply(1L : model$N, function(k)
    which((b[k]<a_0[k,]) & (b[k]>=a_J[k,])) )
  Ey_i <- factor(levels(model$y_i)[Ey_i],levels(model$y_i))
  Ey_i
}


#' @keywords internal
#' @noRd
sort.terms <- function(x) {
  tmp <- decomposeinteractions(x)
  tmp <-lapply(tmp, sort)
  sapply(tmp, paste, collapse=':')
}


#' @keywords internal
#' @noRd
# x- to be sorted, y-pattern
order.as.in<-function(x, y){
  tmp <- data.frame(y=y,id=seq_along(y))
  tmp <- tmp[order(y),]
  tmp2 <- data.frame(x=x,id=seq_along(x))
  tmp2 <- tmp2[order(x),]
  tmp <-cbind(tmp, tmp2)
  tmp[order(tmp[,2]),4]
}
#order.as.in(x=c('aa','cc','bb'),y=c('c','a','b'))


# INTERNAL: Use glm() to get starting parameters
#
# @param object \code{hopit} object.
# @param data data.frame with data used to fit the object.
# @return updated \code{hopit} object.
#' @keywords internal
# @author Maciej J. Danko
#' @noRd
start_glm<-function(object, data){
  g <- as.numeric(object$y_i)
  Y <- sapply(2:max(g), function(k) g<k)
  res <- sapply(seq_len(ncol(Y)),function(yi){
    zdat <- data
    zdat$yi <- Y[,yi]
    weights <- zdat$weights <- object$weights
    f1 <- object$latent.formula
    f1[[2]] <- as.name('yi')
    f1 <- stats::update(f1, paste('.~.+',
                                  paste(unclass(object$thresh.formula))[-1]))
    gl <- stats::glm(f1, data=zdat,
                     family=stats::binomial(link=object$link),
                     weights = weights)
    gl$coef
  })

  st.cn.l.mm <- sort.terms(colnames(object$latent.mm))
  st.cn.t.mm <- sort.terms(colnames(object$thresh.mm))
  st.cn.res  <- sort.terms(rownames(res))
  lind <- which(st.cn.res %in% st.cn.l.mm)
  glm.lambda <- res[which(grepl('Intercept',rownames(res))),]

  glm.latent <- res[lind, ]
  if (object$parcount[1]>1) glm.latent <- - rowMeans(glm.latent) else
    glm.latent <- - mean(glm.latent)
  glm.latent <- glm.latent[order.as.in(st.cn.res[lind], st.cn.l.mm)]
  if (!object$thresh.no.cov) {
    thr.ext.nam <-as.character(interaction(expand.grid(seq_len(object$J-1),
                               colnames(object$thresh.mm))[,2:1],sep=':'))
    indx <- which(st.cn.res %in% st.cn.t.mm)
    glm.gamma <- res[indx,]
    if (object$parcount[3]>object$parcount[2]) glm.gamma <-
      glm.gamma[order.as.in(st.cn.res[indx], st.cn.t.mm),]
    glm.gamma <- as.vector(t(glm.gamma))
    names(glm.gamma) <- thr.ext.nam
  } else {
    glm.gamma <- NULL
  }

  if (any(glm.lambda!=sort(glm.lambda))) stop(call.=NULL, hopit_msg(98))
  object$glm.start <- c(glm.latent,glm.lambda,glm.gamma)
  object$glm.start.ls <-list(latent.params = glm.latent,
                            thresh.lambda = glm.lambda,
                            thresh.gamma = glm.gamma)
  object
}


#' INTERNAL: Get the starting parameters
#'
#' @param object a \code{hopit} object.
#' @param data a data.frame with data used to fit the object.
#' @return an updated \code{hopit} object.
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
get.hopit.start<-function(object, data){
  logSigma <- 0
  object <- suppressWarnings(start_glm(object, data))
  par.ls <- object$glm.start.ls

  if (object$thresh.no.cov){
    z <- glm2hopit_nogamma(par.ls$latent.params,
                           par.ls$thresh.lambda,
                           thresh_1_exp = object$control$thresh.1.exp)
  } else {
    z <- glm2hopit(par.ls$latent.params,
                   par.ls$thresh.lambda,
                   par.ls$thresh.gamma,
                   thresh_1_exp = object$control$thresh.1.exp)
  }

  if (object$hasdisp) {
    object$start <- c(z$coef, logSigma)
  } else object$start <- z$coef

  object
}


#' @keywords internal
#' @noRd
drop.levels.response<-function(formula, data){
  response.name <- colnames(stats::model.frame(
    stats::update(formula,'.~1'), data))
  if (length(response.name) || nchar(response.name)) {
    tmp <- data[[response.name]]
    lv <- levels(tmp)
    tmp <- droplevels(tmp)
    if (length(lv)!=length(levels(tmp)))
      message(paste(hopit_msg(104),'"',response.name,'".',sep=''))
    data[[response.name]] <- tmp
  }
  data
}

#' @keywords internal
#' @noRd
drop.levels.data<-function(formula, data){
  varlist <- attr(stats::terms(formula),"term.labels")
  ind <- which(!grepl(':',varlist))
  varlist <- varlist[ind]
  ind <- which(sapply(varlist, function(k) is.factor(data[[k]])))
  varfactlist <- varlist[ind]
  for (k in seq_along(varfactlist)) {
    tmp <- data[[varfactlist[k]]]
    lv <- levels(tmp)
    tmp <- droplevels(tmp)
    if (length(lv)!=length(levels(tmp)))
      message(paste(hopit_msg(104),'"',varfactlist[k],'".',sep=''))
    data[[varfactlist[k]]] <- tmp
  }
  data
}


#' @keywords internal
#' @noRd
transform.data<-function(formula, data, type = 'min'){
  varlist <- attr(stats::terms(formula),"term.labels")
  ind <- which(!grepl(':',varlist))
  varlist <- varlist[ind]
  ind <- which(sapply(varlist, function(k) is.numeric(data[[k]])))
  varnumlist <- varlist[ind]
  type <- type [1]
  for (k in seq_along(varnumlist)) {
    tmp <- data[[varnumlist[k]]]
    if (tolower(type) == 'min') {
      tmp <- (tmp - min(tmp))
    } else if (tolower(type) == 'scale_01') {
      tmp <- (tmp - min(tmp))/(max(tmp)-min(tmp))
    } else if (tolower(type) == 'standardize' || tolower(type) == 'standardise') {
      tmp <- (tmp - mean(tmp))/stats::sd(tmp)
    } else if (tolower(type) == 'standardize_trunc' || tolower(type) == 'standardise_trunc') {
      tmp <- (tmp - min(tmp))/stats::sd(tmp)
    } else if (tolower(type) != 'none') stop(call.=NULL, hopit_msg(99))
    data[[varnumlist[k]]] <- tmp
  }
  data
}


#' @keywords internal
#' @noRd
analyse.formulas<-function(object, latent.formula, thresh.formula, data){

  thresh.terms <- attr(stats::terms(thresh.formula),'term.labels')
  latent.terms <- attr(stats::terms(latent.formula),'term.labels')
  if(any(thresh.terms %in% latent.terms)){
    tmp <- thresh.terms[thresh.terms %in% latent.terms]
    stop(paste(hopit_msg(91), ' ', paste(tmp,collapse=', '),'. ',
               hopit_msg(92), sep = ''), call. = NULL)
  }

  response.name <- colnames(stats::model.frame(stats::update(latent.formula,
                                                             '.~1'), data))
  if (length(response.name)>1) stop()
  thresh.list <- decomposeinteractions(thresh.terms)
  latent.list <- decomposeinteractions(latent.terms)

  if (response.name %in% unlist(latent.list) ||
      response.name %in% unlist(thresh.list)) stop(hopit_msg(102), call. = NULL)

  thresh.list.L <- sapply(thresh.list, length)
  latent.list.L <- sapply(latent.list, length)
  thresh.main <- unlist(thresh.list[which(thresh.list.L<2)])
  latent.main <- unlist(latent.list[which(latent.list.L<2)])

  if (length(thresh.list) && length(latent.list))
    sapply(thresh.list[which(thresh.list.L>1)], function(k)
      if (!all((k%in%thresh.main) | (k%in%latent.main)))
        stop(paste(hopit_msg(93), paste(k,collapse = ':'),
                   hopit_msg(94)), call. = NULL) else
                     paste(paste(k, collapse = ':'), ': OK', sep = ''))

  if (length(latent.list)) {
    if (length(thresh.list)) {
      sapply(latent.list[which(latent.list.L>1)], function(k)
        if (!all((k%in%latent.main) | (k%in%thresh.main)))
          stop(paste(hopit_msg(93), paste(k, collapse = ':'),
                     hopit_msg(94)), call. = NULL) else
                       paste(paste(k,collapse = ':'), ': OK', sep = ''))
    } else {
      sapply(latent.list[which(latent.list.L>1)], function(k)
        if (!all((k%in%latent.main)))
          stop(paste(hopit_msg(93), paste(k,collapse = ':'),
                     hopit_msg(94)), call. = NULL) else
                       paste(paste(k, collapse = ':'), ': OK', sep = ''))
    }

  cross.inter.latent.list <- latent.list[which(sapply(latent.list, function(k)
    !all(k%in%latent.main)))]

  } else cross.inter.latent.list <- list()

  if (length(cross.inter.latent.list)) {
    cross.inter.latent <- paste(sapply(cross.inter.latent.list, paste,
                                       collapse=':'),collapse=' + ')
    cross.main.latent.list <- unlist(cross.inter.latent.list)
    cross.main.latent.list <- cross.main.latent.list[
      which(cross.main.latent.list%notin%latent.main)]
    cross.main.latent <- paste(cross.main.latent.list, collapse=' + ')
    latent.formulaA <-stats::update(latent.formula,
                                    paste('.~. + ', cross.main.latent))
    m1<-stats::model.matrix(latent.formulaA, data)
    m2<-stats::model.matrix(latent.formula, data)
    cm1<-colnames(m1)
    rco<-cm1[cm1 %notin% colnames(m2)]
    latent.mm <- m1[,which(cm1 %notin% rco)]
  } else {
    latent.mm <- stats::model.matrix(latent.formula, data)
    cross.inter.latent <- NULL
  }

  if (length(thresh.list)) {
    cross.inter.thresh.list <- thresh.list[which(sapply(thresh.list,
                               function(k) !all(k%in%thresh.main)))]
  } else cross.inter.thresh.list <- NULL
  if (length(cross.inter.thresh.list)) {
    cross.inter.thresh <- paste(sapply(cross.inter.thresh.list,
                                       paste, collapse=':'),collapse=' + ')
    cross.main.thresh.list <- unlist(cross.inter.thresh.list)
    cross.main.thresh.list <- cross.main.thresh.list[
      which(cross.main.thresh.list%notin%thresh.main)]
    cross.main.thresh <- paste(cross.main.thresh.list, collapse=' + ')
    thresh.formulaA <-stats::update(thresh.formula,
                                    paste('~. + ',cross.main.thresh))
    m1<-stats::model.matrix(thresh.formulaA, data)
    m2<-stats::model.matrix(thresh.formula, data)
    cm1<-colnames(m1)
    rco<-cm1[cm1 %notin% colnames(m2)]
    thresh.mm <- m1[,which(cm1 %notin% rco)]
  } else {
    thresh.mm <- stats::model.matrix(thresh.formula, data)
    cross.inter.thresh <- NULL
  }

  latent.names <- colnames(latent.mm)
  grpi <- findintercept(latent.names)
  indi <- which(!grpi)
  latent.mm <- as.matrix(latent.mm[,indi])
  latent.names <- latent.names[indi]
  colnames(latent.mm) <- latent.names

  thresh.names <- colnames(thresh.mm)
  grpi <- findintercept(thresh.names)
  indi <- which(!grpi)
  thresh.mm <- as.matrix(thresh.mm[,indi])
  thresh.names <- thresh.names[indi]
  colnames(thresh.mm) <- thresh.names

  object$latent.names <- latent.names
  object$latent.formula <- latent.formula
  object$latent.terms <- latent.terms
  object$latent.mm <- latent.mm
  object$cross.inter.latent <- cross.inter.latent
  object$thresh.names <- thresh.names
  object$thresh.formula <- thresh.formula
  object$thresh.terms <- thresh.terms
  object$thresh.mm <- thresh.mm
  object$cross.inter.thresh <- cross.inter.thresh

  if (any(dim(object$thresh.mm) == 0L)) {
    object$thresh.no.cov <- TRUE
  } else {
    object$thresh.no.cov <- FALSE
  }

  object
}
