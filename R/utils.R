#' Not \%in\% function
#'
#' @param x, y numeric vectors
#' @usage x \%notin\% y
#' @export
'%notin%' <-function(x, y) match(x, y, nomatch = 0L) == 0L

#' Check if one set is a subset of an another subset
#'
#' @param x, y numeric vectors
#' @usage x \%c\% y
#' @export
'%c%' <-function(x, y) all(match(x, y, nomatch = 0L))

#' Not \%notc\% function
#'
#' @param x, y numeric vectors
#' @usage x \%notc\% y
#' @export
'%notc%' <- function(x, y) !all(match(x, y, nomatch = 0L))

#' @noRd
rep_row <- function(mat, times) t(matrix(t(mat), NCOL(mat), NROW(mat) * times))

#' INTERNAL: Calculate special matrices for gradient calaculation
#' @param model fitted model.
#' @return updated model
#' @keywords internal
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

#' INTERNAL: numerical gradient
#' @keywords internal
#' @param fn function
#' @param par parameters
#' @param eps epsilon
#' @param ... other parameters passed to fn
my.grad <- function(fn, par, eps, ...){
  sapply(1L : length(par), function(k){
    epsi <- rep(0L, length(par))
    epsi[k] <- eps
    (fn(par + epsi, ...) - fn(par - epsi, ...))/2/eps
  })
}

#' INTERNAL: Decode link parmeter
#' @param model fitted model.
#' @return updated model
#' @keywords internal
hopit_c_link<-function(model){
  model$link <- tolower(model$link)
  if (model$link %in% c('probit','logit')){
    if (model$link=='probit') link=0 else link=1
  } else stop(paste(hopit_msg(17),model$link),call. = NULL)
  link
}


#' INTERNAL: Converts a vector of an categorical variable into a matrix with dummies in columns
#'
#' @param V a vector of categories.
#' @author Maciej J. Danko
#' @keywords internal
Vector2DummyMat<-function(V) sapply(levels(as.factor(V)), function(k) as.factor(V) == k)*1L


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
formula2classes <- function(formula, data, sep='_', add.var.names = FALSE, return.matrix = FALSE){
  tmp <- model.frame(formula, data)
  mod.mat <- tmp
  lv <- lapply(seq_len(NCOL(tmp)),function (k) levels(as.factor(tmp[,k])))
  names(lv) <-colnames(tmp)
  tmp2 <- expand.grid(lv)
  if (add.var.names) tmp2 <- sapply(seq_len(NCOL(tmp2)), function (k) paste(colnames(tmp2)[k],'[',tmp2[,k],']',sep=''))
  nlv <- levels(interaction(as.data.frame(tmp2),sep=sep))
  if (add.var.names) tmp <- sapply(seq_len(NCOL(tmp)), function (k) paste(colnames(tmp)[k],'[',tmp[,k],']',sep=''))
  tmp <- interaction(as.data.frame(tmp),sep=sep)
  tmp <- factor(tmp, levels=nlv)
  if (return.matrix) list(x = tmp, mat = mod.mat, class.mat = tmp2) else tmp
}


#' INTERNAL: Do cumsum() in each row of a matrix
#'
#' @param mat a matrix.
#' @keywords internal
cumsum_row<-function(mat) t(apply(as.matrix(mat), 1L, cumsum))


#' INTERNAL: Column path in a matrix
#'
#' For each row of a matrix \code{mat} extracts a value corresponding to a column stored in a vector \code{y}.
#' @param mat a matrix
#' @param y a vector with column indices corresponding to each row in the matrix \code{mat}.
#' @author Maciej J. Danko
#' @keywords internal
col_path<-function(mat, y, offset = 0) colpath(mat, y, offset) # RcppEigen


#' Convert individual data to frequency table of unique combination of dependent and independent variables
#'
#' @param formula formula indicating, which variables will be used to construct new database.
#' @param data data.frame including all variables listed in formula.
#' @param FreqNam name of the column with frequencies.
#' @author Maciej J. Danko
#' @export
data2freq<-function(formula, data, FreqNam='Freq'){
  what <- c(deparse(formula[[2]]),attr(terms(formula),"term.labels"))
  tmp <- data[,which(names(data)%in%what)]
  V <- as.numeric(as.factor(apply(tmp,1,paste,collapse='',sep='')))
  tmp <- cbind(tmp,V)
  tmp <- tmp[order(V),]
  V <- V[order(V)]
  V2 <- sapply(unique(V),function(k) sum(V==k))
  newd <- tmp[match(unique(V), V),]
  newd[,NCOL(newd)] <- V2
  colnames(newd)[NCOL(newd)] <- FreqNam
  newd
}

# #' @keywords internal
# unravel <-function(mat, freq)  {
#  mat <- cbind(mat, freq)
#  FreqInd <- NCOL(mat)
#  ls <- apply(mat,1,function(k) ((rep_row(as.matrix(k[-FreqInd]),k[FreqInd]))))
#  do.call('rbind', ls)
#}

#' INTERNAL: Clasify individuals according to the latent.params and calculated thresholds
#' @keywords internal
#' @param model \code{hopit} object.
classify.ind<-function(model){
  p <- hopit_ExtractParameters(model, model$coef)
  a <- hopit_Threshold(thresh.lambda = p$thresh.lambda, thresh.gamma = p$thresh.gamma,
                      model = model)
  b <- hopit_Latent(p$latent.params, model)
  a_0=a[,-1]
  a_J=a[,-ncol(a)]
  Ey_i <- sapply(1L : model$N, function(k) which((b[k]<a_0[k,]) & (b[k]>=a_J[k,])) )
  Ey_i <- factor(levels(model$y_i)[Ey_i],levels(model$y_i))
  Ey_i
}

#' @keywords internal
order.as.in<-function(a,b){  #a = in
  if (!all(a %in% b) || !all(b %in% a)) stop()
  x <- data.frame(x=a, idx = seq_along(a), stringsAsFactors = FALSE)
  x <- x[order(x$x),]
  y <- data.frame(y=b, idy = seq_along(b), stringsAsFactors = FALSE)
  y <- y[order(y$y),]
  z <- data.frame(zx=x$idx, zy=y$idy)
  z <- z[order(z$zx),]
  z$zy
}

#' @keywords internal
greplin<-function(p, x) sapply(p, function(y) any(grepl(y, x, fixed = TRUE)))

#' @keywords internal
extractCoef <- function(COEF, model){
  Lind <- grepl('Intercept',names(COEF),fixed='TRUE')
  xglm.lambda <- sort(COEF[Lind])
  #Rind <- names(COEF) %in% colnames(model$latent.mm)
  Rind <- greplin(names(COEF), colnames(model$latent.mm))
  xglm.latent <-  COEF[Rind]
  oi <- order.as.in(a=colnames(model$latent.mm), b=names(xglm.latent))
  xglm.latent <- - xglm.latent[oi] #remove negative sign in latent
  xglm.gamma <- COEF[!Lind & !Rind]
  thr.ext.nam<-as.character(interaction(expand.grid(seq_len(model$J-1),colnames(model$thresh.mm))[,2:1],sep=':'))
  oi <- order.as.in(a=thr.ext.nam, b=names(xglm.gamma))
  xglm.gamma <- xglm.gamma[oi]
  list(start.ls = list(latent.params = xglm.latent,
                       thresh.lambda = xglm.lambda,
                       thresh.gamma = xglm.gamma),
       start <- c(xglm.latent, xglm.lambda, xglm.gamma))
}

#' INTERNAL: Use vglm to get starting parameters
#' @param model \code{hopit} object.
#' @param data data.frame with data used to fit the model.
#' @return updated model
#' importFrom VGAM vglm
#' @keywords internal
start.vglm <-function(model, data){
  logTheta <- 0
  latent.formula <- model$latent.formula
  thresh.formula <- model$thresh.formula
  if (length(thresh.formula)>2) thresh.formula[[2]] <- NULL
  thrf <- deparse(thresh.formula[[2]])
  Ynam <- deparse(latent.formula[[2]])
  tmr <- model$latent.terms
  tmt <- model$thresh.terms
  model$J <- length(levels(as.factor(data[,Ynam]))) #update J if not calculated
  #ltmt <- sapply(tmt,function(k) length(levels(as.factor(data[,k]))))-1
  ltmt <- ncol(model$thresh.mm)
  model$parcount <- c(length(tmr),(model$J-1),ltmt*(model$J-1)) #update model parcount if not calculated
  if (!length(model$weights)) model$weights <- rep(1,model$N)
  incc <- tmr %in% tmt
  if (any(incc)) {
    message(hopit_msg(18))
    ignored.var<- tmr[incc]
    latent.formula <- update(latent.formula, paste('~ . ',paste(' -', ignored.var,collapse='')))
  } else ignored.var <- NULL
  big.formula <- update(latent.formula, paste('~ ', thrf,' + . + 1'))
  Y <<- Vector2DummyMat(data[,paste(latent.formula[[2]])])
  w <- model$weights
  data$w <- w
  big.formula[[2]] <- as.name('Y')
  small.formula <- formula(paste('FALSE ~', thrf))
  mv2<-switch(model$link,
              probit = VGAM::vglm(big.formula, weights = w, data = data,
                                  family = VGAM::cumulative(parallel = small.formula, link = 'probit')), #direct substitution of link doesn't work
              logit = VGAM::vglm(big.formula, weights = w, data = data,
                                 family = VGAM::cumulative(parallel = small.formula, link = 'logit')))
  rm(Y, envir = .GlobalEnv)
  cmv2 <- coef(mv2)
  model$vglm <- mv2
  model$vglm.LL<-VGAM::logLik(mv2)
  # Lind <- grepl('Intercept',names(cmv2),fixed='TRUE')
  # vglm.lambdas <- sort(cmv2[Lind])
  # #Rind <- names(cmv2) %in% colnames(model$latent.mm)
  # Rind <- greplin(names(cmv2), colnames(model$latent.mm))
  # vglm.latent <-  cmv2[Rind]
  # oi <- order.as.in(a=colnames(model$latent.mm), b=names(vglm.latent))
  # vglm.latent <- - vglm.latent[oi] #remove negative sign in latent
  # vglm.gamma <- cmv2[!Lind & !Rind]
  # thr.ext.nam<-as.character(interaction(expand.grid(seq_len(model$J-1),colnames(model$thresh.mm))[,2:1],sep=':'))
  # oi <- order.as.in(a=thr.ext.nam, b=names(vglm.gamma))
  # vglm.gamma <- vglm.gamma[oi]
  # model$vglm.start.ls = list(latent.params = vglm.latent,
  #                            thresh.lambda = vglm.lambdas,
  #                            thresh.gamma = vglm.gamma)
  # model$vglm.start <- c(vglm.latent, vglm.lambdas, vglm.gamma)
  z <- extractCoef(cmv2, model)
  model$vglm.start.ls <- z$start.ls
  model$vglm.start <- z$start
  if (model$hasdisp) model$vglm.start <- c(model$vglm.start, logTheta)
  parcount <- model$parcount
  parcount[1] <- parcount[1] - length(ignored.var) #check!
  list(vglm.model=model,ignored.latent.var=ignored.var,new.latent.formula=latent.formula)
}

compareStr<-function(s1,s2) sapply(seq_along(s1), function(k) s1[k] == substr(s2[k],1,nchar(s1[k])))

#' INTERNAL: Use glm to get starting parameters
#' @param model \code{hopit} object.
#' @param data data.frame with data used to fit the model.
#' @return updated model
#' importFrom VGAM vglm
#' @keywords internal
start.glm<-function(model, data){
  g <- as.numeric(model$y_i)
  Y <- sapply(2:max(g), function(k) g<k)
  res <- sapply(seq_len(ncol(Y)),function(yi){
    zdat <- data
    zdat$yi <- Y[,yi]
    f1 <- model$latent.formula
    f1[[2]] <- as.name('yi')
    f1 <- update(f1, paste('.~.+',deparse(model$thresh.formula[[-1]])))
    gl <- glm(f1,data=zdat,family=binomial(link=model$link))
    if (!gl$converged) stop(hopit_msg(19), call.=NULL)
    gl$coef
    #check convergence
  })
  glm.lambda <-  res[which(grepl('Intercept',rownames(res))),]
  glm.latent <- res[rownames(res)%in%colnames(model$latent.mm),]
  #glm.latent <- res[unlist(sapply(attr(terms(model$latent.formula),'term.labels'),grep, x=rownames(res))),]
  glm.latent <- - rowMeans(glm.latent)
  if (!model$thresh.no.cov) {
    thr.ext.nam <-as.character(interaction(expand.grid(seq_len(model$J-1),colnames(model$thresh.mm))[,2:1],sep=':'))
    glm.gamma <- as.vector(t(res[rownames(res)%in%colnames(model$thresh.mm),]))
    names(glm.gamma) <- thr.ext.nam
  } else {
    glm.gamma <- NULL
  }

  model$glm.start <- c(glm.latent,glm.lambda,glm.gamma)
  model$glm.start.ls <-list(latent.params = glm.latent,
                            thresh.lambda = glm.lambda,
                            thresh.gamma = glm.gamma)
  model
}


#' INTERNAL: Translate vglm to hopit start parameters
#'
#' @param model \code{hopit} object.
#' @param data data.frame with data used to fit the model.
#' @return updated model
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
get.hopit.start<-function(model, data){
  logTheta <- 0
  #if ((!model$method) && (model$start.method=='glm')) {
    m <- suppressWarnings(start.glm(model, data))
    model <- m
    par.ls <- model$glm.start.ls
  # } else if ((model$method) || (model$start.method=='vglm')) {
  #   m <- suppressWarnings(start.vglm(model, data))
  #   model <- m$vglm.model
  #   par.ls <- model$vglm.start.ls
  # } else stop(hopit_msg(20),call.=NULL)

  #if ((!model$method) || (model$start.method=='glm')) {
    if (model$thresh.no.cov){
      z <- glm2hopit_nogamma(par.ls$latent.params, par.ls$thresh.lambda, thresh_1_exp = model$control$thresh.1.exp)
    } else {
      z <- glm2hopit(par.ls$latent.params, par.ls$thresh.lambda, par.ls$thresh.gamma, thresh_1_exp = model$control$thresh.1.exp)
    }
    if (length(m$ignored.latent.var)){
      ini.mis <- mean(z$latent_params)
      if (any(class(data[,m$ignored.latent.var])!='factor')) stop(hopit_msg(21),call. = NULL)
      npar <- sum(sapply(m$ignored.latent.var, function(k) length(levels(data[,k]))-1))
      model$latent.formula <- update(m$new.latent.formula,paste('~ . + ',paste(m$ignored.latent.var,collapse='',sep=' + ')))
      z$latent_params <- c(z$latent_params, rep(ini.mis, npar))
      z$coef <- c(z$latent_params,z$thresh_lambda,z$thresh_gamma)
    }

    model$start.ls$latent.params <- z$latent_params
    model$start.ls$lambda <- z$thresh_lambda
    model$start.ls$gamma.start <-z$thresh_gamma

    if (model$hasdisp) {
      model$start.ls$logTheta <- logTheta
      model$start <- c(z$coef, logTheta)
    } else model$start <- z$coef

  # } else {
  #   model$start.ls <- model$vglm.start.ls
  #   model$start <- model$vglm.start
  #
  # }
  model
}
