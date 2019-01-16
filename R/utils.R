#' Not \%in\% function
#'
#' @param x,y numeric vectors
#' @usage x \%notin\% y
#' @author Maciej J. Danko
#' @export
'%notin%' <-function(x, y) match(x, y, nomatch = 0L) == 0L


#' Check if one set is a subset of an another subset
#'
#' @param x,y numeric vectors
#' @usage x \%c\% y
#' @author Maciej J. Danko
#' @export
'%c%' <-function(x, y) all(match(x, y, nomatch = 0L))


#' Not \%c\% function
#'
#' @param x,y numeric vectors
#' @usage x \%notc\% y
#' @author Maciej J. Danko
#' @export
'%notc%' <- function(x, y) !all(match(x, y, nomatch = 0L))


#' @noRd
rep_row <- function(mat, times) t(matrix(t(mat), NCOL(mat), NROW(mat) * times))


#' INTERNAL: Calculate special matrices for gradient calaculation
#' @param model fitted model.
#' @return updated model
#' @keywords internal
#' @author Maciej J. Danko
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
#' @author Maciej J. Danko
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
#' @author Maciej J. Danko
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
#' @author Maciej J. Danko
cumsum_row<-function(mat) t(apply(as.matrix(mat), 1L, cumsum))


#' INTERNAL: Clasify individuals according to the latent.params and calculated thresholds
#'
#' @keywords internal
#' @param model \code{hopit} object.
#' @author Maciej J. Danko
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


#' INTERNAL: Use glm to get starting parameters
#' @param model \code{hopit} object.
#' @param data data.frame with data used to fit the model.
#' @return updated model
#' @keywords internal
#' @author Maciej J. Danko
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
    if (!gl$converged) warning(hopit_msg(19), call.=NA)
    gl$coef
    #check convergence
  })
  glm.lambda <-  res[which(grepl('Intercept',rownames(res))),]
  glm.latent <- res[rownames(res)%in%colnames(model$latent.mm),]
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
  model <- suppressWarnings(start.glm(model, data))
  par.ls <- model$glm.start.ls

  if (model$thresh.no.cov){
    z <- glm2hopit_nogamma(par.ls$latent.params, par.ls$thresh.lambda, thresh_1_exp = model$control$thresh.1.exp)
  } else {
    z <- glm2hopit(par.ls$latent.params, par.ls$thresh.lambda, par.ls$thresh.gamma, thresh_1_exp = model$control$thresh.1.exp)
  }

  if (model$hasdisp) {
    model$start <- c(z$coef, logTheta)
  } else model$start <- z$coef

  model
}
