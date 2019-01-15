#' INTERNAL : updating model uses new set of latent variable coeficients
#'
#' @keywords internal
#' @param model a fitted \code{hopit} model.
#' @param newregcoef a new set of latent variable coeficients. Vector of length model$parcount[1].
#' @param data used to fit original model.
#' @param hessian logical indicating if to calculate estfun, hessian and var-cov matrix.
update.latent <-function(model, newregcoef, data, hessian=FALSE){
  coefnames <- names(model$coef)
  thresh.names <- colnames(model$thresh.mm)
  model$coef[seq_len(model$parcount[1])] <- newregcoef
  class(model) <- 'hopit'
  names(model$coef) <- coefnames
  colnames(model$thresh.mm) <- thresh.names
  p <- hopit_ExtractParameters(model)

  model$alpha <- hopit_Threshold(p$thresh.lambda, p$thresh.gamma, model)
  model$y_latent_i <- hopit_Latent(p$reg.params, model)
  model$maxobservedlatentrange <-  range(model$y_latent_i)
  model$Ey_i <- factor(colSums(sapply(1L : model$N, function(k) model$alpha[k,]<model$y_latent_i[k])),levels=1L:model$J)
  levels(model$Ey_i) <- levels(model$y_i)
  model$coef.ls <- p
  model$deviance <- -2 * model$LL

  if(hessian) {
    hes <- my.grad(fn = hopit_derivLL, par = model$coef, model=model, eps = 1e-4, collapse = TRUE, negative=FALSE)
    if (model$hasdisp && remove.theta) {
      hes <- hes[-nrow(hes),-ncol(hes)] #remove theta from vcov
      model$coef <- model$coef[-length(model$coef)] #remove from coef
    }
    model$hessian <- hes

    model$vcov.basic <- try(base::solve(-hes), silent = FALSE)
    if (class(model$vcov) == 'try-error') {
      stop(call. = NULL, 'Model is probably unidentifiable, $vcov (variance-covariance matrix) cannot be computed.')
      model$vcov.basic <- NA
    }
    if (model$hasdisp && remove.theta) COEF <- c(model$coef,model$coef.ls$logTheta) else COEF <- model$coef
    model$estfun <- hopit_derivLL(COEF, model, collapse = FALSE)
    if (remove.theta) model$estfun <- model$estfun[,-ncol(model$estfun)]

    if (length(model$design)) {
      model$vcov <- svy.varcoef.hopit(model$vcov.basic, model$estfun, design)
    } else {
      model$vcov <- model$vcov.basic
    }
  }
  if (!length(model$design)) {
    k <- 2
    model$AIC <- model$deviance + k * (length(model$coef.ls$reg.params)+
                                         length(model$coef.ls$thresh.lambda)+
                                         length(model$coef.ls$thresh.gamma)+model$hasdisp)

  } else model$AIC <- NA

  return(model)
}

#' Bootstraping hopit model
#'
#' @param model a fitted \code{Hopit} model.
#' @param data data used to fit the model
#' @param func function to be bootstrapped of the form func(model, data, ...).
#' @param nboot number of bootstrap replicates.
#' @importFrom MASS mvrnorm
#' @author Maciej J. Danko
#' @export
boot_latent_hopit<-function(model, data, func, nboot=500, ...){
  N <- seq_len(model$parcount[1])
  if (length(model$vcov)<2) stop(call.=NULL, 'No vcov detected.')
  bootsample <- MASS::mvrnorm(nboot, mu=model$coef[N], Sigma=model$vcov[N,N])
  boots <- lapply(seq_len(nboot), function(k) func(model=update.latent(model, bootsample[k,N],data=data),data=data,...))
  boots
}
