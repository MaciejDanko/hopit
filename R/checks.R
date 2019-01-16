#' @noRd
#' @keywords internal
check_decreasing.levels<-function(decreasing.levels, levels_y_i){
  if (length(decreasing.levels) == 0) {
    message(paste(hopit_msg(24),toString(levels_y_i),hopit_msg(25)))
    message(hopit_msg(26))
    stop(hopit_msg(27),call.=FALSE)
  }
}

#' @noRd
#' @keywords internal
check_hopit_method<-function(method){
  method <- tolower(method[1])
  if (method=='vglm') method <- 1 else if (method=='hopit') method <- 0 else stop(call.=NULL, hopit_msg(28))
  method
}

#' @noRd
#' @keywords internal
check_thresh_formula <- function(thresh.formula){
  thresh.formula <- as.formula(thresh.formula)
  if (length(thresh.formula)>2L){
    warning(call. = FALSE, hopit_msg(29))
    thresh.formula[[2]] <- NULL
  }
  thresh.formula <- update.formula(thresh.formula, '~.+1')
  if (any(grepl('offset(',tolower(as.character(thresh.formula[[2]])),fixed=TRUE))) stop(hopit_msg(31), call.=NULL)
  thresh.formula
}

#' @noRd
#' @keywords internal
check_strata_formula <- function(strata.formula){
  strata.formula <- as.formula(strata.formula)
  if (length(strata.formula)>2L){
    warning(call. = FALSE, hopit_msg(30))
    strata.formula[[2]] <- NULL
  }
  strata.formula <- update.formula(strata.formula, '~.+1')
  if (any(grepl('offset(',tolower(as.character(strata.formula[[2]])),fixed=TRUE))) stop(hopit_msg(31), call.=NULL)
  strata.formula
}

#' @noRd
#' @keywords internal
check_latent_formula <- function(latent.formula){
  latent.formula <- as.formula(latent.formula)
  latent.formula <- update.formula(latent.formula, '~.+1')
  if (any(grepl('offset(',as.character(latent.formula[[3]]),fixed=TRUE))) stop(hopit_msg(31), call.=NULL)
  latent.formula
}

#' @noRd
#' @keywords internal
check_vcov<-function(vcov){
  if (class(vcov) == 'try-error') {
    warning(call. = FALSE, hopit_msg(32))
    vcov <- NA
  }
  vcov
}

#' @noRd
#' @keywords internal
check_response<-function(response){
  if (!is.factor(response)) stop(hopit_msg(33), call.=NULL)
  if (length(levels(response))<3L) stop (hopit_msg(34), call.=NULL)
}

#' @noRd
#' @keywords internal
check_design<-function(weights, design, N){
  if (length(weights) && length(design)) stop(hopit_msg(35), call.=NULL)
  if (length(weights) && (length(weights) != N)) {
    stop(hopit_msg(36), call.=NULL)
  }
}
