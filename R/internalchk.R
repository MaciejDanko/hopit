#' @noRd
#' @keywords internal
check_decreasing.levels<-function(decreasing.levels, levels_y_i){
  if (length(decreasing.levels) != 1) {
    message(paste(hopit_msg(24), toString(levels_y_i),hopit_msg(25)))
    message(hopit_msg(26))
    stop(hopit_msg(27),call.=FALSE)
  }
}


#' @noRd
#' @keywords internal
check_thresh_formula <- function(thresh.formula, data){
  thresh.formula <- stats::as.formula(thresh.formula)
  thresh.formula <- stats::update.formula(thresh.formula, '~.+1')
  MF <- stats::model.frame(thresh.formula, data)
  if (length(stats::model.offset(MF))) stop(hopit_msg(31), call.=NULL)
  if (length(stats::model.response(MF))) {
    warning(call. = FALSE, hopit_msg(29))
    thresh.formula[[2]] <- NULL
  }
  LT <- attr(stats::terms(thresh.formula),"term.labels")
  if (any(grepl('I(',LT,fixed=TRUE))) stop(hopit_msg(97), call.=NULL)
  thresh.formula
}


#' @noRd
#' @keywords internal
check_latent_formula <- function(latent.formula, data){
  latent.formula <- stats::as.formula(latent.formula)
  latent.formula <- stats::update.formula(latent.formula, '~.+1')
  MF <- stats::model.frame(latent.formula, data)
  if (!ncol(MF)) {
    stop(paste(hopit_msg(100), hopit_msg(101),sep='\n'), call.=NULL)
  }
  if (length(stats::model.offset(MF))) stop(hopit_msg(31), call.=NULL)
  LT <- attr(stats::terms(latent.formula),"term.labels")
  if (!length(LT)) stop(hopit_msg(100), call.=NULL)
  if (any(grepl('I(',LT,fixed=TRUE))) stop(hopit_msg(97), call.=NULL)
  latent.formula
}


#' @noRd
#' @keywords internal
check_vcov<-function(vcov){
  if ('try-error' %in% class(vcov)) {
    warning(call. = FALSE, hopit_msg(32))
    vcov <- NA
  }
  vcov
}


#' @noRd
#' @keywords internal
check_response<-function(response){
  if (!length(response)) stop(hopit_msg(101), call.=NULL)
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
