#' @noRd
#' @keywords internal
check_decreasing.levels<-function(decreasing.levels, levels_y_i){
  if (length(decreasing.levels) == 0) {
    message(paste(hopit_msg(24), toString(levels_y_i),hopit_msg(25)))
    message(hopit_msg(26))
    stop(hopit_msg(27),call.=FALSE)
  }
}


#' @noRd
#' @keywords internal
check_thresh_formula <- function(thresh.formula){
  thresh.formula <- stats::as.formula(thresh.formula)
  if (length(thresh.formula)>2L){
    warning(call. = FALSE, hopit_msg(29))
    thresh.formula[[2]] <- NULL
  }
  thresh.formula <- stats::update.formula(thresh.formula, '~.+1')
  if (any(grepl('offset(',tolower(as.character(thresh.formula[[2]])),
                fixed=TRUE))) stop(hopit_msg(31), call.=NULL)
  if (any(grepl('I(',attr(stats::terms(thresh.formula),"term.labels"),
                fixed=TRUE))) stop(hopit_msg(97), call.=NULL)
  thresh.formula
}


#' @noRd
#' @keywords internal
check_latent_formula <- function(latent.formula){
  latent.formula <- stats::as.formula(latent.formula)
  latent.formula <- stats::update.formula(latent.formula, '~.+1')
  if (any(grepl('offset(',as.character(latent.formula[[3]]),fixed=TRUE)))
    stop(hopit_msg(31), call.=NULL)
  if (any(grepl('I(',attr(stats::terms(latent.formula),"term.labels"),
                fixed=TRUE))) stop(hopit_msg(97), call.=NULL)
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
