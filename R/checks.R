#' @noRd
#' @keywords internal
check_decreasing.levels<-function(decreasing.levels, levels_y_i){
  if (length(decreasing.levels) == 0) {
    message(paste('Are the levels',toString(levels_y_i),'sorted in incresing order?'))
    message('if not please set revf to TRUE.')
    stop('"revf" must be given',call.=FALSE)
  }
}

#' @noRd
#' @keywords internal
check_hopit_method<-function(method){
  method <- tolower(method[1])
  if (method=='vglm') method <- 1 else if (method=='hopit') method <- 0 else stop('Unknown method')
  method
}

#' @noRd
#' @keywords internal
check_thresh_formula <- function(thresh.formula){
  if (length(thresh.formula)>2L){
    warning(call. = FALSE, 'The treshold formula should be given without dependent variable.')
    thresh.formula[[2]] <- NULL
  }
  thresh.formula <- update.formula(thresh.formula, '~.+1')
  if (any(grepl('offset(',tolower(as.character(thresh.formula[[2]])),fixed=TRUE))) stop('Offset not supprted.', call.=NULL)
  thresh.formula
}

#' @noRd
#' @keywords internal
check_reg_formula <- function(reg.formula){
  reg.formula <- update.formula(reg.formula, '~.+1')
  if (any(grepl('offset(',as.character(reg.formula[[3]]),fixed=TRUE))) stop('Offset not supported.', call.=NULL)
  reg.formula
}

#' @noRd
#' @keywords internal
check_vcov<-function(vcov){
  if (class(vcov) == 'try-error') {
    warning(call. = FALSE, 'Model is probably unidentifiable, $vcov (variance-covariance matrix) cannot be computed.')
    vcov <- NA
  }
  vcov
}

#' @noRd
#' @keywords internal
check_response<-function(response){
  if (!is.factor(response)) stop('Response must be a factor with ordered levels.', call.=NULL)
  if (length(levels(response))<3L) stop ('Response must have 3 or more levels.', call.=NULL)
}

#' @noRd
#' @keywords internal
check_design<-function(weights, design, N){
  if (length(weights) && length(design)) stop('Multiple weights specification detected. Please use either design or weights parameter.', call.=NULL)
  if (length(weights) && (length(weights) != N)) {
    stop('Vector of survey weights must be of the same length as data.', call.=NULL)
  }
}
