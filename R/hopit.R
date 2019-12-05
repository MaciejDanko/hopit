#' INTERNAL: Calculate the model cut-points (alpha)
#'
#' @author Maciej J. Danko
#' @keywords internal
#' @param thresh.lambda,thresh.gamma vectors with model parameters.
#' @param model a \code{hopit} object.
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_Threshold<-function(thresh.lambda, thresh.gamma, model){
  if (model$thresh.no.cov) thresh.gamma = NA
  getThresholds(model$thresh.mm,
                thresh.lambda,
                thresh.gamma,
                model$thresh.no.cov,
                thresh_start = model$control$thresh.start,
                thresh_1_exp = model$control$thresh.1.exp) #RcppEigen
}


#' INTERNAL: Calculate the predicted continuous latent measure (h_i).
#'
#' @param latent.params vectors with model parameters.
#' @param model a \code{hopit} object
#' @author Maciej J. Danko
#' @keywords internal
hopit_Latent <- function(latent.params, model = NULL)
  model$latent.mm %*% (as.matrix(latent.params))


# INTERNAL: Extract model parameters as a list
#
# Extract model parameters as a list.
# @param model a \code{hopit} object.
# @param parameters model parameters (optional). If not delivered then taken from \code{model$coef}.
# @param parcount vector with parameter counts for latent, lambda, and gamma.
# @author Maciej J. Danko
#' @noRd
#' @keywords internal
hopit_ExtractParameters <- function(model,
                                    parameters,
                                    parcount = model$parcount){
  logSigma <- 0
  if (!length(parcount)) stop(call.=NULL, hopit_msg(1))
  if (missing(parameters)) {
    parameters <- model$coef
    if (!length(parameters)) stop(hopit_msg(2))
  }
  if (length(parameters) != sum(parcount) + model$hasdisp)
    stop(call.=NULL,hopit_msg(3))

  latent.params <- parameters[1L : parcount[1L]]
  cpc <- cumsum(parcount)

  if (parcount[2L]) {
    thresh.lambda <- parameters[(cpc[1L] + 1L) : cpc[2L]]
  } else {
    stop(call.=NULL, hopit_msg(4))
  }

  if (parcount[3L]) {
    thresh.gamma <- parameters[(cpc[2L] + 1L) : cpc[3L]]
  } else {
    thresh.gamma <- NULL
  }

  if (model$hasdisp) {
    list(latent.params = latent.params,
         thresh.lambda = thresh.lambda,
         thresh.gamma = thresh.gamma,
         logSigma = parameters[length(parameters)])
  } else {
    list(latent.params = latent.params,
         thresh.lambda = thresh.lambda,
         thresh.gamma = thresh.gamma,
         logSigma = logSigma)
  }
}


#' INTERNAL: The log likelihood function
#'
#' @param parameters model parameters (optional). If the parameters are not delivered, they are taken from \code{model$coef}.
#' @param model a \code{hopit} object.
#' @param collapse a logical indicating whether to sum the individual LL contributions.
#' @param use_weights a logical indicating whether to use model weights.
#' @param negative a logical indicating whether the function should return negative.
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_negLL <- function(parameters = model$coef,
                        model,
                        collapse = TRUE,
                        use_weights,
                        negative = TRUE){
  if (missing(use_weights)) {
    if (length(model$use.weights)) use_weights <- model$use.weights else
      stop(hopit_msg(95))
  }
  link = hopit_c_link(model)
  if (collapse) {
    LL <- LLFunc(parameters,
                 yi = as.numeric(unclass(model$y_i)),
                 latent_mm = model$latent.mm,
                 thresh_mm = model$thresh.mm,
                 parcount = model$parcount,
                 hasdisp = model$hasdisp,
                 link = link,
                 thresh_no_cov = model$thresh.no.cov,
                 negative = negative,
                 thresh_1_exp = model$control$thresh.1.exp,
                 weights = model$weights,
                 use_weights = use_weights,
                 thresh_start = model$control$thresh.start,
                 out_val = model$control$LL_out_val)
  } else {
    LL <- LLFuncIndv(parameters,
                     yi = as.numeric(unclass(model$y_i)),
                     latent_mm = model$latent.mm,
                     thresh_mm = model$thresh.mm,
                     parcount = model$parcount,
                     hasdisp = model$hasdisp,
                     link = link,
                     thresh_no_cov = model$thresh.no.cov,
                     negative = negative,
                     thresh_1_exp = model$control$thresh.1.exp,
                     weights = model$weights,
                     thresh_start = model$control$thresh.start,
                     use_weights = use_weights)
  }
  if (use_weights) {
    LL <- LL / sum(model$weights) * model$N #scale likelihood
  }
  LL
}


#' INTERNAL: The gradient of the log likelihood function
#'
#' @param parameters model parameters (optional). If the parameters not delivered, they are taken from the \code{model$coef}.
#' @param model a \code{hopit} object.
#' @param collapse a logical indicating whether to sum individual LL contributions.
#' @param use_weights a logical indicating whether to use model weights.
#' @param negative  a logical indicating whether the function should return negative LL.
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_derivLL <- function(parameters=model$coef, model,
                          collapse = TRUE, use_weights, negative = FALSE){
  if (missing(use_weights)) {
    if (length(model$use.weights)) use_weights <- model$use.weights else
      stop(hopit_msg(95))
  }

  link <- hopit_c_link(model)

  if (collapse) {
    LLgr <- LLGradFunc(parameters,
                       yi = as.numeric(unclass(model$y_i)),
                       YYY1 = model$YYY1,
                       YYY2 = model$YYY2,
                       YYY3 = model$YYY3[,-model$J],
                       YYY4 = model$YYY3[,-1],
                       hasdisp = model$hasdisp,
                       latent_mm = model$latent.mm,
                       thresh_mm = model$thresh.mm,
                       thresh_extd = model$thresh.extd,
                       parcount = model$parcount,
                       link = link,
                       thresh_no_cov = model$thresh.no.cov,
                       negative = negative,
                       thresh_1_exp = model$control$thresh.1.exp,
                       weights = model$weights,
                       thresh_start = model$control$thresh.start,
                       use_weights = use_weights)
  } else {
    LLgr <- LLGradFuncIndv(parameters,
                           yi = as.numeric(unclass(model$y_i)),
                           YYY1 = model$YYY1,
                           YYY2 = model$YYY2,
                           YYY3 = model$YYY3[,-model$J],
                           YYY4 = model$YYY3[,-1],
                           latent_mm = model$latent.mm,
                           thresh_mm = model$thresh.mm,
                           thresh_extd = model$thresh.extd,
                           parcount = model$parcount,
                           hasdisp = model$hasdisp,
                           link = link,
                           thresh_no_cov = model$thresh.no.cov,
                           negative = negative,
                           thresh_1_exp = model$control$thresh.1.exp,
                           weights = model$weights,
                           thresh_start = model$control$thresh.start,
                           use_weights = use_weights)
  }

  if (use_weights) {
    LLgr <- LLgr / sum(model$weights) * model$N #scale likelihood
  }
  LLgr
}


#' INTERNAL: Fit a \code{hopit} model given the starting parameters
#'
#' Fit the model.
#' @param model a \code{hopit} object.
#' @param start starting parameters.
#' @param use_weights a logical indicating whether to use model weights.
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_fitter <- function(model, start = model$start, use_weights){

  if (missing(use_weights)) {
    if (length(model$use.weights)) use_weights <- model$use.weights else
      stop(hopit_msg(95), call.=NULL)
  }

  #be compatible with older versions
  if ((!length(model$YYY1)) || (!length(model$YYY1))){
    model <- calcYYY(model)
  }
  control <- model$control

  link <- hopit_c_link(model)

  LLgr <- function(par, neg = TRUE){
    LLGradFunc(par,
               yi = as.numeric(unclass(model$y_i)),
               YYY1 = model$YYY1,
               YYY2 = model$YYY2,
               YYY3 = model$YYY3[,-model$J],
               YYY4 = model$YYY3[,-1],
               hasdisp = model$hasdisp,
               latent_mm = model$latent.mm,
               thresh_mm = model$thresh.mm,
               thresh_extd = model$thresh.extd,
               parcount = model$parcount,
               link = link,
               thresh_no_cov = model$thresh.no.cov,
               negative = neg,
               thresh_1_exp = model$control$thresh.1.exp,
               weights = model$weights,
               thresh_start = model$control$thresh.start,
               use_weights = use_weights)
  }

  LLfn <- function(par, neg = TRUE){
    LLFunc(par,
           yi = as.numeric(unclass(model$y_i)),
           latent_mm = model$latent.mm,
           thresh_mm = model$thresh.mm,
           parcount = model$parcount,
           link = link,
           thresh_no_cov=model$thresh.no.cov,
           negative = neg,
           thresh_1_exp = model$control$thresh.1.exp,
           weights = model$weights,
           use_weights = use_weights,
           thresh_start=model$control$thresh.start,
           out_val = model$control$LL_out_val,
           hasdisp = model$hasdisp)
  }
  #Two methods one after each other
  fastgradfit <- function(fit, meto){
    #BFGS and CG method
    #try({
    if (meto == 'BFGS') {
      fit <- stats::optim(par = fit$par, fn = LLfn, gr = LLgr,
                          method = 'BFGS', hessian = FALSE,
                          control = list(maxit = model$control$bgfs.maxit,
                                         reltol = model$control$bgfs.reltol))
    } else if (meto == 'CG'){
      fit <- stats::optim(par = fit$par, fn = LLfn, gr = LLgr,
                          method = 'CG', hessian = FALSE,
                          control = list(maxit = model$control$cg.maxit,
                                         reltol = model$control$cg.reltol))
    }
    #}, silent = TRUE)

    return(fit)
  }

  #z <- try({
  z1 <- z2 <- NULL
  fit <- list(par = start)
  if ('CG' %in% control$fit.methods)
    z1 <- try({fit <- fastgradfit(fit, meto = 'CG')}, silent=TRUE)
  if ('BFGS' %in% control$fit.methods)
    z2 <- try({fit <- fastgradfit(fit, meto = 'BFGS')}, silent=TRUE)
  if (model$control$nlm.fit || !length(fit$par)) {
    if (!length(fit$par)) fit$par <- start #delete?
    if (model$control$trace) cat(hopit_msg(13),hopit_msg(5),sep='')
    fit <- suppressWarnings(stats::nlm(f = LLfn, p=fit$par,
                                       gradtol = model$control$nlm.gradtol,
                                       steptol = model$control$nlm.steptol,
                                       hessian = FALSE,
                                       iterlim=model$control$nlm.maxit))
    fit <- list(par=fit$estimate, value=fit$minimum)
  }
  fit
  #}, silent = TRUE)
  both <- length(z1) & length(z2)
  err1 <- "try-error" %in% class(z1)
  err2 <- "try-error" %in% class(z2)
  if ((both && err1 && err2) ||
      (!both && length(z1) && err1) ||
      (!both && length(z2) && err2) ||
      LLfn(fit$par) == Inf)  stop(call.=NULL, hopit_msg(6))

  model$coef <- fit$par
  model$LL <- unname(-fit$value)

  if (use_weights) {
    model$LL <- model$LL / sum(model$weights) * model$N #scale likelihood
  }
  model
}


#' Auxiliary for controlling the fitting of a \code{hopit} model
#'
#' @description
#' An auxiliary function for controlling the fitting of a \code{hopit} model.
#' Use this function to set the control
#' parameters of the \code{\link{hopit}} and other related functions.
#' @param grad.eps an epsilon parameter ("a very small number") used to calculate the Hessian from the gradient function.
#' @param nlm.fit a logical; if FALSE (default) the \code{nlm} optimization method
#' is omitted and only the BFGS and/or the CG methods are run.
#' @param bgfs.maxit,cg.maxit,nlm.maxit the maximum number of iterations.
#' See \code{\link{optim}} and \code{\link{nlm}} for details.
#' @param bgfs.reltol,cg.reltol the relative convergence tolerances for the BFGS and the CG methods.
#' See \code{\link{optim}} for details.
#' @param nlm.gradtol,nlm.steptol a tolerance at which the scaled gradient is
#' considered close enough to zero and
#' a minimum allowable relative step length for the nlm method. See \code{\link{nlm}}.
#' @param fit.methods "CG", "BFGS", or both. If both, the CG is run first, followed by the BFGS. See \code{\link{optim}}.
#' @param trace a logical for whether to trace the process of model fitting.
#' @param transform.latent,transform.thresh a type of transformation applied to
#' the all of the latent's or all of the threshold's numeric variables. Possible values:
#' \itemize{
#'   \item{"none"} {- no transformation}
#'   \item{"min"} {- subtract the minimum from a variable}
#'   \item{"scale_01"} {- transform the variable to fit the range from 0 to 1}
#'   \item{"standardize" or "standardise"} {- subtract the mean from a variable then divide it by it's standard deviation}
#'   \item{"standardize_trunc" or "standardise_trunc"} {- subtract the minimum from a variable then divide it by it's standard deviation}
#' }
#' @seealso \code{\link{hopit}}
#' @author Maciej J. Danko
#' @export
hopit.control<-function(grad.eps = 3e-5,
                        bgfs.maxit = 1e4,
                        cg.maxit = 1e4,
                        nlm.maxit = 150,
                        bgfs.reltol = 5e-10,
                        cg.reltol = 5e-10,
                        nlm.gradtol = 1e-7,
                        nlm.steptol = 1e-7,
                        fit.methods = 'BFGS',
                        nlm.fit = FALSE,
                        trace = TRUE,
                        transform.latent = 'none',
                        transform.thresh = 'none'){

  if (!length(fit.methods)) stop(hopit_msg(8),call. = NULL) else
    fit.methods <- toupper(fit.methods)
  if (any(fit.methods %notin% c('CG','BFGS'))) stop (hopit_msg(9),call.=NULL)

  list(grad.eps = grad.eps,
       bgfs.maxit = bgfs.maxit,
       cg.maxit = cg.maxit,
       nlm.maxit = nlm.maxit,
       bgfs.reltol = bgfs.reltol,
       cg.reltol = cg.reltol,
       nlm.gradtol = nlm.gradtol,
       nlm.steptol = nlm.steptol,
       fit.methods = fit.methods,
       nlm.fit = nlm.fit,
       trace = trace,
       transform.latent = transform.latent,
       transform.thresh = transform.thresh)
}


#' Generalized hierarchical ordered threshold models.
#'
#' @description
#' The ordered response data classify a measure of interest into ordered categories
#' collected during a survey. For example, if the dependent variable is a happiness
#' rating, a respondent typically answers a question such as: “Taking all things
#' together, would you say you are ... ?" and then selects from response options
#' along the lines of: "very happy", "pretty happy", "not too happy", and "very unhappy"
#' \insertCite{Liao2005}{hopit}. Similarly, if interviewees are asked to evaluate their
#' health in general (e.g., “Would you say your health is ... ?”) they, can typically choose among
#' several categories, such as "very good", "good", "fair", "bad", and "very bad"
#' \insertCite{King2004,Jurges2007,Rebelo2014,OKSUZYAN2019}{hopit}. In political science, a respondent
#' may be asked for an opinion about recent legislation (e.g. “Rate your feelings about
#' the proposed legislation.") and asked to choose among categories like: "strongly
#' oppose", "mildly oppose", "indifferent", "mildly support", and "strongly support"
#' \insertCite{GreeneHensher2010}{hopit}. It is easy to imagine other multi-level ordinal
#' variables that might be used during a survey and to which the methodology described
#' below could be applied.\cr
#'
#' In practice, it is assumed that when responding to a survey question about their general
#' happiness, health, feelings, attitudes or other status, participants are
#' assessing their true value of this unobserved continuous variable, and
#' project it onto the discrete scale provided. The thresholds that individuals
#' use to categorize their true status by selecting a specific response option
#' may be affected by the reference group chosen, their earlier life experiences,
#' and cross-cultural differences in using scales. Thus, the responses of
#' individuals may differ depending on their gender, age, cultural background,
#' education, and personality traits; among other factors
#' \insertCite{King2004,Jurges2007,OKSUZYAN2019}{hopit}.\cr

#' From the perspective of reporting behavior modeling, one of the main tasks
#' researchers face is to compute this continuous estimate of the underlying,
#' latent measures of individuals based on several specific characteristics
#' of the responses considered (e.g., health variables or happiness variables),
#' and to account for variations in reporting across socio-demographic and
#' cultural groups. More specifically, to build a latent, underlying measure,
#' a generalized hierarchical ordered threshold model is fitted that regresses
#' the reported status/attitude/feeling on two sets of independent variables
#' \insertCite{Boes2006,Green2014}{hopit}. When the dependent reported ordered
#' variable is self-rated health status, then the first set of variables –
#' i.e., health variables – assess specific aspects of individuals’ health,
#' such as measures of chronic conditions, mobility, difficulties with a range
#' of daily activities, grip strength, anthropometric characteristics, and
#' lifestyle behaviors. Using the second set of independent variables
#' (threshold variables), the model also adjusts for differences across
#' socio-demographic and cultural groups, such as differences in cultural
#' background, gender, age, and education
#' \insertCite{King2004,Jurges2007,OKSUZYAN2019}{hopit}.\cr
#'
#' Ordered threshold models are used to fit ordered categorical dependent variables.
#' The generalized ordered threshold models \insertCite{Terza1985,Boes2006,Green2014}{hopit}
#' are an extension of the ordered threshold models \insertCite{McKelvey1975}{hopit}.
#' Whereas in the latter models, the thresholds are constant, in the generalized models the
#' thresholds are allowed to be dependent on covariates.
#' \insertCite{GreeneHensher2010,Green2014;textual}{hopit} pointed out that for a
#' model to make sense, the thresholds must also be ordered.
#' This observation motivated Greene and coauthors to call these models *HOPIT*, which stands
#' for hierarchical ordered probit models.
#'
#' The fitted *hopit* model is used to analyze heterogeneity in reporting behavior.
#' See \code{\link{standardizeCoef}}, \code{\link{latentIndex}},
#' \code{\link{getCutPoints}}, \code{\link{getLevels}}, and \code{\link{boot_hopit}}.
#' @details
#' The function fits generalized hierarchical ordered threshold models.\cr
#'
#' \code{latent.formula} models the latent variable.
#' If the response variable is self-rated health, then the latent measure can depend on different health
#' conditions and diseases (latent variables are called health variables).
#' Latent variables are modeled with the parallel regression assumption. According to this assumption, the coefficients
#' that describe the relationship between the lowest response category and all of the higher response categories, are the same as the coefficients
#' that describe the relationship between another (e.g., adjacent) lowest response category and the remaining higher response categories.
#' The predicted latent variable is modeled as a linear function of the health variables and the corresponding coefficients.\cr
#'
#' \code{thresh.formula} models the threshold variable.
#' The thresholds (cut-points, \code{alpha}) are modeled by the threshold variables \code{gamma} and the intercepts \code{lambda}.
#' It is assumed that they model the contextual characteristics of the respondent (e.g., country, gender, and age).
#' The threshold variables are modeled without the parallel regression assumption; thus, each threshold is modeled by
#' a variable independently \insertCite{Boes2006,Green2014}{hopit}.
#' The \code{hopit}() function uses the parameterization of thresholds proposed by \insertCite{Jurges2007;textual}{hopit}.\cr
#'
#' \code{decreasing.levels} it is the logical that determines the ordering of the levels of the categorical response variable.
#' It is always advisable to first check the ordering of the levels before starting (see example 1)\cr
#'
#' It is possible to model the interactions, including interactions between the latent and the threshold variables. The interactions added to the latent formula
#' only model the latent measure, and the interactions modeled in the threshold formula only model the thresholds.
#' The general rule for modeling any kind of interaction is to use "*" to specify interactions within a latent (or threshold) formula and to
#' use ':' to specify interactions between the latent and the threshold variables. In the latter case, the main effects of an interaction must also be specified;
#' i.e., the main latent effects must be specified in the latent formula, and the main threshold effect must be speciffied in the threshold formula.
#' See also \code{Example 3} below.\cr
#'
#' For more details, please see the package vignette, which is also available under this link:
#' \href{https://github.com/MaciejDanko/hopit/blob/master/vignettes/vig_hopit.pdf}{vig_hopit.pdf}
#'
#' @param latent.formula a formula used to model the latent variable. It should not contain any threshold variable.
#' To specify the interactions between the latent and the threshold variables, see details.
#' @param thresh.formula a formula used to model the threshold variable. It should not contain any latent variable.
#' To specify interactions between the latent and the threshold variables, see details.
#' Any dependent variable (left side of "~" in the formula) will be ignored.
#' @param data a data frame that includes all modeled variables.
#' @param decreasing.levels a logical indicating whether self-reported health classes are ordered in decreasing order.
#' @param fit.sigma a logical indicating whether to fit an additional parameter sigma,
#' which models a standard deviation of the error term (e.g., the standard deviation of the cumulative normal distribution in the probit model).
#' @param design an optional survey design. Use the \code{\link[survey]{svydesign}} function to specify the design.
#' The design cannot be specified together with parameter \code{weights}.
#' @param weights optional model weights. Use the design to construct survey weights.
#' @param link a link function. The possible values are \code{"probit"} (default) and \code{"logit"}.
#' @param start a vector with starting coefficient values in the form \code{c(latent_parameters, threshold_lambdas, threshold_gammas)} or
#' \code{c(latent_parameters, threshold_lambdas, threshold_gammas, logSigma)} if the \code{fit.sigma == TRUE}.
#' @param control a list with control parameters. See \code{\link{hopit.control}}.
#' @param na.action a function that indicates what should happen when the \code{data} contain \code{NA}s.
#' The default is \code{\link[stats]{na.fail}},
#' which generates an error if any missing value is found. The alternative is \code{\link[stats]{na.omit}}
#' (or \code{\link[stats]{na.exclude}} equivalently), which removes rows with missing
#' values from the \code{data}. Using \code{\link[stats]{na.pass}} will lead to an error.
#' @importFrom stats na.fail
#' @importFrom stats sigma
#' @importFrom stats model.frame
#' @importFrom stats coef
#' @importFrom Rdpack reprompt
#' @return a \code{hopit} object used by other functions and methods. The object is a list with the following components:
#'  \item{control}{ a list with control parameters. See \code{\link{hopit.control}}.}
#'  \item{link}{ a link function used.}
#'  \item{hasdisp}{ a logical indicating whether fit.sigma was modeled.}
#'  \item{use.weights}{ a logical indicating whether any weights were used.}
#'  \item{weights}{ a vector with model weights.}
#'  \item{frame}{ a model frame.}
#'  \item{latent.formula}{ a latent formula used to fit the model.}
#'  \item{latent.mm}{ a latent model matrix.}
#'  \item{latent.terms}{ latent variables used, and their interactions.}
#'  \item{cross.inter.latent}{ a part of the latent formula used for modeling cross-interactions in the latent model}
#'  \item{thresh.formula}{ a threshold formula used to fit the model.}
#'  \item{thresh.mm}{ a threshold model matrix.}
#'  \item{thresh.extd}{ an extended threshold model matrix.}
#'  \item{thresh.terms}{ threshold variables used, and their interactions.}
#'  \item{cross.inter.thresh}{ a part of the threshold formula used for modeling cross-interactions in the threshold model}
#'  \item{thresh.no.cov}{ a logical indicating whether gamma parameters are present.}
#'  \item{parcount}{ a 3-element vector with a number of parameters for the latent variables (beta),
#'  the threshold intercepts (lambda), and the threshold covariates (gamma).}
#'  \item{coef}{ a vector with model coefficients.}
#'  \item{coef.ls}{ model coefficients as a list.}
#'  \item{start}{ a vector with the starting values of the coefficients.}
#'  \item{alpha}{ estimated individual-specific thresholds.}
#'  \item{y_i}{ a vector with individual responses - the response variable.}
#'  \item{y_latent_i}{ a vector with predicted latent measures for each individual.}
#'  \item{Ey_i}{ a vector with predicted categorical responses for each individual.}
#'  \item{J}{ a number of response levels.}
#'  \item{N}{ a number of observations.}
#'  \item{deviance}{ a deviance.}
#'  \item{LL}{ a log likelihood.}
#'  \item{AIC}{ an AIC for models without a survey design.}
#'  \item{vcov}{ a variance-covariance matrix.}
#'  \item{vcov.basic}{ a variance-covariance matrix that ignores the survey design.}
#'  \item{hessian}{ a Hessian matrix.}
#'  \item{estfun}{ a gradient (a vector of partial derivatives) of the log likelihood function at the estimated coefficient values.}
#'  \item{YYY1,YYY2,YYY3}{ an internal objects used for the calculation of gradient and Hessian functions.}
#' @references \insertAllCited{}
#' @export
#' @author Maciej J. Danko
#' @seealso
#' \code{\link{coef.hopit}},
#' \code{\link{profile.hopit}},
#' \code{\link{hopit.control}},
#' \code{\link{anova.hopit}},
#' \code{\link{vcov.hopit}},
#' \code{\link{logLik.hopit}},
#' \code{\link{AIC.hopit}},
#' \code{\link{summary.hopit}},
#' \code{\link[survey]{svydesign}}, \cr\cr
#' For heterogeneity in reporting behavior analysis see:\cr
#' \code{\link{standardizeCoef}},
#' \code{\link{latentIndex}},
#' \code{\link{getCutPoints}},
#' \code{\link{getLevels}},
#' \code{\link{boot_hopit}},
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # first determine the order of the levels of the dependent variable
#' levels(healthsurvey$health)
#'
#' # the order of response levels decreases from the best health to
#' # the worst health; hence the hopit() parameter decreasing.levels
#' # is set to TRUE
#'
#' # Example 1 ---------------------
#'
#' # fitting the model:
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # summarize the fit:
#' summary(model1)
#'
#' # extract parameters in the form of a list
#' cm1 <- coef(model1, aslist = TRUE)
#'
#' # names of the returned coefficients
#' names(cm1)
#'
#' # extract the latent health coefficients
#' cm1$latent.params
#'
#' # check the fit
#' \donttest{
#' profile(model1)
#' }
#' # Example 2 ---------------------
#'
#' \donttest{
#' # incorporate the survey design
#' design <- svydesign(ids = ~ country + psu, weights = healthsurvey$csw,
#' data = healthsurvey)
#'
#' model2 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                   heart_attack_or_stroke + poor_mobility +
#'                   very_poor_grip + depression + respiratory_problems +
#'                   IADL_problems + obese + diabetes + other_diseases,
#'                 thresh.formula = ~ sex + ageclass + country,
#'                 decreasing.levels = TRUE,
#'                 design = design,
#'                 control = list(trace = FALSE),
#'                 data = healthsurvey)
#'
#' # compare the latent variables
#' cbind('No survey design' = coef(model1, aslist = TRUE)$latent.par,
#' 'Has survey design' = coef(model2, aslist = TRUE)$latent.par)
#' }
#' \donttest{
#' # Example 3 ---------------------
#'
#' # defining the interactions between the threshold and the latent variables
#'
#' # correctly defined interactions:
#' model3 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility * very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases +
#'                 sex : depression + sex : diabetes + ageclass:obese,
#'               thresh.formula = ~ sex * ageclass + country + sex : obese,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#' }
#' \dontrun{
#' # badly defined interactions:
#'
#' # 1) lack of a main effect of "other_diseases" in any formula
#' # it can be solved by adding " + other_diseases" to the latent formula
#' model3a <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases : sex,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # 2) the main effect of sex is present in both formulas.
#' # it can be solved by replacing "*" with ":" in "other_diseases * sex"
#' model3b <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases * sex,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' }
#' # Example 4 ---------------------
#'
#' \donttest{
#' # construct a naive continuous variable:
#' hs <- healthsurvey
#' hs$cont_var <- sample(5000:5020,nrow(hs),replace=TRUE)
#'
#' latent.formula = health ~ hypertension + high_cholesterol +
#'   heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'   depression + respiratory_problems +
#'   IADL_problems + obese + diabetes + other_diseases
#'
#' # in some cases, when continuous variables are used, the hopit:::get.hopit.start() function
#' # do not find starting parameters (R version 3.4.4 (2018-03-15)):
#' \dontrun{
#' model4 <- hopit(latent.formula = latent.formula,
#'                 thresh.formula = ~ sex + cont_var,
#'                 decreasing.levels = TRUE,
#'                 data = hs)
#' }
#' # one of the solutions is to transform one or more continuous variables:
#' hs$cont_var_t <- hs$cont_var-min(hs$cont_var)
#'
#' model4b <- hopit(latent.formula = latent.formula,
#'                  thresh.formula = ~ sex + cont_var_t,
#'                  decreasing.levels = TRUE,
#'                  data = hs)
#'
#' # this can also be done automatically using the the control parameter
#' model4c <- hopit(latent.formula = latent.formula,
#'                  thresh.formula = ~ sex + cont_var,
#'                  decreasing.levels = TRUE,
#'                  control = list(transform.thresh = 'min',
#'                                 transform.latent = 'none'),
#'                  data = hs)
#'
#' model4d <- hopit(latent.formula = latent.formula,
#'                  thresh.formula = ~ sex + cont_var,
#'                  decreasing.levels = TRUE,
#'                  control = list(transform.thresh = 'scale_01',
#'                                 transform.latent = 'none'),
#'                  data = hs)
#'
#' model4e <- hopit(latent.formula = latent.formula,
#'                  thresh.formula = ~ sex + cont_var,
#'                  decreasing.levels = TRUE,
#'                  control = list(transform.thresh = 'standardize',
#'                                 transform.latent = 'none'),
#'                  data = hs)
#'
#' model4f <- hopit(latent.formula = latent.formula,
#'                  thresh.formula = ~ sex + cont_var,
#'                  decreasing.levels = TRUE,
#'                  control = list(transform.thresh = 'standardize_trunc',
#'                                 transform.latent = 'none'),
#'                  data = hs)
#'
#' round(t(rbind(coef(model4b),
#'               coef(model4c),
#'               coef(model4d),
#'               coef(model4e),
#'               coef(model4f))),4)
#'
#' }
hopit<- function(latent.formula,
                 thresh.formula = ~ 1,
                 data,
                 decreasing.levels,
                 start = NULL,
                 fit.sigma = FALSE,
                 design = list(),
                 weights = NULL,
                 link = c('probit', 'logit'),
                 control = list(),
                 na.action = na.fail){

  if (!fit.sigma) remove.sigma = FALSE else remove.sigma = TRUE
  if (missing(data)) data <- environment(latent.formula)
  data <- na.action(data)
  link <- match.arg(link)
  control <- do.call("hopit.control", control)
  control$thresh.start <- control$LL_out_val <- -Inf;
  control$thresh.1.exp <- FALSE;

  model <- NULL
  model$control <- control
  model$link <- link[1]
  model$hasdisp <- fit.sigma
  model$na.action <- na.action

  thresh.formula <- check_thresh_formula(thresh.formula, data)
  latent.formula <- check_latent_formula(latent.formula, data)

  data <- drop.levels.response(latent.formula, data)
  check_response(stats::model.response(
    stats::model.frame(latent.formula, data)))

  if (control$transform.latent != 'none')
    data <- transform.data(latent.formula, data, control$transform.latent)
  if (control$transform.thresh != 'none')
    data <- transform.data(thresh.formula, data, control$transform.thresh)

  data <- drop.levels.data(latent.formula, data)
  data <- drop.levels.data(thresh.formula, data)

  model <- analyse.formulas(model, latent.formula, thresh.formula, data)

  #model$y_i <- stats::model.frame(model$latent.formula,
  #                          data = data)[,all.vars(model$latent.formula[[2]])]
  model$y_i <- stats::model.response(
    stats::model.frame(model$latent.formula, data = data))
  #check_response(model$y_i) #already tested so maybe remove this

  if (missing(decreasing.levels)) decreasing.levels = NULL
  check_decreasing.levels(decreasing.levels, levels(model$y_i))
  if (!decreasing.levels) model$y_i <- factor(model$y_i, rev(levels(model$y_i)))
  model$decreasing.levels <- decreasing.levels
  model$y_latent_i <- NA # latent
  model$Ey_i <- NA # ordinal classified utput
  model$J <- length(levels(model$y_i))
  model$N <- length(model$y_i)

  model$thresh.extd <- matrix(rep_row(model$thresh.mm, model$J-1),model$N,
                              NCOL(model$thresh.mm)*(model$J-1))

  Cr <- dim(model$latent.mm)[2L]
  Ct <- dim(model$thresh.mm)[2L]
  if (model$thresh.no.cov) Ct <- 0L
  model$parcount <- c(Cr, model$J - 1L, Ct*(model$J - 1L))
  model$parcount[3L] <- model$parcount[3L]*(model$thresh.no.cov == FALSE)
  interce <- paste(1L : (model$J - 1L), 2L : (model$J), sep = '|')
  if (model$thresh.no.cov){
    tmp <- NULL
  } else {
    tmp <- as.matrix(expand.grid(interce, model$thresh.names,
                                 KEEP.OUT.ATTRS = FALSE))
    tmp <- tmp[,c(2,1)]
    tmp <- paste('(G)', apply(tmp, 1L, paste, sep = '', collapse = '.'),
                 sep = '.')
  }

  coefnames <-  c(model$latent.names, paste('(L)', interce, sep = '.'), tmp)
  if (model$hasdisp) coefnames <- c(coefnames, 'logSigma')

  model$weights <- NULL
  check_design(weights, design, model$N)
  model$design <- design
  if (length(design)) {
    model$weights <- design$prob
  } else if (length(weights))
    model$weights <- weights
  if (!length(model$weights)) {
    model$weights <- rep(1, model$N)
    model$use.weights <- FALSE
  } else model$use.weights <- TRUE

  model$weights <- as.vector(matrix(model$weights, 1L, model$N))
  #scaling weights
  model$weights <- model$N * model$weights / sum(model$weights)

  model$frame <- cbind.data.frame(stats::model.frame(latent.formula, data),
                                  stats::model.frame(thresh.formula, data))

  #calculate special matrices for gradient calaculation
  model <- calcYYY(model)

  if (!length(start)) {
    if (model$control$trace) cat(hopit_msg(10))
    model <- suppressWarnings(get.hopit.start(model, data))
    if (model$control$trace) cat(hopit_msg(13))
    if (any(model$parcount!=sapply(model$glm.start.ls,length)) ||
        !length(model$glm.start.ls)) {
      stop(paste(hopit_msg(96),hopit_msg(103)),call. = NULL)
    }
  } else {
    model$start <- start
  }

  if ((sum(model$parcount) + model$hasdisp) != length(model$start))
    stop(hopit_msg(96),call. = NULL)

  if (model$control$trace) cat(hopit_msg(11))
  model <- hopit_fitter(model, start = model$start)
  class(model) <- 'hopit'
  #colnames(model$thresh.mm) <- model$thresh.names
  names(model$coef) <- coefnames
  if (model$control$trace) cat(hopit_msg(13),hopit_msg(12),sep='')

  p <- hopit_ExtractParameters(model)
  model$alpha <- hopit_Threshold(p$thresh.lambda, p$thresh.gamma, model)
  model$Ey_i <- classify.ind(model)
  model$y_latent_i <- hopit_Latent(model$coef[seq_len(model$parcount[1])],
                                   model)
  model$coef.ls <- p
  model$deviance <- -2 * model$LL
  if (model$control$trace) cat(hopit_msg(13),hopit_msg(14),sep='')

  hes <- my.grad(fn = hopit_derivLL, par = model$coef, model=model, eps = 1e-4,
                 collapse = TRUE, negative=FALSE)
  if (model$hasdisp && remove.sigma) {
    hes <- hes[-nrow(hes),-ncol(hes)] #remove sigma from vcov
    model$coef <- model$coef[-length(model$coef)] #remove from coef
  }
  model$hessian <- hes

  model$vcov.basic <- try(base::solve(-hes), silent = FALSE)
  model$vcov.basic <- check_vcov(model$vcov.basic)
  if (model$control$trace) cat(hopit_msg(13),hopit_msg(15),sep='')
  if (model$hasdisp && remove.sigma) COEF <-
    c(model$coef,model$coef.ls$logSigma) else COEF <- model$coef
  model$estfun <- hopit_derivLL(COEF, model, collapse = FALSE)
  if (remove.sigma) model$estfun <- model$estfun[,-ncol(model$estfun)]
  if (model$control$trace) cat(hopit_msg(13))

  if (length(model$design)) {
    if (model$control$trace) cat(hopit_msg(16))
    model$vcov <- svy.varcoef_hopit(model$vcov.basic, model$estfun, design)
    model$AIC <- NA
    if (model$control$trace) cat(hopit_msg(13))
  } else {
    k <- 2
    model$vcov <- model$vcov.basic
    model$AIC <- model$deviance + k * (length(model$coef.ls$latent.params) +
                                       length(model$coef.ls$thresh.lambda) +
                                       length(model$coef.ls$thresh.gamma) +
                                         model$hasdisp)
  }
  model$glm.start <- model$glm.start.ls <- NULL
  return(model)
}
