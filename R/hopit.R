#' INTERNAL: Calculate model cut-points (thresholds)
#'
#' @author Maciej J. Danko
#' @keywords internal
#' @param thresh.lambda,thresh.gamma vectors with model parameters.
#' @param model a \code{hopit} object.
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_Threshold<-function(thresh.lambda, thresh.gamma, model){
  if (model$thresh.no.cov) thresh.gamma = NA
  getThresholds(model$thresh.mm, thresh.lambda, thresh.gamma, model$thresh.no.cov,
                thresh_start = model$control$thresh.start, thresh_1_exp = model$control$thresh.1.exp) #RcppEigen
}

#' INTERNAL: Calculate predicted continous latent measure.
#'
#' @param latent.params vectors with model parameters.
#' @param model a \code{hopit} object
#' @author Maciej J. Danko
#' @keywords internal
hopit_Latent <- function(latent.params, model = NULL) model$latent.mm %*% (as.matrix(latent.params))


# #' INTERNAL: Calculate maximum possible latent range
# #' @param model a fitted \code{hopit} model.
# #' @param data a data used to fit the model.
# #' @keywords internal
# hopit_latentrange <- function (model, data) {
#
#   cfm <- model$coef[seq_len(model$parcount[1])]
#
#   L <- apply(model$latent.mm,2,function(x)length(unique(x))) -1
#   FACTORS <- as.logical(L==1)
#   if (any(FACTORS)) {
#     d <- cfm[FACTORS]
#     r <- c(sum(d[d<0]), sum(d[d>0]))
#   }
#   if (any(!FACTORS)) {
#     #message('Continuous covariate detected, model may not give realistic results.',appendLF = FALSE)
#     d <- matrix(cfm[!FACTORS],length(cfm[!FACTORS]),model$N) * t(model$latent.mm[,!FACTORS])
#     r <- c(min(d[d<0],r), max(r,d[d>0]))
#   }
#   r
# }


#' INTERNAL: Extract model parameters as a list
#'
#' Extract model parameters as a list.
#' @param model \code{hopit} object.
#' @param parameters model parameters (optional). If not delivered then taken from \code{model$coef}.
#' @param parcount vector with parameters counts for latent, lambda, and gamma.
#' @author Maciej J. Danko
#' @keywords internal
hopit_ExtractParameters <- function(model, parameters, parcount = model$parcount){
  logTheta <- 0
  if (!length(parcount)) stop(call.=NULL, hopit_msg(1))
  if (missing(parameters)) {
    parameters <- model$coef
    if (!length(parameters)) stop(hopit_msg(2))
  }
  if (length(parameters) != sum(parcount) + model$hasdisp) stop(call.=NULL,hopit_msg(3))

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
    list(latent.params = latent.params, thresh.lambda = thresh.lambda, thresh.gamma = thresh.gamma, logTheta = parameters[length(parameters)])
  } else {
    list(latent.params = latent.params, thresh.lambda = thresh.lambda, thresh.gamma = thresh.gamma, logTheta = logTheta)
  }
}


#' INTERNAL: The log likelihood function
#'
#' @param parameters model parameters (optional). If not delivered then taken from \code{model$coef}.
#' @param model \code{hopit} object.
#' @param collapse a logical indicating if to sum individual LL contributions.
#' @param use_weights a logical indicating if to use model weights.
#' @param negative a logical indicating if the function should return negative LL or not.
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_negLL <- function(parameters=model$coef, model, collapse = TRUE, use_weights = model$use.weights, negative = TRUE){
  link = hopit_c_link(model)
  if (collapse) {
    LL <- LLFunc(parameters, yi=as.numeric(unclass(model$y_i)),latent_mm=model$latent.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
                 hasdisp = model$hasdisp,
                 link=link,thresh_no_cov=model$thresh.no.cov, negative=negative, thresh_1_exp = model$control$thresh.1.exp,
                 weights=model$weights,use_weights = use_weights, thresh_start=model$control$thresh.start, out_val = model$control$LL_out_val)
  } else {
    LL <- LLFuncIndv(parameters, yi=as.numeric(unclass(model$y_i)),latent_mm=model$latent.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
                     hasdisp = model$hasdisp,
                     link=link,thresh_no_cov=model$thresh.no.cov, negative=negative, thresh_1_exp = model$control$thresh.1.exp,
                     weights=model$weights, thresh_start = model$control$thresh.start, use_weights = use_weights)
  }
  if (use_weights) {
    LL <- LL / sum(model$weights) * model$N #scale likelihood
  }
  LL
}


#' INTERNAL: The gradient of the log likelihood function
#'
#' @param parameters model parameters (optional). If not delivered then taken from \code{model$coef}.
#' @param model \code{hopit} object.
#' @param collapse a logical indicating if to sum individual LL contributions.
#' @param use_weights a logical indicating if to use model weights.
#' @param negative  a logical indicating if the function should return negative LL or not.
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_derivLL <- function(parameters=model$coef, model,
                          collapse = TRUE, use_weights = model$use.weights, negative = FALSE){
  link = hopit_c_link(model)

  if (collapse) {
    LLgr <- LLGradFunc(parameters, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2,
                       YYY3=model$YYY3[,-model$J],YYY4=model$YYY3[,-1], hasdisp = model$hasdisp,
                       latent_mm=model$latent.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
                       link=link,thresh_no_cov=model$thresh.no.cov,negative=negative, thresh_1_exp = model$control$thresh.1.exp,
                       weights=model$weights, thresh_start = model$control$thresh.start, use_weights = use_weights)
  } else {
    LLgr <- LLGradFuncIndv(parameters, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2, YYY3=model$YYY3[,-model$J], YYY4=model$YYY3[,-1],
                           latent_mm=model$latent.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
                           hasdisp = model$hasdisp,
                           link=link,thresh_no_cov=model$thresh.no.cov, negative=negative, thresh_1_exp = model$control$thresh.1.exp,
                           weights=model$weights,thresh_start = model$control$thresh.start, use_weights = use_weights)
  }

  if (use_weights) {
    LLgr <- LLgr / sum(model$weights) * model$N #scale likelihood
  }
  LLgr
}


#' INTERNAL: Fit \code{hopit} model given a starting parameters
#'
#' Fit the model.
#' @param model a \code{hopit} object.
#' @param start starting parameters.
#' @param use_weights a logical indicating if to use model weights.
#' @author Maciej J. Danko
#' @keywords internal
#' @useDynLib hopit
#' @importFrom Rcpp evalCpp
hopit_fitter <- function(model, start = model$start, use_weights = model$use.weights){

  #be compatible with older versions
  if ((!length(model$YYY1)) || (!length(model$YYY1))){
    model <- calcYYY(model)
  }
  control <- model$control

  link <- hopit_c_link(model)

  LLgr <- function(par, neg = TRUE) LLGradFunc(par, yi=as.numeric(unclass(model$y_i)), YYY1=model$YYY1, YYY2=model$YYY2, YYY3=model$YYY3[,-model$J],
                                               YYY4=model$YYY3[,-1],hasdisp = model$hasdisp,
                                               latent_mm=model$latent.mm, thresh_mm=model$thresh.mm, thresh_extd=model$thresh.extd, parcount=model$parcount,
                                               link=link,thresh_no_cov=model$thresh.no.cov, negative=neg, thresh_1_exp = model$control$thresh.1.exp,
                                               weights=model$weights, thresh_start=model$control$thresh.start, use_weights = use_weights)

  LLfn <- function(par, neg = TRUE) LLFunc(par, yi=as.numeric(unclass(model$y_i)),latent_mm=model$latent.mm, thresh_mm=model$thresh.mm, parcount=model$parcount,
                                           link=link, thresh_no_cov=model$thresh.no.cov, negative=neg, thresh_1_exp = model$control$thresh.1.exp,
                                           weights=model$weights,use_weights = use_weights, thresh_start=model$control$thresh.start,
                                           out_val = model$control$LL_out_val, hasdisp = model$hasdisp)
  #Two methods one after each other
  fastgradfit <- function(fit, meto){
    #BFGS and CG method
    try({
      if (meto == 'BFGS') {
        fit <- stats::optim(par = fit$par, fn = LLfn, gr = LLgr,
                     method = 'BFGS', hessian = FALSE, control=list(maxit=model$control$bgfs.maxit, reltol=model$control$bgfs.reltol))
      } else if (meto == 'CG'){
        fit <- stats::optim(par = fit$par, fn = LLfn, gr = LLgr,
                     method = 'CG', hessian = FALSE, control=list(maxit=model$control$cg.maxit, reltol=model$control$cg.reltol))
      }
    }, silent = FALSE)

    return(fit)
  }

  z <- try({
    fit <- list(par = start)
    if ('CG' %in% control$fit.methods) fit <- fastgradfit(fit, meto = 'CG')
    if ('BFGS' %in% control$fit.methods) fit <- fastgradfit(fit, meto = 'BFGS')
    if (!model$control$quick.fit || !length(fit$par)) {
      if (!length(fit$par)) fit$par <- start
      if (model$control$trace) cat(hopit_msg(13),hopit_msg(5),sep='')
      fit <- suppressWarnings(stats::nlm(f = LLfn, p=fit$par, gradtol = model$control$nlm.gradtol, steptol = model$control$nlm.steptol, hessian = FALSE, iterlim=model$control$nlm.maxit))
      fit <- list(par=fit$estimate, value=fit$minimum)
    }
  }, silent = FALSE)
  if (class(z) == "try-error") stop(call.=NULL, hopit_msg(6))

  model$coef <- fit$par
  model$LL <- unname(-fit$value)

  if (use_weights) {
    model$LL <- model$LL / sum(model$weights) * model$N #scale likelihood
  }
  model
}


#' Auxiliary for controlling the fitting of \code{hopit} model
#'
#' @description
#' Auxiliary function for controlling the fitting of \code{hopit} model. Use this function to set control
#' parameters of the \code{\link{hopit}} and other related functions.
#' @param grad.eps epsilon for numerical Hessian function.
#' @param quick.fit logical, if TRUE extensive \code{nlm} optimization method is ignored and only BFGS and CG methods are run.
#' @param bgfs.maxit,cg.maxit,nlm.maxit the maximum number of iterations. See \code{\link{optim}} and \code{\link{nlm}} for details.
#' @param bgfs.reltol,cg.reltol relative convergence tolerance. See \code{\link{optim}} for details.
#' @param nlm.gradtol,nlm.steptol tolerance at which the scaled gradient is considered close enough to zero and
#' minimum allowable relative step length. See \code{\link{nlm}}.
#' @param fit.methods either 'CG' or 'BFGS'. See \code{\link{optim}}.
#' @param trace logical, if to trace model fitting.
#' @param LL_out_val,thresh.1.exp,thresh.start internal parameters under development, do not change.
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
                        fit.methods = c('CG','BFGS'),
                        quick.fit = TRUE,
                        trace = TRUE,
                        thresh.start = -Inf,
                        thresh.1.exp = FALSE,
                        LL_out_val = -Inf){

  if (!length(fit.methods)) stop(hopit_msg(8),call. = NULL) else fit.methods <- toupper(fit.methods)
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
       quick.fit = quick.fit,
       trace = trace,
       thresh.start = thresh.start,
       thresh.1.exp = thresh.1.exp,
       LL_out_val = LL_out_val)
}


#' Extract Theta parameter from the \code{hopit} model
#'
#' @param model a fitted \code{hopit} model.
#' @export
#' @author Maciej J. Danko
getTheta <- function(model) unname(exp(model$coef.ls$logTheta))


#' Generalized hierarchical ordered threshold models.
#'
#' @description
#' The ordered response data classifies a measure of interest into ordered categories collected during a survey.
#' If the dependent variable is a happiness then a respondent typically answers a question:
#' “Taking all things together, would you say you are ...? “ and have some response options e.g.
#' "very happy", "pretty happy", "not too happy", "very unhappy" \insertCite{Liao2005}{hopit}.
#' Similarly if interviewees are asked to evaluate their health in general (e.g. “Would you say your health is ...?”)
#' they may choose among several categories, e.g. very good, good, fair, bad, and very bad
#' \insertCite{King2004,Jurges2007,Rebelo2014}{hopit}. In political sciences a respondent may be
#' asked for an opinion about recent legislation (e.g. “Rate your feelings about the proposed legislation")
#' and asked to choose among several categories "strongly oppose", "mildly oppose", "indifferent",
#' "mildly support", "strongly support" \insertCite{GreeneHensher2010}{hopit}.
#' It is easy to imagine other multi-level ordinal variables that might by use during the survey and
#' to which methodology described below could be applied with.\cr
#'
#' Practically, it is assumed that when responding to a survey question about their general happiness,
#' health, feeling, attitude or other status, participants assess their true value of this
#' unobserved continuous variable, and project it to a provided discrete scale.
#' The thresholds that each individual uses to categorize their true status into a specific response
#' option may be affected by the choice of a reference group, earlier life experiences, and cross-cultural
#' differences in using scales, and thus, may differ across individuals depending on their gender,
#' age, cultural background, education, and personality traits, among other factors.\cr
#'
#' From the reporting-styles modeling perspective, one of the main tasks is to compute this continuous estimate of individuals’
#' underlying, latent measure based on several specific characteristics of the considered response
#' (e.g. health variables or happiness variables) and accounting also for variations in reporting across
#' socio-demographic and cultural groups. More specifically, to build the latent, underlying measure
#' a generalized hierarchical ordered threshold model is fitted, which regresses the reported
#' status/attitude/feeling on two sets of independent variables \insertCite{Boes2006,Green2014}{hopit}.
#' When a dependent reported ordered variable is self-rated health status then the first set of variables
#' (health variables) assesses individuals’ specific aspects of health, and might include chronic conditions,
#' mobility level, difficulties with a range of daily activities, performance on grip strength test,
#' anthropometric measures, lifestyle behaviors, etc. Using the second set of independent variables
#' (threshold variables), the model also adjusts for the differences across socio-demographic and cultural
#' groups like cultural background, gender, age, education, etc. \insertCite{King2004,Jurges2007}{hopit}.\cr
#'
#' Ordered threshold models are used to fit ordered categorical dependent variables.
#' The generalized ordered threshold models \insertCite{Terza1985,Boes2006,Green2014}{hopit}
#' are an extension to the ordered threshold models \insertCite{McKelvey1975}{hopit}.
#' In the latter models the thresholds are constant, whereas generalized models allows
#' thresholds to be dependent on covariates.
#' \insertCite{GreeneHensher2010,Green2014;textual}{hopit} pointed out that also
#' thresholds must be ordered so that a model has a sense.
#' This motivated Greene and coauthors to call this models *HOPIT*, which stands for hierarchical ordered probit models.
#'
#' The fitted *hopit* model is used to analyse reporting styles. See \code{\link{standardizeCoef}}, \code{\link{latentIndex}},
#' \code{\link{getCutPoints}}, and \code{\link{getLevels}}.
#' @details
#' The function fits generelaized hierarchical ordered threshold models.\cr
#'
#' \code{latent.formula} models latent variable.
#' if the response variable is self-rated health then latent measure can depend on different health
#' conditions and diseases (latent variables are called health variables).
#' Latent variables are modeled with parallel regression assumption. According to the assumption, coefficients,
#' which describe the relationship between lowest and all higher response categories, are the same as those coefficients,
#' which describe the relationship between another (e.g. adjacent) lowest and the remaining higher response categories.
#' The predicted latent variable is modeled as a linear function of health variables and corresponding coefficients.\cr
#'
#' \code{thresh.formula} models threshold variable.
#' The thresholds (cut points, \code{alpha}) are modeled by threshold variables \code{gamma} and intercepts \code{lambda}.
#' It is assumed that they model contextual characteristics of the respondent (e.g. country, gender, age, etc. ).
#' Threshold variables are modeled without parallel regression assumption, thus each threshold is modeled by
#' a variable independently \insertCite{Boes2006,Green2014}{hopit}.
#' \code{hopit}() function uses parametrization of tresholds proposed by \insertCite{Jurges2007;textual}{hopit}.\cr
#'
#' \code{decreasing.levels} it is the logical that determines the ordering of levels of the categorical response variable.
#' It is always good to check first the ordering of the levels before starting (see example 1)\cr
#'
#' @param latent.formula formula used to model latent variable. It should not contain any threshold variable.
#' To specify interactions between latent and threshold variables see \code{crossinter.formula}
#' @param thresh.formula formula used to model threshold variable. It should not contain any latent variable.
#' To specify interactions between latent and threshold variables see \code{crossinter.formula}
#' Any dependent variable (left side of "~") will be ignored.
#' @param crossinter.formula formula used to model interactions between threshold and latent variables.
#' The formula should be in the form \cr\code{"~ thresh.var.X : latent.var.Y + thresh.var.Z : latent.var.Q + ..." }.\cr
#' All \code{"latent.vars"} must be present in threshold formula and all \code{"thresh.vars"} must be present in threshold formula.
#' Any dependent variable (left side of "~") will be ignored.
#' @param data a data frame including all modeled variables.
#' @param decreasing.levels logical indicating if self-reported health classes are ordered in decreasing order.
#' @param overdispersion logical indicting if to fit additional parameter theta, which models a variance of the error term.
#' @param design an optional survey design. Use \code{\link[survey]{svydesign}} function to specify the design.
#' The design cannot be specified together with parameter \code{weights}.
#' @param weights optional model weights. Use design to construct survey weights.
#' @param link the link function. The possible values are \code{"probit"} (default) and \code{"logit"}.
#' @param start a vector with starting coefficient values in the form \code{c(latent_parameters, threshold_lambdas, threshold_gammas)}.
#' @param control a list with control parameters. See \code{\link{hopit.control}}.
#' @return a \code{hopit} object used by other functions and methods. The object is a list with following components:
#'  \item{control}{ a list with control parameters. See \code{\link{hopit.control}}.}
#'  \item{link}{ the used link funtion.}
#'  \item{hasdisp}{ logical, was overdispersion modeled?}
#'  \item{latent.formula}{ used latent formula.}
#'  \item{latent.mm}{ latent model matrix.}
#'  \item{latent.terms}{ used latent variables.}
#'  \item{thresh.formula}{ used threshold formula.}
#'  \item{thresh.mm}{ threshold model matrix.}
#'  \item{thresh.extd}{ threshold extended model matrix.}
#'  \item{thresh.terms}{ used threshold variables.}
#'  \item{thresh.no.cov}{ logical, are gamma parameters present?}
#'  \item{parcount}{ 3-element vector with number of parmeters for latent latent variable (beta),
#'  threshold intercept (lambda), and threshold covariates (gamma).}
#'  \item{coef}{ a vector with coefficients.}
#'  \item{coef.ls}{ coefficients as a list.}
#'  \item{start}{ vector with starting vlues of coefficients.}
#'  \item{alpha}{ estimated individual-specific thresholds.}
#'  \item{y_i}{ the response variable.}
#'  \item{y_latent_i}{ predicted latent measure.}
#'  \item{Ey_i}{ predicted categorical response.}
#'  \item{J}{ the number of response levels.}
#'  \item{N}{ the number of observations.}
#'  \item{deviance}{ deviance.}
#'  \item{LL}{ log likelihood.}
#'  \item{AIC}{ AIC for models without survey design.}
#'  \item{vcov}{ variance-covariance matrix.}
#'  \item{hessian}{ a Hessian matrix.}
#'  \item{estfun}{ gradient of the log likelihood function at estimated coefficient values.}
#'  \item{YYY1,YYY2,YY3}{ internal objects used for calculation of gradient and Hessian functions.}
#'  \item{use.weights,vcov.basic,glm.start,glm.start.ls}{ other internal objects.}
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
#' For reporting styles analysis see:\cr
#' \code{\link{standardizeCoef}},
#' \code{\link{latentIndex}},
#' \code{\link{getCutPoints}},
#' \code{\link{getLevels}}.
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # first determine the order of the dependent variable
#' levels(healthsurvey$health)
#'
#' # Example 1 ---------------------
#'
#' # the order is decreasing (from the best health to the worst health)
#' # so we set: decreasing.levels = TRUE
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
#' # extract parameters in a form of list
#' cm1 <- coef(model1, aslist = TRUE)
#'
#' # names of returned coefficients
#' names(cm1)
#'
#' # extracting latent health coefficients
#' cm1$latent.params
#'
#' # Checking the fit
#' profile(model1)
#'
#' # Example 2 ---------------------
#'
#' # incorporating survey design
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
#' # compare latent variables
#' cbind('No survey design' = coef(model1, aslist = TRUE)$latent.par,
#' 'Has survey design' = coef(model2, aslist = TRUE)$latent.par)
#'
#' # Example 3 ---------------------
#'
#' # using crossinter.formula
#' model3 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               crossinter.formula = ~ sex : depression + sex : diabetes,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # print the model:
#' coef(model3, aslist = TRUE)$latent.params
#'
#' # compare fit of model1 and model 3
#' # Likelihood Ratio Test
#' anova(model1, model3)
#'
#' # AIC
#' AIC(model1, model3)
#'
#' # Example 4 ---------------------
#'
#' # model with overdispersion
#' model4 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               overdispersion = TRUE,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # estimated variance of the error term:
#' getTheta(model4)
#'
#' # compare fit of model1 and model4
#' # Likelihood Ratio Test
#' print(anova(model1, model4), short = TRUE)
hopit<- function(latent.formula,
                 thresh.formula = ~ 1,
                 crossinter.formula = ~ 1,
                 data,
                 decreasing.levels,
                 start = NULL,
                 overdispersion = FALSE,
                 design = list(),
                 weights = NULL,
                 link = c('probit', 'logit'),
                 control = list()){

  if (!overdispersion) remove.theta = FALSE else remove.theta = TRUE
  if (missing(data)) data <- environment(latent.formula)
  link <- match.arg(link)
  control <- do.call("hopit.control", control)

  model <- NULL
  model$control <- control
  model$link <- link[1]
  model$hasdisp <- overdispersion

  thresh.formula <- check_thresh_formula(thresh.formula)
  latent.formula <- check_latent_formula(latent.formula)

  thresh.terms <- attr(terms(thresh.formula),'term.labels')
  latent.terms <- attr(terms(latent.formula),'term.labels')
  thresh.list <- unlist(decomposeformula(thresh.terms))
  latent.list <- unlist(decomposeformula(latent.terms))
  if(any(thresh.list %in% latent.list) || any(latent.list %in% thresh.list)) stop(hopit_msg(89), call.=NULL)

  crossinter.formula <- check_crossinter_formula(crossinter.formula)
  crossinter <- attr(stats::terms(crossinter.formula),'term.labels')
  if (length(crossinter)) {
    tmp <- length(crossinter)
    crossinter <- crossinter[which(grepl(":", unlist(crossinter), fixed=TRUE))]
    if (length(crossinter) != tmp) warning(hopit_msg(90), call. = NA)
  }
  if (length(crossinter)) {
    # check if one component is in the thresh and second in the latent
    crossinter.list <- decomposeformula(crossinter)
    crossinter.test <- sapply(crossinter.list, function(k) any(k %in% thresh.terms) & any(k %in% latent.terms))
    if (!all(crossinter.test)) {
      ll <- paste(crossinter[which(!crossinter.test)],collapse = ' & ')
      stop(paste(hopit_msg(87),ll, hopit_msg(88)), call.=NULL)
    }
    crossinter.thresh <-unique(sapply(crossinter.list, function(k) k[k %in% thresh.terms]))
    crossinter.f <- paste(crossinter, collapse=' + ')
    latent.formulaA <-update(latent.formula, paste('.~. + ',crossinter.f,' + ',crossinter.thresh))
    latent.formulaB <-update(latent.formulaA, paste('.~. - ',crossinter.thresh))
    m1<-stats::model.matrix(latent.formulaA, data)
    m2<-stats::model.matrix(latent.formulaB, data)
    cm1<-colnames(m1)
    rco<-cm1[!(cm1 %in% colnames(m2))]
    model$latent.mm <- m1[,which(!(cm1 %in% rco))]
    model$latent.formula <- latent.formulaB
    model$crossinter.formula <- crossinter.formula
    model$crossinter <- crossinter
  } else {
    model$latent.formula <- latent.formula
    model$latent.mm <- as.matrix(stats::model.matrix(latent.formula, data = data))
  }

  latent.names <- colnames(model$latent.mm)
  grpi <- findintercept(latent.names)
  model$latent.mm <- as.matrix(model$latent.mm[,!grpi])
  latent.names <- latent.names[!grpi]

  model$thresh.formula <- thresh.formula
  model$thresh.mm <- as.matrix(stats::model.matrix(thresh.formula, data = data))
  thresh.names <- colnames(model$thresh.mm)
  grpi <- findintercept(thresh.names)
  model$thresh.mm <- as.matrix(model$thresh.mm[,!grpi])
  thresh.names <- thresh.names[!grpi]
  if (any(dim(model$thresh.mm) == 0L)) {
    model$thresh.no.cov <- TRUE
  } else {
    model$thresh.no.cov <- FALSE
  }

  model$latent.terms <- attr(stats::terms(stats::as.formula(latent.formula)),"term.labels")
  model$thresh.terms <-attr(stats::terms(stats::as.formula(thresh.formula)),"term.labels")

  model$y_i <- stats::model.frame(latent.formula, data = data)[,all.vars(latent.formula[[2]])]
  check_response(model$y_i)
  if (missing(decreasing.levels)) decreasing.levels = NULL
  check_decreasing.levels(decreasing.levels, levels(model$y_i))
  if (!decreasing.levels) model$y_i <- factor(model$y_i, rev(levels(model$y_i)))
  model$y_latent_i <- NA # latent
  model$Ey_i <- NA # ordinal classified utput
  model$J <- length(levels(model$y_i))
  model$N <- length(model$y_i)

  model$thresh.extd <- matrix(rep_row(model$thresh.mm, model$J-1),model$N, NCOL(model$thresh.mm)*(model$J-1))

  Cr <- dim(model$latent.mm)[2L]
  Ct <- dim(model$thresh.mm)[2L]
  if (model$thresh.no.cov) Ct <- 0L
  model$parcount <- c(Cr, model$J - 1L, Ct*(model$J - 1L))
  model$parcount[3L] <- model$parcount[3L]*(model$thresh.no.cov == FALSE)
  interce <- paste(1L : (model$J - 1L), 2L : (model$J), sep = '|')
  if (model$thresh.no.cov){
    tmp <- NULL
  } else {
    tmp <- as.matrix(expand.grid(interce, thresh.names, KEEP.OUT.ATTRS = FALSE))
    tmp <- tmp[,c(2,1)]
    tmp <- paste('(G)', apply(tmp, 1L, paste, sep = '', collapse = '.'), sep = '.')
  }

  coefnames <-  c(latent.names, paste('(L)', interce, sep = '.'), tmp)
  if (model$hasdisp) coefnames <- c(coefnames, 'logTheta')

  model$weights <- NULL
  check_design(weights, design, model$N)
  if (length(design)) model$weights <- design$prob else if (length(weights)) model$weights <- weights
  if (!length(model$weights)) {
    model$weights <- rep(1, model$N)
    model$use.weights <- FALSE
  } else model$use.weights <- TRUE

  model$weights <- as.vector(matrix(model$weights, 1L, model$N))
  #scaling weights
  model$weights <- model$N * model$weights / sum(model$weights)

  #calculate special matrices for gradient calaculation
  model <- calcYYY(model)

  if (!length(start)) {
    if (model$control$trace) cat(hopit_msg(10))
    model <- suppressWarnings(get.hopit.start(model, data)) #rename this function get.hopit.start
    if (model$control$trace) cat(hopit_msg(13))
  } else {
    model$start <- start
  }
  if (model$control$trace) cat(hopit_msg(11))
  model <- hopit_fitter(model, start = model$start)
  class(model) <- 'hopit'
  colnames(model$thresh.mm) <- thresh.names
  names(model$coef) <- coefnames
  if (model$control$trace) cat(hopit_msg(13),hopit_msg(12),sep='')

  p <- hopit_ExtractParameters(model)
  model$alpha <- hopit_Threshold(p$thresh.lambda, p$thresh.gamma, model)
  model$Ey_i <- classify.ind(model)
  model$y_latent_i <- hopit_Latent(model$coef[seq_len(model$parcount[1])], model)
  model$coef.ls <- p
  model$deviance <- -2 * model$LL
  if (model$control$trace) cat(hopit_msg(13),hopit_msg(14),sep='')

  hes <- my.grad(fn = hopit_derivLL, par = model$coef, model=model, eps = 1e-4, collapse = TRUE, negative=FALSE)
  if (model$hasdisp && remove.theta) {
    hes <- hes[-nrow(hes),-ncol(hes)] #remove theta from vcov
    model$coef <- model$coef[-length(model$coef)] #remove from coef
  }
  model$hessian <- hes

  model$vcov.basic <- try(base::solve(-hes), silent = FALSE)
  model$vcov.basic <- check_vcov(model$vcov.basic)
  if (model$control$trace) cat(hopit_msg(13),hopit_msg(15),sep='')
  if (model$hasdisp && remove.theta) COEF <- c(model$coef,model$coef.ls$logTheta) else COEF <- model$coef
  model$estfun <- hopit_derivLL(COEF, model, collapse = FALSE)
  if (remove.theta) model$estfun <- model$estfun[,-ncol(model$estfun)]
  if (model$control$trace) cat(hopit_msg(13))

  if (length(model$design)) {
    if (model$control$trace) cat(hopit_msg(16))
    model$vcov <- svy.varcoef.hopit(model$vcov.basic, model$estfun, design)
    model$AIC <- NA
    if (model$control$trace) cat(hopit_msg(13))
  } else {
    k <- 2
    model$vcov <- model$vcov.basic
    model$AIC <- model$deviance + k * (length(model$coef.ls$latent.params)+
                                       length(model$coef.ls$thresh.lambda)+
                                       length(model$coef.ls$thresh.gamma)+model$hasdisp)
  }
  return(model)
}
