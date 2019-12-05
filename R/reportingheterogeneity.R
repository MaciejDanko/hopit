#' Calculate the latent index
#' @description
#' Calculate the latent index from the fitted model. The latent index is a standardized latent measure that takes values from 0 to 1, where
#' 0 refers to the worst predicted state (the maximal observed value for the latent measure) and 1 refers
#' to the best predicted state (the minimal observed value for the latent measure).
#' @param model a fitted \code{hopit} model.
#' @param subset an optional vector that specifies a subset of observations.
#' @return a vector with a latent index for each individual.
#' @references
#'  \insertRef{Jurges2007}{hopit}\cr\cr
#'  \insertRef{OKSUZYAN2019}{hopit}
#' @author Maciej J. Danko
#' @export
#' @seealso \code{\link{standardizeCoef}}, \code{\link{getCutPoints}}, \code{\link{getLevels}}, \code{\link{hopit}}.
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # the order of response levels decreases from the best health to
#' # the worst health; hence the hopit() parameter decreasing.levels
#' # is set to TRUE
#' levels(healthsurvey$health)
#'
#' # Example 1 ---------------------
#'
#' # fit a model
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                   heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                   depression + respiratory_problems +
#'                   IADL_problems + obese + diabetes + other_diseases,
#'                 thresh.formula = ~ sex + ageclass + country,
#'                 decreasing.levels = TRUE,
#'                 control = list(trace = FALSE),
#'                 data = healthsurvey)
#'
#' # calculate the health index
#' hi <- latentIndex(model1)
#'
#' summary(hi)
#'
#' # plot a simple histogram of the function output
#' hist(hi, col='deepskyblue3')
#'
#' #plot the reported health status versus the health index.
#' plot(hi, response = "data", ylab = 'Health index',
#'      col='deepskyblue3', main = 'Reported health levels')
#'
#' # plot the model-predicted health levels versus the health index.
#' plot(hi, response = "fitted", ylab = 'Health index',
#'      col='deepskyblue3', main = 'Model-predicted health levels')
latentIndex <- function(model, subset = NULL) {
  if (length(subset) == 0) subset=seq_len(model$N)
  r <- range(model$y_latent_i[subset])
  hi <- (1 - ((model$y_latent_i - r[1]) / diff(r)))[subset]
  mi <- data.frame(y_i=model$y_i[subset], Ey_i=model$Ey_i[subset])
  attr(hi,'model.info')<-mi
  class(hi) <- c('hopitHI','vector')
  hi
}


#' @rdname latentIndex
#' @export
healthIndex <- latentIndex


#' Standardization of the coefficients
#' @description
#' Calculate standardized the coefficients (e.g. disability weights for the health variables) using
#' the predicted latent measure obtained from the model.\cr
#' In the self-rated health example the standardized coefficients are called disability weights \insertCite{Jurges2007;textual}{hopit}
#' and are calculated for each health variable to provide information about the impact of a specific health measure on the latent index
#' (see \code{\link{latentIndex}}). The disability weight for a health variable is equal to the ratio of the corresponding health coefficient
#'  and the difference between the lowest and the highest values of the predicted latent health. In other words, the disability weight reduces
#'  the latent index by some given amount or percentage (i.e., the latent index of every individual is reduced by the same amount if the person had a heart attack or other
#'  heart problems)\insertCite{Jurges2007}{hopit}.
#' @param model a fitted \code{hopit} model.
#' @param namesf a vector of the names of coefficients or one argument function that modifies the names of coefficients.
#' @name standardizeCoef
#' @importFrom stats symnum
#' @return a vector with standardized coefficients.
#' @references
#'  \insertRef{Jurges2007}{hopit}\cr\cr
#'  \insertRef{OKSUZYAN2019}{hopit}
#' @author Maciej J. Danko
#' @export
#' @seealso \code{\link{latentIndex}}, \code{\link{getCutPoints}}, \code{\link{getLevels}}, \code{\link{hopit}}.
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # the order of response levels decreases from the best health to
#' # the worst health; hence the hopit() parameter decreasing.levels
#' # is set to TRUE
#' levels(healthsurvey$health)
#'
#' # Example 1 ---------------------
#'
#' # fit a model
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # a function that modifies the coefficient names.
#' txtfun <- function(x) gsub('_',' ',substr(x,1,nchar(x)-3))
#'
#' # calculate and plot the disability weights
#' sc <- standardizeCoef(model1, namesf = txtfun)
#' sc
#'
#' summary(sc)
#'
#' plot(sc)
standardizeCoef <- function (model,
                             namesf = identity) {

  cfm <- suppressMessages(summary(model)$coef[seq_len(model$parcount[1]),])
  oldnames<-rownames(cfm)
  if (class(namesf)[1]=='function') newnames <- namesf(oldnames) else {
    newnames <- namesf
  }

  r <- range(model$y_latent_i)
  res <- as.matrix((cfm[,1])/diff(r))
  sigstar <- stats::symnum(cfm$`Pr(>|z|)`, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  res2<- data.frame('old.names'=oldnames, 'Pr(>|z|)'=cfm$`Pr(>|z|)`,' '=sigstar,check.names = FALSE, fix.empty.names = FALSE,
                   stringsAsFactors = FALSE)
  rownames(res) <- rownames(res2) <- newnames
  colnames(res) <- 'Std. coef'
  class(res)<-c('hopitDW')
  attr(res, 'Ptab')<- res2
  attr(res, 'legend')<-attr(sigstar, 'legend')
  return(res)
}


#' @rdname standardizeCoef
#' @export
standardiseCoef<-standardizeCoef


#' @rdname standardizeCoef
#' @export
disabilityWeights<-standardizeCoef


#' Calculate the threshold cut-points and individual adjusted responses using Jurges' method
#' @description
#' Calculate the threshold cut-points and individual adjusted responses using Jurges' method
#' @param model a fitted \code{hopit} model.
#' @param decreasing.levels a logical indicating whether self-reported health classes are ordered in increasing order.
#' @param subset an optional vector specifying a subset of observations.
#' @return a list with the following components:
#'  \item{cutpoints}{ cut-points for the adjusted categorical response levels with the corresponding percentiles of the latent index.}
#'  \item{adjusted.levels}{ adjusted categorical response levels for each individual.}
#' @references
#'  \insertRef{Jurges2007}{hopit}\cr\cr
#'  \insertRef{OKSUZYAN2019}{hopit}
#' @author Maciej J. Danko
#' @export
#' @seealso \code{\link{latentIndex}}, \code{\link{standardiseCoef}}, \code{\link{getLevels}}, \code{\link{hopit}}.
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # the order of response levels decreases from the best health to
#' # the worst health; hence the hopit() parameter decreasing.levels
#' # is set to TRUE
#' levels(healthsurvey$health)
#'
#' # Example 1 ---------------------
#'
#' # fit a model
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # calculate the health index cut-points
#' z <- getCutPoints(model = model1)
#' z$cutpoints
#'
#' plot(z)
#'
#' # tabulate the adjusted health levels for individuals (Jurges method):
#' rev(table(z$adjusted.levels))
#'
#' # tabulate the original health levels for individuals
#' table(model1$y_i)
#'
#' # tabulate the predicted health levels
#' table(model1$Ey_i)
getCutPoints <- function(model,
                         decreasing.levels=model$decreasing.levels,
                         subset=NULL){

  if (length(subset) == 0) subset=seq_along(model$y_i)
  Y <- model$y_i[subset]
  if (!length(decreasing.levels)) decreasing.levels=model$decreasing.levels
  if (decreasing.levels) dorev <- rev.default else dorev <- identity
  h.index <- healthIndex(model, subset)
  tY <- table(Y)
  tmp <- dorev(as.vector(tY))
  invcs <- (cumsum(tmp)/sum(tmp))[-length(tmp)]
  R1 <- stats::quantile(h.index, invcs)

  if (length(h.index)){
    CIN <- c(0,R1,1)
    if (anyNA(CIN)){
      ind=seq_along(CIN)[is.na(CIN)]
      for (k in ind) CIN[ind]=CIN[ind-1]
    }
    if (anyDuplicated(CIN)) {
      warning(hopit_msg(85), call.=NA)
      CIN <- CIN + cumsum(duplicated(CIN)*CIN/1e7)
      CIN <- CIN / max(CIN)
    }
    adjusted.levels<- cut(h.index, CIN,labels= dorev(levels(Y)),
                         include.lowest = TRUE)
  } else adjusted.levels <- NA
  res <- list(decreasing.levels=decreasing.levels, cutpoints=R1,
              adjusted.levels=adjusted.levels, h.index=h.index)
  attr(res, 'y_i') <- model$y_i
  class(res) <- c('hopitCP','list')
  res
}


#' Summarize the adjusted and the original self-rated response levels
#' @description
#' Summarize the adjusted and the original self-rated response levels.
#' @param model a fitted \code{hopit} model.
#' @param formula a formula containing the grouping variables. It is by default set to threshold formula.
#' @param data data used to fit the model.
#' @param sep a separator for the level names.
#' @param decreasing.levels a logical indicating whether self-reported health classes are ordered in increasing order.
#' @param sort.flag a logical indicating whether to sort the levels.
#' @return a list with the following components:
#'  \item{original}{ frequencies of original response levels for selected groups/categories.}
#'  \item{adjusted}{ frequencies of adjusted response levels (Jurges 2007 method) for selected groups/categories.}
#'  \item{N.original}{ the number of original response levels for selected groups/categories.}
#'  \item{N.adjusted}{ the number of adjusted response levels (Jurges 2007 method) for selected groups/categories.}
#'  \item{categories}{ selected groups/categories used in summary.}
#'  \item{tab}{ an original vs. an adjusted contingency table.}
#'  \item{mat}{ a matrix with columns: grouping variables, original response levels, adjusted response levels.
#'  Each row corresponds to a single individual from the data used to fit the model.}
#' @references
#'  \insertRef{Jurges2007}{hopit}\cr\cr
#'  \insertRef{OKSUZYAN2019}{hopit}
#' @author Maciej J. Danko
#' @export
#' @seealso \code{\link{getCutPoints}}, \code{\link{latentIndex}}, \code{\link{standardiseCoef}}, \code{\link{hopit}}.
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # the order of response levels decreases from the best health to
#' # the worst health; hence the hopit() parameter decreasing.levels
#' # is set to TRUE
#' levels(healthsurvey$health)
#'
#' # fit a model
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # Example 1 ---------------------
#'
#' # calculate a summary by country
#' hl <- getLevels(model=model1, formula=~ country, sep=' ')
#' plot(hl, las=1, mar = c(3,2,1.5,0.5))
#'
#' # differences between frequencies of original and adjusted health levels
#' round(100*(hl$original - hl$adjusted),2)
#'
#' # extract good and bad health levels (combined levels)
#' Org <- cbind(bad = rowSums(hl$original[,1:2]),
#'              good = rowSums(hl$original[,4:5]))
#' Adj <- cbind(bad = rowSums(hl$adjusted[,1:2]),
#'              good = rowSums(hl$adjusted[,4:5]))
#' round(100*(Org - Adj),2)
#'
#' # plot the differences
#' barplot(t(Org - Adj), beside = TRUE, density = 20, angle = c(-45, 45),
#'         col = c('pink4', 'green2'),
#'         ylab = 'Original - adjusted reported health frequencies')
#' abline(h = 0); box()
#' legend('top', c('Bad health','Good health'),
#'        density = 20, angle = c(-45, 45),
#'        fill = c('pink4', 'green2'), bty = 'n', cex = 1.2)
#'
#' # in country X, bad health seems to be over-reported while good health
#' # is under-reported; in country Z, good health is highly over-reported.
#'
#' # Example 2 ---------------------
#'
#' # summary by gender and age
#' hl <- getLevels(model = model1, formula=~ sex + ageclass, sep=' ')
#' plot(hl)
#'
#' # differences between frequencies of original and adjusted health levels
#' round(100*(hl$original - hl$adjusted),2)
#'
#' # extract good health levels (combined "Very good" and "Excellent" levels)
#' Org <- rowSums(hl$original[,4:5])
#' Adj <- rowSums(hl$adjusted[,4:5])
#' round(100*(Org - Adj),2)
#'
#' pmar <- par('mar'); par(mar = c(9.5, pmar[2:4]))
#' barplot(Org-Adj,
#'         ylab = 'Original - adjusted reported good health frequencies',
#'         las = 3,
#'         density = 20, angle = c(45, -45), col = c('blue', 'orange'))
#' abline(h = 0); box(); par(mar = pmar)
#' legend('top', c('Man','Woman'), density = 20, angle = c(-45, 45),
#'        fill = c('blue', 'orange'), bty = 'n', cex = 1.2)
#'
#' # results show that women in general tend to over-report good health,
#' # while men aged 50-59 greatly under-report good health.
#'
#' # more examples can be found in the description of the boot_hopit() function.
getLevels<-function(model,
                    formula=model$thresh.formula,
                    data = model$frame,
                    sep = '_',
                    decreasing.levels = model$decreasing.levels,
                    sort.flag = FALSE){

  data <- model$na.action(data)
  sort.flag <- sort.flag[1]
  if (model$control$transform.latent != 'none')
    data <- transform.data(model$latent.formula, data,
                           model$control$transform.latent)
  if (model$control$transform.thresh != 'none')
    data <- transform.data(model$thresh.formula, data,
                           model$control$transform.thresh)

  if (class(formula)[1]=='formula')
    inte_ <- formula2classes(formula, data, sep = sep,
                             return.matrix = TRUE) else
      stop(call.=NULL, hopit_msg(86))
  inte <- inte_$x
  namind <- inte_$class.mat
  nam <- levels(inte)

  cpall <- getCutPoints(model, decreasing.levels = decreasing.levels)
  TAB1 <- round(table(original = model$y_i,
                      adjusted = cpall$adjusted.levels)*100/length(model$y_i),2)
  tmp <- untable(t(table(factor(model$y_i,
                            levels=levels(cpall$adjusted.levels)), inte)))
  N1 <- tmp
  tmp <-tmp/rowSums(tmp)
  if (sort.flag) {
    oD1 <- order(tmp[,NCOL(tmp)]+tmp[,NCOL(tmp)-1])
    tmp <- tmp[oD1,]
    orignalind <- namind[oD1,]
  } else orignalind <- namind

  tmp2 <- untable(t(table(cpall$adjusted.levels, inte)))
  N2 <- tmp2
  tmp2 <- tmp2/rowSums(tmp2)
  if (sort.flag) {
    oD2 <- order(tmp2[,NCOL(tmp2)]+tmp2[,NCOL(tmp2)-1])
    tmp2 <- tmp2[oD2,]
    adjustedind <- namind[oD2,]
  } else adjustedind <- namind

  res <- list(original= tmp,
              adjusted= tmp2,
              tab= TAB1,
              N.original= N1,
              N.adjusted= N2,
              categories = orignalind, # = adjustedind
              mat=as.data.frame(cbind(group=inte_$mat,
                        original= model$y_i,
                        adjusted= cpall$adjusted.levels)))
  class(res) <- c('hopitLV','list')
  res
}
