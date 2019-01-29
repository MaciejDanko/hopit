#' Calculate latent index
#' @description
#' Calculate latent index from the latent variable. It takes values from 0 to 1, where
#' zero is prescribed to the worse predicted state (maximal observed value for the latent variable) and 1 is prescribed
#' to the best predicted health (minimal observed value fo the latent variable).
#' @param model a fitted \code{hopit} model.
#' @param decreasing.levels logical indicating if self-reported (e.g. health) classes are ordered in decreasing order.
#' @param subset an optional vector specifying a subset of observations.
#' @param plotf logical indicating if to plot summary figure.
#' @param response X axis plotting option, choose \code{'data'} for raw responses and \code{'fitted'} for model reclassified responses
#' @param ylab a label of y axis.
#' @param ... further parameters passed to the \code{\link{plot}} function.
#' @return a vector with latent index for each individual.
#' @references \insertRef{Jurges2007}{hopit}
#' @author Maciej J. Danko
#' @export
#' @seealso \code{\link{standardizeCoef}}, \code{\link{getCutPoints}}, \code{\link{getLevels}}, \code{\link{hopit}}.
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # the order of response levels is decreasing (from the best health to the worst health)
#' levels(healthsurvey$health)
#'
#' # Example 1 ---------------------
#'
#' # fitting a model
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # calculate health index and plotting reported health status
#' # vs. health index.
#' hi <- latentIndex(model1, plotf = TRUE, response = "data",
#'                   ylab = 'Health index', col='deepskyblue3')
#'
#' # a simple histogram of the function output
#' hist(hi)
#'
#' # calculate health index and plotting adjusted health status (Jurges 2007)
#' # vs. health index.
#' latentIndex(model1, plotf = TRUE, response = "Jurges",
#'                  ylab = 'Health index', col='deepskyblue3')
#'
#' # calculate health index and plotting predicted health status
#' # vs. health index.
#' latentIndex(model1, plotf = TRUE, response = "fitted",
#'                  ylab = 'Health index', col='deepskyblue3')
latentIndex <- function(model, decreasing.levels = TRUE,
                        subset = NULL, plotf = FALSE,
                        response = c('data','fitted','Jurges'),
                        ylab = 'Latent index', ...) {
  if (length(subset) == 0) subset=seq_len(model$N)
  r <- range(model$y_latent_i[subset])
  hi <- (1 - ((model$y_latent_i - r[1]) / diff(r)))[subset]
  if (plotf) {
    response <- tolower(match.arg(response))
    if (response=='data') YY <- model$y_i else if (response=='fitted') YY <- model$Ey_i else if (response=='jurges') {
      z <- getCutPoints(model=model, decreasing.levels = decreasing.levels, plotf = FALSE)
      YY <- factor(z$adjused.levels,levels(model$y_i))
    } else stop(hopit_msg(83),call.=NULL)
    graphics::plot(YY[subset], hi,las=3, ylab=ylab, ...)
  }
  if (plotf) invisible(hi) else return(hi)
}


#' @rdname latentIndex
#' @export
healthIndex <- latentIndex


#' Standardization of coefficients
#' @description
#' Calculate standardized coefficients - disability weights
#' computed as the latent coefficients from the generalised ordered probit model divided by the maximum possible range
#' of its linear prediction. The range is calculated as difference between maximum and minimum possible value of the latent variable
#' given estimated parameters.
#' @param model a fitted \code{hopit} model.
#' @param ordered logical indicating if to order the disability weights.
#' @param plotf logical indicating if to plot results.
#' @param plotpval logical indicating if to plot p-values.
#' @param mar,oma see \code{\link{par}}.
#' @param YLab,YLab.cex label and size of the label for y axis.
#' @param namesf  a vector of names of coefficients or one argument function that modifies names of coefficients.
#' @param ... arguments passed to \code{\link{boxplot}}.
#' @name standardizeCoef
#' @return a vector with standardized coefficients.
#' @references \insertRef{Jurges2007}{hopit}
#' @author Maciej J. Danko
#' @export
#' @seealso \code{\link{latentIndex}}, \code{\link{getCutPoints}}, \code{\link{getLevels}}, \code{\link{hopit}}.
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # the order of response levels is decreasing (from the best health to the worst health)
#' levels(healthsurvey$health)
#'
#' # Example 1 ---------------------
#'
#' # fitting a model
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # A function that modifies coefficient names.
#' txtfun <- function(x) gsub('_',' ',substr(x,1,nchar(x)-3))
#'
#' # Calculate and plot disability weights
#' sc <- standardizeCoef(model1, plotf = TRUE, namesf = txtfun)
#' sc
standardizeCoef <- function (model,
                             ordered = TRUE,
                             plotf = FALSE,
                             plotpval = FALSE,
                             mar = c(15, 4, 1, 1), oma = c(0, 0, 0, 0),
                             YLab = "Disability weight",
                             YLab.cex = 1.1,
                             namesf = identity, ...) {

  z <- (model$coef)[seq_len(model$parcount[1])]
  if (class(namesf)=='function') names(z) <- namesf(names(z)) else names(z) <- namesf

  if (ordered) {
    oz <- order(z, decreasing = TRUE)
    cfm <- z[oz]
  } else cfm <- z
  r <- range(model$y_latent_i)
  res <- as.matrix((cfm)/diff(r))

  if (plotf) {
    opar <- graphics::par(c("mar", "oma"))
    graphics::par(mar = mar, oma = oma)
    rr <- graphics::barplot(t(res), las = 3, ...)
    if(plotpval) {
      y <- summary(model)
      pval <- format(round(y$coef$`Pr(>|z|)`[seq_len(model$parcount[1])],4),digits=4,scientific=FALSE)[oz]
      yr <- res/2
      ind <- yr < max(res)*0.1
      yr[ind] <- (res+max(res)*0.1)[ind]
      graphics::text(rr,yr,paste('P =',pval),srt=90,col=c('white','black')[1+ind])
    }
    graphics::mtext(YLab, 2, cex = YLab.cex, line = 2.5)
    suppressWarnings(graphics::par(opar))
  }
  if (plotf)
    invisible(res)
  else return(res)
}


#' @rdname standardizeCoef
#' @export
standardiseCoef<-standardizeCoef


#' @rdname standardizeCoef
#' @export
disabilityWeights<-standardizeCoef


#' Calculate threshold cut-points using Jurges' method
#' @description
#' Calculate threshold cut-points using Jurges' method.
#' @param model a fitted \code{hopit} model.
#' @param decreasing.levels logical indicating if self-reported health classes are ordered in decreasing order.
#' @param subset an optional vector specifying a subset of observations.
#' @param plotf a logical indicating if to plot the results.
#' @param XLab,XLab.cex label and size of the label for x axis.
#' @param YLab,YLab.cex label and size of the label for y axis.
#' @param mar,oma see \code{\link{par}}.
#' @param group.labels.type position of the legend. One of \code{middel}, \code{border}, or \code{none}.
#' @return a list with following components:
#'  \item{cutpoints}{ cutpoints for adjusted categorical response levels with corresponding percentiles of latent index.}
#'  \item{adjused.levels}{ adjusted categorical response levels for each individual.}
#' @references \insertRef{Jurges2007}{hopit}
#' @author Maciej J. Danko
#' @export
#' @seealso \code{\link{latentIndex}}, \code{\link{standardiseCoef}}, \code{\link{getLevels}}, \code{\link{hopit}}.
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # the order of response levels is decreasing (from the best health to the worst health)
#' levels(healthsurvey$health)
#'
#' # Example 1 ---------------------
#'
#' # fitting a model
#' model1 <- hopit(latent.formula = health ~ hypertension + high_cholesterol +
#'                 heart_attack_or_stroke + poor_mobility + very_poor_grip +
#'                 depression + respiratory_problems +
#'                 IADL_problems + obese + diabetes + other_diseases,
#'               thresh.formula = ~ sex + ageclass + country,
#'               decreasing.levels = TRUE,
#'               control = list(trace = FALSE),
#'               data = healthsurvey)
#'
#' # health index cut-points
#' z <- getCutPoints(model = model1)
#' z$cutpoints
#'
#' # adjusted health levels for individuals: Jurges method
#' rev(table(z$adjused.levels))
#'
#' # original health levels for individuals
#' table(model1$y_i)
#'
#' # adjusted health levels for individuals: Estimated model thresholds
#' table(model1$Ey_i)
getCutPoints <- function(model, subset=NULL, plotf = TRUE, mar=c(4,4,1,1),oma=c(0,0,0,0),
                         XLab='Health index', XLab.cex=1.1, YLab='Counts', YLab.cex=1.1,
                         decreasing.levels=TRUE, group.labels.type=c('middle','border','none')){
  if (length(subset) == 0) subset=seq_along(model$y_i)
  Y <- model$y_i[subset]
  if (decreasing.levels) dorev <- rev else dorev <- identity
  h.index <- healthIndex(model, subset, plotf = FALSE)
  tY <- table(Y)
  tmp <- dorev(as.vector(tY))
  invcs <- (cumsum(tmp)/sum(tmp))[-length(tmp)]
  R1 <- stats::quantile(h.index, invcs)

  lv <- dorev(as.character(levels(model$y_i)))
  Nm <- paste(lv[-length(lv)],lv[-1],sep=' | ')
  Nm <- sapply(Nm, function(k) bquote(bold(.(k))))
  if (plotf) {
    group.labels.type<-tolower(group.labels.type[1])
    if(group.labels.type %notin%  c('middle','border','none')) stop (hopit_msg(84),call.=NULL)
    opar <- graphics::par(c('mar','oma'))
    graphics::par(mar=mar, oma=oma)
    z<-graphics::hist(h.index, 200,xlab='',ylab='' ,
            main='', yaxs='i', col=grDevices::grey(0.8, alpha = 0.5),border=grDevices::grey(0.4, alpha = 0.5))
    if (group.labels.type == 'border') {
      for (j in seq_along(Nm)) graphics::text(x=R1[j],y=(1.1*max(z$counts))/2,labels=Nm[[j]],
                                    srt=90,pos=2,offset=0.67,col=2)
    } else if (group.labels.type == 'middle'){
      R11=-diff(c(0,R1,1))/2+c(R1,1)+graphics::strheight('S',units='figure')/2
      for (j in seq_along(lv)) graphics::text(x=R11[j],y=(3*1.1*max(z$counts))/4,labels=lv[j],
                                    srt=90,pos=3,offset=0.67,col=2)
    }
    graphics::box()
    graphics::abline(v=R1,lwd=2,col=2)
    graphics::mtext(XLab, 1, cex=XLab.cex, line = 2.5)
    graphics::mtext(YLab, 2, cex=YLab.cex, line = 2.5)
    suppressWarnings(graphics::par(opar))
  }
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
    adjused.levels<- cut(h.index, CIN,labels= dorev(levels(Y)))
  } else adjused.levels <- NA
  res <- list(cutpoints=R1, adjused.levels=(adjused.levels))
  if (plotf) invisible(res) else return(res)
}


#' Summarize adjusted and original self-rated response levels
#' @description
#' Summarize adjusted and original self-rated response levels.
#' @param model a fitted \code{hopit} model.
#' @param formula a formula containing the grouping variables. It is by default set to threshold formula.
#' @param data data used to fit the model.
#' @param plotf a logical indicating if to plot the results.
#' @param sep separator for levels names.
#' @param decreasing.levels logical indicating if self-reported health classes are ordered in increasing order.
#' @param sort.flag logical indicating if to sort the levels.
#' @param mar,oma see \code{\link{par}}.
#' @param YLab,YLab.cex label and size of the label for y axis.
#' @param legbg legend background color. See \code{bg} parameter in \code{\link{legend}}.
#' @param legbty legend box type. See \code{bty} parameter in \code{\link{legend}}.
#' @return a list with following components:
#'  \item{original}{ frequencies of original response levels for selected groups/categories.}
#'  \item{adjusted}{ frequencies of adjusted response levels (Jurges 2007 method) for selected groups/categories.}
#'  \item{N.original}{ numbers of original response levels for selected groups/categories.}
#'  \item{N.adjusted}{ numbers of adjusted response levels (Jurges 2007 method) for selected groups/categories.}
#'  \item{categories}{ selected groups/categories used in summary.}
#'  \item{tab}{ original vs. adjusted contingency table.}
#'  \item{mat}{ a matrix with columns: grouping variables, original response levels, adjusted response levels.
#'  Each row corresponds to a single individual from the data used to fit the model.}
#' @references \insertRef{Jurges2007}{hopit}
#' @author Maciej J. Danko
#' @export
#' @seealso \code{\link{getCutPoints}}, \code{\link{latentIndex}}, \code{\link{standardiseCoef}}, \code{\link{hopit}}.
#' @examples
#' # DATA
#' data(healthsurvey)
#'
#' # the order of response levels is decreasing (from the best health to the worst health)
#' levels(healthsurvey$health)
#'
#' # fitting a model
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
#' # summary by country
#' hl <- getLevels(model=model1, formula=~ country,
#'                 data = healthsurvey,
#'                 sep=' ', plotf=TRUE)
#'
#' # differences in frequencies between original and adjusted health levels
#' round(100*(hl$original - hl$adjusted),2)
#'
#' # extract good and bad health (combine levels)
#' Org <- cbind(bad = rowSums(hl$original[,1:2]), good = rowSums(hl$original[,4:5]))
#' Adj <- cbind(bad = rowSums(hl$adjusted[,1:2]), good = rowSums(hl$adjusted[,4:5]))
#' round(100*(Org - Adj),2)
#'
#' # plot the differences
#' barplot(t(Org - Adj), beside = TRUE, density = 20, angle = c(-45, 45),
#'         col = c('pink4', 'green2'), ylab = 'Original - adjusted reported health frequencies')
#' abline(h = 0); box()
#' legend('top', c('Bad health','Good health'), density = 20, angle = c(-45, 45),
#'        fill = c('pink4', 'green2'), bty = 'n', cex = 1.2)
#'
#' # in country X the bad health seems to be over-reported and good health under reported,
#' # in country Z the good health is highly over-reported.
#'
#' # Example 2 ---------------------
#'
#' # summary by gender and age
#' hl <- getLevels(model = model1, formula=~ sex + ageclass,
#'                 data = healthsurvey,
#'                 sep=' ', plotf=TRUE)
#'
#' # differences in frequencies between original and adjusted health levels
#' round(100*(hl$original - hl$adjusted),2)
#'
#' # extract good health levels (combined "Very good" and "Excelent" levels)
#' Org <- rowSums(hl$original[,4:5])
#' Adj <- rowSums(hl$adjusted[,4:5])
#' round(100*(Org - Adj),2)
#'
#' pmar <- par('mar'); par(mar = c(9.5, pmar[2:4]))
#' barplot(Org-Adj, ylab = 'Original - adjusted reported good health frequencies', las = 3,
#'         density = 20, angle = c(45, -45), col = c('blue', 'orange'))
#' abline(h = 0); box(); par(mar = pmar)
#' legend('top', c('Man','Woman'), density = 20, angle = c(-45, 45),
#'        fill = c('blue', 'orange'), bty = 'n', cex = 1.2)
#'
#' # the results show that women in general tends to over-report good health.
#' # men in ages 50-59 greatly under-report good health.
#'
#' # more examples can be found in the description of boot.hopit() function.
getLevels<-function(model,
                    formula=model$thresh.formula,
                    data = environment(model$thresh.formula),
                    decreasing.levels = TRUE,
                    sort.flag = FALSE,
                    plotf = TRUE, sep='_',
                    mar = c(7,2,1.5,0.5), oma = c(0,3,0,0),
                    YLab = 'Fraction [%]',
                    YLab.cex = 1.1,
                    legbg = grDevices::adjustcolor('white',alpha.f=0.4),
                    legbty = 'o'
){
  if (class(formula)=='formula') inte_ <- formula2classes(formula, data, sep=sep, return.matrix = TRUE) else
    stop(call.=NULL, hopit_msg(86))
  inte <- inte_$x
  namind <- inte_$class.mat
  nam <- levels(inte)
  cpall <- getCutPoints(model, plotf = FALSE, decreasing.levels = decreasing.levels)
  TAB1 <- round(table(original=model$y_i, adjusted=cpall$adjused.levels)*100/length(model$y_i),2)
  tmp <- untable(t(table(factor(model$y_i,levels=levels(cpall$adjused.levels)), inte)))
  N1 <- tmp
  tmp <-tmp/rowSums(tmp)
  if (sort.flag) {
    oD1 <- order(tmp[,NCOL(tmp)]+tmp[,NCOL(tmp)-1])
    tmp <- tmp[oD1,]
    orignalind <- namind[oD1,]
  } else orignalind <- namind

  tmp2 <- untable(t(table(cpall$adjused.levels, inte)))
  N2 <- tmp2
  tmp2 <- tmp2/rowSums(tmp2)
  if (sort.flag) {
    oD2 <- order(tmp2[,NCOL(tmp2)]+tmp2[,NCOL(tmp2)-1])
    tmp2 <- tmp2[oD2,]
    adjustedind <- namind[oD2,]
  } else adjustedind <- namind

  if (plotf) {
    opar <- graphics::par(c('mar','oma','mfrow'))
    graphics::par(mfrow=c(1,2))
    graphics::par(mar=mar,oma=oma)
    graphics::barplot(t(tmp),las=3,main='Original')
    graphics::barplot(t(tmp2),las=3,main='Adjusted', legend.text=TRUE,
            args.legend = list(x='center', box.col=NA,
                               bg=legbg, bty=legbty))
    graphics::par(mfrow=c(1,1))
    graphics::par(mar=mar,oma=rep(0,4))
    graphics::mtext(YLab,2,cex=YLab.cex)
    suppressWarnings(graphics::par(opar))
  }
  res <- list(original= tmp,
              adjusted= tmp2,
              tab= TAB1,
              N.original= N1,
              N.adjusted= N2,
              categories = orignalind, # = adjustedind
              mat=cbind(inte_$mat,
                        original= model$y_i,
                        adjusted= cpall$adjused.levels))
  class(res) <- 'healthlevels'
  if (plotf) invisible(res) else return(res)
}
