## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
unlink('vignettes/vignette_cache', recursive = TRUE)

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  library(devtools)
#  install_github("maciejdanko/hopit")

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  library(hopit)

## ---- echo=FALSE,  results='hide', eval=TRUE, include=FALSE--------------
g <- capture.output(library(hopit))

## ---- echo=TRUE----------------------------------------------------------
# load *healthsurvey* dataset
data(healthsurvey)

# horizontal view on the dataset (omitting ID)
print(t(healthsurvey[1:6,-1]), quote=FALSE, na.print='NA', right=TRUE)

## ---- echo=TRUE----------------------------------------------------------
# extract parameters in a form of list
cm1 <- coef(model1, aslist = TRUE)

# names of returned coefficients
names(cm1)

# latent health variables
cm1$latent.params

## ---- echo=TRUE----------------------------------------------------------
AIC(model2, model1)

## ---- echo=TRUE----------------------------------------------------------
anova(model2, model1)

## ---- echo=TRUE----------------------------------------------------------
model3<- hopit(latent.formula = health ~ hypertenssion * high_cholesterol + 
                             heart_atack_or_stroke + poor_mobility + 
                             very_poor_grip + depression + respiratory_problems + 
                             IADL_problems + obese + diabetes + other_diseases, 
               thresh.formula = ~ sex * ageclass + country,
               decreasing.levels = TRUE,
               control=list(trace=FALSE),
               data = healthsurvey)

print(anova(model3,model2), short=TRUE)


## ---- echo=TRUE----------------------------------------------------------
cbind('No survey design'=coef(model2,aslist=TRUE)$latent.par,
      'Has survey design'=coef(model2s,aslist=TRUE)$latent.par)

## ---- echo=TRUE----------------------------------------------------------
profile(model3)

## ---- echo=TRUE----------------------------------------------------------
model3$coef.ls$latent.params

## ---- echo=TRUE, fig.height = 5, fig.width = 5, fig.align = "center"-----
# A function that modifies coefficient names.  
txtfun <- function(x) gsub('_',' ',substr(x,1,nchar(x)-3))

# Calcualte and plot disability weights
sc <- standardizeCoef(model3, plotf = TRUE, namesf = txtfun)
sc

## ---- echo=TRUE, fig.height = 4, fig.width = 5, fig.align = "center"-----
hi <- latentIndex(model3, plotf = TRUE, response = "data", 
                  ylab = 'Health index', col='deepskyblue3')

## ---- echo=TRUE, fig.height = 4, fig.width = 5, fig.align = "center"-----
hi <- latentIndex(model3, plotf = TRUE, response = "fitted", 
                  ylab = 'Health index', col='deepskyblue3')

## ---- echo=TRUE, fig.height = 4, fig.width = 5, fig.align = "center"-----
hi <- latentIndex(model3, plotf = TRUE, response = "Jurges", 
                  ylab = 'Health index', col='deepskyblue3')

## ---- echo=TRUE, fig.height = 4, fig.width = 5, fig.align = "center"-----
z=getCutPoints(model=model3)

# Health index cut-points
z$cutpoints

#Adjusted health levels for individuals: Jurges method
table(z$adjused.health.levels)

#Adjusted health levels for individuals: Estimated model thresholds
table(model3$Ey_i)

#Original health levels for individuals
table(model3$y_i)


## ---- echo=TRUE----------------------------------------------------------
round(100*(hl$original - hl$adjusted),2)


