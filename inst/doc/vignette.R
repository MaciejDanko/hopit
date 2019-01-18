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

## ---- echo=TRUE, cache=TRUE----------------------------------------------
# first determine the order of the dependent variable
levels(healthsurvey$health)

# the order is decreasing (from best health to the worst health)
# so we set: decreasing.levels = TRUE
model1<- hopit(latent.formula = health ~ hypertenssion + high_cholesterol + 
                             heart_atack_or_stroke + poor_mobility + very_poor_grip + 
                             depression + respiratory_problems + 
                             IADL_problems + obese + diabetes + other_diseases, 
               thresh.formula = ~ sex + ageclass,
               decreasing.levels = TRUE,
               control=list(trace=FALSE),
               data = healthsurvey)

summary(model1)

## ---- echo=TRUE----------------------------------------------------------
# extract parameters in a form of list
cm1 <- coef(model1, aslist = TRUE)

# names of returned coefficients
names(cm1)

# extracting latent health coefficients
cm1$latent.params

## ---- echo=TRUE, cache=TRUE----------------------------------------------
model2<- hopit(latent.formula = health ~ hypertenssion + high_cholesterol + 
                      heart_atack_or_stroke + poor_mobility + 
                      very_poor_grip + depression + respiratory_problems + 
                      IADL_problems + obese + diabetes + other_diseases, 
               thresh.formula = ~ sex + ageclass + country,
               decreasing.levels = TRUE,
               control=list(trace=FALSE),
               data = healthsurvey)

## ---- echo=TRUE----------------------------------------------------------
AIC(model2, model1)

## ---- echo=TRUE----------------------------------------------------------
anova(model2, model1)

## ---- echo=TRUE----------------------------------------------------------
model3<- hopit(latent.formula = health ~ hypertenssion + high_cholesterol + 
                      heart_atack_or_stroke + poor_mobility + 
                      very_poor_grip + depression + respiratory_problems + 
                      IADL_problems + obese + diabetes + other_diseases, 
               thresh.formula = ~ sex * ageclass + country,
               decreasing.levels = TRUE,
               control=list(trace=FALSE),
               data = healthsurvey)

print(anova(model3,model2), short=TRUE)


## ---- echo=TRUE, cache=TRUE----------------------------------------------
design <- svydesign(ids = ~ country + psu, weights = healthsurvey$csw, 
                    data = healthsurvey)

model2s<- hopit(latent.formula = health ~ hypertenssion + high_cholesterol + 
                       heart_atack_or_stroke + poor_mobility + 
                       very_poor_grip + depression + respiratory_problems + 
                       IADL_problems + obese + diabetes + other_diseases, 
                thresh.formula = ~ sex + ageclass + country,
                decreasing.levels = TRUE,
                design = design,
                control=list(trace=FALSE),
                data = healthsurvey)

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


## ---- echo=TRUE, cache=TRUE, fig.height = 4, fig.width = 6, fig.align = "center"----
# Health levels for combination of age and gender, and pooled country of origin.
hl <- getLevels(model=model3, formula=~ sex + ageclass, data = healthsurvey, 
                      sep=' ', plotf=TRUE, legbty = 'n')

## ---- echo=TRUE----------------------------------------------------------
round(100*(hl$original - hl$adjusted),2)


## ---- echo=TRUE, fig.height = 6, fig.width = 5, fig.align = "center", cache=TRUE----
# Function to be bootstraped
diff_BadHealth <- function(model, data) {
  hl <- getLevels(model=model, formula=~ sex + ageclass, data = data, 
                  sep=' ', plotf=FALSE)
  hl$original[,1] + hl$original[,2] - hl$adjusted[,1]- hl$adjusted[,2]
}

# Estimate of the difference
est.org <- diff_BadHealth(model = model3, data = healthsurvey)

# Perform the bootstrap
B <- boot_hopit(model = model3, data = healthsurvey, 
                func = diff_BadHealth, nboot = 100)

# Calcualte lower and upper bounds using percentile method
est.CI <- boot_hopit_CI(B)

# Plotting the difference and its (assymetrical) confidence intervals
pmar <- par('mar'); par(mar = c(9.5,pmar[2:4]))
m <- max(abs(est.CI))
pos <- barplot(est.org, names.arg = names(est.org), las = 3, ylab = 'Orginal - Adjusted', 
               ylim=c(-m, m), density = 20, angle = c(45, -45), col = c('blue', 'orange'))
for (k in seq_along(pos)) lines(c(pos[k,1],pos[k,1]), est.CI[,k], lwd = 2, col = 2)
abline(h = 0); box(); par(mar = pmar)


