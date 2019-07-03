rm(list = ls())
graphics.off()

library(testthat)
library(hopit)

# Test function -------------------------------

test_hopit <- function (object, data) {
  skip_on_cran()
  N <- object$N
  cat('class ')
  expect_s3_class(object, "hopit")
  cat('print ')
  expect_output(print(object))
  cat('print.summary ')
  expect_output(print(summary(object)))
  cat('non-empty tests ')
  expect_true(length(object$estfun)>0)
  expect_true(length(object$hessian)>0)
  expect_true(length(object$vcov)>0)
  expect_true(length(object$vcov.basic)>0)
  expect_equal(length(object$y_i), N)
  expect_equal(length(object$y_latent_i), N)
  expect_equal(length(object$Ey_i), N)
  expect_equal(nrow(data), N)
  cat('LL ')
  expect_gt(object$LL, -Inf)
  expect_false(is.na(object$LL))
  cat('design ')
  if (length(object$design)) expect_true(is.na(object$AIC))
  if (!length(object$design)) expect_identical(object$vcov, object$vcov.basic)
  # tmp <- profile(object)
  # expect_output(print(tmp, plotf = FALSE))
  # expect_equal(sum(is.na(tmp)),0)
  # expect_equal(capture_output(print(tmp, plotf = FALSE)),
  #             "All parameters seem to be at arg.max (at optimum).")
  cat('HI ')
  expect_true(length(healthIndex(object))>0)
  expect_equal(length(healthIndex(object)), N)
  cat('D ')
  expect_true(length(disabilityWeights(object))>0)
  expect_equal(length(disabilityWeights(object)),
               length(coef(object,aslist = T)$latent.params))
  cat('coef ')
  expect_equal(length(coef(object)), sum(object$parcount))
  cat('sigma ')
  if (!object$hasdisp) {
    expect_equal(round(sigma(object),8), 1)
  }
  cat('levels ')
  hl <- getLevels(model=object, formula=~ sex + ageclass, data = healthsurvey,
                  sep=' ')
  expect_equal(dim(hl$original), dim(hl$adjusted))
  expect_equal(length(hl), 7)
  cat('cut-points ')
  cp <- getCutPoints(object)
  expect_equal(length(cp$cutpoints), object$J-1)
  expect_equal(length(cp$adjusted.levels), N)
  # print(which(is.na(cp$adjusted.levels)))
  expect_equal(sum(is.na(cp$adjusted.levels)),0)
  cat('end \n')
  invisible()
}

# Initialize ----------------------------------

# extended data
newhealthsurvey <- healthsurvey
newhealthsurvey2 <- healthsurvey
newhealthsurvey2$sex[1000:1005] <- NA
newhealthsurvey2$diabetes[1002:1008] <- NA
newhealthsurvey2$health[1005:1012] <- NA

newhealthsurvey$cont_var <- sample(5000:5015,nrow(newhealthsurvey),replace=TRUE)
newhealthsurvey$age <- 50+10*runif(1e4)*as.numeric(newhealthsurvey$ageclass)
newhealthsurvey$health2 <- as.numeric(newhealthsurvey$health)
newhealthsurvey$age2 <- 5*(newhealthsurvey$age)
newhealthsurvey$var1 <- newhealthsurvey$health2*newhealthsurvey$age
newhealthsurvey$var2 <- newhealthsurvey$health2-1
newhealthsurvey$var3 <- newhealthsurvey$sex
newhealthsurvey$health3 <- newhealthsurvey$health
levels(newhealthsurvey$health3)<-c(levels(newhealthsurvey$health3), 'very poor')
newhealthsurvey$sex3 <- newhealthsurvey$sex
levels(newhealthsurvey$sex3)<-c(levels(newhealthsurvey$sex), 'unknown')
newhealthsurvey$diabetes3 <- newhealthsurvey$diabetes
levels(newhealthsurvey$diabetes3)<-c(levels(newhealthsurvey$diabetes), 'unknown')
newhealthsurvey$depression3 <- newhealthsurvey$depression
levels(newhealthsurvey$depression3) <- c(levels(newhealthsurvey$depression), 'unknown', 'very unknown')

# formulas to be tested
latent.formula.1 <- health ~ hypertension + high_cholesterol +
  heart_attack_or_stroke + poor_mobility + very_poor_grip +
  depression + respiratory_problems +
  IADL_problems + obese + diabetes + other_diseases
latent.formula.2 <- health ~ hypertension
latent.formula.3 <- health ~ age
latent.formula.4 <- health ~ 1
latent.formula.5 <- ~ 1
latent.formula.6 <- health ~ hypertension + high_cholesterol +
  heart_attack_or_stroke + poor_mobility + very_poor_grip +
  depression + respiratory_problems +
  IADL_problems + obese + diabetes : sex + other_diseases
latent.formula.7 <- health ~ hypertension + high_cholesterol +
  heart_attack_or_stroke + poor_mobility + very_poor_grip +
  depression + respiratory_problems +
  IADL_problems + obese + diabetes * sex + other_diseases
latent.formula.8 <- health ~ hypertension + high_cholesterol +
  heart_attack_or_stroke + poor_mobility + very_poor_grip +
  depression + respiratory_problems +
  IADL_problems + obese + diabetes + other_diseases + diabetes : sex
latent.formula.9 <- health ~ age * heart_attack_or_stroke
latent.formula.A <- health ~ hypertension + offset(age)
latent.formula.B <- health ~ 1
latent.formula.C <- ~ hypertension + high_cholesterol
latent.formula.D <- ~ offset(age)
latent.formula.E <- health ~ hypertension + I(age^2)
latent.formula.F <- health2 ~ hypertension + high_cholesterol +
  heart_attack_or_stroke + poor_mobility + very_poor_grip +
  depression + respiratory_problems +
  IADL_problems + obese + diabetes + other_diseases
latent.formula.G <- health ~ hypertension + high_cholesterol +
  heart_attack_or_stroke + poor_mobility + very_poor_grip +
  depression + respiratory_problems +
  IADL_problems + obese + diabetes + other_diseases + var3
latent.formula.H <- health ~ hypertension + high_cholesterol +
  heart_attack_or_stroke + poor_mobility + very_poor_grip +
  depression + respiratory_problems +
  IADL_problems + obese + other_diseases
latent.formula.I <- health3 ~ hypertension + high_cholesterol +
  heart_attack_or_stroke + poor_mobility + very_poor_grip +
  depression3 + respiratory_problems +
  IADL_problems + obese + diabetes3 + other_diseases

thresh.formula.1 <-  ~ sex + ageclass + country
thresh.formula.2 <-  health ~ sex + ageclass + country
thresh.formula.3 <-  ~ country
thresh.formula.4 <-  ~ 1
thresh.formula.5 <-  ~ cont_var #check also transformed data
thresh.formula.6 <-  ~ cont_var + country #check also transformed data
thresh.formula.7 <-  ~ cont_var * country #check also transformed data
thresh.formula.8 <-  ~ sex : country
thresh.formula.9 <-  ~ sex * country
thresh.formula.A <-  ~ sex + country : heart_attack_or_stroke
thresh.formula.B <-  ~ sex + country * heart_attack_or_stroke
thresh.formula.C <-  ~ sex + country + country : heart_attack_or_stroke
thresh.formula.D <-  ~ offset(age)
thresh.formula.E <-  ~ I(age^2)
thresh.formula.F <-  ~ 2
#thresh.formula.G <-  ~ -1
#thresh.formula.H <-  ~ 0
thresh.formula.I <-  ~ sex + health
thresh.formula.J <-  ~ sex * ageclass * country
thresh.formula.K <-  ~ sex + country + age2
thresh.formula.L <-  ~ sex + country + var1
thresh.formula.M <-  ~ sex + country + var2
thresh.formula.N <-  ~ sex * ageclass + country
thresh.formula.O <-  ~ sex3 + ageclass + country


# Messages

expect_message(hopit(latent.formula = latent.formula.I,
                     thresh.formula = thresh.formula.O,
                     data = newhealthsurvey,
                     decreasing.levels = TRUE))

# Warnings ------------------------

# non-finite value supplied by optim
expect_warning(hopit(latent.formula = latent.formula.1,
                     thresh.formula = thresh.formula.2,
                     data = newhealthsurvey,
                     decreasing.levels = TRUE))


# Errors --------------------------

#Error in na.fail.default(data) : missing values in object
expect_error(hopit(latent.formula = latent.formula.1,
                   thresh.formula = thresh.formula.1,
                   data = newhealthsurvey2,
                   decreasing.levels = TRUE))
# hopit(latent.formula = latent.formula.1,
#       thresh.formula = thresh.formula.1,
#       data = newhealthsurvey2,
#       na.action = na.omit,
#       decreasing.levels = TRUE)

# initial value in 'vmmin' is not finite
expect_error(hopit(latent.formula = latent.formula.G,
                   thresh.formula = thresh.formula.1,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# initial value in 'vmmin' is not finite
expect_error(hopit(latent.formula = latent.formula.3,
                   thresh.formula = thresh.formula.M,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# initial value in 'vmmin' is not finite
expect_error(hopit(latent.formula = latent.formula.3,
                   thresh.formula = thresh.formula.L,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# non-finite value supplied by optim
expect_error(hopit(latent.formula = latent.formula.3,
                   thresh.formula = thresh.formula.K,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# No latent variales detected in the latent formula.
expect_error(hopit(latent.formula = latent.formula.B,
                   thresh.formula = thresh.formula.1,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# Response must be a factor with ordered levels.
expect_error(hopit(latent.formula = latent.formula.F,
                   thresh.formula = thresh.formula.1,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# No main effects found for the country:heart_attack_or_stroke interaction.
expect_error(hopit(latent.formula = latent.formula.1,
                   thresh.formula = thresh.formula.A,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# No main effects found for the country:heart_attack_or_stroke interaction.
expect_error(hopit(latent.formula = latent.formula.2,
                   thresh.formula = thresh.formula.A,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# The same terms found in latent and threshold model: heart_attack_or_stroke.
expect_error(hopit(latent.formula = latent.formula.1,
                   thresh.formula = thresh.formula.B,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# No response detected in the latent formula.
expect_error(hopit(latent.formula = latent.formula.C,
                   thresh.formula = thresh.formula.9,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# No main effects found for the sex:country interaction.
expect_error(hopit(latent.formula = latent.formula.1,
                   thresh.formula = thresh.formula.8,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# # Cannot find starting parameters. Please consider scaling, tran....
# expect_error(hopit(latent.formula = latent.formula.1,
#                    thresh.formula = thresh.formula.5,
#                    data = newhealthsurvey,
#                    decreasing.levels = TRUE))
#
# # Cannot find starting parameters. Please consider scaling, tran....
# expect_error(hopit(latent.formula = latent.formula.1,
#                    thresh.formula = thresh.formula.6,
#                    data = newhealthsurvey,
#                    decreasing.levels = TRUE))

# Offset not supported.
expect_error(hopit(latent.formula = latent.formula.D,
                   thresh.formula = thresh.formula.3,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# Offset not supported.
expect_error(hopit(latent.formula = latent.formula.A,
                   thresh.formula = thresh.formula.3,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# Offset not supported.
expect_error(hopit(latent.formula = latent.formula.1,
                   thresh.formula = thresh.formula.D,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# I(...) not implemented yet.
expect_error(hopit(latent.formula = latent.formula.1,
                   thresh.formula = thresh.formula.E,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# I(...) not implemented yet.
expect_error(hopit(latent.formula = latent.formula.E,
                   thresh.formula = thresh.formula.3,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# invalid model formula in ExtractVars
expect_error(hopit(latent.formula = latent.formula.1,
                   thresh.formula = thresh.formula.F,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# Response found in latent or threshold formulas.
expect_error(hopit(latent.formula = latent.formula.1,
                   thresh.formula = thresh.formula.I,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# No latent variales detected in the latent formula.
expect_error(hopit(latent.formula = latent.formula.4,
                   thresh.formula = thresh.formula.4,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# No latent variales detected in the latent formula. No response detected in
# the latent formula.
expect_error(hopit(latent.formula = latent.formula.5,
                   thresh.formula = thresh.formula.4,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# No main effects found for the diabetes:sex interaction.
expect_error(hopit(latent.formula = latent.formula.6,
                   thresh.formula = thresh.formula.4,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# No main effects found for the diabetes:sex interaction.
expect_error(hopit(latent.formula = latent.formula.6,
                   thresh.formula = thresh.formula.1,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# The same terms found in latent and threshold model: sex. ...
expect_error(hopit(latent.formula = latent.formula.7,
                   thresh.formula = thresh.formula.1,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# No main effects found for the diabetes:sex interaction.
expect_error(hopit(latent.formula = latent.formula.8,
                   thresh.formula = thresh.formula.3,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# No main effects found for the country:heart_attack_or_stroke interaction.
expect_error(hopit(latent.formula = latent.formula.2,
                   thresh.formula = thresh.formula.C,
                   data = newhealthsurvey,
                   decreasing.levels = TRUE))

# Multiple weights specification detected.
# Please use either design or weights parameter.
expect_error(hopit(latent.formula = latent.formula.1,
                   thresh.formula = thresh.formula.1,
                   data = newhealthsurvey,
                   weights = runif(10000),
                   design = svydesign(ids = ~ country + psu,
                                      weights = healthsurvey$csw,
                                      data = newhealthsurvey),
                   fit.sigma = TRUE,
                   decreasing.levels = TRUE))

# # Cannot find starting parameters. Please
# expect_warning(hopit(latent.formula = latent.formula.1,
#                    thresh.formula = thresh.formula.7,
#                    data = newhealthsurvey,
#                    decreasing.levels = TRUE))


# Passes -------------------------

test_that("Hopit works",{
  skip_on_cran()
  print('m0')
  m0  <<- hopit(latent.formula = latent.formula.8,
                thresh.formula = thresh.formula.1,
                data = newhealthsurvey,
                decreasing.levels = TRUE)
  print('m1')
  m1  <- hopit(latent.formula = latent.formula.8,
               thresh.formula = thresh.formula.1,
               data = newhealthsurvey,
               decreasing.levels = TRUE)
  print('m2')
  m2  <- hopit(latent.formula = latent.formula.7,
               thresh.formula = thresh.formula.3,
               data = newhealthsurvey,
               decreasing.levels = TRUE)

  # m3  <- hopit(latent.formula = latent.formula.1,
  #              thresh.formula = thresh.formula.G,
  #              data = newhealthsurvey,
  #              decreasing.levels = TRUE)
  #
  # m4  <- hopit(latent.formula = latent.formula.1,
  #              thresh.formula = thresh.formula.H,
  #              data = newhealthsurvey,
  #              decreasing.levels = TRUE)
  print('m5')
  m5  <- hopit(latent.formula = latent.formula.1,
               thresh.formula = thresh.formula.9,
               data = newhealthsurvey,
               decreasing.levels = TRUE)
  print('m6')
  m6  <- hopit(latent.formula = latent.formula.9,
               thresh.formula = thresh.formula.9,
               data = newhealthsurvey,
               control = list(transform.latent = 'min'),
               decreasing.levels = TRUE)
  print('m7')
  m7  <- hopit(latent.formula = latent.formula.1,
               thresh.formula = thresh.formula.5,
               data = newhealthsurvey,
               control = list(transform.thresh = 'standardize'),
               decreasing.levels = TRUE)
  print('m8')
  m8  <<- hopit(latent.formula = latent.formula.1,
                thresh.formula = thresh.formula.5,
                data = newhealthsurvey,
                control = list(transform.thresh = 'min'),
                decreasing.levels = TRUE)
  print('m9')
  m9  <- hopit(latent.formula = latent.formula.3,
               thresh.formula = thresh.formula.5,
               data = newhealthsurvey,
               control = list(transform.thresh = 'min'),
               decreasing.levels = TRUE)
  print('mA')
  mA  <- hopit(latent.formula = latent.formula.3,
               thresh.formula = thresh.formula.1,
               data = newhealthsurvey,
               control = list(transform.thresh = 'min'),
               decreasing.levels = TRUE)
  print('mB')
  mB  <- hopit(latent.formula = latent.formula.1,
               thresh.formula = thresh.formula.7,
               data = newhealthsurvey,
               control = list(transform.thresh = 'standardize'),
               decreasing.levels = TRUE)
  print('mC')
  mC  <- hopit(latent.formula = latent.formula.2,
               thresh.formula = thresh.formula.B,
               data = newhealthsurvey,
               decreasing.levels = TRUE)
  print('mD')
  mD  <- hopit(latent.formula = latent.formula.1,
               thresh.formula = thresh.formula.C,
               data = newhealthsurvey,
               decreasing.levels = TRUE)
  print('mE')
  mE  <<- hopit(latent.formula = latent.formula.9,
                thresh.formula = thresh.formula.1,
                data = newhealthsurvey,
                decreasing.levels = TRUE)
  print('mF')
  mF  <- hopit(latent.formula = latent.formula.1,
               thresh.formula = thresh.formula.J,
               data = newhealthsurvey,
               control = list(bgfs.maxit = 1e5,
                              bgfs.reltol = 1e-12),
               decreasing.levels = TRUE)
  print('mG')
  mG  <- hopit(latent.formula = latent.formula.1,
               thresh.formula = thresh.formula.N,
               data = newhealthsurvey,
               decreasing.levels = TRUE)
  print('mH')
  mH  <<- hopit(latent.formula = latent.formula.1,
                thresh.formula = thresh.formula.1,
                data = newhealthsurvey,
                fit.sigma = TRUE,
                decreasing.levels = TRUE)
  print('mI')
  mI  <<- hopit(latent.formula = latent.formula.1,
                thresh.formula = thresh.formula.1,
                data = newhealthsurvey,
                design = svydesign(ids = ~ country + psu,
                                   weights = healthsurvey$csw,
                                   data = newhealthsurvey),
                fit.sigma = TRUE,
                decreasing.levels = TRUE)
  print('mJ')
  mJ  <- hopit(latent.formula = latent.formula.1,
               thresh.formula = thresh.formula.7,
               data = newhealthsurvey,
               control = list(transform.latent = 'standardize',
                              transform.thresh = 'standardise'),
               decreasing.levels = TRUE)
  print('mK')
  mK  <- hopit(latent.formula = latent.formula.1,
               thresh.formula = thresh.formula.1,
               data = newhealthsurvey,
               weights = runif(10000),
               fit.sigma = FALSE,
               decreasing.levels = TRUE)
  print('mL')
  mL  <<- hopit(latent.formula = latent.formula.1,
                thresh.formula = thresh.formula.1,
                data = newhealthsurvey,
                start = c(coef(mK)),
                fit.sigma = FALSE,
                decreasing.levels = TRUE)
  print('mM')
  mM  <- hopit(latent.formula = latent.formula.1,
               thresh.formula = thresh.formula.1,
               data = newhealthsurvey,
               start = c(coef(mH), logSigma = log(sigma(mH))),
               fit.sigma = TRUE,
               decreasing.levels = TRUE)
  print('mN')
  mN  <- hopit(latent.formula = latent.formula.3,
               thresh.formula = thresh.formula.5,
               data = newhealthsurvey,
               control = list(transform.latent = 'standardize',
                              transform.thresh = 'standardise'),
               decreasing.levels = TRUE)

  print('mO')
  mO  <<- hopit(latent.formula = latent.formula.H,
                thresh.formula = thresh.formula.1,
                data = newhealthsurvey,
                fit.sigma = FALSE,
                decreasing.levels = TRUE)

})

test_that("Hopit works2",{
  skip_on_cran()
  # test hopit object -------------------------
  print('test obj')
  test_hopit(object=mH, data = newhealthsurvey)
  test_hopit(object=mE, data = newhealthsurvey)
  test_hopit(object=m8, data = newhealthsurvey)
  test_hopit(object=mI, data = newhealthsurvey)

  # Anova -------------------------
  print('test anova')
  an1 <- anova(m0, mL)
  expect_equal(length(an1),5)
  expect_message(anova(mL, mH))
  expect_true(all(diff(c(length(coef(mO)),length(coef(mL)),length(coef(m0))))>0))
  an2 <- anova(mO, mL, m0, method='sequential', direction = 'increasing')
  expect_output(print(an2))
  expect_equal(dim(an2$table),c(2,3))
  an3 <- anova(mO, mL, m0, method='with.most.complex', direction = 'increasing')
  expect_equal(dim(an3$table),c(2,3))
  an4 <- anova(mO, mL, m0, method='with.least.complex', direction = 'increasing')
  expect_equal(dim(an4$table),c(2,3))

  # other errors --------------------------
  print('other')
  expect_error(anova(mO, mL, m0, method='with.least.complex', direction = 'decreasing'))

  expect_error(hopit(latent.formula = latent.formula.1,
                     thresh.formula = thresh.formula.1,
                     data = newhealthsurvey,
                     start = coef(mH),
                     fit.sigma = TRUE,
                     decreasing.levels = TRUE))

  expect_error(hopit(latent.formula = latent.formula.1,
                     thresh.formula = thresh.formula.1,
                     data = newhealthsurvey,
                     start = c(coef(mH), logSigma = log(sigma(mH))),
                     fit.sigma = FALSE,
                     decreasing.levels = TRUE))
})
