#' Artificially generated health survey data
#'
#' A dataset containing artificialy generated survey data
#'
#' @format A data frame with 10000 rows and 11 variables:
#' \describe{
#'   \item{ID}{personal indentification number}
#'   \item{health}{reported health, 5 levels}
#'   \item{diabetes}{has diabetes? Yes or No}
#'   \item{obese}{has obese? Yes or No}
#'   \item{IADL_problems}{problems in Instrumental Activities of Daily Living? Yes or No}
#'   \item{hypertenssion}{has hypertenssion? Yes or No}
#'   \item{high_cholesterol}{has high cholesterol? Yes or No}
#'   \item{respiratory_problems}{has respirator problems? Yes or No}
#'   \item{heart_atack_or_stroke}{had stroke or heart attack? Yes or No}
#'   \item{poor_mobility}{has poor mobility? Yes or No}
#'   \item{very_poor_grip}{Cannot perform grip strength? Yes or No}
#'   \item{depression}{has depression? Yes or No}
#'   \item{other_diseases}{has other diseases? Yes or No}
#'   \item{sex}{sex/gender: woman or man}
#'   \item{ageclass}{categorized age: [50,60), [60,70), [70,80), [80,120)}
#'   \item{education}{two levels of education: primary or lower and secondary or higher}
#'   \item{contHM}{a continous health measure}
#'   \item{country}{country: X, Y, or Z}
#'   \item{csw}{cross-sectional survey weights}
#'   \item{psu}{primary statistical unit}
#'   \item{ssu}{secondary statistical unit}
#' }
#'
#' @source Data was randomly generated using probabilities of occurence of particular
#' combiantion of diseases, coditions, sex, age, education, reported health, etc.
#' The structure of the data and some probabilities were inspired by WAVE1 SHARE database
#' (DOIs: 10.6103/SHARE.w1.600), see Börsch-Supan et al for methodological details (Börsch-Supan et al. 2013).
#'
#' None of the records represent any true record/individual of SHARE database
#'
#' @references Börsch-Supan A, Brandt M, Hunkler C, et al (2013) Data resource profile: The survey of health, ageing and retirement in europe (share). Int J Epidemiol 42:992–1001. doi: 10.1093/ije/dyt088
#' @name healthsurvey
"healthsurvey"
