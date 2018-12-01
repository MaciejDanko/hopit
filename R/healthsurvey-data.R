#' Artificially generated health survey data
#'
#' A dataset containing artificial survey data
#'
#' @format A data frame with 10000 rows and 10 variables:
#' \describe{
#'   \item{ID}{Personal number}
#'   \item{health}{reported health, 5 levels}
#'   \item{hypertenssion}{has hypertenssion? Yes or No}
#'   \item{cholesterol}{has high cholesterol? Yes or No}
#'   \item{heart_stroke_respiratory}{has respiratory problems or had stroke or heart attack? Yes or No}
#'   \item{poor_mobility_or_grip}{has poor mobility or cannot perform grip strength? Yes or No}
#'   \item{depression}{has depression? Yes or No}
#'   \item{sex}{sex: woman or man}
#'   \item{ageclass}{categorized age: [50,60), [60,70), [70,80), [80,120)}
#' }
#'
#' @source Data was generated using probabilities of occurence of particular
#' combiantion of diseases, coditions, sex, age, and reported health of SHARE database.
#' @name healthsurvey
"healthsurvey"
