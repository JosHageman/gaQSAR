#' Compute Q2 (cross-validated R-squared)
#'
#' Compute the $Q^2$ statistic from observed and predicted values:
#' \deqn{Q^2 = 1 - \frac{\sum (y_{obs} - y_{pred})^2}{\sum (y_{obs} - \bar{y}_{obs})^2}.}
#' This function does not perform cross-validation; it only evaluates the formula
#' for supplied predictions (e.g., from LOOCV or an external test set).
#'
#' @param y Numeric vector of observed response values.
#' @param yhat Numeric vector of predicted response values, aligned with `y`.
#'
#' @details The vectors `y` and `yhat` must have the same length. The denominator
#' is the total sum of squares of `y`; if `y` is constant this denominator is 0,
#' making $Q^2$ undefined (the result will be `NaN` or `Inf`). Missing values are
#' not handled specially; if present, the result may be `NA`.
#'
#' @return A numeric scalar $Q^2$ value. Values can be negative when predictive
#' performance is poor and approach 1 for perfect predictions.
#'
#' @seealso [singleCV()], [createQ2Plot()], [predictOOBObjects()]
#'
#' @export
Q2 <- function(y, yhat) {

  if (length(y) != length(yhat)) {
    stop("Observed and predicted vectors must have the same length.", call. = FALSE)
  }

  return(1 - (sum((y - yhat)^2) / sum((y - mean(y))^2)))
}
