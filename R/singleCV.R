#' LOOCV Q2 fitness function for small datasets
#'
#' Compute a leave-one-out cross-validated $Q^2$ for a candidate set of
#' predictors. Designed for use as a fitness function inside a genetic
#' algorithm. Duplicate predictors are penalized with a large negative value,
#' and candidate sets with variance inflation factor (VIF) above a threshold are
#' rejected (fitness 0).
#'
#' @param predictors Integer vector of 1-based predictor indices referring to
#'   columns of `x`.
#' @param x Matrix or data.frame of predictors (rows = observations,
#'   columns = candidate variables).
#' @param y Numeric response vector aligned with rows of `x`.
#' @param vifThreshold Numeric; maximum allowed VIF (default 5).
#'
#' @details An intercept is added internally. VIF values are computed by
#' regressing each selected predictor on the remaining selected predictors; if
#' any VIF exceeds `vifThreshold`, the fitness is 0. LOOCV is performed via
#' explicit refits for each left-out observation.
#'
#' @return A numeric scalar: $Q^2$ if valid; 0 if the VIF constraint fails; or
#'   -100 if duplicate predictors are present.
#'
#' @export
singleCV <- function(predictors, x, y, vifThreshold = 5) {

  if (length(unique(predictors)) != length(predictors)) {
    return(-100)
  }

  x <- cbind(1, x[, predictors, drop = FALSE])

  # VIF check
  vif <- vapply(2:ncol(x), function(i) {
    xOther <- x[, -i, drop = FALSE]
    xTarget <- x[, i]

    qrObj <- qr(xOther)
    if (qrObj$rank < ncol(xOther)) return(Inf)

    beta <- qr.coef(qrObj, xTarget)
    resid <- xTarget - xOther %*% beta

    sse <- sum(resid * resid)
    tss <- sum((xTarget - mean(xTarget))^2)

    if (tss == 0) Inf else 1 / (sse / tss)
  }, numeric(1))

  if (any(vif > vifThreshold)) {
    return(0)
  }


  # LOOCV
  yHat <- vapply(seq_len(nrow(x)), function(i) {
    xTrain <- x[-i, , drop = FALSE]
    yTrain <- y[-i]

    qrObj <- qr(xTrain)
    if (qrObj$rank < ncol(xTrain)) {
      return(NA_real_)
    }

    beta <- qr.coef(qrObj, yTrain)
    sum(x[i, ] * beta)
  }, numeric(1))

  if (anyNA(yHat)) {
    return(0)
  }

  Q2(y, yHat)
}
