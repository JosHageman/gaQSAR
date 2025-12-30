#' Predict out-of-bag objects and compute external Q2
#'
#' Compute external predictions for test objects using stored coefficient
#' vectors (per selected model) and derive the external $Q^2$ metric. For each
#' element in `output`, predictor columns in `x` are matched by name to the
#' coefficient names (excluding the intercept), an intercept term is prepended,
#' and predictions are obtained via a matrix product. The resulting external
#' $Q^2$ and per-object predictions/residuals are attached to each element.
#'
#' @param output A list-like GA selection result where each element contains a
#'   numeric coefficient vector `model` for the selected predictors, with the
#'   intercept in the first position and remaining names matching columns in `x`.
#' @param x A data.frame or matrix with predictors in columns and observations in
#'   rows. Column names must cover the predictor names used in `model` (excluding
#'   the intercept). Row order must align with `y`.
#' @param y Numeric vector of observed response values aligned with the rows of `x`.
#' @param verbose Logical; if `TRUE`, prints the computed external $Q^2$ for each
#'   model element.
#'
#' @details Predictor order does not matter: columns are matched by name to the
#' coefficient vector (excluding the intercept). An intercept column of ones is
#' automatically added before prediction. Missing matches in `x` will yield `NA`
#' indices and may cause subsetting errors.
#'
#' @return The input `output` object, with each element augmented by:
#'   - `Q2Ext`: numeric external $Q^2$ computed by `Q2(y, yHat)`.
#'   - `yExt`: data.frame with columns `y`, `yHat`, and `residual` for test data.
#'
#' @seealso [Q2()], [createQ2Plot()], [createWilliamsPlot()]
#'
#' @export
predictOOBObjects <- function(output, x, y, verbose=FALSE) {

  #find corresponding predicor columns
  for (i in seq_along(output)) {
    idx <- match(names(output[[i]]$model)[-1], colnames(x))
    #get predicted values
    yHat <- output[[i]]$model %*% t(cbind(1, x[, idx, drop = FALSE]))
    Q2Ext <- Q2(y, yHat)

    if (verbose) message(sprintf("External Q2: %.4f", Q2Ext))

    output[[i]]$Q2Ext <- Q2Ext
    output[[i]]$yExt <- data.frame(y = y, yHat = t(yHat), residual = (y - t(yHat)))
  }

  #return new object with Q2 and predictions for ext validation set
  return(output)
}
