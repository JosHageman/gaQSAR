#' Plot method for gaQSAR objects
#'
#' Plot diagnostics for a `gaQSAR` result.
#'
#' @param x An object of class "gaQSAR" (returned from `gaVariableSelection()`).
#' @param type Character scalar specifying which plot to create or print.
#'   Supported values are `"all"` (default), `"fitness"`, `"predObs"` and `"williams"`.
#' @param ... Additional arguments (currently unused).
#'
#' @details `type = "all"` prints available plots separately, one after another.
#'  The fitness plot is always created with [createBestFitnessPlot()]. The 
#' observed-versus-predicted plot is created from `x$yTrain` and, when available,
#' `x$yExt`. The Williams plot use stored plot objects created by [createWilliamsPlot()].
#'
#' @return Invisibly returns the created or printed `ggplot2` object for a
#'   single plot type, or a named list of plot objects for `type = "all"`.
#'
#' @seealso [gaVariableSelection()], [createWilliamsPlot()]
#'
#' @export
plot.gaQSAR <- function(x, type = "all", ...) {

  if (!inherits(x, "gaQSAR")) {
    stop("`x` must be a gaQSAR object.", call. = FALSE)
  }

  type <- match.arg(type, c("all", "fitness", "predObs", "williams"))

  if (type == "fitness") {
    p <- createBestFitnessPlot(x)
    print(p)
    return(invisible(p))
  }

  if (type == "predObs") {
    p <- createPredObsPlot(x)
    print(p)
    return(invisible(p))
  }

  if (type == "williams") {
    p <- getStoredPlot(
      x,
      plotName = "williamsPlot",
      missingMessage = "Williams plot is not available. Run createWilliamsPlot() first."
    )
    print(p)
    return(invisible(p))
  }

  plotsList <- list()

  plotsList$fitness <- createBestFitnessPlot(x)
  print(plotsList$fitness)

  if (hasUsablePredictionData(x)) {
    plotsList$predObs <- createPredObsPlot(x)
    print(plotsList$predObs)
  } else {
    message(
      "Observed-versus-predicted plot is not available because neither ",
      "`x$yTrain` nor `x$yExt` contains usable prediction data."
    )
  }

  if (inherits(x$williamsPlot, "ggplot")) {
    plotsList$williams <- x$williamsPlot
    print(plotsList$williams)
  } else {
    message("Williams plot is not available. Run createWilliamsPlot() first.")
  }

  invisible(plotsList)
}


#' Create observed-versus-predicted plot
#'
#' @keywords internal
#' @noRd
createPredObsPlot <- function(x) {
  plotData <- getPredObsData(x)

  ggplot2::ggplot(
    plotData,
    ggplot2::aes(x = observed, y = predicted, color = dataset, shape = dataset)
  ) +
    ggplot2::geom_point(size = 2.5, alpha = 0.8) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      color = "red",
      linetype = "dashed",
      linewidth = 0.7
    ) +
    ggplot2::labs(
      title = "Observed vs Predicted",
      x = "Observed",
      y = "Predicted",
      color = "Dataset",
      shape = "Dataset"
    ) +
    ggplot2::theme_minimal(base_size = 14)
}


#' Get observed-versus-predicted data
#'
#' @keywords internal
#' @noRd
getPredObsData <- function(x) {
  dataList <- list()

  if (hasPredictionColumns(x$yTrain, observedName = "y", predictedName = "yhat")) {
    dataList$training <- data.frame(
      observed = x$yTrain$y,
      predicted = x$yTrain$yhat,
      dataset = "Training",
      stringsAsFactors = FALSE
    )
  }

  if (hasPredictionColumns(x$yExt, observedName = "y", predictedName = "yHat")) {
    dataList$external <- data.frame(
      observed = x$yExt$y,
      predicted = x$yExt$yHat,
      dataset = "External",
      stringsAsFactors = FALSE
    )
  }

  if (length(dataList) == 0L) {
    stop(
      "Observed-versus-predicted plot requires `x$yTrain` with columns `y` ",
      "and `yhat`, or `x$yExt` with columns `y` and `yHat`.",
      call. = FALSE
    )
  }

  plotData <- do.call(rbind, dataList)
  plotData <- plotData[stats::complete.cases(plotData[, c("observed", "predicted")]), , drop = FALSE]

  if (nrow(plotData) == 0L) {
    stop(
      "Observed-versus-predicted plot requires at least one complete observed ",
      "and predicted value pair.",
      call. = FALSE
    )
  }

  plotData$dataset <- factor(plotData$dataset, levels = c("Training", "External"))
  plotData
}


#' Check for usable prediction data
#'
#' @keywords internal
#' @noRd
hasUsablePredictionData <- function(x) {
  hasPredictionColumns(x$yTrain, observedName = "y", predictedName = "yhat") ||
    hasPredictionColumns(x$yExt, observedName = "y", predictedName = "yHat")
}


#' Check prediction columns
#'
#' @keywords internal
#' @noRd
hasPredictionColumns <- function(data, observedName, predictedName) {
  is.data.frame(data) &&
    all(c(observedName, predictedName) %in% names(data)) &&
    any(stats::complete.cases(data[, c(observedName, predictedName), drop = FALSE]))
}


#' Get stored ggplot object
#'
#' @keywords internal
#' @noRd
getStoredPlot <- function(x, plotName, missingMessage) {
  if (!inherits(x[[plotName]], "ggplot")) {
    stop(missingMessage, call. = FALSE)
  }

  x[[plotName]]
}
