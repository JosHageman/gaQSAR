#' Plot best fitness per generation
#'
#' Builds a simple line plot of the GA best fitness value per generation.
#' Works with both `gaQSAR` and `gaQSAR_dcv` objects. For nested CV objects
#' (gaQSAR_dcv), one series is shown per outer fold (each inner GA run).
#'
#' @param object A `gaQSAR` object or a `gaQSAR_dcv` object.
#' @param title Optional character plot title. If missing, a default title is
#'   chosen based on the object type.
#'
#' @return A `ggplot2` object.
#'
#' @export
createBestFitnessPlot <- function(object, title = NULL) {

  buildDf <- function(runLabel, bestFitness) {
    if (is.null(bestFitness) || length(bestFitness) == 0) return(NULL)
    data.frame(
      generation = seq_along(bestFitness),
      bestFitness = as.numeric(bestFitness),
      run = runLabel,
      stringsAsFactors = FALSE
    )
  }

  if (inherits(object, "gaQSAR")) {
    df <- buildDf("GA run", object$bestFitnessPerGeneration)
    plotTitle <- if (is.null(title)) "Best fitness per generation" else title
  } else if (inherits(object, "gaQSAR_dcv")) {
    fm <- object$foldModels
    if (is.null(fm) || length(fm) == 0) {
      stop("`object` does not contain any fold models to plot.", call. = FALSE)
    }

    dfList <- lapply(seq_along(fm), function(idx) {
      model <- fm[[idx]]
      if (is.null(model) || is.null(model$bestFitnessPerGeneration)) return(NULL)
      buildDf(paste0("Fold ", idx), model$bestFitnessPerGeneration)
    })
    df <- do.call(rbind, dfList)
    plotTitle <- if (is.null(title)) "Best fitness per generation (inner GA, all folds)" else title
  } else {
    stop("`object` must be a gaQSAR or gaQSAR_dcv object.", call. = FALSE)
  }

  if (is.null(df) || nrow(df) == 0) {
    stop("No best-fitness data available to plot.", call. = FALSE)
  }

  df <- df[is.finite(df$bestFitness), , drop = FALSE]
  if (nrow(df) == 0) {
    stop("Best-fitness data contains no finite values to plot.", call. = FALSE)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = generation, y = bestFitness, color = run)) +
    ggplot2::geom_line(alpha = 0.85) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::labs(
      title = plotTitle,
      x = "Generation",
      y = "Best fitness",
      color = "Run"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.4)
    )

  p
}
