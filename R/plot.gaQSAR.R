#' Plot method for gaQSAR objects
#'
#' Plot all ggplot objects contained in a gaQSAR result. This includes the Q2
#' plot (if created) and Williams plots (if created) for each model size.
#'
#' @param x An object of class "gaQSAR" (returned from `gaVariableSelection()`).
#' @param ... Additional arguments (currently unused).
#'
#' @details Extracts all ggplot2 objects from the result list and prints them
#' sequentially. The function searches for objects ending in "Plot" (e.g., 
#' `q2Plot`, `williamsPlot`). Each plot is displayed separately.
#'
#' @return Invisibly returns a list of plot objects, or `NULL` if no plots are found.
#'
#' @seealso [gaVariableSelection()], [createQ2Plot()], [createWilliamsPlot()]
#'
#' @export
plot.gaQSAR <- function(x, ...) {

  if (!inherits(x, "gaQSAR")) {
    stop("`x` must be a gaQSAR object.", call. = FALSE)
  }

  model <- x

  # Collect all plot objects from the result
  plots_list <- list()

  # Extract all components that are ggplot objects
  for (name in names(model)) {
    if (grepl("Plot$", name)) {
      if (inherits(model[[name]], "ggplot")) {
        plots_list[[name]] <- model[[name]]
      }
    }
  }

  # If no plots found, inform the user
  if (length(plots_list) == 0) {
    message("No plots found in the gaQSAR object. ")
    message("Run createQ2Plot() and/or createWilliamsPlot() first.")
    return(invisible(NULL))
  }

  # Print each plot separately in sequence
  if (length(plots_list) > 0) {
    for (p in plots_list) {
      print(p)
    }
    return(invisible(plots_list))
  }
}
