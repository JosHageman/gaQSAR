#' Plot Q2 versus number of predictors
#'
#' Create a diagnostic plot of cross-validated performance against model size.
#' The figure shows $Q^2$ for the training set (LOOCV) and, when available,
#' the external validation set as a function of the number of selected
#' predictors. Use this for visual inspection of GA-based variable selection
#' results and to compare model parsimony versus predictive ability.
#'
#' @param output A list-like object where each element corresponds to a model
#'   of a given size and contains named numeric entries `numVar`, `Q2Loocv`,
#'   and optionally `Q2Ext`. Typical inputs are elements from GA variable
#'   selection runs. If `Q2Ext` is absent, only the training curve is drawn.
#' @param label Character string used as the plot title (set to `""` for no title).
#'
#' @details
#' Lines are drawn only for groups that contain two or more points; otherwise
#' individual points are shown. Point shape and color encode whether values
#' originate from the training (`Q2Loocv`) or external (`Q2Ext`) set. The
#' x-axis shows the number of predictors; tick labels are limited to integers.
#'
#' @return Invisibly returns the created `ggplot2` object for further styling
#'   or saving.
#'
#' @seealso [gaVariableSelection()], [Q2()], [createWilliamsPlot()]
#'
#' @export
createQ2Plot <- function(output, label = "") {

  #plot Q2 curve with ggplot
  tmp <- vapply(seq_along(output), function(i) {
    c(numVar = output[[i]]$numVar, Q2Loocv = output[[i]]$Q2Loocv, Q2Ext = output[[i]]$Q2Ext)
  }, numeric(3))

  sdat2 <- reshape2::melt(as.data.frame(t(tmp)), id.vars = "numVar")
  sdat2_line <- subset(sdat2, stats::ave(numVar, variable, FUN = length) > 1)

  p0 <- ggplot2::ggplot(sdat2, ggplot2::aes(x = numVar, y = value)) +
    ggplot2::geom_line(
      data = sdat2_line,
      ggplot2::aes(col = variable, linetype = variable),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(pch = variable, col = variable),
      size = 4
    ) +
    ggplot2::labs(
      x = "Number of predictors",
      y = bquote("Q"^2),
      title = label
    ) +
    ggplot2::scale_shape_discrete(
      name = "Set", breaks = c("Q2Loocv", "Q2Ext"),
      labels = c(bquote("Q"^2 ~ " Training set"), bquote("Q"^2 ~ " Validation set"))
    ) +
    ggplot2::scale_colour_discrete(
      name = "Set", breaks = c("Q2Loocv", "Q2Ext"),
      labels = c(bquote("Q"^2 ~ " Training set"), bquote("Q"^2 ~ " Validation set"))
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(),
      labels = function(x) ifelse(x %% 1 == 0, x, "")
    )

  output[[1]]$q2Plot <- p0
  return(output)
}
