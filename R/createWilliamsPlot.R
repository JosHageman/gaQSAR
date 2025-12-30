#' Create a Williams plot for QSAR model diagnostics
#'
#' Generate a Williams plot (standardized residuals versus leverage) for each
#' selected model in a GA variable selection result. Training leverages and
#' standardized residuals are computed from an OLS fit on `xtrain`/`ytrain`.
#' External leverages are obtained from the test design matrix using the
#' training inverse information, and external standardized residuals are derived
#' using the training residual standard deviation.
#'
#' @param output A list-like object whose elements describe selected models.
#'   Each element must contain `importantPredictors` (names or indices of the
#'   predictors to use). The function augments each element with `williamsData`
#'   (a data.frame) and `williamsPlot` (a `ggplot2` object).
#' @param xtrain A data.frame or matrix of training predictors. Row names are used
#'   as molecule identifiers; columns must cover `importantPredictors`.
#' @param ytrain Numeric vector of training responses aligned with `xtrain` rows.
#' @param xtest A data.frame or matrix of test predictors with the same columns as
#'   `xtrain` (or a superset). Row names are used as molecule identifiers.
#' @param ytest Numeric vector of test responses aligned with `xtest` rows.
#' @param label Optional character string appended to the plot title.
#' @param residualThreshold Numeric threshold for standardized residuals.
#'
#' @details The plot displays training and external points with distinct colors.
#' Horizontal reference lines are drawn at 0 and ±3 standardized residuals. The
#' vertical leverage threshold is set to $h^* = 3p/n$, where $p$ is the number of
#' model coefficients (including the intercept) and $n$ is the number of training
#' observations.
#'
#' @return The input `output` object, with each element augmented by:
#'   - `williamsData`: data.frame with columns `leverages`, `residuals`, and `type`.
#'   - `williamsPlot`: a `ggplot2` object containing the Williams plot.
#'
#' @seealso [createQ2Plot()], [Q2()], [predictOOBObjects()]
#'
#' @export
createWilliamsPlot <- function(output, xtrain, ytrain, xtest, ytest,
                               label = "", residualThreshold = 2.5) {

  for (i in seq_len(length(output))) {

    xtrainSubset <- xtrain[, output[[i]]$importantPredictors]
    xtestSubset <- xtest[, output[[i]]$importantPredictors]

    #fit on training data
    myModel <- stats::lm(ytrain ~ ., data = data.frame(ytrain, xtrainSubset))

    #leverages training set
    leverages_train <- stats::hatvalues(myModel)

    #standardized residuals (training)
    residuals_train <- stats::rstandard(myModel)

    Xtrain <- stats::model.matrix(myModel)
    Xtest  <- stats::model.matrix(~ ., data = as.data.frame(xtestSubset))

    # Leverages test set
    # Calculate (X'X)^(-1) based on training data
    XtX_inv <- MASS::ginv(t(Xtrain) %*% Xtrain)

    # Leverages for test set
    leverages_test <- diag(Xtest %*% XtX_inv %*% t(Xtest))

    # External residuals for test set
    ytest_pred <- stats::predict(myModel, newdata = data.frame(xtestSubset))
    residuals_test_raw <- ytest - ytest_pred

    # Externally standardized test residuals
    sigma_train <- summary(myModel)$sigma
    residuals_test <- residuals_test_raw / (sigma_train * sqrt(1 - leverages_test))


    plotTitle <- paste0(
      "Williams plot", label, "\n -",
      paste(colnames(xtrainSubset), collapse = "\n -")
    )

    williamsData <- data.frame(leverages=c(leverages_train, leverages_test),
                               residuals=c(residuals_train, residuals_test),
                               type=c(rep("Training", length(leverages_train)),
                                      rep("Validation", length(leverages_test)))   )

    leverageThreshold <- (3 * length(stats::coef(myModel))) / nrow(xtrainSubset)

    p1 <- ggplot2::ggplot(williamsData,
                          ggplot2::aes(x = leverages, y = residuals,
                                       color = type, label = rownames(williamsData)
      )) +

      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggrepel::geom_text_repel(size = 3.5, max.overlaps = Inf) +

      # Horizontal residualThreshold lines
      ggplot2::geom_hline( yintercept = c(-residualThreshold, 0, residualThreshold), linetype = c("dashed", "solid", "dashed"),
        color = c("red", "black", "red"), linewidth = c(0.6, 0.4, 0.6)
      ) +

      # Vertical leverage threshold line
      ggplot2::geom_vline( xintercept = leverageThreshold,
        linetype = "dashed", color = "blue", linewidth = 0.7 ) +

      ggplot2::scale_color_manual(
        values = c("Training" = "#1f77b4", "Validation" = "#ff7f0e") ) +

      ggplot2::labs(title = plotTitle, x = "Leverage (h)", y = "Standardized residuals",
        color = "Dataset") +
      ggplot2::theme_minimal(base_size = 14)

      # Store results
      output[[i]]$williamsData <- williamsData
      output[[i]]$williamsPlot <- p1
  }

  return(output)
}
