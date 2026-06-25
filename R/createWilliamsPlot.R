#' Create Williams plots for QSAR model diagnostics
#'
#' Creates Williams plots for one `gaQSAR` object or for a list of `gaQSAR`
#' objects. A Williams plot shows standardized residuals against leverage and is
#' commonly used to inspect outliers and the applicability domain of a QSAR
#' model.
#'
#' For each model, the selected predictors are taken from `importantPredictors`.
#' The model is refitted on `xtrain` and `ytrain`. Training leverages are computed
#' from the fitted design matrix. Validation leverages are computed from the
#' validation design matrix using the inverse information matrix from the training
#' data.
#'
#' @param output A `gaQSAR` object or a list of `gaQSAR` objects. Each object must
#'   contain `importantPredictors`, either as predictor names or as column indices
#'   referring to `xtrain`.
#' @param xtrain A data.frame or matrix with training predictors. Row names are
#'   used as object identifiers when available.
#' @param ytrain Numeric vector with training responses, aligned with the rows of
#'   `xtrain`.
#' @param xtest A data.frame or matrix with validation predictors. It must contain
#'   the selected predictor columns. Row names are used as object identifiers when
#'   available.
#' @param ytest Numeric vector with validation responses, aligned with the rows of
#'   `xtest`.
#' @param label Optional character string appended to the plot title.
#' @param residualThreshold Numeric threshold for standardized residuals. The
#'   default is 2.5.
#' @param labelPoints Character scalar controlling point labels: `"outliers"`
#'   labels points with high standardized residuals and/or high leverage,
#'   `"all"` labels all points, and `"none"` suppresses point labels. The default
#'   is `"outliers"`.
#'
#' @details Training standardized residuals are computed as
#'   `residual / (sigma * sqrt(1 - h))`. Validation residuals are standardized as
#'   prediction residuals using `residual / (sigma * sqrt(1 + h))`, where `h` is
#'   the leverage of the validation object relative to the training design. This
#'   avoids undefined residuals when validation leverages are larger than one.
#'
#'   Horizontal reference lines are drawn at 0 and at
#'   `+/- residualThreshold`. The leverage threshold is `3 * p / n`, where `p` is
#'   the number of model coefficients including the intercept and `n` is the
#'   number of training observations.
#'
#'   Points are treated as Williams plot outliers when they have an absolute
#'   standardized residual larger than `residualThreshold` and/or a leverage
#'   larger than the leverage threshold.
#'
#' @return The input object, with each `gaQSAR` object augmented by:
#'     - `williamsData`: data.frame with object identifiers, leverages, residuals and dataset type.
#'     - `williamsPlot`: a `ggplot2` Williams plot.
#'     - `williamsOutliers`: counts of high-residual and high-leverage objects
#'
#' @seealso [Q2()], [predictOOBObjects()]
#'
#' @export
createWilliamsPlot <- function(output, xtrain, ytrain, xtest, ytest,
                               label = "", residualThreshold = 2.5,
                               labelPoints = c("outliers", "all", "none")) {

  singleObject <- inherits(output, "gaQSAR")
  labelPoints <- match.arg(labelPoints)

  if (singleObject) {
    output <- list(output)
  }

  isGaQsars <- vapply(output, inherits, logical(1), what = "gaQSAR")

  if (!all(isGaQsars)) {
    stop("`output` must be a gaQSAR object or a list of gaQSAR objects.", call. = FALSE)
  }

  xtrain <- as.data.frame(xtrain, check.names = FALSE)
  xtest <- as.data.frame(xtest, check.names = FALSE)
  ytrain <- as.numeric(ytrain)
  ytest <- as.numeric(ytest)

  if (nrow(xtrain) != length(ytrain)) {
    stop("`ytrain` must have the same length as the number of rows in `xtrain`.", call. = FALSE)
  }

  if (nrow(xtest) != length(ytest)) {
    stop("`ytest` must have the same length as the number of rows in `xtest`.", call. = FALSE)
  }

  if (!is.numeric(residualThreshold) || length(residualThreshold) != 1L ||
      is.na(residualThreshold) || residualThreshold <= 0) {
    stop("`residualThreshold` must be a single positive number.", call. = FALSE)
  }

  if (anyNA(ytrain) || anyNA(ytest)) {
    stop("`ytrain` and `ytest` must not contain missing values.", call. = FALSE)
  }

  trainingObjectId <- rownames(xtrain)
  validationObjectId <- rownames(xtest)

  if (is.null(trainingObjectId)) {
    trainingObjectId <- as.character(seq_len(nrow(xtrain)))
  }

  if (is.null(validationObjectId)) {
    validationObjectId <- as.character(seq_len(nrow(xtest)))
  }

  for (i in seq_along(output)) {

    predictorSpec <- output[[i]]$importantPredictors

    if (is.null(predictorSpec) || length(predictorSpec) == 0L) {
      stop("Each `gaQSAR` object must contain non-empty `importantPredictors`.", call. = FALSE)
    }

    if (is.numeric(predictorSpec)) {
      if (anyNA(predictorSpec) || any(predictorSpec < 1) || any(predictorSpec > ncol(xtrain))) {
        stop("`importantPredictors` contains invalid indices for `xtrain`.", call. = FALSE)
      }
      predictorNames <- colnames(xtrain)[as.integer(predictorSpec)]
    } else {
      predictorNames <- as.character(predictorSpec)
    }

    idxTrain <- match(predictorNames, colnames(xtrain))

    if (anyNA(idxTrain)) {
      missingNames <- predictorNames[is.na(idxTrain)]
      stop(
        sprintf(
          "Missing predictor columns in `xtrain`: %s",
          paste(missingNames, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    idxTest <- match(predictorNames, colnames(xtest))

    if (anyNA(idxTest)) {
      missingNames <- predictorNames[is.na(idxTest)]
      stop(
        sprintf(
          "Missing predictor columns in `xtest`: %s",
          paste(missingNames, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    xtrainSubset <- xtrain[, idxTrain, drop = FALSE]
    xtestSubset <- xtest[, idxTest, drop = FALSE]

    if (anyNA(xtrainSubset) || anyNA(xtestSubset)) {
      stop("Selected predictor columns must not contain missing values.", call. = FALSE)
    }

    xtrainMatrix <- data.matrix(xtrainSubset)
    xtestMatrix <- data.matrix(xtestSubset)

    if (!all(is.finite(xtrainMatrix)) || !all(is.finite(xtestMatrix))) {
      stop("Selected predictor columns must be numeric and finite.", call. = FALSE)
    }

    xtrainDesign <- cbind(Intercept = 1, xtrainMatrix)
    xtestDesign <- cbind(Intercept = 1, xtestMatrix)

    fit <- stats::lm.fit(x = xtrainDesign, y = ytrain)
    coefficients <- fit$coefficients

    if (anyNA(coefficients)) {
      stop(
        "The selected predictors produce a rank-deficient linear model.",
        call. = FALSE
      )
    }

    nTrain <- nrow(xtrainDesign)
    numberOfCoefficients <- ncol(xtrainDesign)
    residualDf <- nTrain - numberOfCoefficients

    if (residualDf <= 0) {
      stop(
        "There are not enough training observations for the selected number of predictors.",
        call. = FALSE
      )
    }

    fittedTrain <- as.vector(xtrainDesign %*% coefficients)
    predictedTest <- as.vector(xtestDesign %*% coefficients)

    residualsTrainRaw <- ytrain - fittedTrain
    residualsTestRaw <- ytest - predictedTest

    sigmaTrain <- sqrt(sum(residualsTrainRaw^2) / residualDf)

    if (!is.finite(sigmaTrain) || sigmaTrain <= 0) {
      stop("The fitted model has zero or non-finite residual standard deviation.", call. = FALSE)
    }

    qrTrain <- qr(xtrainDesign)
    qTrain <- qr.Q(qrTrain)
    leveragesTrain <- rowSums(qTrain^2)

    xtxInverse <- chol2inv(qr.R(qrTrain))
    leveragesTest <- rowSums((xtestDesign %*% xtxInverse) * xtestDesign)

    residualsTrain <- residualsTrainRaw / (sigmaTrain * sqrt(pmax(1 - leveragesTrain, 0)))
    residualsTest <- residualsTestRaw / (sigmaTrain * sqrt(1 + leveragesTest))

    objectId <- c(trainingObjectId, validationObjectId)

    williamsData <- data.frame(
      objectId = objectId,
      leverages = c(leveragesTrain, leveragesTest),
      residuals = c(residualsTrain, residualsTest),
      type = c(
        rep("Training", length(leveragesTrain)),
        rep("Validation", length(leveragesTest))
      ),
      stringsAsFactors = FALSE
    )

    leverageThreshold <- (3 * numberOfCoefficients) / nTrain

    trainingMask <- williamsData$type == "Training"
    validationMask <- williamsData$type == "Validation"
    williamsData$highResidual <- abs(williamsData$residuals) > residualThreshold
    williamsData$highLeverage <- williamsData$leverages > leverageThreshold

    trainingHighResiduals <- sum(
      trainingMask & williamsData$highResidual,
      na.rm = TRUE
    )
    trainingHighLeverage <- sum(
      trainingMask & williamsData$highLeverage,
      na.rm = TRUE
    )
    validationHighResiduals <- sum(
      validationMask & williamsData$highResidual,
      na.rm = TRUE
    )
    validationHighLeverage <- sum(
      validationMask & williamsData$highLeverage,
      na.rm = TRUE
    )

    labelData <- switch(
      labelPoints,
      all = williamsData,
      outliers = williamsData[
        williamsData$highResidual | williamsData$highLeverage,
        , drop = FALSE
      ],
      none = williamsData[0, , drop = FALSE]
    )

    plotTitle <- paste0(
      "Williams plot",
      if (nzchar(label)) paste0(" - ", label) else "",
      "\n - ",
      paste(predictorNames, collapse = "\n - ")
    )

    williamsPlot <- ggplot2::ggplot(
      williamsData,
      ggplot2::aes(
        x = leverages,
        y = residuals,
        color = type
      )
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::geom_hline(
        yintercept = c(-residualThreshold, residualThreshold),
        linetype = "dashed",
        linewidth = 0.6
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "solid",
        linewidth = 0.4
      ) +
      ggplot2::geom_vline(
        xintercept = leverageThreshold,
        linetype = "dashed",
        linewidth = 0.7
      ) +
      ggplot2::scale_color_manual(
        values = c("Training" = "#1f77b4", "Validation" = "#ff7f0e")
      ) +
      ggplot2::labs(
        title = plotTitle,
        x = "Leverage (h)",
        y = "Standardized residuals",
        color = "Dataset"
      ) +
      ggplot2::theme_minimal(base_size = 14)

    if (nrow(labelData) > 0L) {
      williamsPlot <- williamsPlot +
        ggrepel::geom_text_repel(
          data = labelData,
          ggplot2::aes(label = objectId),
          size = 3.5,
          max.overlaps = Inf
        )
    }

    output[[i]]$williamsData <- williamsData
    output[[i]]$williamsPlot <- williamsPlot
    output[[i]]$williamsOutliers <- list(
      training = list(
        high_residuals = trainingHighResiduals,
        high_leverage = trainingHighLeverage
      ),
      validation = list(
        high_residuals = validationHighResiduals,
        high_leverage = validationHighLeverage
      )
    )
  }

  if (singleObject) {
    return(output[[1]])
  }

  output
}
