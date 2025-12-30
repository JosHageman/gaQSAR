#' Split data into training and test sets
#'
#' Split observations into training and test sets using either a Kennard–Stone
#' split ("KS") based on sample representativeness or a simple random split.
#' The function returns row indices for both sets.
#'
#' @param data Data frame or matrix with observations in rows (predictors and/or
#'   response may be present; only row indices are returned).
#' @param method Character; one of `"KS"` (Kennard–Stone) or `"random"`.
#' @param trainPercentage Numeric in (0, 1); fraction of observations assigned to
#'   the training set.
#' @param pc Numeric; proportion of variance to retain in the PCA step used by
#'   `prospectr::kenStone()` when `method = "KS"`.
#' @param verbose Logical; if `TRUE`, prints the chosen training and test IDs (row names).
#'
#' @details For `method = "KS"`, the function requires the `prospectr` package
#' and delegates to `prospectr::kenStone()` with `k = round(trainPercentage * n)`
#' on a scaled version of `data`. For `method = "random"`, rows are permuted and
#' the first `k` indices are assigned to training.
#'
#' @return A list with integer vectors `model` (training indices) and `test` (test indices).
#'
#' @export
splitUp <- function(data, method = c("KS", "random"), trainPercentage = 0.80,
    pc = 0.95, verbose = FALSE) {

  method <- match.arg(method)

  stopifnot(is.numeric(trainPercentage), length(trainPercentage) == 1, is.finite(trainPercentage))
  stopifnot(trainPercentage > 0, trainPercentage < 1)
  stopifnot(is.numeric(pc), length(pc) == 1, is.finite(pc))
  stopifnot(is.logical(verbose), length(verbose) == 1)

  n <- nrow(data)
  if (is.null(n) || n < 2) {
    stop("`data` must have at least 2 rows.", call. = FALSE)
  }

  nTrain <- max(1, min(n - 1, as.integer(round(trainPercentage * n))))

  if (method == "KS") {
    if (!requireNamespace("prospectr", quietly = TRUE)) {
      stop("Package 'prospectr' is required for method = 'KS'.", call. = FALSE)
    }

    if (isTRUE(verbose)) {
      message(sprintf("KS split (%.2f, %.2f)", trainPercentage, pc))
    }

    # prospectr::kenStone returns a list with $model and $test indices.
    split <- prospectr::kenStone(data, k = nTrain, pc = pc, .scale = TRUE)

    # Drop any extra fields to keep the return object stable.
    split$pc <- NULL
  } else {
    if (isTRUE(verbose)) {
      message(sprintf("Random split (%.2f)", trainPercentage))
    }

    perm <- sample.int(n, size = n, replace = FALSE)
    split <- list(
      model = perm[seq_len(nTrain)],
      test = perm[-seq_len(nTrain)]
    )
  }

  if (isTRUE(verbose)) {
    ids <- rownames(data)
    if (is.null(ids)) {
      ids <- as.character(seq_len(n))
    }

    message(sprintf("Training set IDs: %s", paste(ids[split$model], collapse = ", ")))
    message(sprintf("Test set IDs: %s", paste(ids[split$test], collapse = ", ")))
  }

  split
}
