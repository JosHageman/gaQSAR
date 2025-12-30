#' Integer-valued GA mutation operator
#'
#' Mutate an integer-encoded chromosome by independently replacing each gene,
#' with probability `object@pmutation`, by a new value drawn uniformly from
#' `1:object@upper[i]`. Intended for use as the `mutation` function in `GA::ga()`.
#'
#' @param object A `GA` object from the `GA` package. The operator reads
#'   `object@upper` (gene-wise upper bounds), `object@pmutation` (per-gene
#'   mutation probability), and `object@population` (matrix of chromosomes).
#' @param parent Integer; row index of the parent chromosome in `object@population`.
#' @param ... Ignored; for compatibility with the `GA` package interface.
#'
#' @details Mutation decisions are drawn independently per gene using
#' `stats::runif()`. When mutation occurs for gene `i`, a new integer in the
#' range `[1, upper[i]]` is sampled with replacement via `sample.int()`.
#'
#' @return A numeric vector containing the mutated chromosome.
#'
#' @seealso GA::ga
#'
#' @export
gaintegerMutation <- function(object, parent, ...) {
  ups <- object@upper
  mutated <- as.vector(object@population[parent, ])
  numGenes <- length(mutated)

  doMutate <- stats::runif(numGenes) < object@pmutation
  for (i in seq_len(numGenes)) {
    if (isTRUE(doMutate[i])) {
      mutated[i] <- sample.int(ups[i], size = 1L, replace = TRUE)
    }
  }

  mutated
}
