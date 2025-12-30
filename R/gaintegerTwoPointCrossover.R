#' Integer-valued two-point crossover operator
#'
#' Perform a two-point crossover on two parent chromosomes. Two crossover points
#' are sampled uniformly, ordered as `p1 < p2`, and the outer segments
#' (positions `1..p1` and `p2..end`) are swapped between parents. Intended for
#' use as the `crossover` function in `GA::ga()`.
#'
#' @param object A `GA` object from the `GA` package; uses `object@population`
#'   as the matrix of chromosomes.
#' @param parents Integer vector of length 2 giving the row indices of the two
#'   parents in `object@population`.
#' @param ... Ignored; included for compatibility with the `GA` package interface.
#'
#' @details Returns two children with the same integer encoding and length as
#' the parents. Fitness values are not computed here and are returned as
#' `NA_real_` placeholders.
#'
#' @return A list with elements `children` (2 x nGenes integer matrix) and
#'   `fitness` (numeric vector of length 2 filled with `NA`).
#'
#' @seealso GA::ga
#'
#' @export
gaintegerTwoPointCrossover <- function(object, parents, ...) {
  stopifnot(length(parents) == 2L)

  parentMatrix <- object@population[parents, , drop = FALSE]
  nGenes <- ncol(parentMatrix)

  children <- parentMatrix

  if (nGenes > 1L) {
    crossoverPoints <- sort(sample.int(nGenes, size = 2L))
    p1 <- crossoverPoints[1]
    p2 <- crossoverPoints[2]

    swapIndex <- c(seq_len(p1), seq.int(p2, nGenes))

    children[1, swapIndex] <- parentMatrix[2, swapIndex]
    children[2, swapIndex] <- parentMatrix[1, swapIndex]
  }

  list(children = children, fitness = rep(NA_real_, 2L))
}
