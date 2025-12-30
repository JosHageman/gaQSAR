#' Integer-valued one-point crossover operator
#'
#' Perform one-point crossover on two parent chromosomes. A single crossover
#' point is sampled uniformly over gene positions; genes to the left of the
#' point (inclusive) are swapped between parents. Intended for use as the
#' `crossover` function in `GA::ga()`.
#'
#' @param object A `GA` object from the `GA` package; uses `object@population` as
#'   the matrix of chromosomes.
#' @param parents Integer vector of length 2 giving the row indices of the two
#'   parents in `object@population`.
#' @param ... Ignored; included for compatibility with the `GA` package interface.
#'
#' @details The function returns two children with identical lengths and integer
#' encoding as the parents. The fitness values are not evaluated here and are
#' returned as `NA_real_` placeholders, to be computed by the GA engine.
#'
#' @return A list with elements `children` (2 x nGenes integer matrix) and
#'   `fitness` (numeric vector of length 2 filled with `NA`).
#'
#' @seealso GA::ga
#'
#' @export
gaintegerOnePointCrossover <- function(object, parents, ...) {
  stopifnot(length(parents) == 2L)

  parentMatrix <- object@population[parents, , drop = FALSE]
  nGenes <- ncol(parentMatrix)

  children <- parentMatrix

  if (nGenes >= 1L) {
    crossoverPoint <- sample.int(nGenes, size = 1L)
    swapIndex <- seq_len(crossoverPoint)

    children[1, swapIndex] <- parentMatrix[2, swapIndex]
    children[2, swapIndex] <- parentMatrix[1, swapIndex]
  }

  list(children = children, fitness = rep(NA_real_, 2L))
}
