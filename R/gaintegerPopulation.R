#' Integer-valued GA population initializer
#'
#' Create an initial population for integer-encoded chromosomes. For each gene
#' position `j`, values are drawn independently and uniformly from the range
#' `1:object@upper[j]`. Intended for use as the `population` function in
#' `GA::ga()`.
#'
#' @param object A `GA` object from the `GA` package; reads `object@upper`
#'   (gene-wise upper bounds) and `object@popSize` (number of individuals).
#' @param ... Ignored; present for compatibility with the `GA` package interface.
#'
#' @details The returned matrix has one row per individual and one column per
#' gene. Sampling uses `sample.int()` with replacement.
#'
#' @return A numeric matrix with dimensions `popSize x length(upper)` containing
#'   the initial population.
#'
#' @seealso GA::ga
#'
#' @export
gaintegerPopulation <- function(object, ...) {
  ups <- object@upper
  population <- vapply(ups, function(u) {
    sample.int(u, size = object@popSize, replace = TRUE)
  }, integer(object@popSize))

  population
}
