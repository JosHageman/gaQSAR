#' GA monitor function for QSAR variable selection
#'
#' Create a monitor callback for `GA::ga()` that prints concise progress lines
#' at fixed intervals: generation number, best fitness value, number of selected
#' variables, and their indices.
#'
#' @param interval Integer; print progress every `interval` generations.
#'
#' @details The returned function closes over `interval` and is intended for the
#' `monitor` argument of `GA::ga()`. It inspects the current GA object to obtain
#' the best individual and reports its fitness and selected variable indices.
#' Output is written to standard output via `cat()` and nothing is returned.
#'
#' @return A function suitable for the `monitor` argument of `GA::ga()`.
#'
#' @seealso GA::ga
#'
#' @export
QSARMonitorFactory <- function(interval = 50) {

  force(interval)

  function(object) {

    if (object@iter %% interval != 0) {
      return(invisible(NULL))
    }

    bestIdx <- which.max(object@fitness)
    bestFitness <- object@fitness[bestIdx]
    selectedVars <- sort(object@population[bestIdx, ])

    msg <- sprintf(
      "Iter %4d | Best Q2: %.3f | #vars: %d | vars: %s\n",
      object@iter,
      bestFitness,
      length(selectedVars),
      paste(selectedVars, collapse = " ")
    )

    cat(msg)

    return(invisible(NULL))
  }
}
