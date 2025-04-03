#' @title Plot the resolution results.
#'
#' @description
#' Plots the resolution results (if they can be found in the class). The x-axis
#' reflects the different resolutions and the y axis the modularity observed
#' with that resolution.
#'
#' @param object The class, either `bixverse::rbh_graph` or `bixverse::bulk_coexp`.
#' @param print_head Boolean. Print the Top5 resolution parameters and their
#' meta data. Only applicable for `bulk_coexp` objects.
#' @param ... Additional arguments to parse to the functions.
#'
#' @return Plots the result, if the results were found in the class. Otherwise,
#' throws a warning and returns NULL.
#'
#' @export
plot_resolution_res <- S7::new_generic(
  name = 'plot_resolution_res',
  dispatch_args = 'object',
  fun = function(object, print_head = TRUE, ...) {
    S7::S7_dispatch()
  }
)
