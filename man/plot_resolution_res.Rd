% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_generics.R
\name{plot_resolution_res}
\alias{plot_resolution_res}
\title{Plot the resolution results.}
\usage{
plot_resolution_res(object, print_head = TRUE, ...)
}
\arguments{
\item{object}{The class, either \code{bixverse::rbh_graph} or \code{bixverse::bulk_coexp}.}

\item{print_head}{Boolean. Print the Top5 resolution parameters and their
meta data. Only applicable for \code{bulk_coexp} objects.}

\item{...}{Additional arguments to parse to the functions.}
}
\value{
Plots the result, if the results were found in the class. Otherwise,
throws a warning and returns NULL.
}
\description{
Plots the resolution results (if they can be found in the class). The x-axis
reflects the different resolutions and the y axis the modularity observed
with that resolution.
}
