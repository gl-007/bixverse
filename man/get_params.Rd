% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_class.R
\name{get_params}
\alias{get_params}
\title{Get the parameters that were used.}
\usage{
get_params(object, to_json = FALSE, pretty_json = FALSE)
}
\arguments{
\item{object}{A class within bixverse that inherits from \code{\link[=bixverse_base_class]{bixverse_base_class()}}.}

\item{to_json}{Shall the params be returned as a JSON string.}

\item{pretty_json}{Shall the params be returned as a pretty JSON string.}
}
\value{
Depending on parameters either the R list or a (pretty) JSON string.
}
\description{
Extracts parameters from the \code{bixverse_base_class} class (or child
classes) and has options to return (pretty) JSONs.
}
