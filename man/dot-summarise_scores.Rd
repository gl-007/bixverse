% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_graphs.R
\name{.summarise_scores}
\alias{.summarise_scores}
\title{Summarise gene scores if they are duplicates.}
\usage{
.summarise_scores(x, summarisation = c("max", "mean", "harmonic_sum"))
}
\arguments{
\item{x}{Named numeric.}

\item{summarisation}{String. Which summary function to use.}
}
\value{
Named numeric.
}
\description{
Summarise gene scores if they are duplicates.
}
