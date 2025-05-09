% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_bulk.R
\name{run_limma_voom}
\alias{run_limma_voom}
\title{Wrapper for a Limma Voom analysis}
\usage{
run_limma_voom(
  meta_data,
  main_contrast,
  dge_list,
  co_variates = NULL,
  ...,
  .verbose = TRUE
)
}
\arguments{
\item{meta_data}{data.table. The meta information about the experiment in
which the contrast info (and potential co-variates) can be found.}

\item{main_contrast}{String. Which column contains the main groups you want
to test differential gene expression with the Limma-Voom workflow for.}

\item{dge_list}{DGEList, see \code{\link[edgeR:DGEList]{edgeR::DGEList()}}.}

\item{co_variates}{String or NULL. Optional co-variates you wish to consider
during model fitting.}

\item{...}{Additional parameters to forward to \code{\link[limma:ebayes]{limma::eBayes()}} or
\code{\link[limma:voom]{limma::voom()}}.}

\item{.verbose}{Boolean. Controls verbosity of the function.}
}
\value{
A data.table with all the DGE results from \code{\link[limma:toptable]{limma::topTable()}} for
the identified contrast pairs.
}
\description{
Wrapper function to run Limma Voom workflows.
}
