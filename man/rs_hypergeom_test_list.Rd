% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extendr-wrappers.R
\name{rs_hypergeom_test_list}
\alias{rs_hypergeom_test_list}
\title{Run a hypergeometric test over a list of target genes}
\usage{
rs_hypergeom_test_list(target_genes_list, gene_sets, gene_universe)
}
\arguments{
\item{target_genes_list}{A character vector representing the target gene set.}

\item{gene_sets}{A list of strings that represent the gene sets to test
against.}

\item{gene_universe}{A character vector representing the gene universe from
which the target genes and gene sets are sampled from.}
}
\value{
A list containing:
\itemize{
\item pvals - The p-values from the hypergeometric test
\item odds ratios - The calculated odds ratios
\item overlap - The size of the overlap
\item gene_set_lengths - The length of the gene sets.
}
}
\description{
Given a list of target gene sets, this function will test for
each of the individual target genes the hypergeoemetric enrichment against
the specified gene sets. WARNING! Incorrect use can cause kernel crashes.
Wrapper around the Rust functions with type checks are provided in the
package.
}
