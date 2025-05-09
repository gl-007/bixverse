% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_bulk.R
\name{.check_pca_grp_differences}
\alias{.check_pca_grp_differences}
\title{Helper function to check if PC1/2 distinguish groups}
\usage{
.check_pca_grp_differences(pc1, pc2, grps)
}
\arguments{
\item{pc1}{Numeric vector. The values of PC1.}

\item{pc2}{Numeric vector. The values of PC2.}

\item{grps}{Factor or character vector. The group vector.}
}
\value{
A data.table with two columns
}
\description{
Runs an ANOVA for the group variable vs PC1 and PC2 and checks if either
principal component can separate the groups.
}
