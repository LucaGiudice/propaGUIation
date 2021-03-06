% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_genes_annotation.R
\name{download_genes_annotation}
\alias{download_genes_annotation}
\title{This function downloads the genes annotations for human and mouse}
\usage{
download_genes_annotation(organism = "human", gene_names)
}
\arguments{
\item{organism}{,string representing the orgnanism for which downloading the
genes annotations}

\item{gene_names, }{vector of strings, gene names to search annotation for}
}
\value{
list with the gene id, chromosome, start and end positions.
}
\description{
This function downloads the genes annotations for human and mouse
}
