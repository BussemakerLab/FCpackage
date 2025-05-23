% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FCfunctions.R
\name{scoreMatrixSeed}
\alias{scoreMatrixSeed}
\title{Score Matrix Seed}
\usage{
scoreMatrixSeed(matrix, seed, weight = c())
}
\arguments{
\item{matrix}{motif matrix in frequency measurements}

\item{seed}{seed matrix describing binding motif.}

\item{weight}{numeric vector of weights of each motif position of the seed}
}
\value{
Minimum distance score
}
\description{
Score a motif matrix with a seed, the matrix and seed needs to be in the same length
}
