% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FCfunctions.R
\name{plotPredTrue}
\alias{plotPredTrue}
\title{Plot Prediction vs. True}
\usage{
plotPredTrue(
  predTrue,
  xlab = "Experimental -ΔΔG/RT",
  ylab = "Predicted -ΔΔG/RT",
  main = "prediction",
  xlim = c(1.05 * min(predTrue$true), 1.05 * max(predTrue$true)),
  ylim = c(1.05 * min(predTrue$pred), 1.05 * max(predTrue$pred)),
  throughZero = T
)
}
\arguments{
\item{predTrue}{data frame with 2 columns true and pred, result from getBaseLineAccuracy(), closestSeqPred(), or SVDregression.Iterative.CV()}

\item{xlab}{x axis label}

\item{ylab}{y axis label}

\item{main}{plot title}

\item{xlim}{x axis limits}

\item{ylim}{y axis limits}

\item{throughZero}{If true the regression line passes through zero}
}
\value{
Scatter plot showing prediction accuracy
}
\description{
Make scatter plot between predicted and true value of motif prediction at a position
}
