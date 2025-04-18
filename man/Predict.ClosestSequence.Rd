% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FCfunctions.R
\name{Predict.ClosestSequence}
\alias{Predict.ClosestSequence}
\title{Predict.ClosestSequence}
\usage{
Predict.ClosestSequence(mono_motifs, Alignment, testSeq)
}
\arguments{
\item{mono_motifs}{List of motifs that are input as frequency values (0-1 with highest of 1)}

\item{Alignment}{The alignment resulted from concatAli}

\item{testSeq}{The sequence to be predicted}
}
\value{
motif model of closest sequence prediction
}
\description{
Predicting binding motif with ClosesSeqeunce method
}
