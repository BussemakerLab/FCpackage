% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FCfunctions.R
\name{SRpred}
\alias{SRpred}
\title{Similarity Regression Prediction}
\usage{
SRpred(mono_motifs, Alignment, pos, weightfile)
}
\arguments{
\item{mono_motifs}{Motif samples to incldue in the model, result from loadMono_motifs() or filterMotifList().}

\item{Alignment}{The alignment resulted from concatAli, ideally after matchAliMotif().}

\item{pos}{positions to show in the PSAM}

\item{weightfile}{json file that contains pretrained weights from Lamber et. al.,}
}
\value{
predTrue data frame
}
\description{
Predict motif model with Similarity regression
}
