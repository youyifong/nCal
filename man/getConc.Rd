\name{getConc}
\alias{getConc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Concentration Estimation
}
\description{
Estimate analyte concentrations based on observed outcome and a fitted curve.
}
\usage{
getConc(fit, y, n.replicate = 1, check.out.of.range = 1, x.range = NULL, 
    y.range = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
drc object or bcrm object}
  \item{y}{
numeric vector. Observed outcome in the samples of interest. Each element corresponds to one sample}
  \item{n.replicate}{
integer. Number of replicates that are averaged to generate y}
  \item{check.out.of.range}{
integer.}
  \item{x.range}{
vector of 2 numbers. The minimal and maximal concentration of the standard samples
}
  \item{y.range}{
vector of 2 numbers. The minimal and maximal observed response of the standard samples
}
  \item{verbose}{
Boolean.
}
}
\details{
Vectorized for y.
}
