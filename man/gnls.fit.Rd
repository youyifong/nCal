\name{gnls.fit}
\alias{gnls.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit with gnls function
}
\description{
Fit with gnls function
}
\usage{
gnls.fit(formula, data, fit.4pl = FALSE, startVal = NULL, varFun = nlme::varPower(),
    verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{fit.4pl}{
%%     ~~Describe \code{fit.4pl} here~~
}
  \item{startVal}{
%%     ~~Describe \code{startVal} here~~
}
  \item{varFun}{
%%     ~~Describe \code{varFun} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
When startVal is given, varPower(value) does not seem to change numerical values of the fit
}
