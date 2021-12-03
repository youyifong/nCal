\name{plot5PL}
\alias{plot5PL}
\alias{lines5PL}
\title{ Plot a 5PL function}
\description{ 
  \code{plot5PL} plots a 5PL function
}

\usage{
plot5PL (param, xlim, ylim=NULL, col=NULL, lty=NULL, lwd=1, plot.legend=FALSE,
    add=FALSE, legend=NULL, main=NULL, xlab=NULL, ylab=NULL, xaxt="s", 
    yaxis.log.scale=FALSE, expy=FALSE, logy=FALSE, ...) 
lines5PL (param, xlim, ...) 

}

\arguments{

  \item{param}{Can be a vector or a matrix. If matrix, each row represents one curve. First check for classical parameterization, then gh, then ED50 parameterization. e/b/f can be provided in log form.}
  \item{xlim}{range of t, i.e. log (concentration)}
  \item{ylim}{ylim}
  \item{col}{col}
  \item{lty}{lty}
  \item{lwd}{lwd}
  \item{plot.legend}{Boolean. Indicates whether to plot legend.}
  \item{add}{Boolean. Indicates whether to add to an existing plot or create a new plot.}
  \item{legend}{vector. Legend.}
  \item{main}{string. Title.}
  \item{xlab}{xlab}
  \item{ylab}{ylab}
  \item{xaxt}{xaxt}
  \item{yaxis.log.scale}{Boolean. Controls whether to draw y axis on the log scale}
  \item{expy}{Boolean. Controls whether to exponentiate 5PL function values}
  \item{logy}{Boolean. Controls whether to log 5PL function values}
  \item{...}{additional arguments}
}

\author{ Youyi Fong \email{yfong@fhcrc.org} }
\details{ 
x axis is always drawn in the log scale.

}

\keyword{ distribution }
