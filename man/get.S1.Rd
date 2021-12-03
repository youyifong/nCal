\name{get.abc}
\alias{get.S1}
\alias{get.S2}
\alias{get.abc}
\alias{get.abs.dev}
\title{
Compute a measure of distance between two curves. 
}
\description{
abc criterion is area between curves divivded by the width of t.range. 
S1 criterion between two curves is defined as the integrated squared distance divivded by the width of t.range.
S2 is relative bias divivded by the width of t.range.
}
\usage{
get.abc(p1, p2, t.range)
get.S1(p1, p2, t.range)
get.S2(p1, p2, t.range)
get.abs.dev(p1, p2, t.range, y.range)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
a vector of coefficients specify the first curve
}
  \item{p2}{
a vector of coefficients specify the second curve
}
  \item{t.range}{
a vector two real numbers. The range of log standard concentrations. The first one is the minimum t and the second one is the maximum t
}
  \item{y.range}{
%%     ~~Describe \code{y.range} here~~
}
}
\details{
p1 and p2 can be in either classical or gh parameterization.
}
\value{
a real number
}
\examples{

get.abc(p.eotaxin[1,], p.eotaxin[2,], t.range=log(c(0.51,1e4)))
get.S1(p.eotaxin[1,], p.eotaxin[2,], t.range=log(c(0.51,1e4)))
get.S2(p.eotaxin[1,], p.eotaxin[2,], t.range=log(c(0.51,1e4)))
get.abs.dev(p.eotaxin[1,], p.eotaxin[2,], t.range=log(c(0.51,1e4)), y.range=c(5,6))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
