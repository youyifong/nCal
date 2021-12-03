\name{ED5PL}
\alias{ED5PL}
\alias{FivePL.t}
\alias{FivePL.t.func}
\alias{FivePL.x}
\alias{FivePL.x.inv}
\alias{FivePL.x.inv.func}
\alias{FivePL.t.inv}
\alias{FivePL.t.inv.func}
\alias{FourPL.x.inv}
\alias{FourPL.x}
\alias{FourPL.t.func}
\alias{cla2gh}
\alias{gh2cla}
\alias{cla2ed50}
\alias{cla2ed50b}
\alias{ed502cla}
\alias{ed50b2cla}
\alias{get.curve.param.list}
\alias{simulate1curve}
\alias{vpl1.deriv}
\alias{vpl1.deriv.func}
\alias{vpl2.deriv}
\alias{vpl2.deriv.func}
\alias{vpl3.deriv}
\alias{vpl3.deriv.func}

\title{
Functions Related to 5PL Function 
}
\description{
5PL functions.
}
\usage{
FivePL.t (t,param)
FivePL.t.func (param) 
FivePL.x (x,param)
FivePL.x.inv (y,param)
FivePL.x.inv.func (param)
FivePL.t.inv (y,param) 
FivePL.t.inv.func (param) 
FourPL.x.inv (y, param)  
FourPL.x (x, param)
FourPL.t.func (param)
cla2gh (param)
gh2cla (param)
cla2ed50 (param)
cla2ed50b (param)
ed502cla (param)
ed50b2cla (param)
get.curve.param.list (param)
simulate1curve (param, t, sd.e=0.1, expy=TRUE, gamma=0) 
vpl1.deriv (x,param)
vpl1.deriv.func (param)
vpl2.deriv (x,param)
vpl2.deriv.func (param)
vpl3.deriv (x,param)
vpl3.deriv.func (param)

ED5PL (param, tao)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{vector of numbers. Parameters of the 5PL curve.}
  \item{t}{numeric vector. Log concentrations.}
  \item{x}{numeric vector. Concentrations.}
  \item{y}{numeric vector. Response.}
  \item{sd.e}{sd.e}
  \item{tao}{vector of numbers or single number. Effective doses. Vectorized for this argument.}
  \item{expy}{Boolean. Controls whether to exponentiate y}
  \item{gamma}{power variance function parameter}
}
\details{
FivePL.t and other functions are vectorized for the x and the y arguments.

Four parameterizations are implemented. 
Classical: b,c,d,e,f
gh: c,d,f,g,h
ED50: c,d,f,b,tao(ED50)
ED50b: c,d,f,h,tao(ED50). This is called tao-h in Cumberland et al. (2014)


}
\author{ Youyi Fong \email{yfong@fhcrc.org}, Xuesong Yu, William N. Cumberland }
\examples{

FivePL.t(5:6, p.eotaxin[1,])
FivePL.t.func(p.eotaxin[1,])
FivePL.x.inv(c(4,5,11), p.eotaxin[1,])
FivePL.t.inv.func(p.eotaxin[1,])



}
