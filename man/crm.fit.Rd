\name{crm.fit}
\alias{crm.fit}
\alias{deviance.crm}
\alias{lines.crm}
\alias{coef.crm}
\alias{print.crm}
\title{
Fit Concentration Response Model
}
\description{
crm.fit can fit a constant or power variance function or log transform both sides.
}
\usage{
crm.fit (formula, data, fit.4pl=FALSE, var.model=c("constant","power"),
    robust="mean", method=c("gls-pl","gnls","mle"), max.iter=50, reltol=1e-3, 
    gof.threshold=0.2, log.both.sides=FALSE, verbose=FALSE) 


\method{deviance}{crm}(object, ...)
\method{print}{crm}(x, ..., digits=3)
\method{lines}{crm}(x, ...)
\method{coef}{crm}(object, parameterization=c("cla","gh","ed50b","edb50"), ...)

}

\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{fit.4pl}{
Boolean}
  \item{var.model}{
string}
  \item{robust}{
string}
  \item{method}{
string}
  \item{max.iter}{
number}
  \item{digits}{
number}
  \item{reltol}{
numeric}
  \item{gof.threshold}{
numeric}
  \item{verbose}{Boolean}
  \item{log.both.sides}{Boolean, log transform both sides}
  \item{object, x}{crm object}
  \item{parameterization}{string, output parameterization}
  \item{...}{additional argument}


}
\details{
crm.fit implements an iterative method for estimating a model with power variance.
method: gls-pl means GLS-PL (see reference)
log.both.sides: transform both sides (see reference)
}
\value{
An object of crm and drm type.
  \item{var.power}{estimated power parameter in the power variance function}
}
\references{
Fong, Y., Yu, X. (2014) Transformation Model Choice in Nonlinear Regression Analysis of Serial Dilution Assays, submitted
}
\examples{

\dontrun{
dat.std=dat.QIL3[dat.QIL3$assay_id=="LMX001",]

# run 3 iter to save time for examples
fit.1=crm.fit(fi~expected_conc, dat.std, var.model="power", verbose=TRUE, max.iter=2)
fit.2=crm.fit(log(fi)~expected_conc, dat.std, verbose=TRUE)
fit.3=crm.fit(log(fi)~expected_conc, dat.std, var.model="power", verbose=TRUE, max.iter=2)

sapply(list(fit.1, fit.2, fit.3), coef)
fit.1$var.power
fit.2$var.power
fit.3$var.power

plot(fit.1, log="xy", type="all", lwd=3, pch="*")
lines(fit.2, expy=TRUE, col=2, lwd=3)
lines(fit.3, expy=TRUE, col=4, lty=2, lwd=3)
    
}


}
