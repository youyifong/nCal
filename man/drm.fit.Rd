\name{drm.fit}
\alias{drm.fit}
\alias{getVarComponent.drc}
\title{Fit drm}
\description{ 
  \code{drm.fit} fit concentration-response curves using drm function from drc package. 

}
    

\usage{

drm.fit (formula, data, robust="mean", fit.4pl=FALSE, w=NULL, gof.threshold=.2, 
    verbose=FALSE, bcVal = NULL, bcAdd = 0)

\method{getVarComponent}{drc} (object, ...) 

}

\arguments{

  \item{formula}{
a formula object.}
  \item{data}{
a data frame object.}
  \item{robust}{
a string. Passed to drm. See ?drm for more details.}
  \item{fit.4pl}{
boolean. If TRUE, 4PL model is fitted. If FALSE, 5PL model is fitted.}
  \item{gof.threshold}{a threshold to determine when to try more self start functions}
  \item{w}{weights}
  \item{bcVal}{numeric, passed to drm}
  \item{bcAdd}{numeric, passed to drm}
  \item{...}{...}
  \item{verbose}{
Boolean. If TRUE, verbose messages are printed.
}
  \item{object}{
a drm object.
}

}

\details{
drm.fit differs from drc::drm in several aspects. 

(1) It tries several self start functions in order to get better fits. 

(2) It uses gof.threshold to report lack of fit.

(3) It tried to determine whether the standard deviation of the parameter estimates can be estimated.

}

\value{
An object of type drm. 
}

\examples{

# simulate a dataset
set.seed(1)
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
fi=simulate1curve (p.eotaxin[1,], rep(log.conc,each=n.replicate), sd.e=0.2)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="Test", 
    assay_id="Run 1", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate),
    replicate=rep(1:n.replicate, 10))

fit = drm.fit(log(fi) ~ expected_conc, dat = dat.std)
plot(fit, log="xy")
fit

}

\keyword{ distribution }
