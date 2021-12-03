\name{bcrm}
\alias{bcrm}
\alias{getVarComponent.bcrm}
\alias{vcov.bcrm}
\alias{coef.bcrm}
\alias{plot.bcrm}
\alias{print.bcrm}
\alias{get.single.fit}
\title{ Bayesian Concentration-Response Model}
\description{ 
  \code{bcrm} fit concentration-response curves with a Bayesian random effects model using JAGS
}
    

\usage{

bcrm (formula, data, 
    parameterization=c("gh","classical"), 
    error.model=c("norm","t4","mixnorm","mix2","replicate_re","tar1","lar1"), 
    prior=c("cytokine","BAMA","RT-PCR","ELISA","default"), 
    prior.sensitivity=c("none","1","2","3","4"),
    mean.model=c("5PL","4PL"),
    n.iter=1e5, jags.seed=1, n.thin=NULL, T.init=NULL, 
    keep.jags.samples=FALSE, standards.only=TRUE, n.adapt=1e3,
    t.unk.truth=NULL, params.true=NULL, # for simulation study use
    verbose=FALSE
)

\method{plot}{bcrm} (x, 
    assay_id=NULL, fit.2=NULL, fit.3=NULL, 
    points.only=FALSE, all.lines.only=FALSE, 
    same.ylim=FALSE, lty3=NULL, lcol2=NULL, lcol3=NULL, 
    lcol=1, lwd=.1, lty=1, # for lines
    t=NULL, log="x", col.outliers=TRUE, pch.outliers=TRUE, 
    use.dif.pch.for.replicate=FALSE, main=NULL,
    additional.plot.func=NULL, add=FALSE, ...
)


\method{print}{bcrm} (x, ...)
\method{coef}{bcrm} (object, type="gh", ...)
\method{vcov}{bcrm} (object, type="gh", ...)
\method{getVarComponent}{bcrm} (object, ...) 
get.single.fit (fit, assay_id) 

}

\arguments{

  \item{formula}{formula. Gives the response column and concentration column.}
  \item{data}{ a data frame. Each row represents the measurement from one well/bead_type. See details }
  \item{parameterization}{ string.}
  \item{error.model}{string.}
  \item{prior}{string.}
  \item{mean.model}{mean model}
  \item{pch.outliers}{pch for outliers}
  \item{n.iter}{ a number indicating the number of iterations to run.}
  \item{jags.seed}{ a number to seed the random number generator used within jags.}
  \item{n.thin}{ a number specifying the thinning factor applied to the jags samples.}
  \item{keep.jags.samples}{boolean. If TRUE, the fit object being returned has an element named "jags.samples". coef samples will always be saved in "coef.samples".}
  \item{t.unk.truth}{True unknown concentrations, for simulation study use only.}
  \item{params.true}{True curve parameters, for simulation study use only.}
  \item{T.init}{a integer vector. Initial values for mixture indicators.}
  \item{prior.sensitivity}{integer. A number between 1 and 4. Change priors.}
  \item{standards.only}{boolean. If TRUE, data is subset to standard samples only.}
  \item{n.adapt}{integer. Passed to jags.model. If 0, then no adaptation happens and reproducible results can be obtained from jags.model.}
  \item{verbose}{boolean. If TRUE, debug messages are printed.}
  
  \item{x}{bcrm fit object.}
  \item{object}{bcrm fit object.}
  \item{type}{string. 5PL parameterization.}
  \item{fit}{bcrm fit object.}
  \item{...}{...}
  \item{assay_id}{string. Label for the assay run.}


  \item{add}{Boolean. If TRUE, adding to an existing plot.}
  \item{lcol}{integer. Line color.}
  \item{fit.2}{a bcrm object. A second fit object to be plotted together with x.}
  \item{lwd}{numeric. Line width.}
  \item{points.only}{Boolean. If TRUE, only the data points are plotted and not the fitted curves}
  \item{all.lines.only}{Boolean. If TRUE, only the fitted curves are plotted.}
  \item{t}{a numeric vector. The log concentrations.}
  \item{same.ylim}{Boolean. If TRUE, all fitted curves are plotted with the same ylim.}
  \item{lty3}{integer. lty for plotting fit.3.}
  \item{fit.3}{bcrm object. A third optional bcrm object to be plotted.}
  \item{lcol2}{integer. Line color for plotting fit.2.}
  \item{lcol3}{integer. Line color for plotting fit.3.}
  \item{col.outliers}{Boolean. If TRUE, outliers are colored differently.}
  \item{main}{string.}
  \item{lty}{integer. Line type for plotting x.}
  \item{log}{string. If it is "x", then the x axis is labeled on the scale of concentration; otherwise, it is labeled on the scale of log concentration.}
  \item{additional.plot.func}{function, to be called after plotting each fit}
  \item{use.dif.pch.for.replicate}{boolean}

}

\details{

data is expected to contain one to many plates with the same analyte. 
    \itemize{
        \item well_role     Defines the role of a well. This should be from {Standard, Unknown, ...}. Standard wells are used to generate standard curves, and concentrations of the substance in the Unknown well will be estimated
        \item assay_id      Identifies an assay, which is defined to be a collection of Standard and non-Standard wells. Measured fi from the Standard wells are used to create a set of standard curves, one of each bead type.        
                            Based on the standard curves and the measured fi from the non-Standard wells, concentrations of the substance in the non-Standard wells will be estimated. 
                            An assay can be a plate, if every plate has Standard wells; or it can be multiple plates run by one technician at one time, if there are only one set of Standard wells on these plates
        \item dilution      Standard samples are often prepared by starting with one sample and doing serial dilutions. 
                            Unknown samples may be measured at several dilutions so that one of the dilutions may fall into the more reliable region of the standard curve
        \item replicate     Index of technical replicates for a sample. Typical values are 1 or 2. May be used in plotting. Optional
        \item expected_conc Standard sapmles have expected concentrations. If this column is present, the dilution and starting_conc are optional and will not be used. This column does not apply to non-Standard samples
    }
    
Main error.model supported: drc, classical_norm, classical_t4, classical_mixnorm, classical_lar1, gh_norm, gh_mixnorm, gh_lar1 
Also support: classical_replicate_re, gh_replicate_re, gh_tar1

Only two replicates are supported for now for the correlated noise models.

Sometimes jags.model fails with one .RNG.seed. The function will increase it by 1 and try again. The function tries three times before giving up.


}

\value{
An object of type bcrm. 
}

\author{ Youyi Fong \email{yfong@fhcrc.org} }
\references{
Fong, Y., Wakefield, J., DeRosa, S., Frahm, N. (2012) A robust Bayesian random effects model for nonlinear calibration problems, Biometrics, 68:1103-1112.
}
\examples{


\dontrun{

set.seed(1)
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
fi=simulate1curve (p.eotaxin[1,], rep(log.conc,each=n.replicate), sd.e=0.3)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="test", 
    assay_id="assay1", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate))
# add unknown
dat.unk=rbind(
      data.frame(fi=exp(6.75), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=1, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(6.70), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=2, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(3), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=3, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(4.4), expected_conc=NA, analyte="test", assay_id="assay1", 
        sample_id=4, well_role="Unknown", dilution=10)
)
dat=rbind(dat.std, dat.unk)
# second plate
fi=simulate1curve (p.eotaxin[2,], rep(log.conc,each=n.replicate), sd.e=0.3)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="test", 
   assay_id="assay2", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate))
# add unknown
dat.unk=rbind(
      data.frame(fi=exp(6.75), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=1, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(6.70), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=2, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(3), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=3, well_role="Unknown", dilution=1)
    , data.frame(fi=exp(4.4), expected_conc=NA, analyte="test", assay_id="assay2", 
        sample_id=4, well_role="Unknown", dilution=10)
)
dat=rbind(dat, dat.std, dat.unk)

fits = bcrm(log(fi)~expected_conc, dat, parameterization="gh", error.model="norm", 
    prior="cytokine", n.iter=6e3)
par(mfrow=c(1,2))
plot(fits)


# takes longer

# Example from Fong et al. (2012) 
fits.t4 = bcrm (log(fi)~expected_conc, dat.QIL3, parameterization="gh", error.model="t4", 
    prior="cytokine")
par(mfrow=c(2,3))
plot(fits.t4)

fits.norm = bcrm (log(fi)~expected_conc, dat.QIL3, parameterization="gh", error.model="norm",
    prior="cytokine")
par(mfrow=c(2,3))
plot(fits.norm)


    
}

}
\keyword{ distribution }
