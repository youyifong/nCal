\name{ncal}
\alias{nCal}
\alias{ncal}
\alias{rumi}
\alias{ncal.formula}
\alias{ncal.character}
\alias{ncalGUI}
\title{Main function for the nCal package}
\description{ 
  \code{ncal} fits standard curves and estimates unknown sample concentrations. rumi exists for backwards compatibility. 
}

\usage{

\method{ncal}{formula} (formula, data, 
    bcrm.fit=FALSE, bcrm.model="t4", robust="mean", bcrm.prior="default",
    force.fit=TRUE, fit.4pl=FALSE, return.fits=FALSE, 
    plot=TRUE, auto.layout=TRUE, plot.se.profile=TRUE, plot.log="x", plot.unknown=TRUE,
    test.LOD=FALSE, find.LOD=FALSE, find.LOQ=FALSE, grid.len=50, lod.ci=95,
    unk.replicate=NULL, find.best.dilution=FALSE, unk.median=FALSE,
    control.jags=list(n.iter=1e5, jags.seed=1, n.thin=NULL, 
        keep.jags.samples=FALSE, n.adapt=1e3), 
    cex=.5, additional.plot.func=NULL, check.out.of.range=1, xlab=NULL, ylab=NULL, 
    main=NULL,
    var.model=c("constant","power"), log.both.sides=FALSE,
    control.crm.fit=list(max.iter=20),
    verbose=FALSE, 
...)


\method{ncal}{character} (file, is.luminex.xls, formula, bcrm.fit, verbose=FALSE, ...)

rumi(data, ...)

ncalGUI (verbose=FALSE)

}

\arguments{

  \item{formula}{ a formula. If the first argument is formula, bcrm.formula is called. }
  \item{data}{ a data frame. If the first argument is dat, bcrm.data.frame is called. Each row of the data frame represents the measurement from one well/bead_type }
  
  \item{bcrm.fit}{Boolean. If TRUE, a Bayesian random effects model is fitted, else a drm model is fitted.}
  \item{bcrm.model}{string. Noise model used in bcrm. See error.model in the help file for bcrm for more information.}
  \item{bcrm.prior}{string.}
  \item{var.model}{constant: constant variance. power: power variance function}
  \item{robust}{string. Passed to drm.}
  
  \item{plot}{ a Boolean, default FALSE. This controls whether or not to make plots of the fitted curves and the standard error profiles }
  \item{plot.unknown}{ a Boolean. This controls whether or not to make plots of the fitted curves with unknown sample concentrations layered on top }
  \item{auto.layout}{ a Boolean, default TRUE This controls whether or not to let bcrm controls the layout of the plots }
  \item{test.LOD}{ a Boolean, default TRUE This controls whether or not to test for limits of detection }
  \item{plot.se.profile}{ a Boolean, default FALSE. This controls whether or not to make plots of standard error profile, this doubles as find.LOQ }
  \item{force.fit}{ a Boolean, default FALSE. If FALSE, return NULL when the goodness of fit is bad; otherwise, always return the fit }
  \item{find.LOD}{ Boolean. Controls whether or not to compute limits of detection.}
  \item{find.LOQ}{ Boolean. Controls whether or not to compute limits of quantification.}
  \item{find.best.dilution}{ a Boolean, default FALSE. When there are more than one dilutions for a given non-standard sample, this controls whether or not to find the best dilution }
  \item{return.fits}{ a Boolean, default FALSE. This controls whether or not to return as an attribute the curve fits }
  \item{grid.len}{ an integer, default 500 This determines the resolution of the standard error profile plot and limits of quantification }
  \item{unk.replicate}{ an integer, default NULL. This is the number of replicates for an unknown sample (at a single dilution if multiple dilutions are available) }
  \item{unk.median}{ a Boolean, default FALSE. If TRUE, median of unknown replicates are used to estimate conc, instead of mean}
  \item{fit.4pl}{ a Boolean, default FALSE }
  \item{lod.ci}{ an integer, default 95. If default, one sided 50 percent CI is used to determine LODi. }
  \item{plot.log}{ a string, default "x". Used to populate the log argument of the plot function.}
  \item{control.jags}{list. Parameters controlling the behavior of posterior sampling by JAGS.}
  \item{is.luminex.xls}{Boolean. Indicates whether the file is in the Luminex Excel file format outputted by the Bio-Plex software.}
  \item{cex}{numeric. Passed to plot functions.}
  \item{main}{plotting parameter}
  \item{xlab}{plotting parameter}
  \item{ylab}{plotting parameter}
  \item{additional.plot.func}{function. Will be called after the first panel is drawn.}
  \item{check.out.of.range}{integer. If 1, an estimated concentration is set to half of the smallest standard concentration, or the largest standard if it is outside the range of standard concentration. If 2, do this only when the estimated concentration is 0 or Inf}
  \item{verbose}{ a Boolean, default FALSE. This controls whether or not to print detailed messages }
  \item{log.both.sides}{Boolean, whether to log transform both sides of the formula}
  \item{control.crm.fit}{list. Parameters controlling the behavior of crm.fit}

  \item{file}{string. Name of the data file.}
  
  
  \item{...}{ passed to the function plotting fitted curves }

}

\details{

    Certain columns are expected of the input data frame:
    \itemize{
        \item analyte		Identifies a type of bead. Beads can be named after the substance that is recognized by the antibody that coats the bead, or they can be named by the protein that coats the bead. 
                            This is also known as analyte or antigen in different contexts. To be retro-compatible, this column can also be named analyte
        \item well_role     Defines the role of a well. This should be from {Standard, Unknown, ...}. Standard wells are used to generate standard curves, and concentrations of the substance in the Unknown well will be estimated
        \item assay_id      Identifies an assay, which is defined to be a collection of Standard and non-Standard wells. Measured fi from the Standard wells are used to create a set of standard curves, one of each bead type.        
                            Based on the standard curves and the measured fi from the non-Standard wells, concentrations of the substance in the non-Standard wells will be estimated. 
                            An assay can be a plate, if every plate has Standard wells; or it can be multiple plates run by one technician at one time, if there are only one set of Standard wells on these plates
        \item dilution      Standard samples are often prepared by starting with one sample and doing serial dilutions. 
                            Unknown samples may be measured at several dilutions so that one of the dilutions may fall into the more reliable region of the standard curve
        \item replicate     Index of technical replicates for a sample. Typical values are 1 or 2. May be used in plotting. Optional
        \item sample_id     Identifies a non-Standard sample. One sample may be measured in replicates and/or in several dilutions. Should be defined meaningfully for Unknowns
    }
    
    If no formula is specified, we also expect 
    \itemize{
        \item fi            Fluorescence intensity
        \item starting_conc Standard samples are often prepared by starting with one sample and doing serial dilutions. This column does not apply to non-Standard samples
        \item expected_conc Standard sapmles have expected concentrations. If this column is present, the dilution and starting_conc are optional and will not be used. This column does not apply to non-Standard samples
    }
    
    The program processes each assay separately. For each assay, each bead_type is processed separately. First fit.drc() is called to fit dose-response curves using information from the Standard wells. A plot of the fitted curve is generated.
    Then for each non-Standard sample, the number of dilutions is determined. For each dilution, log(fi) from replicate wells are averaged and used to compute the estimated concentration and the standard error of the estimate.
    
    Bad fits can happen for three reasons. 1) An error occurs in drm. 2) Estimated variance for some curve parameter is negative. 3) A goodness of fit statistics is above a pre-set threshold. 
    When a bad fit happens, by default we plot the observed fi, not the fitted curve, and do not proceed to estimation of concentration. 
    However, when force.fit is set to TRUE, in the latter two cases of bad fits, we will proceed to plot the fitted curve and use it to estimate concentrations. 
    In the second case, the standard errors of the estimated concentrations are set to NA. In the third case, standard errors will be computed.
    The default value for force.fit is set to FALSE to promote caution. 
    
    An important factor affecting the success of the fitting procedure is the choice of self start function. In addition to trying the default self start function in the most current drc package, we also try a self start function that
    is based on four parameter log-logistic model, and the self start function in version 1.5.2 of the drc package.
    
    Standard errors convey the uncertainty we have about estimated concentrations. Standard error profiles show the relationship between standard errors and estimated concentrations. 
    A sequence of hypothetical fi between the expected fi for the smallest and largest concentrations on the Standard curve are generated. The number of hypothetical fi is controlled by the variable grid.len. 
    For each hypothetical fi, the estimated concentration and associated standard error are compuated. 
    There are two sources of uncertainty in an estimated concentration. One part of the uncertainty, we call replication-sensitive, comes from the fact that the measured fi has measurement error in it. It can be reduced by doing 
    replicates of the non-Standard samples. The number of replicates used in computing standard error profile is controlled by the variable rep.se.profile. 
    The other part, we call replication-insensitive, comes from the fact that there are uncertainty about the Standard curve as well. 
    The total uncertainty is plotted in black, the replication-sensitive in red, and the replication-insensitive in blue. 
    The two parts of uncertainty add up to the total not on the standard error scale, but on the variance scale, which is the square of standard errors. 
    We choose to plot on the standard error scale because this scale is more comparable to the estimated concentration.
    
    When unk.replicate is NULL (default), the program sets it to the number of replicates of the first non-Standard sample it encounters; if there is no non-Standard sample, it is set to 1. This value is only used in the computation of LOD. It does not affect se of unknown sample conc estimate.
    
    When find.best.dilution is false, the estimated concentrations for all dilutions (after adjusting for dilutions) are returned; otherwise, the best dilution, determined as the dilution having the smallest standard error, is returned.
    
    When plot is FALSE, no plot is made. When plot is TRUE, but plot.rep.profile is false, only fitted curves are plotted.
    
    When auto.layout is TRUE, bcrm will choose a layout that makes sense for showing one analyte per page. For example, if plot.se.profile is TRUE, it will be 2x2.
    
    When test.LOD is true, estimated concentration will be tested against the null hypothesis that it is not different from the extremem data points of the standard samples
    
    The difference between Inf and NA for se: Inf is if the fi is outside certain ranges, NA is if the s.e. is bad for some reason.
    
    find.LOD, return an attribute "LOD", that is the lowest and highest concentration detectable, defined in the sense that ...
    plot.se.profile also leads to the return of an attribute "LOQ", that is the lowest and highest concentration at which percent cv is at 20
    
    drm requires that weights are evaluated in the global env. drm.weights (drm.fit.R) is our answer to that. 
    
    
    

}

\value{
A data frame, each row contains one estimated concentration. All columns of the input data frame are preserved, in addition two new columns are added: est.log.conc and se. est.log.conc contains estimated log concentration,
 while se is the standard error of the est.log.conc. 

When return.fits is TRUE, the returned data frame has an attribute "fits", that is a list of the fitted curves.

Standard curves are plotted. When plot.se.profile is TRUE, error profiles are also plotted. Two error profile are plotted, one is the standard error of the estimated log concentration versus estimated log concentration. 
The other is 100 x (standard error of estimated concentration / estimated concentration) versus log estimated concentration. 
}

\author{ Youyi Fong \email{yfong@fhcrc.org}, Xuesong Yu \email{xyu@scharp.org}.}

\references{
Fong, Y., Yu, X., Sebestyen, K., Gilbert, P. and Self, S. (2013) nCal: a R package for nonlinear calibration. Bioinformatics
}

\examples{

#begin=Sys.time()

# basic example

# simulate a dataset
set.seed(1)
log.conc=log(1e4)-log(3)*9:0
n.replicate=2
fi=simulate1curve (p.eotaxin[1,], rep(log.conc,each=n.replicate), sd.e=0.2)
dat.std=data.frame(fi, expected_conc=exp(rep(log.conc,each=n.replicate)), analyte="Test", 
    assay_id="Run 1", sample_id=NA, well_role="Standard", dilution=rep(3**(9:0), each=n.replicate),
    replicate=rep(1:n.replicate, 10))
# add unknown
dat.unk=rbind(
      data.frame(fi=exp(6.75), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=1, 
        well_role="Unknown", dilution=1, replicate=1)
    , data.frame(fi=exp(6.70), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=2, 
        well_role="Unknown", dilution=1, replicate=1)
    , data.frame(fi=exp(3), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=3, 
        well_role="Unknown", dilution=1, replicate=1)
    , data.frame(fi=exp(4.4), expected_conc=NA, analyte="Test", assay_id="Run 1", sample_id=4, 
        well_role="Unknown", dilution=10, replicate=1)
)
dat=rbind(dat.std, dat.unk)


\dontrun{
# to save some time


# does drm fit
out = ncal(log(fi)~expected_conc, dat, return.fits = TRUE, plot.se.profile=TRUE) 
fit=attr(out, "fits")[[1]]

# does jags fit and collect 1e5 posterior samples, it may be better to set n.iter higher
out.norm = ncal(log(fi)~expected_conc, dat, bcrm.fit=TRUE, bcrm.model="norm",
    return.fits = TRUE, plot.se.profile=TRUE,
    control.jags=list(n.iter=1e4), verbose=FALSE)
fit.norm=attr(out.norm, "fits")

# compare drm fit with bcrm fit
rbind(out, out.norm)
rbind(cla2gh(coef(fit)), coef(fit.norm))
rbind(sqrt(diag(vcov(fit))), sqrt(diag(vcov(fit.norm, type="classical"))) )
sd.est = c(summary(fit)$rseMat[1], getVarComponent.bcrm(fit.norm)^.5)
sd.est



# do ncal with imported data from a raw Luminex output file
# the importing step can step a litte time
out = ncal (paste(system.file(package="nCal")[1],'/misc/02-14A22-IgA-Biotin-tiny.xls',sep=""), 
    is.luminex.xls=TRUE, formula=log(fi)~expected_conc, bcrm.fit=FALSE, return.fits=TRUE)
fit=attr(out, "fits")[[1]]
getConc(fit, c(5,6))

}

#end=Sys.time();print(end-begin)


}
\keyword{ distribution }
