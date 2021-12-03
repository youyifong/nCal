\name{dat.QIL3}
\alias{dat.QIL3}
\docType{data}
\title{
An example for hierarchical modeling used in Fong, Wakefield, De Rosa, Frahm (2012)
}
\description{
An example for hierarchical modeling used in Fong, Wakefield, De Rosa, Frahm (2012)
}
\format{
  A data frame with 120 observations on the following 13 variables.
  \describe{
    \item{\code{well}}{a character vector. Well identifier on a microplate.}
    \item{\code{assay_id}}{a character vector. Assay identifier.}
    \item{\code{analyte}}{a character vector. Substance to be measured.}
    \item{\code{well_role}}{a character vector. Defines the role of a well. See ?ncal for more information.}
    \item{\code{beadct}}{a numeric vector. Number of beads in the well, specific to multiplex bead array assay.}
    \item{\code{dilution}}{a numeric vector. The dilution factor of a sample.}
    \item{\code{expected_conc}}{a numeric vector. The expected concentrations of a sample.}
    \item{\code{fi}}{a numeric vector. Fluoresence intensity readout.}
    \item{\code{ptid}}{a integer vector. Participant ID.}
    \item{\code{sample_id}}{a Boolean vector. All NA.}
    \item{\code{thawdt}}{a character vector. All empty string}
    \item{\code{visit}}{a integer vector. All NA.}
    \item{\code{replicate}}{a numeric vector. Techinical replicate identifier.}
  }
}
\references{
Fong, Y., Wakefield, J., DeRosa, S., Frahm, N. (2012) A robust Bayesian random effects model for nonlinear calibration problems, Biometrics, 68:1103-1112.
}
