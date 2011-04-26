\name{ssGAM2Bugs}
\alias{ssGAM2Bugs}
\title{Convert samples from a model fitted with spikeSlabGAM into a bugs-object...}
\usage{ssGAM2Bugs(m, rm=c("alpha", "ksi", "gamma"))
}
\description{Convert samples from a model fitted with \code{spikeSlabGAM} into a \code{bugs}-object}
\details{Use \code{plot} and \code{print} on the returned \code{\link[R2WinBUGS]{bugs}} object for convergence diagnostics.
Lack of convergence for \eqn{\alpha} or \eqn{\xi} does not necessarily mean that the sampler has not converged for \eqn{\beta}.}
\value{a \code{\link[R2WinBUGS]{bugs}} object which has convenience functions for assessing convergence based on parallel MCMC runs}
\author{Fabian Scheipl}
\arguments{\item{m}{a model fitted with \code{spikeSlabGAM}}
\item{rm}{which parameter blocks to omit from the conversion. Defaults to omission of \eqn{\alpha, \xi}, and \eqn{\gamma}.  
Set to NULL to keep all MCMC samples.}
}
\examples{#see help for spikeSlabGAM}

