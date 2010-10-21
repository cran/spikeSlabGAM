\name{print.summary.spikeSlabGAM}
\alias{print.summary.spikeSlabGAM}
\title{Print summary for posterior of a spikeSlabGAM fit...}
\usage{\method{print}{summary.spikeSlabGAM}(x, digits=3, printPGamma=TRUE, printModels=TRUE, showModels,
    ...)
}
\description{Print summary for posterior of a \code{spikeSlabGAM} fit}
\details{The model table shows at least the 10 and at most the 20 most probable models as found by SSVS, or enough models to
account for 90\% of the probability mass in the model space as found by SSVS by default.}
\seealso{\code{\link{summary.spikeSlabGAM}}}
\value{invisibly returns \code{x}}
\author{Fabian Scheipl}
\arguments{\item{x}{an object of class \code{summary.spikeSlabGAM}}
\item{printPGamma}{print marginal inclusion probabilities and norm of the linear predictor for each term? Defaults to TRUE}
\item{printModels}{print most probable models visited by the SSVS? Defaults to TRUE}
\item{showModels}{how many of th most probable models to display. See details for default.}
\item{digits}{see \code{\link[base]{options}}, defaults to 3}
\item{...}{arguments based from or to other methods (not used)}
}

