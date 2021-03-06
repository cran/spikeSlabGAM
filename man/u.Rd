% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/terms.R
\name{u}
\alias{u}
\title{Generate design for an always included covariate}
\usage{
u(x, ...)
}
\arguments{
\item{x}{covariate}

\item{...}{arguments passed to \code{\link{model.matrix}}}
}
\value{
a design matrix for x
}
\description{
Basically a wrapper for \code{model.matrix(~ x, ...)}.
}
\author{
Fabian Scheipl
}
