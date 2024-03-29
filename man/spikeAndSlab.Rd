% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spikeAndSlab.R
\name{spikeAndSlab}
\alias{spikeAndSlab}
\title{Set up and sample a spike-and-slab prior model.}
\usage{
spikeAndSlab(
  y,
  X,
  family = c("gaussian", "binomial", "poisson"),
  hyperparameters = list(),
  model = list(),
  mcmc = list(),
  start = list()
)
}
\arguments{
\item{y}{response}

\item{X}{design matrix}

\item{family}{(character) the family of the response, defaults to
normal/Gaussian response}

\item{hyperparameters}{a list of hyperparameters controlling the priors (see
details)}

\item{model}{a list with information about the grouping structure of the
model (see details)}

\item{mcmc}{(optional) list setting arguments for the sampler (see details)}

\item{start}{(optional) list containing the starting values for \eqn{\beta,
\gamma, \tau^2, \sigma^2, w} and, optionally, the random seed}
}
\value{
a list with components:\describe{ \item{\code{formula}}{see
  arguments} \item{\code{data}}{see arguments} \item{\code{family}}{see
  arguments} \item{\code{y}}{see arguments} \item{\code{X}}{see arguments}
  \item{\code{hyperparameters}}{see arguments} \item{\code{model}}{see
  arguments} \item{\code{mcmc}}{see arguments} \item{\code{start}}{see
  arguments} \item{\code{posteriorPred}}{a list with entries \code{mu} and
  \code{y} containing samples of the expected values and realizations of the
  response from the posterior predictive} \item{\code{postMeans}}{a list
  containing the posterior means of the parameters: \describe{
  \item{\code{beta}}{the regression coefficients} \item{\code{alpha}}{}
  \item{\code{ksi}}{} \item{\code{tau}}{hypervariances of the penalized model
  terms} \item{\code{gamma}}{inclusion indicator variables of the model
  terms} \item{\code{pV1}}{\eqn{P(\gamma = 1)}} \item{\code{w}}{hyperparameter
  for \code{gamma}} \item{\code{sigma2}}{error variance (for Gaussian data)}
  \item{\code{logLik}}{log likelihood} \item{\code{logPost}}{log of
  (unnormalized) posterior}}} \item{\code{samples}}{a list containing the
  posterior samples of the parameters, see above for explanation of the
  entries} \item{\code{DIC}}{a vector with \eqn{DIC, pD, \bar{D},\hat{D}}.
  Usually doesn't make much sense for this kind of model because of the
  posterior's multimodality.} \item{\code{fitted}}{a matrix with the
  posterior mean of the linear predictor in the first column and the
  posterior mean of the expected response in the second.}
  \item{\code{runTime}}{of the sampler, in seconds}}
}
\description{
This function sets up a spike-and-slab model for variable selection and model
choice in generalized additive models and samples its posterior. It uses a
blockwise Metropolis-within-Gibbs sampler and the redundant multiplicative
parameter expansion described in the reference. This routine is not meant to
be called directly by the user -- \code{\link{spikeSlabGAM}} provides a
formula-based interface for specifying models and takes care of (most of) the
housekeeping. Sampling of the chains is done in parallel using package
\code{parallel}. A "SOCK" cluster is set up under Windows to do so (and
closed after computations are done, I try to clean up after myself), see
\code{\link[parallel]{makeCluster}} etc. Use \code{options(mc.cores =<foo>)}
to set the (maximal) number of processes forked by the parallelization. If
\code{options()$mc.cores} is unspecified, it is set to 2.
}
\details{
Details for model specification: \describe{
\item{\code{hyperparameters}}{\describe{
\item{\code{w}}{hyperparameters for the \eqn{Beta}-prior for \eqn{w};
defaults to \code{c(alphaW = 1, betaW = 1)}, i.e. a uniform distribution.}
\item{\code{tau2}}{hyperparameters for the \eqn{\Gamma^{-1}}-prior of the
hypervariances \eqn{\tau^2}; defaults to \code{c(a1 = 5, a2 = 25)}}
\item{\code{gamma}}{sets \eqn{v_0}, the ratio between the spike and slab
variances, defaults to \code{c(v0 = 0.00025)}}
\item{\code{sigma2}}{hyperparameters for \eqn{\Gamma^{-1}}-prior for error
variance; defaults to \code{c(b1 = 1e-4, b2 = 1e-4)}. Only relevant for Gaussian
response.} \item{\code{varKsi}}{variance for prior of \eqn{\xi}, defaults to
1} \item{\code{ksiDF}}{defaults to 0 for a gaussian prior for \eqn{\xi}, else
induces a t-prior for \eqn{\xi}} with \code{ksiDF} degrees of freedom.}}
\item{\code{model}}{\describe{
\item{\code{groupIndicators}}{a factor that maps the columns of X to the
different model terms} \item{\code{H}}{a matrix containing the hierarchy of
the penalized model terms} \item{\code{n}}{number of observations}
\item{\code{q}}{length of \eqn{\beta}} \item{\code{scale}}{scale/weights of
the response, defaults to \code{rep(1, n)}, use this to specify number of
trials for binomial data} \item{\code{offset}}{defaults to \code{rep(0,
n)}}}} 
\item{\code{mcmc}}{\describe{
\item{\code{nChains}}{how many parallel chains to run: defaults to 3}
\item{\code{chainLength}}{how many samples should be generated per chain,
defaults to 500} \item{\code{burnin}}{how many initial iterations should be
discarded, defaults to 100} \item{\code{thin}}{save only every \code{thin}-th
iteration, defaults to 5} \item{\code{verbose}}{verbose output and  report
progress? defaults to TRUE} \item{\code{returnSamples}}{defaults to TRUE}
\item{\code{sampleY}}{generate samples of y and its conditional expectation
from posterior predictive? defaults to FALSE}
\item{\code{useRandomStart}}{use random draw or ridge estimate for beta as
starting value? defaults to TRUE, i.e. random starting values.}
\item{\code{blocksize}}{approx. blocksizes of the updates for \eqn{\alpha,
\xi}. Defaults to 50 for gaussian responses and 5/15 for non-gaussian
responses.} \item{\code{scalemode}}{how to do term-wise rescaling of
subvectors of \eqn{\xi} in each iteration:  0 means no rescaling, 1 means
rescaling s.t. each mean\eqn{(|\xi_g|) = 1},  2 means rescaling s.t. each
max\eqn{(|\xi_g|) = 1}} \item{\code{modeSwitching}}{probability to do P-IWLS
with the mode of the proposal set to the current value, which is useful if
the chain gets stuck. Defaults to \eqn{0.05}. Increase this if acceptance rates
are too low.} \item{\code{reduceRet}}{don't return data and samples for
\eqn{\alpha, \xi, \tau^2}? defaults to FALSE}}}
\item{\code{start}}{\describe{\item{\code{beta}}{starting
values for \eqn{\beta}. Defaults to a modified approximate ridge-penalized ML
estimate. See vignette for details on default specification.}
\item{\code{gamma}}{starting values for \eqn{\gamma}. Defaults to a vector of
1's if \code{mcmc$useRandomStart} is \code{FALSE}, otherwise drawn from the
prior.} \item{\code{tau2}}{starting values for \eqn{\tau^2}. Defaults to the
mode of the prior if \code{mcmc$useRandomStart} is \code{FALSE}, otherwise
drawn from the prior.} \item{\code{sigma2}}{starting values for
\eqn{\sigma^2}. Only relevant for Gaussian response. Defaults to the variance
of the response divided by the number of covariates if
\code{mcmc$useRandomStart} is \code{FALSE}, otherwise drawn from the prior.}
\item{\code{w}}{starting value for \eqn{w}. Defaults to the mean of the prior
if \code{mcmc$useRandomStart} is \code{FALSE}, otherwise drawn from the
prior.} \item{\code{seed}}{Sets RNG seed for reproducible results. Parallel
chains are seeded with this seed incremented by the number of the chain.}}}}
}
\references{
Scheipl, F. (2010) Normal-Mixture-of-Inverse-Gamma Priors for
  Bayesian Regularization and Model Selection in Structured Additive
  Regression Models. \emph{LMU Munich, Department of Statistics}: Technical
  Reports, No.84 (\url{https://epub.ub.uni-muenchen.de/11785/})
}
\author{
Fabian Scheipl, Daniel Sabanes Bove
}
