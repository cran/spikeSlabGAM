\name{plot.spikeSlabGAM}
\alias{plot.spikeSlabGAM}
\title{Generates graphical summaries of a fitted model}
\usage{
  \method{plot}{spikeSlabGAM} (x, labels = NULL, cumulative
  = TRUE, commonEtaScale = FALSE, aggregate = mean,
  quantiles = c(0.1, 0.9), gridlength = 20, base_size = 12,
  ggElems = list(), rows = NULL, cols = NULL, widths =
  NULL, heights = NULL, align = NULL, maxPlotsPerPage = 6,
  ...)
}
\arguments{
  \item{x}{a fitted \code{spikeSlabGAM} model}

  \item{labels}{a character vector of names of model terms
  to be plotted}

  \item{cumulative}{Defaults to TRUE, in which case all
  lower order terms that are involved in an interaction are
  cumulated and then plotted (e.g, if a model contains 2
  smooth effects and their interaction, ONLY the sum of the
  marginal smooth and linear terms and all their
  interactions are plotted.) If FALSE, a separate plot for
  every term in the model is produced.}

  \item{commonEtaScale}{use the same limits for all
  vertical axes of the different panels? Defaults to FALSE.
  Can be useful to compare effect sizes more effectively
  between panels, but tends to mess up the scales.}

  \item{aggregate}{(function) which summary statistic to
  use for the posterior of the model terms. Defaults to the
  mean.}

  \item{quantiles}{which quantiles to use for the borders
  of credible regions. Defaults to 10\% and 90\%
  percentiles. Set to NULL to omit plotting credible
  regions.}

  \item{gridlength}{length of the (univariate) grids for
  the covariates on which to evaluate the posterior.
  Defaults to 20.}

  \item{base_size}{default base font size for plot (see
  e.g. \code{\link[ggplot2]{theme_gray}})}

  \item{ggElems}{a list of plot elements to give to
  \code{ggplot}. Use this to supply custom themes or
  colors, fortify the plot(s) with partial residuals etc.
  Defaults to an empty list. Unless specified differently
  here, the default ggplot-theme
  (\code{\link[ggplot2]{theme_gray}}) is changed to a white
  background with major gridlines in gray
  (\code{'grey95'}), no minor grid lines, and smaller text
  for the legends.}

  \item{rows}{a list controlling the positioning of the
  various plots on the device. See
  \code{\link{arrange.ggplots}}}

  \item{cols}{a list controlling the positioning of the
  various plots on the device. See
  \code{\link{arrange.ggplots}}}

  \item{widths}{a list controlling the positioning of the
  various plots on the device. See
  \code{\link{arrange.ggplots}}}

  \item{heights}{a list controlling the positioning of the
  various plots on the device. See
  \code{\link{arrange.ggplots}}}

  \item{align}{a vector controlling the positioning of the
  various plots on the device. See
  \code{\link{arrange.ggplots}}. Defaults to no alignment.}

  \item{maxPlotsPerPage}{maximum number of plots rendered
  on a single page/device. Defaults to 9.}

  \item{...}{arguments passed from or to other methods (not
  used)}
}
\value{
  a list of \code{\link[ggplot2]{ggplot}}-objects
  (invisible)
}
\description{
  This function plots the estimated linear predictors of
  the terms in a model on a grid of values. By default
  displays all 3-way, 2-way interactions and main effects
  present in the model.
}
\note{
  Note that \code{cumulative=TRUE} will only find all
  relevant terms to accumulate if, for all numeric
  covariates that have a smooth term, the smooth term is
  specified \emph{after} the linear term in the formula.
}
\examples{
#see ?spikeSlabGAM
}
\author{
  Fabian Scheipl
}
\seealso{
  \code{\link{plotTerm}} for more details on the specific
  plots
}

