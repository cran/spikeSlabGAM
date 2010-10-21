\name{arrange.ggplots}
\alias{arrange.ggplots}
\title{Arrange multiple ggplots on the same device, trying to align axes etc.}
\usage{arrange.ggplots(..., rows=as.list(rep(1, length(dots))),
    cols=as.list(1:length(dots)), widths=rep(1,
    max(unlist(rows))), heights=rep(1, max(unlist(cols))),
    align=rep(TRUE, length(dots)))
}
\description{Arrange multiple ggplots on the same device, trying to align axes etc.}
\details{This function takes a list of \code{ggplot}-objects and arranges them in a grid, trying to align the axes.
Horizontal alignment for plots with titles will not be perfect.}
\author{Harish V., based on \code{align.r} in package \code{ggextra}; wrapper function by Fabian Scheipl}
\arguments{\item{...}{ggplot objects}
\item{rows}{a list with as many entries as there are elements of \code{...}, each list entry can be a scalar
or a vector of integers specifying which rows of the grid layout the respective plot is to use. Defaults to a list of \code{1}s,
i.e. all the subplots given in \code{...} are rendered side by side.}
\item{cols}{a list with as many entries as there are elements of \code{...}, each list entry can be a scalar
or a vector of integers specifying which columns of the grid layout the respective plot is to use. Defaults to a list of \code{1} to
\code{length(...)}, i.e. all the subplots given in \code{...} are rendered side by side.}
\item{widths}{a vector of relative widths for the columns. Defaults to the same width for all columns.}
\item{heights}{a vector of relative heights for the rows. Defaults to the same height for all rows.}
\item{align}{a vector of logicals the same length \code{...}. Defaults to TRUE for all subplots. If FALSE, no attempt is made to
align the axes of the respective plot with the surroundding subplots.}
}

