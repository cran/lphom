% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.lphom}
\alias{plot.lphom}
\title{Graphical representation of a RxC ecological inference (vote transfer) matrix}
\usage{
\method{plot}{lphom}(
  x,
  complete = FALSE,
  margins = TRUE,
  digits = 2,
  row.names = NULL,
  col.names = NULL,
  size.numbers = 6,
  size.labels = 4,
  size.margins = 4,
  colour.cells = "deeppink3",
  colour.grid = "blanchedalmond",
  alpha = 0.5,
  which = NULL,
  ...,
  type = "w",
  show.plot = TRUE
)
}
\arguments{
\item{x}{An object output of a \code{lphom}  family algorithm.}

\item{complete}{A \code{TRUE/FALSE} argument informing if the complete matrix should be displayed. In \code{raw}, \code{regular}, \code{ordinary} and \code{enriched} scenarios the plot includes the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1\% in all units. Default, \code{FALSE}.}

\item{margins}{A \verb{TRUE/FALSE`` argument informing if the margins of the matrix should be displayed. Default }TRUE`.}

\item{digits}{Integer indicating the number of decimal places to be shown. Default, 2.}

\item{row.names}{Names to be used for the rows of the matrix.}

\item{col.names}{Names to be used for the columns of the matrix.}

\item{size.numbers}{A reference number indicating the average font size to be used for the transfer numbers. Default, 6.}

\item{size.labels}{A number indicating the font size to be used for labels. Default, 4.}

\item{size.margins}{A number indicating the font size to be used for margin numbers. Default, 4.}

\item{colour.cells}{Background base colour for cells.}

\item{colour.grid}{Colour to be used for grid lines.}

\item{alpha}{A [0,1] number of colour transparency.}

\item{which}{A vector of integers informing the units for which the aggregate transfer matrix should be plotted. Default, \code{NULL}, the global matrix is shown.}

\item{...}{Other arguments passed on to methods. Not currently used.}

\item{type}{A character string indicating the solution (transfer matrix) to be plotted. Only valid for \strong{_dual} algorithms. \code{type = "w"} stands for the weighted solution and \code{type = "a"} for the simple average solution. Default \code{w}.}

\item{show.plot}{A \code{TRUE/FALSE} argument indicating if the plot should be displayed as a side-effect. By default, \code{TRUE}.}
}
\value{
Invisibly returns the (ggplot) description of the plot, which is a list with components that contain the plot itself, the data, information about the scales, panels etc.
}
\description{
Plot method for objects obtained with an algorithm of the lphom-family (lphom, tslphom, nslphom, tslphom_dual, nslphom_joint, ....).
}
\note{
ggplot2 is needed to be installed for this function to work.
}
\examples{
mt.ns <- nslphom(France2017P[, 1:8] , France2017P[, 9:12], new_and_exit_voters= "raw")
p <- plot(mt.ns, show.plot = FALSE)
p

}
\author{
Jose M. Pavia, \email{pavia@uv.es}
}
