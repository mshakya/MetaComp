% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_merged_assignment.R
\name{plot_merged_assignment}
\alias{plot_merged_assignment}
\title{Generates a single column ggplot for a taxonomic assignment table.}
\usage{
plot_merged_assignment(assignment, taxonomy_level,
  sorting_order = "abundance", row_limit = 60, plot_title, filename)
}
\arguments{
\item{assignment}{The gottcha-like merged assignment table.}

\item{taxonomy_level}{The level which need to be plotted.}

\item{sorting_order}{the order in which rows shall be sorted, "abundance" is defult,
"alphabetical" is an alternative.}

\item{row_limit}{the max amount of rows to plot (default is 60).}

\item{plot_title}{The plot title.}

\item{filename}{The output file mask, PDF and SVG files will be produced with Cairo device.}
}
\description{
This implementation...
}

