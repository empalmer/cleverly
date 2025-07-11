% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleverly_visualizations.R
\name{plot_one_cluster_old}
\alias{plot_one_cluster_old}
\title{Plot fits and data for a single cluster}
\usage{
plot_one_cluster_old(
  res,
  cluster_val,
  response_names = NULL,
  curve_type = "baseline",
  nrow = 3,
  scales = "fixed",
  Y_counts = NULL
)
}
\arguments{
\item{res}{cleverly result}

\item{cluster_val}{numeric which cluster to plot}

\item{response_names}{vector of response naes}

\item{curve_type}{"baseline" "refit" or "slope"}

\item{nrow}{number of rows for the facet plot to include}

\item{scales}{"fixed" or "free_y"}

\item{Y_counts}{needed if curve_type = "slope" and continuous Z}
}
\value{
ggplot object for one cluster
}
\description{
Plot fits and data for a single cluster
}
