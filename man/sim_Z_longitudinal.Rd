% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_helpers.R
\name{sim_Z_longitudinal}
\alias{sim_Z_longitudinal}
\title{Use Chenyangs setup to simulate count data wtih 3 clusters}
\usage{
sim_Z_longitudinal(
  n = 20,
  range_start = 5000,
  range_end = 20000,
  nknots = 3,
  K = 12,
  order = 3,
  user_var = 1000,
  cor_str,
  al = 0.4,
  miss_p = 0.6,
  slope_base = "cluster_base_alldiff_slope"
)
}
\arguments{
\item{slope_base}{}
}
\value{
data Matrix with columns time, individual, capture number, totaln, counts
}
\description{
No external variables.
}
