% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StructuralDecompose.R
\name{MeanCleaning}
\alias{MeanCleaning}
\title{Mean level checks}
\usage{
MeanCleaning(timeseries, mean_level = 0.5, breaks, frequency = 52)
}
\arguments{
\item{timeseries}{Given time series}

\item{mean_level}{Mean distance between two levels}

\item{breaks}{breakpoints returned}

\item{frequency}{Timeseries frequency, defaults to 12 points}
}
\value{
The series cleaned with the mean check
}
\description{
Mean level checks
}
\examples{
MeanCleaning(timeseries = StructuralDecompose::Nile_dataset[,1], breaks = c(1,4,5), frequency = 1)

MeanCleaning(timeseries = runif(n = 50, min = 1, max = 10), breaks = c(1,4,5), frequency = 12)
}
