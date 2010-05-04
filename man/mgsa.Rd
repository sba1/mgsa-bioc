\name{mgsa}
\alias{mgsa}
\title{mgsa}
\description{Performs an Mgsa analysis.}
\details{Currently, the Mgsa problem is solved using an MCMC sampling algorithm.}
\references{GOing Bayesian: model-based gene set analysis of genome-scale data.}
\arguments{\item{o}{specifies the items which are observed.}
\item{sets}{an instance of class MgsaSets. Alternatively, a list of sets can be specified.
Each set is a vector that contains associated items. The vector can be of any data type,
for instance, integers or characters.}
\item{population}{}
\item{alpha}{}
\item{beta}{}
\item{p}{}
\item{steps}{defines the number of the MCMC sampler. A recommended value is 1e6 or greater.}
\item{restarts}{defines the number of MCMC restarts. Must be greater or equal to 1.}
\item{threads}{defines the number of threads that should be used. A value of 0 means to use all available cores. Default to 0.}}
