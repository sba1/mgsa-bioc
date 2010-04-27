\name{mgsa.main}
\alias{mgsa.main}
\title{mgsa.main}
\usage{mgsa.main(o, sets, population, alpha=NA, beta=NA, p=NA, steps=1e+06,
    restarts=1, threads=0, as=integer(0), debug=0)}
\description{The main function that treats the case of character and integer inputs.}
\arguments{\item{o}{}
\item{sets}{}
\item{population}{}
\item{alpha}{}
\item{beta}{}
\item{p}{}
\item{steps}{number of MCMC steps to be performed}
\item{restarts}{defines the number of restarts.}
\item{threads}{defines number of threads to be used. Defaults to 0 which means that it
corresponds to the number of available cores.}
\item{as}{define the active sets. No MCMC is run in this case. Just the log likelihood is returned.}}
