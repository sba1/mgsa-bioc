\name{mgsa.main}
\alias{mgsa.main}
\title{mgsa.main}
\usage{mgsa.main(o, sets, population, alpha=NA, beta=NA, p=NA, steps=1e+06,
    restarts=1, threads=0, as=integer(0), debug=0)}
\description{The main function that treats the case of character and integer inputs.}
\keyword{internal}
\arguments{\item{o}{a vector that defines the items that are observed. Items can be anything that
occurs in sets.}
\item{sets}{a list of sets. Each set is a vector that contains all the items. The vector
can be of any data type, for instance, integers or characters.}
\item{population}{defines the set of item that should be considered for this calculation.}
\item{alpha}{}
\item{beta}{}
\item{p}{}
\item{steps}{defines the number of MCMC steps to be performed.}
\item{restarts}{defines the number of restarts.}
\item{threads}{defines number of threads to be used. Defaults to 0 which means that it
corresponds to the number of available cores.}
\item{as}{define the active sets. No MCMC is run in this case. Just the log likelihood is returned.}}
