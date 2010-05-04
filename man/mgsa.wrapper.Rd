\name{mgsa.wrapper}
\alias{mgsa.wrapper}
\title{mgsa.wrapper}
\usage{mgsa.wrapper(o, sets, n, alpha=NA, beta=NA, p=NA, steps=1e+06,
    restarts=1, threads=0, as=integer(0), debug=0)}
\description{Wraps the mgsa.trampoline function. It is dumb with respect to input parameter
but returns a processed result.}
\value{an object of class MgsaMcmcResults.}
\keyword{internal}
\arguments{\item{o}{a vector, which defines the items that are observed. Items are expressed using
integer values between 1 and \code{n}.}
\item{sets}{list of sets. It is expected that each set contains integers in range between
1 and \code{n} that correspond to items to which the set is associated.}
\item{n}{defines the total number of items.}
\item{alpha}{}
\item{beta}{}
\item{p}{}
\item{steps}{defines the number of MCMC steps to be performed}
\item{restarts}{defines the number of restarts.}
\item{threads}{defines number of threads to be used. Defaults to 0 which means that it
corresponds to the number of available cores.}
\item{specifies}{the debug level. Mainly for internal use.}}
