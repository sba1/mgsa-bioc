\name{mcmcSummary}
\alias{mcmcSummary}
\title{mcmcSummary}
\usage{mcmcSummary(x)}
\description{Wraps the mgsa.trampoline function. It is dumb with respect to input parameter
but returns a processed result.}
\value{an object of class MgsaMcmcResults.}
\arguments{\item{o}{which items are observed. Items are represented by integers}
\item{sets}{list of sets. It is expected that each set contains integers
that correspond to item to which the set is associated}
\item{n}{total number of items.}
\item{alpha}{}
\item{beta}{}
\item{p}{}
\item{steps}{defines the number of MCMC steps to be performed}
\item{restarts}{defines the number of restarts.}
\item{threads}{defines number of threads to be used. Defaults to 0 which means that it
corresponds to the number of available cores.}
\item{specifies}{the debug level. Mainly for internal use.}}
