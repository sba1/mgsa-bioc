\name{MgsaMcmcResults}
\alias{MgsaMcmcResults}
\alias{MgsaMcmcResults-class}
\title{Instances of this class are used to hold the additional information
that was provided by running (possibly multiple times) an MCMC algorithm.}
\description{
  Instances of this class are used to hold the additional
  information that was provided by running (possibly
  multiple times) an MCMC algorithm.
}
\section{Slots}{

\itemize{
  \item{nsamples}{how many samples collected per MCMC run}

  \item{steps}{how many steps per MCMC run}

  \item{restarts}{how many MCMC runs}

  \item{alphaMcmcPost}{posterior estimates for each MCMC
  run of the parameter alpha}

  \item{betaMcmcPost}{posterior estimates for each MCMC run
  of the parameter beta}

  \item{pMcmcPost}{posterior estimates for each MCMC run of
  the parameter p}

  \item{setsMcmcPost}{posterior estimates for each MCMC run
  of the sets marginal posterior probabilities

  The columns of the matrices \code{alphaMcmcPost},
  \code{betaMcmcPost}, \code{pMcmcPost} and setsMcmcPost
  stores the posterior estimates for each individual MCMC
  run. The row order matches the one of the slot
  \code{alphaPost}, \code{betaPost}, \code{pPots}, and
  \code{setsResults} respectively.

  Accessor methods exist for each slot.}
}}
\seealso{
  \code{\link{mgsa}}
}

