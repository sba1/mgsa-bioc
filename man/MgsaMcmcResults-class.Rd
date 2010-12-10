\name{MgsaMcmcResults-class}
\alias{MgsaMcmcResults-class}
\alias{MgsaMcmcResults}
\title{MgsaMcmcResults-class}
\description{Instances of this class are used to hold the additional information 
that was provided by running (possibly multiple times) an MCMC algorithm.}
\seealso{\code{\link{mgsa}}}
\section{Slots}{\describe{\item{\code{steps}:}{(\code{\link{numeric}}) how many steps per MCMC run}\item{\code{restarts}:}{(\code{\link{numeric}}) how many MCMC runs}\item{\code{alphaMcmcPost}:}{(\code{\link{matrix}}) posterior estimates for each MCMC run of the parameter alpha}\item{\code{betaMcmcPost}:}{(\code{\link{matrix}}) posterior estimates for each MCMC run of the parameter beta}\item{\code{pMcmcPost}:}{(\code{\link{matrix}}) posterior estimates for each MCMC run of the parameter p}\item{\code{setsMcmcPost}:}{(\code{\link{matrix}}) posterior estimates for each MCMC run of the sets marginal posterior probabilities

The columns of the matrices \code{alphaMcmcPost}, \code{betaMcmcPost}, \code{pMcmcPost} and setsMcmcPost stores the posterior estimates for each individual MCMC run.
The row order matches the one of the slot \code{alphaPost}, \code{betaPost}, \code{pPots}, and \code{setsResults} respectively.

Accessor methods exist for each slot.}}}
\section{Extends}{\code{\linkS4class{MgsaResults}}}
\section{Methods}{\describe{\item{\code{\link[=steps,MgsaMcmcResults-method]{steps}}}{}\item{\code{\link[=restarts,MgsaMcmcResults-method]{restarts}}}{}\item{\code{\link[=alphaMcmcPost,MgsaMcmcResults-method]{alphaMcmcPost}}}{}\item{\code{\link[=betaMcmcPost,MgsaMcmcResults-method]{betaMcmcPost}}}{}\item{\code{\link[=pMcmcPost,MgsaMcmcResults-method]{pMcmcPost}}}{}\item{\code{\link[=setsMcmcPost,MgsaMcmcResults-method]{setsMcmcPost}}}{}}}

