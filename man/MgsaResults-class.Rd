\name{MgsaResults-class}
\alias{MgsaResults-class}
\alias{MgsaResults}
\title{Results of an MGSA analysis}
\description{The results of an MGSA analysis.}
\seealso{\code{\link{mgsa}}}
\section{Slots}{\describe{\item{\code{populationSize}:}{(\code{\link{numeric}}) The number of items in the population.}\item{\code{studySetSizeInPopulation}:}{(\code{\link{numeric}}) The number of items both in the study set and in the population.}\item{\code{alphaPost}:}{(\code{\link{data.frame}}) with columns \code{value}, \code{estimate} and \code{std.error}.}\item{\code{betaPost}:}{(\code{\link{data.frame}}) with columns \code{value}, \code{estimate} and \code{std.error}.}\item{\code{pPost}:}{(\code{\link{data.frame}}) with columns \code{value}, \code{estimate} and \code{std.error}.}\item{\code{setsResults}:}{(\code{\link{data.frame}}) with columns \code{inPopulation}, \code{inStudySet}, \code{estimate} and \code{std.error}. 

The columns of the slots \code{alphaPost}, \code{betaPost}, and \code{pPost} contains a realization value, its posterior estimate and standard error for the parameters alpha, beta and p respectively.

The columns of the slot \code{setsResults} contains the number of items of the set in the population, the number of items of the set in the study set, the estimate of its marginal posterior probability and its standard error.
The \code{\link{rownames}} are the names of the sets if available.

Because an \code{MgsaResults} is the outcome of an MGSA analysis (see \code{\link{mgsa}}), accessors but no replacement methods exist for each slot.}}}
\section{Methods}{\describe{\item{\code{\link[=populationSize,MgsaResults-method]{populationSize}}}{}\item{\code{\link[=studySetSizeInPopulation,MgsaResults-method]{studySetSizeInPopulation}}}{}\item{\code{\link[=alphaPost,MgsaResults-method]{alphaPost}}}{}\item{\code{\link[=betaPost,MgsaResults-method]{betaPost}}}{}\item{\code{\link[=pPost,MgsaResults-method]{pPost}}}{}\item{\code{\link[=setsResults,MgsaResults-method]{setsResults}}}{}\item{\code{\link[=show,MgsaResults-method]{show}}}{\code{signature( = "MgsaResults")}: Show an \code{\linkS4class{MgsaResults}}.}\item{\code{\link[=plot,MgsaResults-method]{plot}}}{\code{signature( = "MgsaResults")}: Plot method for MgsaResults objects}}}

