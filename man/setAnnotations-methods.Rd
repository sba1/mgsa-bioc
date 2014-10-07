\name{setAnnotations}
\alias{setAnnotations}
\alias{setAnnotations-methods}
\alias{setAnnotations,MgsaSets,character-method}
\alias{setAnnotations,MgsaSets,missing-method}
\title{Set annotations of a MgsaSets}
\usage{
  setAnnotations(sets, names)

  \S4method{setAnnotations}{MgsaSets,missing}(sets, names)

  \S4method{setAnnotations}{MgsaSets,character}(sets,
    names)
}
\arguments{
  \item{sets}{an instance of class
  \code{\linkS4class{MgsaSets}}.}

  \item{names}{\code{character} an optional vector
  specifying the names of interest.}
}
\value{
  a \code{data.frame}: the set annotations.
}
\description{
  Set annotations of a \code{\linkS4class{MgsaSets}}.
}

