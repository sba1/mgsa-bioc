\name{itemAnnotations}
\alias{itemAnnotations}
\alias{itemAnnotations-methods}
\alias{itemAnnotations,MgsaSets,character-method}
\alias{itemAnnotations,MgsaSets,missing-method}
\title{Item annotations of a MgsaSets}
\usage{
  itemAnnotations(sets, items)

  \S4method{itemAnnotations}{MgsaSets,missing}(sets, items)

  \S4method{itemAnnotations}{MgsaSets,character}(sets,
    items)
}
\arguments{
  \item{sets}{an instance of class
  \code{\linkS4class{MgsaSets}}.}

  \item{items}{\code{character} an optional vector
  specifying the items of interest.}
}
\value{
  a \code{data.frame}: the item annotations.
}
\description{
  Item annotations of a \code{\linkS4class{MgsaSets}}.
}

