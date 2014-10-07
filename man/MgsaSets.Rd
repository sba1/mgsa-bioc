\name{MgsaSets}
\alias{MgsaSets}
\alias{MgsaSets-class}
\title{Sets of items and their annotations}
\description{
  This class describes sets, items and their annotations.
}
\details{
  Internally, the method \code{\link{mgsa}} indexes all
  elements of the sets before fitting the model. In case
  \code{\link{mgsa}} must be run on several observations
  with the same gene sets, computations can be speeded up
  by performing this indexing once for all. This can be
  achieved by building a \code{\linkS4class{MgsaSets}}. In
  order to ensure consistency of the indexing, no replace
  method for any slot is provided. Accessors are available.

  The data frames \code{setAnnotations} and
  \code{itemAnnotations} allow to store annotations. No
  constraint is imposed on the number and names of their
  columns.
}
\section{Slots}{

\itemize{
  \item{sets}{A list whose elements are vector of item
  indices.}

  \item{itemName2ItemIndex}{The mapping of item names to
  index.}

  \item{numberOfItems}{How many items?}

  \item{setAnnotations}{Annotations of the sets. The
  \code{\link{rownames}} are set names.}

  \item{itemAnnotations}{Annotations of the items. The
  \code{\link{rownames}} are item names.}
}}
\examples{
new("MgsaSets", sets=list(set1=c("a", "b"), set2=c("b", "c")))
}
\seealso{
  \code{\linkS4class{MgsaGoSets}}, \code{\link{readGAF}},
  \code{\link{mgsa}}
}

