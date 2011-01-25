\name{MgsaSets-class}
\alias{MgsaSets-class}
\alias{MgsaSets}
\title{Sets of items and their annotations}
\description{This class describes sets, items and their annotations.}
\details{Internally, the method \code{\link{mgsa}} indexes all elements of the sets before fitting the model.
In case \code{\link{mgsa}} must be run on several observations with the same gene sets, computations can be speeded up by performing this indexing once for all.
This can be achieved by building a \code{\linkS4class{MgsaSets}}.
In order to ensure consistency of the indexing, no replace method for any slot is provided. Accessors are available.

The data frames \code{setAnnotations} and \code{itemAnnotations} allow to store annotations. No constraint is imposed on the number and names of their columns.}
\seealso{\code{\linkS4class{MgsaGoSets}}, \code{\link{readGAF}}, \code{\link{mgsa}}}
\examples{new("MgsaSets", sets=list(set1=c("a", "b"), set2=c("b", "c")))}
\section{Slots}{\describe{\item{\code{sets}:}{(\code{\link{list}}) A list whose elements are vector of item indices.}\item{\code{itemName2ItemIndex}:}{(\code{\link{integer}}) The mapping of item names to index.}\item{\code{numberOfItems}:}{(\code{\link{integer}}) How many items?}\item{\code{setAnnotations}:}{(\code{\link{data.frame}}) Annotations of the sets. The \code{\link{rownames}} are set names.}\item{\code{itemAnnotations}:}{(\code{\link{data.frame}}) Annotations of the items. The \code{\link{rownames}} are item names.}}}
\section{Methods}{\describe{\item{\code{initialize}}{\code{signature( = "MgsaSets")}: Initializes the mapping when some parameters are not specified}\item{\code{\link[=itemAnnotations,MgsaSets,missing-method]{itemAnnotations}}}{}\item{\code{\link[=itemAnnotations,MgsaSets,character-method]{itemAnnotations}}}{}\item{\code{\link[=setAnnotations,MgsaSets,missing-method]{setAnnotations}}}{}\item{\code{\link[=setAnnotations,MgsaSets,character-method]{setAnnotations}}}{}\item{\code{\link[=length,MgsaSets-method]{length}}}{\code{signature( = "MgsaSets")}: Length (number of sets) of \code{\linkS4class{MgsaSets}}.}\item{\code{\link[=itemIndices,MgsaSets,character-method]{itemIndices}}}{}\item{\code{\link[=itemIndices,MgsaSets,numeric-method]{itemIndices}}}{}\item{\code{\link[=show,MgsaSets-method]{show}}}{\code{signature( = "MgsaSets")}: Show an \code{\linkS4class{MgsaSets}}.}\item{\code{\link[=subMgsaSets,MgsaSets,numeric-method]{subMgsaSets}}}{}\item{\code{\link[=mgsa,character,MgsaSets-method]{mgsa}}}{}}}

