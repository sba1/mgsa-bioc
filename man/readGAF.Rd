\name{readGAF}
\alias{readGAF}
\title{Read a Gene Ontology annotation file}
\usage{readGAF(filename, evidence=NULL, aspect=c("P", "F", "C"))}
\description{Creates a MgsaGoSets using gene ontology annotations provided by a file in GAF 1.0 or 2.0
format.}
\details{The function extracts from the annotation file all direct gene annotations and infers from the Gene Ontology all the indirect annotations (due to term relationships).
This is done using the package \code{Go.db} which provides the ontology as a database and \code{RSQLite} for querying the database.}
\value{An \code{\linkS4class{MgsaGoSets}} object.}
\seealso{\code{\linkS4class{MgsaGoSets}}, \code{\link{mgsa}}}
\references{The Gene Ontology Consortium. Gene Ontology: tool for the unification of biology. Nature Genetics, 2000.
The GAF file format: \url{http://www.geneontology.org/GO.format.annotation.shtml}
GO evidence codes: \url{http://www.geneontology.org/GO.evidence.shtml}}
\arguments{\item{filename}{The name of the Gene Ontology annotation file. It must be in the GAF 1.0 or 2.0 format. It may be gzip-compressed.}
\item{evidence}{\code{character} or \code{NULL}. Only annotations with evidence code in \code{evidence} are returned. If \code{NULL} (default), annotations of all evidence codes are returned.}
\item{aspect}{\code{character} with values in P, C or F. Only annotations of the listed GO namespaces P (biological process), F (molecular function) or C (cellular component) are returned. By default, annotations of the three namespaces are returned.}
}
\examples{## parsing provided example file (yeast)
gofile = system.file("example_files/gene_association_head.sgd", package="mgsa")
readGAF(gofile)
## only annoations infered from experiment or a direct assay
readGAF(gofile, evidence=c("EXP", "IDA"))}

