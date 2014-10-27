#
# In this file, we define the MgsaGoSets class and some functions to
# create usable objects from it.
# 
#' @include MgsaSets-class.R
NULL

######## Class definitions
#### MgsaGoSets

#' This class represents gene ontology annotations.
#' 
#' For now, it is identical to the parental class \code{\linkS4class{MgsaSets}}.
#' @title Gene Ontology annotations
#' @seealso \code{\link{readGAF}}
#' @exportClass MgsaGoSets
setClass(
		"MgsaGoSets",
		contains = c("MgsaSets"),
		representation = representation(
		)
)


#' This functions takes a 1:1 mapping of go.ids to items and returns
#' a full MgsaGOSets instance. The structure of GO is gathered from GO.db. It
#' is sufficient to specify just the directly asserted mapping (or annotation), 
#' i.e., the most specific ones. The true path rule is taken account, that is, if an
#' item is annotated to a term then it will be also annotated to more general
#' terms (some people prefer to say that just the transitive closure is calculated).
#'
#' @param go.ids a character vector of GO ids (GO:00001234)
#' @param items a vector of identifiers that are annotated to the term
#'   in the corresponding position of the go.ids vector.
#' @export readGAF 

createMgsaGoSets<-function(go.ids,items)
{
	if (length(go.ids) != length(items))
	{
		stop("Arguments go.ids and items differ in length.")
	}
	
	require(GO.db)
	require(RSQLite)
	require(DBI)
	
	# Prepare the data base stuff
	drv <- DBI::dbDriver("SQLite")
	annotation.file <- tempfile()
	annotation.con <- DBI::dbConnect(drv, dbname = annotation.file)
	DBI::dbWriteTable(annotation.con,"ga",data.frame(go_id=go.ids,items=items),row.names=0)
	
	# We now attach the GO Database
	attachSQL = paste("ATTACH '", GO.db::GO_dbfile(), "' AS goDB;", sep = "")
	DBI::dbGetQuery(annotation.con, attachSQL)
	
	# We now make our call
	# Basically, we query terms of which annotated term is a offspring.
	# All those terms are also annotated to the gene. We union all sub ontologies.
	# We also need to consider the direct annotations (hence the 4th SELECT statement). 
	all<-DBI::dbGetQuery(annotation.con, paste("SELECT DISTINCT ga.items AS items,go_bp_offspring._id AS id, go2.go_id AS go_id",
					"FROM ga,goDB.go_term,goDB.go_bp_offspring,goDB.go_term as go2",
					"WHERE ga.go_id = goDB.go_term.go_id AND goDB.go_term._id = goDB.go_bp_offspring._offspring_id AND goDB.go_bp_offspring._id = go2._id",
					"UNION SELECT DISTINCT ga.items AS items,go_cc_offspring._id AS id, go2.go_id AS go_id",
					"FROM ga,goDB.go_term,goDB.go_cc_offspring,goDB.go_term as go2",
					"WHERE ga.go_id = goDB.go_term.go_id AND goDB.go_term._id = goDB.go_cc_offspring._offspring_id AND goDB.go_cc_offspring._id = go2._id",
					"UNION SELECT DISTINCT ga.items AS items,go_mf_offspring._id AS id, go2.go_id AS go_id",
					"FROM ga,goDB.go_term,goDB.go_mf_offspring,goDB.go_term as go2",
					"WHERE ga.go_id = goDB.go_term.go_id AND goDB.go_term._id = goDB.go_mf_offspring._offspring_id AND goDB.go_mf_offspring._id = go2._id",
					"UNION SELECT DISTINCT ga.items AS items,goDB.go_term._id,goDB.go_term.go_id AS go_id FROM ga,goDB.go_term WHERE ga.go_id = goDB.go_term.go_id"))
	
	
	term.anno.db<-DBI::dbGetQuery(annotation.con,"SELECT go_id,term,definition FROM goDB.go_term");
	
	# Cleanup
	DBI::dbGetQuery(annotation.con, "DETACH goDB" )
	DBI::dbDisconnect(annotation.con)
	unlink(annotation.file)
	
	# Map to unique gene ids
	all.items<-factor(all$items)
	all.items.names<-as.vector(levels(all.items))
	levels(all.items)<-1:length(unique(all.items))
	
	sets<-split(as.integer(all.items),all$go_id)
	itemName2ItemIndex<-1:length(all.items.names)
	names(itemName2ItemIndex)<-all.items.names
	
	term.anno<-data.frame(row.names=term.anno.db[,1],term.anno.db[,-1])
	term.anno<-term.anno[names(sets),]
	
	mapping<-new("MgsaGoSets",sets=sets,itemName2ItemIndex=itemName2ItemIndex,setAnnotations=term.anno);
	
	return(mapping)
}

#' Creates a MgsaGoSets using gene ontology annotations provided by a file in GAF 1.0 or 2.0
#' format.
#' 
#' The function extracts from the annotation file all direct gene annotations and infers from the Gene Ontology all the indirect annotations (due to term relationships).
#' This is done using the package \code{Go.db} which provides the ontology as a database and \code{RSQLite} for querying the database. 
#' @title Read a Gene Ontology annotation file
#' @usage readGAF(filename, evidence=NULL, aspect=c("P", "F", "C"))
#' @param filename The name of the Gene Ontology annotation file. It must be in the GAF 1.0 or 2.0 format. It may be gzip-compressed.
#' @param evidence \code{character} or \code{NULL}. Only annotations with evidence code in \code{evidence} are returned. If \code{NULL} (default), annotations of all evidence codes are returned.   
#' @param aspect \code{character} with values in P, C or F. Only annotations of the listed GO namespaces P (biological process), F (molecular function) or C (cellular component) are returned. By default, annotations of the three namespaces are returned.   
#' @return An \code{\linkS4class{MgsaGoSets}} object.
#' @seealso \code{\linkS4class{MgsaGoSets}}, \code{\link{mgsa}}
#' @references The Gene Ontology Consortium. Gene Ontology: tool for the unification of biology. Nature Genetics, 2000.
#' The GAF file format: \url{http://www.geneontology.org/GO.format.annotation.shtml}
#' GO evidence codes: \url{http://www.geneontology.org/GO.evidence.shtml}
#' @examples ## parsing provided example file (yeast)
#' gofile = system.file("example_files/gene_association_head.sgd", package="mgsa")
#' readGAF(gofile)
#' ## only annoations infered from experiment or a direct assay
#' readGAF(gofile, evidence=c("EXP", "IDA"))
#' @export readGAF 

readGAF = function(filename, evidence=NULL, aspect=c("P", "F", "C")){
	
	## the column IDs of interest according to GAF 1.0 and 2.0
	gene.id.col = 2
	symbol.col = 3
	go.id.col = 5
	evidence.col = 7
	aspect.col = 9
	name.col = 10
	
	## validity of parameters
	if( !( is.null(evidence) | is.character(evidence) ) )
		stop("evidence must be NULL or a character vector.")
	
	if(!all(aspect %in% c("P", "F", "C")))
		stop("aspect must be a character vector with all entries in c(\"P\", \"F\", \"C\")")
	
	## reading the file
	goa = read.delim(gzfile(filename), na.strings = "", header=FALSE, comment.char = "!", sep="\t")
	
	goa = na.omit( 
			data.frame ( 
					go.ids = goa[,go.id.col],
					gene.ids = goa[, gene.id.col],
					evidence.code = goa[, evidence.col],
					symbol = goa[, symbol.col],
					name = goa[, name.col],
					aspect.code = goa[,aspect.col]
			)
	)
	
	goa <- goa[goa$aspect.code %in% aspect, ]
	
	if(!is.null(evidence))
		goa <- goa[goa$evidence.code %in% evidence, ]
	
	if(nrow(goa)==0){
		warning("No genes with annotations. Are the evidence codes too restrictive?")
		return( new("MgsaGoSets") )
	}
	
	sets = createMgsaGoSets(go.ids=goa$go.ids, items=goa$gene.ids)
	item.annot <- unique(goa[, c("gene.ids", "symbol","name")])
	
	if (any(duplicated(item.annot$gene.ids))) stop("At least one DB object ID has multiple DB object symbols or names in gene ontology annotation file.")
	
	sets@itemAnnotations <- data.frame(
			row.names = item.annot[, "gene.ids"],
			symbol = as.character(item.annot[, "symbol"]),
			name = as.character(item.annot[, "name"])
	)
	sets@itemAnnotations <- sets@itemAnnotations[names(sets@itemName2ItemIndex),,drop=FALSE]
	
	return(sets)
}
