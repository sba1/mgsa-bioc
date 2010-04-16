#'
#' In this file, we define the MgsaGoSets class and some functions to
#' create usable objects from it.
#' 

#' Set class for GO terms.
#' 
#' Nothing new yet, it comes later
#' 
setClass(
		"MgsaGoSets",
		representation = representation(
		),
		contains = "MgsaSets"
)



#'
#' This functions takes a 1:1 mapping of go.ids to items and returns
#' a full MgsaGOSets instance. The structure of GO is gathered from GO.db. It
#' is sufficient to specify just the directly asserted mapping (or annotation), 
#' i.e., the most specific ones. The true path rule is taken account, that is, if an
#' item is annotated to a term then it will be also annotated to more general
#' terms (some people prefer to say that just the transitive colsure is calculated).
#'
#' @param go.ids a vector of go ids (GO:00001234)
#' @param items a vector of identfiers that are annotated to the term
#'   in the correspondig positon of the go.ids vector.
#'
CreateMgsaGoSets<-function(go.ids,items)
{
	if (length(go.ids) != length(items))
	{
		stop("Arguments go.ids and items differ in length.")
	}
	
	require(GO.db)
	require(RSQLite)
	
	# Prepare the data base stuff
	drv <- dbDriver("SQLite")
	annotation.file <- tempfile()
	annotation.con <- dbConnect(drv, dbname = annotation.file)
	dbWriteTable(annotation.con,"ga",data.frame(go.id=go.ids,items=items),row.names=0)
	
	# We now attach the GO Database
	attachSQL = paste("ATTACH '", GO_dbfile(), "' AS goDB;", sep = "")
	dbGetQuery(annotation.con, attachSQL)
	
	# We now make our call
	# Basically, we query terms of which annotated term is a offspring.
	# All those terms are also annotated to the gene. We union all sub ontologies.
	# We also need to consider the direct annotations (hence the 4th SELECT statement). 
	all<-dbGetQuery(annotation.con, paste("SELECT DISTINCT ga.items AS items,go_bp_offspring._id AS id, go2.go_id AS go_id",
					"FROM ga,goDB.go_term,goDB.go_bp_offspring,goDB.go_term as go2",
					"WHERE ga.go_id = goDB.go_term.go_id AND goDB.go_term._id = goDB.go_bp_offspring._offspring_id AND goDB.go_bp_offspring._id = go2._id",
					"UNION SELECT DISTINCT ga.items AS items,go_cc_offspring._id AS id, go2.go_id AS go_id",
					"FROM ga,goDB.go_term,goDB.go_cc_offspring,goDB.go_term as go2",
					"WHERE ga.go_id = goDB.go_term.go_id AND goDB.go_term._id = goDB.go_cc_offspring._offspring_id AND goDB.go_cc_offspring._id = go2._id",
					"UNION SELECT DISTINCT ga.items AS items,go_mf_offspring._id AS id, go2.go_id AS go_id",
					"FROM ga,goDB.go_term,goDB.go_mf_offspring,goDB.go_term as go2",
					"WHERE ga.go_id = goDB.go_term.go_id AND goDB.go_term._id = goDB.go_mf_offspring._offspring_id AND goDB.go_mf_offspring._id = go2._id",
					"UNION SELECT DISTINCT ga.items AS items,goDB.go_term._id,goDB.go_term.go_id AS go_id FROM ga,goDB.go_term WHERE ga.go_id = goDB.go_term.go_id"))
	
	
	term.anno.db<-dbGetQuery(annotation.con,"SELECT go_id,term,definition FROM goDB.go_term");
	
	# Cleanup
	dbGetQuery(annotation.con, "DETACH goDB" )
	dbDisconnect(annotation.con)
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

#'
#' Creates a MgsaGoSets using annotation provided by a file in GAF 1.0
#' format. The file specified by the filename may be gzip-compressed.
#'
#' TODO:  provide support for evidence codes, choose a better name
#'
CreateMgsaGoSetsFromGAF<-function(filename, gene.id.col = 3, go.id.col = 5, evidence.col =  7, name.col = 10)
{
	goa = read.delim(gzfile(filename), na.strings = "", header=F, comment.char = "!", sep="\t")
	
	goa = na.omit( 
			data.frame ( 
					go.ids = goa[,go.id.col],
					gene.ids = goa[, gene.id.col],
					evidence = goa[, evidence.col],
					name = goa[, name.col]
			)
	)
	
	sets<-CreateMgsaGoSets(go.ids=goa$go.ids,items=goa$gene.ids)
	item.annot <- unique(goa[, c("gene.ids", "name")])
	
	if (any(duplicated(item.annot$gene.ids))) stop("At least one DB object ID has multiple DB object names in gene ontology annotation file.")
	
	sets@itemAnnotations <- data.frame(row.names = item.annot[, "gene.ids"], name=as.character(item.annot[, "name"]))
	sets@itemAnnotations <- sets@itemAnnotations[names(sets@itemName2ItemIndex),,drop=FALSE]
	
	return(sets)
}
