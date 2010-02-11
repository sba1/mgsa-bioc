#'
#' This functions takes a mapping of go.ids and items and returns
#' a list containing a list ofset -> item mappings suitable for mgsa
#' and a mapping function which can be used to get the index of a item
#' name
#'
mgsa.make.go.mapping<-function(go.ids,items)
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

	# Cleanup
	dbGetQuery(annotation.con, "DETACH goDB" )
	dbDisconnect(annotation.con)

	# Map to unique gene ids
	all.items<-factor(all$items)
	all.items.names<-as.vector(levels(all.items))
	levels(all.items)<-1:length(unique(all.items))

	sets<-split(as.integer(all.items),all$id)
	map.vec<-1:length(all.items.names)
	names(map.vec)<-all.items.names

	mapping<-new("MgsaGoMapping");
	mapping@sets<-sets
	mapping@item.idx.map<-map.vec
	
	return(mapping)
}

#'
#' Makes a mapping using a given goa file (as can be downloaded from Gene Ontology).
#' The file may be gzip-compressed.
#'
#' TODO:  provide support for evidence codes
#'
mgsa.make.go.mapping.from.goa<-function(filename, gene.id.col = 3, go.id.col = 5, evidence.col =  7)
{
    goa = read.delim(gzfile(filename), na.strings = "", header=F, comment.char = "!", sep="\t")
    
    goa = na.omit( 
            data.frame ( 
                    go.ids = goa[,go.id.col],
                    gene.ids = goa[, gene.id.col],
                    evidence = goa[, evidence.col]							
            )
    )

	return(mgsa.make.go.mapping(go.ids=goa$go.ids,items=goa$gene.ids))
}
