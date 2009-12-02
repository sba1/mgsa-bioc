getGeneID_SGD = function (x){
    unlist( lapply ( strsplit(as.character(x[,11]), "\\|"), function (y) y[1] ) )
}


make_GOA_matrices = function(file, getGeneId = function (x) x[,3], goIdCol = 5, catCol = 9 , evidenceCol =  7){
    require(GO.db)
    require(Matrix)
    
    ## goterm: gene ontology data frame (terms themselves without annotations)
    goterm = data.frame(GOID=unlist(eapply(GOTERM, function(x) x@GOID)), Term=unlist(eapply(GOTERM, function(x) x@Term)), Ont=unlist(eapply(GOTERM, function(x) x@Ontology))) 
    goterm = na.omit(goterm)
    
    ## goa: gene ontology annotation data frame
    goa = read.delim(file, na.strings = "", header=F, comment.char = "!", sep="\t")
    
    goa = na.omit( 
            data.frame ( 
                    GOID = goa[,goIdCol],
                    GeneID = getGeneId( goa ),
                    GOCAT = goa[,catCol],
                    evidence = goa[, evidenceCol]							
            )
    )
    
    ## all genes with at least one annotation
    genes = levels( goa$GeneID )
    names(genes) = levels( goa$GeneID )
    
    ## initalize result list
    res = list()
    for(gocat in c("MF","BP","CC")) {
        if(gocat=="MF") {
            go_offspr_list = as.list(GOMFOFFSPRING)
            sub_go_org = goa[goa$GOCAT=="F",] 
        }
        if(gocat=="BP") {
            go_offspr_list = as.list(GOBPOFFSPRING)
            sub_go_org = goa[goa$GOCAT=="P",] 
        }
        if(gocat=="CC") {
            go_offspr_list = as.list(GOCCOFFSPRING)
            sub_go_org = goa[goa$GOCAT=="C",] 
        }
        
        # clean-up step for the list
        go_offspr_list = lapply(go_offspr_list, function(x) unlist( as.vector(x) ) )
        
        # include list component (GOID) names in corresponding (GOID) vectors
        go_offspr_list_temp = lapply(names(go_offspr_list), function(x) c(x, go_offspr_list[[x]]) ) 
        names(go_offspr_list_temp) = names(go_offspr_list) 
        go_offspr_list = go_offspr_list_temp
        
        # remove NAs in vectors
        go_offspr_list = lapply( go_offspr_list, function(x) x[!is.na(x)] ) 
        
        # retrieve gene IDs for GOID vectors
        gene_annot_list = lapply( go_offspr_list, function(x) unique(as.vector( sub_go_org$GeneID[ sub_go_org$GOID %in% x ] ) ) ) 
        
        ## turn those list as matrix: gene 2 term matrix
        ## the matrix is sparse
        g2t = sparseMatrix(
                i = match( unlist( gene_annot_list ), genes), 
                j = rep( 1:length(gene_annot_list), times = sapply( gene_annot_list, length ) ),
                x = 1,
                dims = c( length(genes), length( gene_annot_list) )
        )
        
        rownames( g2t ) = names( genes )
        colnames( g2t ) = names( gene_annot_list )
        
        ## remove empty columns and empty rows
        g2t = g2t[rowSums(g2t)>0, colSums(g2t)>0]
        
        ## stores everything in results
        res[[gocat]] = list()
        res[[gocat]]$g2t = g2t
        
        res[[gocat]]$gene = genes[ rownames(g2t) ]
        
        res[[gocat]]$term = as.character( goterm$Term[ match( colnames(g2t), goterm$GOID) ] )
        names(res[[gocat]]$term) = colnames(g2t)
        
        
    }
    res
}
