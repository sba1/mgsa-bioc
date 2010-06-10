dottedString <- function(x, width=20){
	nc = sapply(x, nchar)
	ifelse( nc<=width, x, paste(substr(x,1,width), "...", sep="") )
}

dottedTable <- function(tab, nrows=5 , ncols=2, width=20 ){
	nrowShow <- min(nrows, nrow(tab))
	ncolShow <- min(ncols, ncol(tab))
	x <- tab[1:nrowShow,1:ncolShow,drop=FALSE]
	for(i in 1:ncol(x)) x[[i]] <- dottedString( as.character(x[[i]]), width=width )
	if(ncol(tab)>ncols){
		x[[i+1]] <- rep("",nrowShow)
		colnames(x)[i+1] <- "..."
	}
	x
}
