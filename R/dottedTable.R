#' @nord
dottedString <- function(x, width=30){
	nc = sapply(x, nchar)
	ifelse( nc<=width, x, paste(substr(x,1,width), "...", sep="") )
}

#' @nord
dottedTable <- function(tab, nrows=5 , ncols=2, width=30 ){
	x <- tab
	if(nrows <  nrow(tab) ) x <- x[1:nrows,,drop=FALSE]
	if(ncols <  ncol(tab) ) x <- x[,1:ncols,drop=FALSE]
	
	for(i in seq(along=x)){
		if(is.character(x[[i]]) | is.factor(x[[i]]))
			x[[i]] <- dottedString( as.character(x[[i]]), width=width )
	}
	
	if(ncol(tab)>ncols){
		x[[i+1]] <- rep("",nrow(x))
		colnames(x)[i+1] <- "..."
	}
	x
}
