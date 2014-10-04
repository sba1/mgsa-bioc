#' @include dottedTable.R
NULL

#'
#' This class describes sets, items and their annotations.
#' 
#' Internally, the method \code{\link{mgsa}} indexes all elements of the sets before fitting the model.
#' In case \code{\link{mgsa}} must be run on several observations with the same gene sets, computations can be speeded up by performing this indexing once for all.
#' This can be achieved by building a \code{\linkS4class{MgsaSets}}.
#' In order to ensure consistency of the indexing, no replace method for any slot is provided. Accessors are available.
#' 
#' The data frames \code{setAnnotations} and \code{itemAnnotations} allow to store annotations. No constraint is imposed on the number and names of their columns. 
#' 
#' @slot sets A list whose elements are vector of item indices.
#' @slot itemName2ItemIndex The mapping of item names to index.
#' @slot numberOfItems How many items?
#' @slot setAnnotations Annotations of the sets. The \code{\link{rownames}} are set names. 
#' @slot itemAnnotations Annotations of the items. The \code{\link{rownames}} are item names.
#' @title Sets of items and their annotations
#' @examples new("MgsaSets", sets=list(set1=c("a", "b"), set2=c("b", "c")))
#' @seealso \code{\linkS4class{MgsaGoSets}}, \code{\link{readGAF}}, \code{\link{mgsa}}
#' @exportClass MgsaSets

setClass(
		Class="MgsaSets",
		representation = representation(
				# Sets also have attribute names(). Elements are
				# just vectors of item indices.
				sets = "list",
				
				# maps any item name to an item index
				itemName2ItemIndex = "integer",
				
				# total number of items 
				numberOfItems = "integer",
				
				# Same order as sets
				setAnnotations = "data.frame",
				
				# Rows are ordered according to item indices and contain the item annotations
				itemAnnotations = "data.frame"
		)
)


#' Initializes the mapping when some parameters are not specified
#' 
#' @keywords internal
#' @noRd

setMethod(f = "initialize",
		signature = c("MgsaSets"),
		def = function(.Object, ...) {
			.Object <- callNextMethod()
			
			# Needs numberOfItems to be initialized?
			if (length(.Object@numberOfItems) == 0)
			{
				# Needs itemName2ItemIndex to be initialized?
				if (length(.Object@itemName2ItemIndex) == 0)
				{
					# But do this only, if there is actually a set
					if (length(.Object@sets) != 0)
					{
						sets<-.Object@sets
						v<-factor(unlist(sets))
						
						# extract names
						names<-as.vector(levels(v))
						
						# start constructing the map vector
						itemName2ItemIndex <- 1:length(names)
						
						# rename
						levels(v)<-itemName2ItemIndex
						
						# relist using old sets as skeleton  
						sets<-relist(as.integer(v),sets)
						
						# create map
						names(itemName2ItemIndex) <- names
						
						# Assign
						.Object@sets <- sets
						.Object@itemName2ItemIndex<-itemName2ItemIndex
						
						if (length(.Object@itemAnnotations)==0)
						{
							.Object@itemAnnotations<-data.frame(row.names=names(itemName2ItemIndex))
						}
					}
				}
				
				if (length(.Object@itemName2ItemIndex) != 0)
				{
					.Object@numberOfItems <- max(.Object@itemName2ItemIndex)
				}
			}
			
			if (length(.Object@numberOfItems) != 0)
			{
				if (any(.Object@itemName2ItemIndex != (1:.Object@numberOfItems)))
				{
					stop("Provided itemName2ItemIndex should be equal to 1:numberOfItems for now.");
				}
			}
			
			if (any(duplicated(names(.Object@itemName2ItemIndex))))
			{
				stop("Provided itemName2ItemIndex should not contain duplicated names.");
			}
			
			if (length(.Object@setAnnotations)==0)
			{
				.Object@setAnnotations<-data.frame(row.names=names(.Object@sets))
			}
			
			if (length(.Object@itemAnnotations)==0)
			{
				.Object@itemAnnotations<-data.frame(row.names=names(.Object@itemName2ItemIndex))
			}
			
			.Object
		})

#' Item annotations of a \code{\linkS4class{MgsaSets}}.
#' 
#' @title Item annotations of a MgsaSets
#' @param sets an instance of class \code{\linkS4class{MgsaSets}}.
#' @param items \code{character} an optional vector specifying the items of interest. 
#' @return a \code{data.frame}: the item annotations.
#' @exportMethod itemAnnotations
#' @rdname itemAnnotations-methods
setGeneric( "itemAnnotations", function(sets,items) standardGeneric( "itemAnnotations" ) )

#' @rdname itemAnnotations-methods
setMethod(
		f="itemAnnotations",
		signature = c( "MgsaSets", "missing" ),
		function( sets, items ) sets@itemAnnotations
)

#' @rdname itemAnnotations-methods
setMethod(
		f="itemAnnotations",
		signature=c( "MgsaSets","character" ),
		function( sets, items ) sets@itemAnnotations[match(items,row.names(sets@itemAnnotations)),,drop=FALSE]
)

#' Set annotations of a \code{\linkS4class{MgsaSets}}.
#' 
#' @title Set annotations of a MgsaSets
#' @param sets an instance of class \code{\linkS4class{MgsaSets}}.
#' @param names \code{character} an optional vector specifying the names of interest. 
#' @return a \code{data.frame}: the set annotations.
#' @rdname setAnnotations-methods
#' @exportMethod setAnnotations

setGeneric( "setAnnotations", function(sets,names) standardGeneric( "setAnnotations" ) )

#' @rdname setAnnotations-methods
setMethod(
		f="setAnnotations",
		signature=c( "MgsaSets", "missing" ),
		function( sets, names ) sets@setAnnotations
)

#' @rdname setAnnotations-methods
setMethod(
		f="setAnnotations",
		signature=c( "MgsaSets", "character" ),
		function( sets, names ) sets@setAnnotations[match(names,row.names(sets@setAnnotations)),,drop=FALSE]
)

#' Length (number of sets) of \code{\linkS4class{MgsaSets}}.
#' 
#' @title Length of a MgsaSets. 
#' @param x an instance of class \code{\linkS4class{MgsaSets}}.
#' @return \code{integer} vector. 
#' @rdname length-methods
#' @exportMethod length

setMethod(
		f="length",
		signature=c("MgsaSets"),
		function(x) length(x@sets)
)


#' Returns the indices corresponding to the items
#'
#' @title Item indices of a MgsaSets
#' @param sets an instance of class \code{\linkS4class{MgsaSets}}.
#' @param items \code{character} or \code{numeric} the items of interest. 
#' @return a \code{integer}: the item indices.
#' @rdname itemIndices-methods
#' @exportMethod itemIndices

setGeneric("itemIndices", function(sets, items) standardGeneric("itemIndices"))

#' @rdname itemIndices-methods
setMethod(
		f="itemIndices",
		signature=c( "MgsaSets","character" ),
		function( sets, items ) sets@itemName2ItemIndex[items]
)

#' @rdname itemIndices-methods
setMethod(
		f="itemIndices",
		signature=c( "MgsaSets", "numeric" ),
		function( sets, items ) sets@itemName2ItemIndex[items]
)


#' Show an \code{\linkS4class{MgsaSets}}.
#'
#' @title Show an MgsaSets
#' @param object an instance of class \code{\linkS4class{MgsaSets}}.
#' @return an invisible \code{NULL}
#' @exportMethod show

setMethod(
		f="show",
		signature=c( "MgsaSets" ),
		function( object ){
			cat(
					"Object of class ",
					class( object ),
					"\n",
					length(object@sets),
					" sets over ",
					object@numberOfItems,
					" unique items.\n",
					sep = ""
			)
			
			cat("\nSet annotations:\n")
			nrowShow <- min (5 , nrow(object@setAnnotations) )
			print( dottedTable(object@setAnnotations, nrows=nrowShow) )
			
			if(nrowShow < nrow(object@setAnnotations) ){
				cat("... and ", nrow(object@setAnnotations) - nrowShow, " other sets.\n" )
			}
			
			cat("\nItem annotations:\n")
			nrowShow <- min (5 , nrow(object@itemAnnotations) )
			print( dottedTable(object@itemAnnotations, nrows=nrowShow) )
			
			if(nrowShow < nrow(object@itemAnnotations) ){
				cat("... and ", nrow(object@itemAnnotations) - nrowShow, " other items.\n" )
			}
			
		}
)


#' Returns a subset of an \code{\linkS4class{MgsaSets}} that contains
#' only the specified items. Empty sets are removed.
#' 
#' @title Subset of an MgsaSets
#' @param sets an \code{\linkS4class{MgsaSets}}.
#' @param items \code{character}. The items to restrict on.
#' @return an \code{\linkS4class{MgsaSets}}.
#' @rdname subMgsaSets-methods
#' @exportMethod subMgsaSets

setGeneric("subMgsaSets", function(sets, items) standardGeneric("subMgsaSets"))

#' @rdname subMgsaSets-methods
setMethod(
		f="subMgsaSets",
		signature=c( "MgsaSets", "character" ),
		function( sets, items )
		{
			# copy the original 'sets' object
			rv <- sets
			# subset the item indices map, keeping the old item indices
			oldindices <- sets@itemName2ItemIndex[items]
			oldindices <- oldindices[!is.na(oldindices)]
			# generate the new indices map
			rv@itemName2ItemIndex <- seq_along(oldindices)
			names(rv@itemName2ItemIndex) <- names(oldindices)
			# subset existing item annotations
			rv@itemAnnotations <- rv@itemAnnotations[rownames(rv@itemAnnotations) %in% items,, drop=FALSE]
			# update items count
			rv@numberOfItems <- length(oldindices)
			# intersect each set with the new 'items' and update the indices of set items
			encode <- function(x){ match(intersect(x, oldindices), oldindices) }
			rv@sets <- lapply(rv@sets, encode)
			# discard empty sets
			rv@sets <- rv@sets[sapply(rv@sets,length)>0]
			# subset set annotations
			rv@setAnnotations <- rv@setAnnotations[rownames(rv@setAnnotations) %in% names(rv@sets),, drop=FALSE]
			return (rv)
		})
