
#' Class that describes sets and their associations
#' The attributes of this class are completely private.
#' 
#' TODO: Add proper accessor functions and find out how to hide
#' the raw attributes.
setClass(
        "MgsaMapping",
        representation = representation(
							# Sets also have attribute names(). Elements are
							# just vectors of item indices.
							sets = "list",

							# maps any item name to an item index
							itemName2ItemIndex = "integer",

							# total number of items 
							numberOfItems = "integer"
			)
)
#' Initialization
#' 
#' This intializes the mapping when some of the parameter are not specifed

setMethod("initialize", "MgsaMapping",
					function(.Object, ...) {
						.Object <- callNextMethod()
						if (length(.Object@numberOfItems) == 0)
						{
							if (length(.Object@itemName2ItemIndex) == 0)
							{
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
								}
							}
							
							if (length(.Object@itemName2ItemIndex) != 0)
							{
								.Object@numberOfItems <- max(.Object@itemName2ItemIndex)
							}
						}
						
						.Object
					})

#' Constructor

setGeneric("MgsaMapping", function(sets) standardGeneric("MgsaMapping"))

setMethod(
		"MgsaMapping",
		signature("list"),
		function(sets) new("MgsaMapping",sets=sets)
)


#
#' Returns the indices corresponding to the items
#'
setGeneric("getItemsIndices", function(mapping, items) standardGeneric("getItemsIndices"))

setMethod(
		"getItemsIndices",
		signature( "MgsaMapping","character" ),
		function( mapping, items ) mapping@itemName2ItemIndex[items]
)

#' The basic case, we just add a validity check on top of this
setMethod(
		"getItemsIndices",
		signature( "MgsaMapping", "numeric" ),
		function( mapping, items ) mapping@itemName2ItemIndex[items])
