
#' Class that describes sets and their associations
#' The attributes of this class are completely private.
#' 
#' TODO: Add proper accessor functions and find out how to hide
#' the raw attributes.
setClass(
        "MgsaSets",
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

setMethod("initialize", "MgsaSets",
					function(.Object, ...) {
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
								}
							}
							
							if (length(.Object@itemName2ItemIndex) != 0)
							{
								.Object@numberOfItems <- max(.Object@itemName2ItemIndex)
							}
						}
						
						.Object
					})

#'
#' Constructor
#' 
setGeneric("MgsaSets", function(sets) standardGeneric("MgsaSets"))

setMethod(
		"MgsaSets",
		signature("list"),
		function(sets) new("MgsaSets",sets=sets)
)


#'
#' Returns the indices corresponding to the items
#'
setGeneric("getItemsIndices", function(mapping, items) standardGeneric("getItemsIndices"))

setMethod(
		"getItemsIndices",
		signature( "MgsaSets","character" ),
		function( mapping, items ) mapping@itemName2ItemIndex[items]
)

setMethod(
		"getItemsIndices",
		signature( "MgsaSets", "numeric" ),
		function( mapping, items ) mapping@itemName2ItemIndex[items])


#'
#' Returns a subset of mapping, i.e., a mapping, in which only the given
#' items are assinged to sets. The number of sets of this mapping may also differ.
#' 
setGeneric("getSubMapping", function(mapping, items) standardGeneric("getSubMapping"))

setMethod(
		"getSubMapping",
		signature( "MgsaSets", "numeric" ),
		function( mapping, items )
		{
			sets<-mapping@sets
			subset.contains<-items

			if (T)
			{
				encode <- function(x){ match( intersect(x, items), items) }
				subsets <- lapply(sets, encode)
				subsets <- subset(subsets,lapply(subsets,length)>0)
			} else
			{
				#
				# We create the subset mapping by the following juggling
				#
				
				# We assume that each set has name
				if (is.null(names(sets))) names(sets)<-1:length(sets)
				
				# First, construct a item->set mapping
				# we assume that each item has at least one set
				set.names<-rep(names(sets),lapply(sets,length))
				set.items<-unlist(sets,use.names=F)
				items<-split(set.names,set.items)
				
				# Take the subset
				items.subset<-items[subset.contains]
				
				# Create a new set->item mapping based on the item subset
				subitem.names<-rep(names(items.subset),lapply(items.subset,length))
				subitem.names.f<-factor(subitem.names)
				levels(subitem.names.f)<-1:length(levels(subitem.names.f))		# relabel the genes
				subitem.sets<-unlist(items.subset,use.names=F)
				subsets<-split(as.vector(subitem.names.f),subitem.sets)
			}

			mapping<-new("MgsaSets",sets=subsets)
			
			return(mapping)
		})
