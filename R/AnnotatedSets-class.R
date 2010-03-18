## the class of AnnotatedSets 
## mgsa() works on sets. Annoatted Sets allows to have annotaion on sets, such as full name, decription, etc.
## 

######## Class definitions
#### AnnotatedSets


setClass(
		"AnnotatedSets",
		representation = representation(
				sets = "list",
				setAnnotation = "data.frame",
				elementAnnotation = "data.frame"
		)
)

######## show
setMethod(
		"show",
		signature( "AnnotatedSets" ),
		function( object ){
			cat(
					"Object of class ",
					class( object ),
					"\n",
					nrow(object@setAnnotation),
					" annotated sets.\n",
					nrow(object@elementAnnotation),
					" annotated elements.\n",
					sep = ""
			)
			
			nsetShow <- min (5 , length(object@sets) )
			print(str(object@sets[1:nsetShow]))
			if(nsetShow < length(object@sets) ){
				cat("... and ", length(object@sets) - nsetShow, " other sets.\n" )
			}
			
			cat("\nSet Annotation:\n")
			nrowShow <- min (5 , nrow(object@setAnnotation) )
			
			print( object@setAnnotation[1:nrowShow, ] )
			if(nrowShow < nrow(object@setAnnotation) ){
				cat("... and ", nrow(object@setAnnotation) - nrowShow, " other sets.\n" )
			}
			
			nrowShow <- min (5 , nrow(object@elementAnnotation) )
			
			cat("\nElement Annotation:\n")
			print( object@elementAnnotation[1:nrowShow, ] )
			if(nrowShow < nrow(object@elementAnnotation) ){
				cat("... and ", nrow(object@elementAnnotation) - nrowShow, " other elements.\n" )
			}
			
		}
)

