# TODO: Add comment
# 
# Author: gagneur
###############################################################################


## S4 implementation
setGeneric(
		"test",
		function( x, y, a=0 ){
			standardGeneric("test")
		}
)

#o, sets, population=NULL, alpha=NA, beta=NA, p=NA, steps=1e6
setMethod(
		"test",
		signature = c(x="numeric", y="numeric"),
		function( x, y, a ) {
			x+y+a
		}
)

