
#' Class that describes sets and their associations
#' The attributes of this class are completely private.
#' TODO: Add proper accessor functions
setClass(
        "MgsaMapping",
        representation = representation(
							sets = "list",
							item.idx.map = "integer"
			)
)

