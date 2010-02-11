
#' Class that describes sets and their associations
#' TODO: Add proper accessor functions
setClass(
        "MgsaMapping",
        representation = representation(
							sets = "list",
							item.idx.map = "integer"
			)
)
