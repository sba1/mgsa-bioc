.First.lib <- function(lib, pkg) {
  print(paste(lib,pkg))

  library.dynam("mgsa", pkg, lib)
#  vers <- paste(sessionInfo()$otherPkg$mgsa$Version,".",sep="")

  print(paste("Package mgsa initialized"))
}

.Last.lib <- function() {
  library.dynam.unload("mgsa", paste(.libPaths(),"/mgsa",sep=""))
}
