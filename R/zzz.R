#.packageName <- 'GSCA'

.First.lib <- function(lib, pkg){
  #library.dynam("dEuc2perm", pkg, lib)
  library.dynam(pkg, pkg, lib)
}

