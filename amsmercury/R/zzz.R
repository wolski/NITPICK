# .packageName <- "Mercury";
# 
# .First.lib <- function(libname, pkgname) {
# 	#cat("ams.mercury - Copyright (c) 2006 Marc Kirchner\nBased on 'emass' (Haimi P and Rockwood AL, 2006)\n");
# 	library.dynam("Mercury", pkgname, libname);
# }
# 
# .Last.lib <- function(libpath) {
# 	library.dynam.unload("Mercury");
# }

.onAttach <- function(lib, pkg){
  if(interactive()){
    version <- packageVersion('Mercury')
    packageStartupMessage("Package 'Mercury' version ", version)
    invisible()
  }
}