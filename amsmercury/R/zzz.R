.packageName <- "amsmercury";

.First.lib <- function(libname, pkgname) {
	#cat("ams.mercury - Copyright (c) 2006 Marc Kirchner\nBased on 'emass' (Haimi P and Rockwood AL, 2006)\n");
	library.dynam("amsmercury", pkgname, libname);
}

.Last.lib <- function(libpath) {
	library.dynam.unload("amsmercury");
}
