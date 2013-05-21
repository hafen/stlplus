.First.lib <- function(lib, pkg)
{
   library.dynam("stl2", pkg, lib)
   cat("stl2 1.0 loaded\n")
}

