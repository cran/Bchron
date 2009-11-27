.onAttach <- function(libname, pkgname) {
    Bchronver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                      fields="Version")
    message(paste(pkgname, Bchronver))
    message("Type help(Bchron) for installation instructions or Bchronmenu() to get started.")
}
