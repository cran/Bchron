.onAttach <- function(libname, pkgname) {
    Bchronver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                      fields="Version")
    packageStartupMessage(paste(pkgname, Bchronver))
    packageStartupMessage("Type help(Bchron) for installation instructions\n")
    packageStartupMessage("See http://mathsci.ucd.ie/~parnell_a/ for further details and updates.")
}