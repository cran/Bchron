newgraphwindow <- function(width,height) {
# This function aims to create graphics windows on different OS's where possible
graphoutcome <- try(windows(width,height),silent=TRUE)
if(length(graphoutcome) > 0) graphoutcome <- try(quartz(width,height),silent=TRUE)
if(length(graphoutcome) > 0) graphoutcome <- try(x11(width,height),silent=TRUE)
if(length(graphoutcome) > 0) print("Graphics unavailable on this workstation")
}
