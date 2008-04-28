Bchronplotevent <- function(Bchrondata) {

cat(paste("Now plotting date estimates for each of the depths in ",
    Bchrondata$name, "eventages.txt \n", sep = ""))

if(length(Bchrondata$eventnames)>1) {
  choices <- c(Bchrondata$eventnames,"All")
  title <- "Which event would you like to plot?"
  choose <- menu(choices, title = title)
} else {
  choose <- 1
}

if(choose == length(Bchrondata$eventnames)+1) {
  Bchrondata$hdreventoutput <- rep("",length(Bchrondata$eventnames))
  for(i in 1:length(Bchrondata$eventnames)) {
    Bchrondata$hdreventoutput[i] <- paste(Bchrondata$path,"/Output/",Bchrondata$name,
        "EventAges",Bchrondata$eventnames[i],"HDRs.txt", sep = "")
      Bchronplotdens(Bchrondata$eventagefile[i],Bchrondata$fullname,Bchrondata$hdreventoutput[i],Bchrondata$eventfullnames[i],Bchrondata$version)
  }
} else {
  Bchrondata$hdreventoutput[choose] <- paste(Bchrondata$path,"/Output/",Bchrondata$name,
      "EventAges",Bchrondata$eventnames[choose],"HDRs.txt", sep = "")
  Bchronplotdens(Bchrondata$eventagefile[choose],Bchrondata$fullname,Bchrondata$hdreventoutput[choose],Bchrondata$eventfullnames[choose],Bchrondata$version)
}

cat("\n")
cat("Done! \n")
cat("Press <Enter> to continue...")
readline()
invisible()

return(Bchrondata)

}