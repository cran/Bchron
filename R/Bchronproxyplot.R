Bchronproxyplot <- function (Bchrondata, proxy, title = NULL, xlabel = "Age (k cal yrs BP)",  ylabel = "Proxy", num = 1000, col = "blue",meancol="red", transp=0.1,...) {

# Read in all chronologies  
Bchrondata$chrons <- as.matrix(read.table(Bchrondata$chronsfile))
Bchrondata$ranges <- read.table(Bchrondata$rangesfile,header=TRUE)

# Now produce plot of lines
xlimits <- c(min(Bchrondata$ranges[, 2]), max(Bchrondata$ranges[,4]))
ylimits <- c(min(proxy), max(proxy))
dev.new(...)
plot(1, 1, type = "n", xlim = xlimits, ylim = ylimits, main = title, xlab = xlabel, ylab = ylabel, las = 1)
grid()

# Draw lines for all values in suitable see-through colours
tmp <- col2rgb(col)
mycol <- rgb(tmp[1, 1]/255, tmp[2, 1]/255, tmp[3, 1]/255)
mycol2 <- paste(mycol, as.character(as.hexmode(round(transp*255, 0))), sep = "")
myrows <- sample(seq(1, nrow(Bchrondata$chrons)), num, replace = FALSE)
for(i in 1:num) lines(Bchrondata$chrons[myrows[i], ],proxy, col = mycol2)

lines(colMeans(Bchrondata$chrons),proxy,col=meancol)


}
