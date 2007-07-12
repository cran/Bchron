`hdr2.den` <-
function (x = NULL, prob = c(50, 95, 99), den = NULL, h = NULL, ylab = "", ...) 
{
    library(hdrcde)

    if (missing(den)) 
        den <- density(x, bw=h)
    else if (missing(x)) 
        x <- sample(den$x, 500, replace = TRUE, prob = den$y)
    hdr <- hdr(x, prob, den, h)
    maxden <- max(den$y)
    plot(den, type = "l", ...)
    cols <- rep(c(0, 4, 2, 3), 3)
    nregions <- nrow(hdr$hdr)
    for (i in 1:nregions) 
    {
        #lines(range(den$x), rep(hdr$falpha[i], 2), col = 5)
        for (j in 1:length(hdr$hdr[i, ])) 
        {
            lines(rep(hdr$hdr[i,j], 2), c((0.01 + (i - 1) * 0.02) * maxden, hdr$falpha[i]),col = 6)
        }
    }
    for (i in 1:nrow(hdr$hdr)) 
    {
        for(j in 1:(length(hdr$hdr[i,])/2))
        {
            lines(c(hdr$hdr[i,2*j-1],hdr$hdr[i,2*j]),rep((0.01 + (i - 1) * 0.02) * maxden,2),lwd=7,lend=1,col=cols[i+1])
        }
        
        #lines(hdr$hdr[i, ], rep((0.01 + (i - 1) * 0.02) * maxden,2), 0.02 * maxden, col = cols[i + 1], horiz = TRUE)
        #lines(hdr$hdr[i, ], rep((0.01 + (i - 1) * 0.02) * maxden,length(hdr$hdr[i,])), lwd=7,lend=1,  col = cols[i + 1])
    }
    return(hdr)
}

