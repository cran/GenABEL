"plot.scan.gwaa" <- 
function (x, y, ..., df=1) {
    if (!(any(c(1,2,"all")==df))) stop ("df parameter must be 1, 2, or \"all\"")
    if (class(x) != "scan.gwaa") stop("Plot must be done on an object returned by scan.gwaa()")
    if (length(x$map) != length(x$P1df)) stop("length of map and scan points not equal!")
    if (df==1) plot(x$map,-log10(x$P1df),main=x$formula,xlab="Map",ylab="-log10(P-value)", ...)
    if (df==2) plot(x$map,-log10(x$P2df),main=x$formula,xlab="Map",ylab="-log10(P-value)", ...)
    if (df=="all") {
    	max1y <- max(-log10(x$P1df),na.rm=TRUE)
	max2y <- max(-log10(x$P2df),na.rm=TRUE)
    	maxy <- max(max1y,max2y,na.rm=TRUE)
    	plot(x$map,-log10(x$P1df),main=x$formula,xlab="Map",ylab="-log10(P-value)", ylim=c(0,maxy), ...)
	clr <- "red"
   	if (par()$col == "red") clr <- "green"
	points(x$map,-log10(x$P2df),col=clr)
    }
}

