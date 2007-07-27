"plot.scan.gwaa" <- 
function (x, y, ..., df=1) {
    if (!(any(c(1,2,"all")==df))) stop ("df parameter must be 1, 2, or \"all\"")
    if (class(x) != "scan.gwaa") stop("Plot must be done on an object returned by scan.gwaa()")
    if (length(x$map) != length(x$P1df)) stop("length of map and scan points not equal!")
    if (any(names(x) == "Pgw1df")) if (!is.null(x$Pgw1df)) {x$P1df <- x$Pgw1df; x$P2df <- x$Pgw2df}
    if (df==1) Pv <- x$P1df else if (df==2) Pv <- x$P2df
    if (df==1 || df==2) {
	if (dim(table(x$chromosome))>1) {
		plot(x$map,-log10(Pv),main=(x$formula),xlab="Chromosome",ylab=expression(-log[10](P-value)), axes=FALSE, ...)
		mxlog <- floor(max(-log10(Pv)))
		if (mxlog==0) mxlog <- max(-log10(Pv))
		if (mxlog<1) 
			axis(2,at=c(0,mxlog))
		else
			axis(2,at=c(0:mxlog))
		chrs <- levels(x$chromosome)
		chpos <- rep(NA,length(chrs))
		for (i in 1:length(chrs)) chpos[i] <- mean(x$map[x$chromosome==chrs[i]])
		axis(1,at=chpos,labels=chrs)
	} else {
		plot(x$map,-log10(Pv),main=(x$formula),xlab="Map position",ylab=expression(-log[10](P-value)), ...)
	}
    }
    if (df=="all") {
    	max1y <- max(-log10(x$P1df),na.rm=TRUE)
	max2y <- max(-log10(x$P2df),na.rm=TRUE)
    	maxy <- max(max1y,max2y,na.rm=TRUE)
	if (dim(table(x$chromosome))>1) {
		plot(x$map,-log10(x$P1df),main=(x$formula),xlab="Chromosome",ylab=expression(-log[10](P-value)), axes=F, ylim=c(0,maxy), ...)
		mxlog <- floor(maxy)
		if (mxlog==0) mxlog <- maxy
		if (mxlog<1) 
			axis(2,at=c(0,mxlog))
		else
			axis(2,at=c(0:mxlog))
		chrs <- levels(x$chromosome)
		chpos <- rep(NA,length(chrs))
		for (i in 1:length(chrs)) chpos[i] <- mean(x$map[x$chromosome==chrs[i]])
		axis(1,at=chpos,labels=chrs)
		clr <- "red"
	   	if (par()$col == "red") clr <- "green"
		points(x$map,-log10(x$P2df),col=clr,cex=0.5)
	} else {
		plot(x$map,-log10(x$P1df),main=(x$formula),xlab="Map position",ylab=expression(-log[10](P-value)), ylim=c(0,maxy), ...)
		clr <- "red"
	   	if (par()$col == "red") clr <- "green"
		points(x$map,-log10(x$P2df),col=clr)
	}
    }
}

