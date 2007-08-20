"plot.scan.gwaa" <- 
function (x, y, ..., df=1) {
    df <- as.character(df)
    if (!(any(c("1","2","Pc1df","all")==df))) stop ("df parameter must be 1, 2, \"Pc1df\" or \"all\"")
    if (class(x) != "scan.gwaa") stop("Plot must be done on an object of scan.gwaa-class")
    if (length(x$map) != length(x$P1df)) stop("length of map and scan points not equal!")
    if (any(names(x) == "Pgw1df")) if (!is.null(x$Pgw1df)) {x$P1df <- x$Pgw1df; x$P2df <- x$Pgw2df}
    Pv <- x$P1df
    if (df=="2") {Pv <- x$P2df;}
    else if (df=="Pc1df") {Pv <- x$Pc1df;}
    Pv <- replace(Pv,(Pv<=0),1.e-16)
    if (df=="all") {
        x$P2df <- replace(x$P2df,(x$P2df<=0),1.e-16)
        x$Pc1df <- replace(x$Pc1df,(x$Pc1df<=0),1.e-16)
    	max1y <- max(-log10(Pv),na.rm=TRUE)
	max2y <- max(-log10(x$P2df),na.rm=TRUE)
	max3y <- max(-log10(x$Pc1df),na.rm=TRUE)
    	maxy <- max(c(max1y,max2y,max3y),na.rm=TRUE)
    } else {
    	maxy <- max(-log10(Pv),na.rm=TRUE)
    }
	if (dim(table(x$chromosome))>1) {
		if (df=="all") 
			plot(x$map,-log10(Pv),main=(x$formula),xlab="Chromosome",ylab=expression(-log[10](P-value)), axes=FALSE, ylim=c(0,maxy), ...)
		else
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
		if (df=="all") 
			plot(x$map,-log10(Pv),main=(x$formula),xlab="Map position",ylab=expression(-log[10](P-value)), ylim=c(0,maxy), ...)
		else
			plot(x$map,-log10(Pv),main=(x$formula),xlab="Map position",ylab=expression(-log[10](P-value)), ...)
	}
    if (df=="all") {
	points(x$map,-log10(x$P2df),col="blue")
	points(x$map,-log10(x$Pc1df),col="green")
    }
}

