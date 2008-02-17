"plot.scan.gwaa" <- 
function (x, y, ..., df=1, ylim) {
    df <- as.character(df)
    if (!(any(c("1","2","Pc1df","Pc2df","all")==df))) stop ("df parameter must be 1, 2, \"Pc1df\", \"Pc2df\" or \"all\"")
    if (class(x) != "scan.gwaa") stop("Plot must be done on an object of scan.gwaa-class")
    if (length(x$map) != length(x$P1df)) stop("length of map and scan points not equal!")
    if (any(names(x) == "Pgw1df")) if (!is.null(x$Pgw1df)) {x$P1df <- x$Pgw1df; x$P2df <- x$Pgw2df}
    Pv <- x$P1df
    if (df=="2") {Pv <- x$P2df;}
    else if (df=="Pc1df") {Pv <- x$Pc1df;}
    else if (df=="Pc2df") {Pv <- x$Pc2df;}
    Pv <- replace(Pv,(Pv<=0),1.e-16)
    if (df=="all") {
        x$P2df <- replace(x$P2df,(x$P2df<=0),1.e-16)
        x$Pc1df <- replace(x$Pc1df,(x$Pc1df<=0),1.e-16)
        x$Pc2df <- replace(x$Pc2df,(x$Pc2df<=0),1.e-16)
    	max1y <- max(-log10(Pv),na.rm=TRUE)
	max2y <- max(-log10(x$P2df),na.rm=TRUE)
	max3y <- max(-log10(x$Pc1df),na.rm=TRUE)
	max4y <- max(-log10(x$Pc2df),na.rm=TRUE)
    	maxy <- max(c(max1y,max2y,max3y,max4y),na.rm=TRUE)
    } else {
    	maxy <- max(-log10(Pv),na.rm=TRUE)
    }
	cargnams <- names(match.call())
	
	if (!missing(ylim)) {
		ylim <- ylim
	} else {
		ylim <- c(0,maxy)
	}
	if (dim(table(x$chromosome))>1) {
		if (df=="all") 
			plot(x$map,-log10(Pv),main=(x$formula),xlab="Chromosome",ylab=expression(-log[10](P-value)), axes=FALSE, ylim=ylim, ...)
		else
			plot(x$map,-log10(Pv),main=(x$formula),xlab="Chromosome",ylab=expression(-log[10](P-value)), axes=FALSE, ylim=ylim, ...)
		mxlog <- floor(max(-log10(Pv),na.rm=T))
		if (mxlog==0) mxlog <- max(-log10(Pv),na.rm=T)
		if (mxlog<1) 
			axis(2,at=ylim)
		else
			axis(2,at=c(ylim[1]:ylim[2]))
		chrs <- levels(x$chromosome)
		chpos <- rep(NA,length(chrs))
		for (i in 1:length(chrs)) chpos[i] <- mean(x$map[x$chromosome==chrs[i]])
		axis(1,at=chpos,labels=chrs)
	} else {
		if (df=="all") 
			plot(x$map,-log10(Pv),main=(x$formula),xlab="Map position",ylab=expression(-log[10](P-value)), ylim=ylim, ...)
		else
			plot(x$map,-log10(Pv),main=(x$formula),xlab="Map position",ylab=expression(-log[10](P-value)), ylim=ylim, ...)
	}
    if (df=="all") {
	points(x$map,-log10(x$P2df),col="blue")
	points(x$map,-log10(x$Pc1df),col="green")
	points(x$map,-log10(x$Pc2df),col="red")
    }
}

