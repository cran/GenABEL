"add.plot" <- 
function (x, ..., df=1, col=c("lightgreen","lightblue"), sort=TRUE, delta = 1) {
    df <- as.character(df)
#    colmatch <- pmatch(names(match.call()),c("color","colour"))
#    if (any(!is.na(colmatch))) {
#	    carg <- names(match.call())[[(colmatch[!is.na(colmatch)])[1]]]
#	    col <- carg
#    }
    if (!(any(c("1","2","Pc1df","Pc2df")==df))) stop ("df parameter must be 1, 2, \"Pc1df\", or \"Pc2df\"")
    if (class(x) != "scan.gwaa" && class(x) != "scan.gwaa.2D") stop("Plot must be done on an object of class scan.gwaa or scan.gwaa.2D")
    if (length(x$map) != length(x$P1df)) stop("length of map and scan points not equal!")

    Pv <- x$P1df
    if (df=="2") {Pv <- x$P2df;}
    else if (df=="Pc1df") {Pv <- x$Pc1df;}
    else if (df=="Pc2df") {Pv <- x$Pc2df;}
    Pv <- replace(Pv,(Pv<=0),1.e-16)

	if (sort) {
		newmap <- sortmap.internal(x$chromosome,x$map,delta=delta)
		x$map <- newmap$cummap
		x$chromosome <- x$chromosome[newmap$ix]
		Pv <- Pv[newmap$ix]
		newchnum <- newmap$chnum
		rm(newmap)
		gc()
	}

    if (class(x) == "scan.gwaa") {
		if (dim(table(x$chromosome))>1) {
			chind <- chrom.char2num(x$chromosome) %% length(col)
			idxCH <- which(chind==0)
      	points(x$map,-log10(Pv), col=col[length(col)], ...)
			for (colidx in c(1:(length(col)-1))) {
				idxCH <- which(chind==colidx)
				points(x$map[idxCH],-log10(Pv[idxCH]),col=col[colidx],...)
			}
		} else {
      	points(x$map,-log10(Pv), col=col[1], ...)
		}
    } else {
		image(x=x$map,y=x$map,z=t(-log10(Pv)),add=TRUE,...)
    }
}

