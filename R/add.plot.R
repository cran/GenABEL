"add.plot" <- 
function (x, ..., df=1) {
    df <- as.character(df)
    if (!(any(c("1","2","Pc1df","Pc2df")==df))) stop ("df parameter must be 1, 2, \"Pc1df\", or \"Pc2df\"")
    if (class(x) != "scan.gwaa" && class(x) != "scan.gwaa.2D") stop("Plot must be done on an object of class scan.gwaa or scan.gwaa.2D")
    if (length(x$map) != length(x$P1df)) stop("length of map and scan points not equal!")
    Pv <- x$P1df
    if (df=="2") {Pv <- x$P2df;}
    else if (df=="Pc1df") {Pv <- x$Pc1df;}
    else if (df=="Pc2df") {Pv <- x$Pc2df;}
    Pv <- replace(Pv,(Pv<=0),1.e-16)
    if (class(x) == "scan.gwaa") {
       points(x$map,-log10(Pv), ...)
    } else {
	image(x=x$map,y=x$map,z=t(-log10(Pv)),add=TRUE,...)
    }
}

