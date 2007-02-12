"add.plot.scan.gwaa" <- 
function (x, ..., df=1) {
    if (class(x) != "scan.gwaa") stop("Plot must be done on an object returned by scan.gwaa()")
    if (length(x$map) != length(x$P1df)) stop("length of map and scan points not equal!")
    if (df==1)
     points(x$map,-log10(x$P1df), ...)
    else 
     points(x$map,-log10(x$P2df), ...)
}

