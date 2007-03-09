"add.plot" <- 
function (x, ..., df=1) {
    if (class(x) != "scan.gwaa" && class(x) != "scan.gwaa.2D") stop("Plot must be done on an object of class scan.gwaa or scan.gwaa.2D")
    if (class(x) == "scan.gwaa") {
      if (length(x$map) != length(x$P1df)) stop("length of map and scan points not equal!")
      if (df==1)
       points(x$map,-log10(x$P1df), ...)
      else 
       points(x$map,-log10(x$P2df), ...)
    } else {
      if (df==1)
	image(x=x$map,y=x$map,z=t(-log10(x$P1df)),add=TRUE,...)
      else
	image(x=x$map,y=x$map,z=t(-log10(x$P2df)),add=TRUE,...)
    }
}

