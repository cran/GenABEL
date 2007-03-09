"perid.summary" <- 
function (data,snpsubset,idsubset) {
	if (class(data)=="gwaa.data") data <- data@gtdata
	if (class(data)!="snp.data") stop("The data argument must be of snp.data-class or gwaa.data-class")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	out <- .C("hom",as.raw(data@gtps),as.integer(data@nids),as.integer(data@nsnps),as.integer(0),sout = double(2*data@nids), PACKAGE="GenABEL")$sout
	out <- matrix(out,ncol=2)
	nmeas <- out[,1]
	hom <- out[,2]
	hom[nmeas>0] <- hom[nmeas>0]/nmeas[nmeas>0]
	hom[nmeas<=0] <- NA
	out <- data.frame(NoMeasured=nmeas,CallPP=nmeas/data@nsnps,Het=(1.-hom))
	rownames(out)=data@idnames
	out
}
