"merge.gwaa.data" <-
function(x, y, ... ) {
	if (class(x) != "gwaa.data-class" || class(y) != "gwaa.data-class")
		stop("x and y should have gwaa.data-class")
	snpdata <- try(merge(x@gtdata,y@gtdata, ... ))
	if (class(snpdata) == "try-error")
		stop("error occured in merging gtdata slots of x and y")
	snpdata <- snpdata$data
	phdata <- merge(x@phdata,y@phdata,by="id",all=T)
	if (any(!(snpdata@idnames %in% phdata$id)))
		stop("some ids present in gtdata are missing from phdata")
	phdata <- phdata[snpdata@idnames,]
	out <- new("gwaa.data",phdata=phdata,gtdata=snpdata)
	out
}
