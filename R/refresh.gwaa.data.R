"refresh.gwaa.data" <- 
function(data) {
	if (class(data) != "gwaa.data") stop("data argument should be of gwaa.data-class")
	err <- try(length(data@gtdata@coding),silent=T)
	if (class(err) != "try-error") stop("data already appear to be of new format")
	data@gtdata@coding <- new("snp.coding",as.raw(rep(1,length(pos))))
	data@gtdata@strand <- new("snp.strand",as.raw(rep(0,length(pos))))
	rownames(data@phdata) <- data@phdata$id
	data
}
