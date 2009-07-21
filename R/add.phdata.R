"add.phdata" <- 
function(data,phdata) {
	if (!is(data,"gwaa.data")) stop("data argument must be of gwaa.data-class")
	if (!is(phdata,"data.frame")) stop("phdata argument must be of data.frame-class")
	nidcols <- sum(names(phdata) %in% "id")
	if (nidcols==0) stop("can not find \"id\" column in phdata")
	if (nidcols>1) stop("more than one \"id\" column in phdata")
	if (length(unique(phdata$id)) != length(phdata$id)) stop("duplicated id names in phdata")
	if (class(phdata$id) != "character") {
		warning("phdata id variable does not have character class. Converted")
		phdata$id <- as.character(phdata$id)
	}

	oldph <- data@phdata
	if (length(unique(oldph$id)) != length(oldph$id)) stop("duplicated id names in oldph!!!")
	if (class(oldph$id) != "character") {
		stop("oldph id variable does not have character class!!!")
	}
	nadded <- sum((oldph$id %in% phdata$id),na.rm=T)
	if (nadded < length(oldph$id)) 
		warning(paste("Of",length(oldph$id),"IDs present in data, only",nadded,"found in phdata"))
	newph <- merge(oldph,phdata,by="id",all.x=T,all.y=F)
	if (dim(newph)[1] != dim(oldph)[1]) stop("something terribly wrong in merge...")
	newph <- newph[match(oldph$id,newph$id),]
	rownames(newph) <- newph$id
	
	if (is.na(match("sex",names(newph))))
	if (!is.na(match("sex.x",names(newph)))) {
		nvar <- match("sex.x",names(newph))
		names(newph)[nvar] <- "sex"
	}

	data@phdata <- newph
	data
}
