"rntransform" <-
function(formula,data,family=gaussian) {
	var <- ztransform(formula,data,family)
	out <- rank(var) - 0.5
	out[is.na(var)] <- NA
	mP <- .5/max(out,na.rm=T)
	out <- out/(max(out,na.rm=T)+.5)
	out <- qnorm(out)
	out
}


