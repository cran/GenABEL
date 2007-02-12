"summary.check.marker" <-
function(object,...) {
	out <- rep(NA,16)
	dim(out) <- c(4,4)
	out[1,2] <- cross(object$nocall,object$nofreq)	
	out[1,3] <- cross(object$nocall,object$nohwe)	
	out[1,4] <- cross(object$nocall,object$redundant)	
	out[1,1] <- length(object$nocall) - out[1,2] - out[1,3] - out[1,4]
	out[2,3] <- cross(object$nofreq,object$nohwe)	
	out[2,4] <- cross(object$nofreq,object$redundant)	
	out[2,2] <- length(object$nofreq) - out[1,2] - out[2,3] - out[2,4]
	out[3,4] <- cross(object$nohwe,object$redundant)	
	out[3,3] <- length(object$nohwe) - out[1,3] - out[2,3] - out[3,4]
	out[4,4] <- length(object$redundant) - out[1,4] - out[2,4] - out[3,4]
	rownames(out) <- c("NoCall","NoMAF","NoHWE","Redundant")
	colnames(out) <- c("NoCall","NoMAF","NoHWE","Redundant")
	print(out)
}

