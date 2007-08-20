"summary.check.marker" <-
function(object,...) {
	out1 <- rep(NA,25)
	dim(out1) <- c(5,5)
	out1[1,2] <- cross(object$nocall,object$nofreq)	
	out1[1,3] <- cross(object$nocall,object$nohwe)	
	out1[1,4] <- cross(object$nocall,object$redundant)	
	out1[1,5] <- cross(object$nocall,object$Xmrkfail)	
	out1[1,1] <- length(object$nocall) - out1[1,2] - out1[1,3] - out1[1,4] - out1[1,5]
	out1[2,3] <- cross(object$nofreq,object$nohwe)	
	out1[2,4] <- cross(object$nofreq,object$redundant)	
	out1[2,5] <- cross(object$nofreq,object$Xmrkfail)	
	out1[2,2] <- length(object$nofreq) - out1[1,2] - out1[2,3] - out1[2,4] -out1[2,5]
	out1[3,4] <- cross(object$nohwe,object$redundant)	
	out1[3,5] <- cross(object$nohwe,object$Xmrkfail)	
	out1[3,3] <- length(object$nohwe) - out1[1,3] - out1[2,3] - out1[3,4] - out1[3,5]
	out1[4,5] <- cross(object$redundant,object$Xmrkfail)	
	out1[4,4] <- length(object$redundant) - out1[1,4] - out1[2,4] - out1[3,4] - out1[4,5]
	out1[5,5] <- length(object$Xmrkfail) - out1[1,5] - out1[2,5] - out1[3,5] - out1[4,5]
	rownames(out1) <- c("NoCall","NoMAF","NoHWE","Redundant","Xsnpfail")
	colnames(out1) <- rownames(out1)
	out2 <- rep(NA,25)
	dim(out2) <- c(5,5)
	out2[1,2] <- cross(object$idnocall,object$hetfail)
	out2[1,3] <- cross(object$idnocall,object$ibsfail)
	out2[1,4] <- cross(object$idnocall,object$isfemale)
	out2[1,5] <- cross(object$idnocall,object$ismale)
	out2[1,1] <- length(object$idnocall) - out2[1,2] - out2[1,3] - out2[1,4] - out2[1,4] - out2[1,5]
	out2[2,3] <- cross(object$hetfail,object$ibsfail)
	out2[2,4] <- cross(object$hetfail,object$isfemale)
	out2[2,5] <- cross(object$hetfail,object$ismale)
	out2[2,2] <- length(object$hetfail) - out2[1,2] - out2[2,3] - out2[2,4] - out2[2,4] - out2[2,5]
	out2[3,4] <- cross(object$ibsfail,object$isfemale)
	out2[3,5] <- cross(object$ibsfail,object$ismale)
	out2[3,3] <- length(object$ibsfail) - out2[1,3] - out2[2,3] - out2[3,4] - out2[3,4] - out2[3,5]
	out2[4,5] <- cross(object$isfemale,object$ismale)
	out2[4,4] <- length(object$isfemale) - out2[1,4] - out2[2,4] - out2[3,4] - out2[4,5]
	out2[5,5] <- length(object$ismale) - out2[1,5] - out2[2,5] - out2[3,5] - out2[4,5]
#	rownames(out2) <- c("IDnoCall","HetFail","IBSFail","Xidfail","isfemale","ismale")
	rownames(out2) <- c("IDnoCall","HetFail","IBSFail","isfemale","ismale")
	colnames(out2) <- rownames(out2)
	out <- list()
	out$"Per-SNP fails statistics" <- out1
	out$"Per-person fails statistics" <- out2
	out
}

