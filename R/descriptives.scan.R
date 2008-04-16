"descriptives.scan" <-
function(data,file,top=10,sortby="P1df",digits=10,sep="\t") {
	if (class(data) != "scan.gwaa") stop("data argument must be of scan.gwaa-class")

	sbargs <- c("no","P1df","P2df","Pc1df")

	if (!(sortby %in% sbargs)) {cat("sortby argument must be one of the following:",sbargs);stop("");}

	data$P1df[is.na(data$P1df)] <- 9.99
	data$P2df[is.na(data$P2df)] <- 9.99
	data$Pc1df[is.na(data$Pc1df)] <- 9.99

	if (sortby=="P1df") ix <- sort.int(data$P1df,index.return=TRUE)$ix[1:top]
	else if (sortby=="P2df") ix <- sort.int(data$P2df,index.return=TRUE)$ix[1:top]
	else if (sortby=="Pc1df") ix <- sort.int(data$Pc1df,index.return=TRUE)$ix[1:top]
	else if (sortby=="no") {top <- length(data$P1df); ix <- c(1:top);}
	else {cat("sortby argument must be one of the following:",sbargs);stop("");}

	out <- matrix(rep(NA,top*7),ncol=7)
	out[c(1:top),1] <- data$N[ix]
	out[c(1:top),2] <- data$effB[ix]
	out[c(1:top),3] <- data$P1df[ix]
	out[c(1:top),4] <- data$Pc1df[ix]
	out[c(1:top),5] <- data$effAB[ix]
	out[c(1:top),6] <- data$effBB[ix]
	out[c(1:top),7] <- data$P2df[ix]
	out <- round(out,digits=digits)
	out <- data.frame(out)
	out <- cbind(data$map[ix],out)
	out <- cbind(data$chromosome[ix],out)
	rownames(out) <- data$snpnames[ix]
	colnames(out) <- c("Chromosome","Position","N","effB","P1df","Pc1df","effAB","effBB","P2df")
	if (!missing(file)) {
		cat(sep,file=file,sep="")
		cat(colnames(out),file=file,sep=sep,append=TRUE)
		cat("\n",file=file,sep="",append=TRUE)
		write.table(out,file=file,sep=sep,append=T,col.names=FALSE)
	}
	out
}
