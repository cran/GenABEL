"descriptives.scan" <-
function(data,file,top=10,sortby="P1df",digits=6) {
	if (class(data) != "scan.gwaa") stop("data argument must be of scan.gwaa-class")

	sbargs <- c("P1df","P2df","Pc1df")
	if (!(sortby %in% sbargs)) {cat("sortby argument must be one of the following:",sbargs);stop("");}

	if (sortby=="P1df") ix <- sort.int(data$P1df,index.return=TRUE)$ix[1:top]
	else if (sortby=="P2df") ix <- sort.int(data$P2df,index.return=TRUE)$ix[1:top]
	else if (sortby=="Pc1df") ix <- sort.int(data$Pc1df,index.return=TRUE)$ix[1:top]
	else  {cat("sortby argument must be one of the following:",sbargs);stop("");}

	out <- matrix(rep(NA,top*6),ncol=6)
	j <- 1
	for (i in ix) {
		out[j,1] <- data$effB[i]
		out[j,2] <- data$P1df[i]
		out[j,3] <- data$Pc1df[i]
		out[j,4] <- data$effAB[i]
		out[j,5] <- data$effBB[i]
		out[j,6] <- data$P2df[i]
		j <- j+1
	}
	out <- round(out,digits=digits)
	out <- data.frame(out)
	out <- cbind(data$map[ix],out)
	out <- cbind(data$chromosome[ix],out)
	rownames(out) <- data$snpnames[ix]
	colnames(out) <- c("Chromosome","Position","effB","P1df","Pc1df","effAB","effBB","P2df")
	if (!missing(file)) {
		cat("\t",file=file,sep="")
		cat(colnames(out),file=file,sep="\t",append=TRUE)
		cat("\n",file=file,sep="",append=TRUE)
		write.table(out,file=file,sep="\t",append=T,col.names=FALSE)
	}
	out
}
