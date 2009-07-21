"descriptives.trait" <-
function (data,subset,file,by.var=NULL,digits=3) {
	if (missing(data)) stop("data argument must be provided")
	if (class(data) == "gwaa.data") data <- data@phdata
	if (class(data) != "data.frame") stop("data argument must be of gwaa.data or data.frame-class")
	len <- dim(data)[1]
	ntra <- dim(data)[2]
	if (!missing(subset)) data <- data[subset,]
	if (!is.null(by.var)) {
		svar <- by.var
		if (length(levels(factor(svar)))!=2) stop("The by.var argument should contain a binary variable")
		out <- matrix(data=NA,nrow=ntra,ncol=9)
		lvls <- levels(factor(svar))
		for (i in (1:ntra)) {
			ctrao <- data[,i]
			if (!is.numeric(ctrao) | all(ctrao==svar)) {
				ctra <- ctrao[svar == lvls[1]]
				out[i,1] <- length(ctra) - sum(is.na(ctra))
				ctra <- ctrao[svar == lvls[2]]
				out[i,4] <- length(ctra) - sum(is.na(ctra))
			} else {
				ctra <- ctrao[svar == lvls[1]]
				out[i,1] <- length(ctra) - sum(is.na(ctra))
				out[i,2] <- mean(ctra,na.rm=TRUE)
				out[i,3] <- sd(ctra,na.rm=TRUE)
				ctra <- ctrao[svar == lvls[2]]
				out[i,4] <- length(ctra) - sum(is.na(ctra))
				out[i,5] <- mean(ctra,na.rm=TRUE)
				out[i,6] <- sd(ctra,na.rm=TRUE)
				out[i,7] <- t.test(ctrao ~ svar)$p.value
				out[i,8] <- kruskal.test(ctrao ~ svar)$p.value
				clv <- length(unique(ctrao))
				if (clv>1 & clv<5) out[i,9] <- fisher.test(ctrao,svar)$p.value
			}
		}
	} else {
		out <- matrix(data=NA,nrow=ntra,ncol=3)
		for (i in (1:ntra)) {
			ctra <- data[,i]
			out[i,1] <- length(ctra) - sum(is.na(ctra))
			if (!is.numeric(ctra)) next;
			out[i,2] <- mean(ctra,na.rm=TRUE)
			out[i,3] <- sd(ctra,na.rm=TRUE)
		}
	}
	out <- round(out,digits=digits)
	out <- data.frame(out)
	rownames(out) <- colnames(data)
	if (is.null(by.var)) 
		colnames(out) <- c("No","Mean","SD")
	else 
		colnames(out) <- c(paste("No(by.var=",lvls[1],")",sep=""),"Mean","SD",paste("No(by.var=",lvls[2],")",sep=""),"Mean","SD","Ptt","Pkw","Pexact")
	if (!missing(file)) {
		cat("\t",file=file,sep="")
		cat(colnames(out),file=file,sep="\t",append=TRUE)
		cat("\n",file=file,sep="",append=TRUE)
		write.table(out,file=file,sep="\t",append=T,col.names=FALSE)
	}
	out
}
