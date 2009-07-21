"export.merlin" <- function(data,pedfile="merlin.ped",datafile="merlin.dat",mapfile="merlin.map",format="merlin",fixstrand="no",extendedmap=TRUE,traits=1) {
	if (class(data) != "gwaa.data") stop("Data argumet should be of gwaa.data-class")
	formats <- c("merlin")
	if (!(match(format,formats,nomatch=0)>0)) {
		out <- paste("fromat argument should be one of",formats,"\n")
		stop(out)
	}
	fixes <- c("no","+","-")
	if (!(match(fixstrand,fixes,nomatch=0)>0)) {
		out <- paste("fixstrand argument should be one of",fixes,"\n")
		stop(out)
	}
	if (fixstrand != "no") {
###### as.raw(1) == "+"
###### as.raw(2) == "-"
		if (fixstrand == "-") {tos <- as.raw(1);fos<-as.raw(2);} 
			else if (fixstrand == "+") {tos <- as.raw(2);fos <- as.raw(1)}
			else {stop("strange strand...")}
		stf <- which(as.raw(data@gtdata@strand) == tos) 
		if (length(stf)>0) {
		  revs <- alleleID.revstrand()
  		  savnam <- names(data@gtdata@strand)
#
#		data@gtdata@coding[stf] <- new("snp.coding",as.raw(revs[as.character(data@gtdata@coding[stf])]))
#
#		ugly fix -- 
#
		  data@gtdata@coding[stf] <- new("snp.coding",as.raw(revs[as.character(as.raw(data@gtdata@coding[stf]))]))
		  names(data@gtdata@coding) <- savnam
		  data@gtdata@strand[stf] <- new("snp.strand",rep(fos,length(stf)))
		  names(data@gtdata@strand) <- savnam
		}
	}
	bstp <- 100
	if (data@gtdata@nids>(bstp*1.5)) {
		steps <- seq(from=0,to=data@gtdata@nids,by=bstp)
		if (data@gtdata@nids != steps[length(steps)]) steps[length(steps)+1] <- data@gtdata@nids
		dump.piece(data=data,from=1,to=steps[2],traits=traits,pedfile=pedfile,append=F)
		for (jjj in c(2:(length(steps)-1))) {
			dump.piece(data=data,from=(steps[jjj]+1),to=(steps[jjj+1]),traits=traits,pedfile=pedfile,append=T)
		}
	} else {
		dump.piece(data=data,from=1,to=data@gtdata@nids,traits=traits,pedfile=pedfile,append=F)
	}
	snps <- data@gtdata@snpnames
	inf <- data.frame(a=rep("M",length(snps)),b=snps)
	if (traits>0) {
		tran <- paste("faket",seq(1:traits),sep="")
		inf0 <- data.frame(a=rep("T",traits),b=tran)
		inf <- rbind(inf0,inf)
	}
	write.table(inf,file=datafile,col.n=FALSE,row.n=FALSE,quote=FALSE)
	map <- data.frame(chromosome=as.character(data@gtdata@chromosome),markername=data@gtdata@snpnames,position=data@gtdata@map)
	write.table(map,file=mapfile,col.n=TRUE,row.n=FALSE,quote=FALSE)
	if (extendedmap) {
		map <- data.frame(chromosome=as.character(data@gtdata@chromosome),markername=data@gtdata@snpnames,position=data@gtdata@map,strand=as.character(data@gtdata@strand),coding=as.character(data@gtdata@coding))
		write.table(map,file=paste(mapfile,".ext",sep=""),col.n=TRUE,row.n=FALSE,quote=FALSE)
	}
}

dump.piece <- function(data,fromid,toid,traits,pedfile,append) {
	if (toid < fromid) stop("toid<fromid")
	x <- as.character(data@gtdata[c(fromid:toid),])
	x <- replace(x,(x==""),"0/0")
	x <- replace(x,is.na(x),"0/0")
	ids <- rownames(x)
	nids <- length(ids)
	sx <- data@phdata$sex[c(fromid:toid)]
	sx <- replace(sx,(sx==0),2)
	if (traits > 0) {
		tr <- matrix(rep(0,nids*traits),ncol=traits)
		x <- data.frame(seq(fromid,toid),ids,rep(0,nids),rep(0,nids),sx,tr,x)
	} else {
		x <- data.frame(seq(fromid,toid),ids,rep(0,nids),rep(0,nids),sx,x)
	}
	write.table(x,file=pedfile,col.n=FALSE,row.n=FALSE,quote=FALSE,append=append)
}
