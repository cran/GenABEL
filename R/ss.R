#
#defining snp.mx -- internal class
#
#setClass("genotype",contains="data.frame",package="GenABEL")
#setClass("hsgeno",contains="data.frame",package="GenABEL")

setClass("snps.cell",contains="raw",package="GenABEL")
setClass("snp.mx",contains="matrix",package="GenABEL")
setMethod("[","snp.mx",
	function (x, i, j , drop = FALSE) {
		if (missing(j)) j=c(1:ncol(x));
		if (missing(drop)) drop = FALSE;
		if (is.logical(j)) j=which(j)
		k <- c(1:nrow(x))
		x <- x@.Data[k, j , drop = FALSE]
		if (missing(i)) i=c(1:(nrow(x)*4));
		if (is.logical(i)) i=which(i)
		x <- sset(as.raw(x),length(j),nrow(x)*4,i)
		dim(x) <- c(ceiling(length(i)/4),length(j))
		if (is.matrix(x)) 
			x <- new("snp.mx",x) 
		else 
			x <- new("snps.cell",x)
		x
	})

setMethod("show","snp.mx",
	function(object) {
		nr <- dim(object)[1]
		for (i in (1:nr)) {
			cat(as.raw(object[i,]),"\n")
		}
	})


setClass("snp.coding",contains="raw",package="GenABEL")
setMethod("[",signature(x="snp.coding",i="ANY",j="missing",drop="missing"),
	function (x,i,j,...,drop) {
		nams <- names(x)
		x <- as.raw(x)
		names(x) <- nams
		out <- new("snp.coding",x[i])
		names(out) <- nams[i]
		out
	})
setMethod("show","snp.coding",
	function (object) {
		snam <- names(object)
#		x <- alleleID.raw2char()
#		out <- x[as.character(object)]
		out <- as.raw(object)
		names(out) <- snam
		print(out)
#		cat(out,"\n");
	})


setClass("snp.strand",contains="raw",package="GenABEL")
setMethod("[",signature(x="snp.strand",i="ANY",j="missing",drop="missing"),
	function (x,i,j,...,drop) {
		nams <- names(x)
		x <- as.raw(x)
		names(x) <- nams
		out <- new("snp.strand",x[i])
		names(out) <- nams[i]
		out
	})
setMethod("show","snp.strand",
	function (object) {
		snam <- names(object)
#		tmpo <- c("u","+","-")[as.numeric(object)+1]
		tmpo <- as.raw(object)
		names(tmpo) <- snam
		print(tmpo)
#		cat(tmpo,"\n");
	})


setClass("snp.data",representation(nbytes="numeric",nids="numeric",nsnps="numeric",
		idnames="character",snpnames="character",chromosome="factor",map="numeric",
		coding="snp.coding",strand="snp.strand",
		male="numeric",gtps="snp.mx"),package="GenABEL")
snp.data <- function (nids,rawdata,
		idnames=as.character(c(1:nids)),
		snpnames=as.character(c(1:(length(rawdata)/ceiling(nids/4)))),
		chromosome=as.factor(rep(1,(length(rawdata)/ceiling(nids/4)))),
		map=as.double(seq(1,(length(rawdata)/ceiling(nids/4)))),
		coding=as.raw(rep(1,length(rawdata)/ceiling(nids/4))),
		strand=as.raw(rep(0,length(rawdata)/ceiling(nids/4))),
		male=rep(0,nids)
		) {
	nbytes <- ceiling(nids/4)
	nsnps <- length(rawdata)/nbytes
	if (nsnps != round(nsnps)) stop("wrong number of ids")
	names(map) <- snpnames
	names(chromosome) <- snpnames
	names(coding) <- snpnames
	names(strand) <- snpnames
	names(male) <- idnames
	t <- as.raw(rawdata)
	dim(t) <- c(nbytes,nsnps)
	t <- new("snp.mx",t)
	a <- new("snp.data",nbytes=nbytes,nids=nids,nsnps=nsnps,gtps=t,
			idnames=idnames,snpnames=snpnames,
			chromosome=chromosome,map=map,
			coding=coding,
			strand=strand,
			male=male)
	a
}

setMethod("[","snp.data",
	function(x, i, j, drop) {
		if (missing(j)) j=c(1:x@nsnps);
		if (missing(i)) i=c(1:x@nids);
		if (missing(drop)) drop = FALSE;
		if (is.logical(i)) i=which(i)
		if (is.logical(j)) j=which(j)
		if (is.character(i)) {
			tmp <- i
			i=match(i,x@idnames)
			if (any(is.na(i))) stop(paste("following IDs were not found:",tmp[which(is.na(i))],"\n"))
		}
		if (is.character(j)) {
			tmp <- j
			j=match(j,x@snpnames)
			if (any(is.na(j))) stop(paste("following SNPs were not found:",tmp[which(is.na(j))],"\n"))
		}
		if (length(i) > x@nids  || max(i) > x@nids ) stop("i out of range")
		if (length(j) > x@nsnps || max(j) > x@nsnps) stop("j out of range")
		if (length(i) <= 0) stop("i out of range (== 0)")
		if (length(j) <= 0) stop("j out of range (== 0)")
		a <- new("snp.data")
		a@nids <- length(i)
		a@nsnps <- length(j)
		a@nbytes <- ceiling(a@nids/4.)
		a@snpnames <- x@snpnames[j]
		a@map <- x@map[j]
		a@chromosome <- as.factor(as.character(x@chromosome[j]))
		names(a@chromosome) <- a@snpnames
		a@coding <- new("snp.coding",x@coding[j])
		names(a@coding) <- a@snpnames
		a@strand <- new("snp.strand",x@strand[j])
		names(a@strand) <- a@snpnames
		a@idnames <- x@idnames[i]
		a@male <- x@male[i]
		names(a@male) <- a@idnames
		a@gtps <- x@gtps[i,j]
		a
	})
#setMethod("summary","snp.data",
#	function(object,...) {
#)

setMethod("show","snp.data",
	function(object) {
		cat("@nids =",object@nids,"\n")
		cat("@nsnps =",object@nsnps,"\n")
		cat("@nbytes =",object@nbytes,"\n")
		cat("@idnames =",object@idnames,"\n")
		cat("@snpnames =",object@snpnames,"\n")
		cat("@chromosome =",object@chromosome,"\n")
		cat("@coding = ",as.raw(object@coding),"\n")
		cat("@strand = ",as.raw(object@strand),"\n")
		cat("@map =",object@map,"\n")
		cat("@male =",object@male,"\n")
		cat("@gtps = \n");
		show(object@gtps)
	})

#setAs("snp.data","double",function(from) {to <- as.double.snp.data(from);to})
#
#setAs("snp.data","numeric",function(from) {to <- as.double.snp.data(from);to})
#
#setAs("snp.data","numeric",
#	function(from) {
##		tnsnps <- from@nsnps
#		to <- .C("get_snps_many",as.raw(from@gtps), as.integer(tnids), as.integer(tnsnps), idata = integer(tnids*tnsnps), PACKAGE="GenABEL")$idata
#		to <- replace(to,(to==0),NA)
#		to <- to - 1
#		dim(to) <- c(tnids,tnsnps)
#		colnames(to) <- from@snpnames
#		rownames(to) <- from@idnames
#		to
#	})
#setAs("snp.data","character",
#	function(from) {
#		a <- as(from,"numeric")
#		dm <- dim(a)
#		to <- ifelse(is.na(a),"",c("A/A","A/B","B/B")[a+1])
#		dim(to) <- dm
#		to
#	})
#setAs("snp.data","raw",
#	function(from) {
#		to <- as.raw(from@gtps)
##	})
#setAs("snp.data","genotype",
#	function(from) {
#	for (i in 2:from@nsnps) {
#		gdta <- cbind(gdta,genotype(as(from[,i],"character")))
#	}
#	colnames(gdta) <- from@snpnames
##	class(gdta) <- "genotype"
#	gdta
#	})
#setAs("snp.data","hsgeno",
#	function(from) {
#	g1 <- as(from[,1],"numeric")
#	a1 <- rep(NA,length(g1))
#	a2 <- rep(NA,length(g1))
#	a1 <- replace(a1,(g1==0 | g1==1),1)
#	a1 <- replace(a1,(g1==2),2)
#	a2 <- replace(a2,(g1==1 | g1==2),2)
#	gdta <- data.frame(a1,a2)
#	for (i in 2:from@nsnps) {
#		g1 <- as(from[,i],"numeric")
#		a1 <- rep(NA,length(g1))
#		a2 <- rep(NA,length(g1))
#		a1 <- replace(a1,(g1==0 | g1==1),1)
#		a1 <- replace(a1,(g1==2),2)
#		a2 <- replace(a2,(g1==0),1)
#		a2 <- replace(a2,(g1==1 | g1==2),2)
#		gdta <- cbind(gdta,a1)
#		gdta <- cbind(gdta,a2)
#	}
#	nams <- c()
#	for (i in 1:from@nsnps) nams <- c(nams,paste(from@snpnames[i],".a1",sep=""),paste(from@snpnames[i],".a2",sep=""))
#	colnames(gdta) <- nams
#	rownames(gdta) <- from@idnames
##	class(gdta) <- "hsgeno"
#	gdta
#	})
#
# end-user class
#
setClass("gwaa.data",representation(phdata="data.frame",gtdata="snp.data"),package="GenABEL")
setMethod("show","gwaa.data",
	function (object) {
		show(object@phdata)
		show(object@gtdata)
	})
#setMethod("summary","gwaa.data",
#	function(object) {
#		ret <- list(phdata = summary(object@phdata),
#			   gtdata = summary(object@gtdata))
#		ret
#	})
setMethod("[","gwaa.data",
	function( x, i, j, drop) {
		if (missing(j)) j=c(1:x@gtdata@nsnps);
		if (missing(i)) i=c(1:x@gtdata@nids);
		if (missing(drop)) drop = FALSE;
		if (is.logical(i)) i=which(i)
		if (is.logical(j)) j=which(j)
		if (is.character(i)) {
			tmp <- i
			i=match(i,x@gtdata@idnames)
			if (any(is.na(i))) stop(paste("following IDs were not found:",tmp[which(is.na(i))],"\n"))
		}
		if (is.character(j)) {
			tmp <- j
			j=match(j,x@gtdata@snpnames)
			if (any(is.na(j))) stop(paste("following SNPs were not found:",tmp[which(is.na(j))],"\n"))
		}
		if (length(i) > x@gtdata@nids  || max(i) > x@gtdata@nids ) stop("i out of range")
		if (length(j) > x@gtdata@nsnps || max(j) > x@gtdata@nsnps) stop("j out of range")
		if (length(i) <= 0) stop("i out of range (== 0)")
		if (length(j) <= 0) stop("j out of range (== 0)")
		a <- x@gtdata[i,j]
		b <- x@phdata[i,]
		out <- new("gwaa.data",phdata=b,gtdata=a)
		out
	})

setClass("scan.gwaa",contains="list",package="GenABEL")
setClass("scan.gwaa.2D",contains="list",package="GenABEL")
