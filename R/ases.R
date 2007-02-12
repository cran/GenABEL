#"as.double" <- 
#function(x, ...) UseMethod("as.double")

"as.double.snp.data" <- 
function(x, ...) {
		tnids <- x@nids
		tnsnps <- x@nsnps
		to <- .C("get_snps_many",as.raw(x@gtps), as.integer(tnids), as.integer(tnsnps), idata = integer(tnids*tnsnps), PACKAGE="GenABEL")$idata
		to <- replace(to,(to==0),NA)
		to <- to - 1
		dim(to) <- c(tnids,tnsnps)
		colnames(to) <- x@snpnames
		rownames(to) <- x@idnames
		to
}

"as.character.snp.data" <- 
function(x, ...) {
	a <- as.double(x)
	dm <- dim(a)
	to <- ifelse(is.na(a),"",c("A/A","A/B","B/B")[a+1])
	dim(to) <- dm
	colnames(to) <- x@snpnames
	rownames(to) <- x@idnames
	to
}

"as.genotype" <- 
function(x, ...) UseMethod("as.genotype")

"as.genotype.snp.data" <- 
function(x, ...) {
	gdta <- data.frame(genotype(as.character(x[,1])))
	if (x@nsnps>1) for (i in 2:x@nsnps) {
		gdta <- cbind(gdta,genotype(as.character(x[,i])))
	}
	colnames(gdta) <- x@snpnames
#	class(gdta) <- "genotype"
	gdta
}

"as.hsgeno" <- 
function(x, ...) UseMethod("as.hsgeno")

"as.hsgeno.snp.data" <-
function(x, ...) {
	g1 <- as.double(x[,1])
	a1 <- rep(NA,length(g1))
	a2 <- rep(NA,length(g1))
	a1 <- replace(a1,(g1==0 | g1==1),1)
	a1 <- replace(a1,(g1==2),2)
	a2 <- replace(a2,(g1==0),1)
	a2 <- replace(a2,(g1==1 | g1==2),2)
	gdta <- data.frame(a1,a2)
	if (x@nsnps>1) for (i in 2:x@nsnps) {
		g1 <- as.double(x[,i])
		a1 <- rep(NA,length(g1))
		a2 <- rep(NA,length(g1))
		a1 <- replace(a1,(g1==0 | g1==1),1)
		a1 <- replace(a1,(g1==2),2)
		a2 <- replace(a2,(g1==0),1)
		a2 <- replace(a2,(g1==1 | g1==2),2)
		gdta <- cbind(gdta,a1)
		gdta <- cbind(gdta,a2)
	}
	nams <- c()
	for (i in 1:x@nsnps) nams <- c(nams,paste(x@snpnames[i],".a1",sep=""),paste(x@snpnames[i],".a2",sep=""))
	colnames(gdta) <- nams
	rownames(gdta) <- x@idnames
#	class(gdta) <- "hsgeno"
	gdta
}

"as.raw.snp.data" <- 
function (x) {
	to <- as.raw(from@gtps)
	to
}

"as.raw.snp.mx" <-
function(x) {
	to <- unclass(x)
	to <- as.raw(to)
	to
}
