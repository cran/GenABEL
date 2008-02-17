alleleID.alleles <- function() {
	a <- list();
	a[[1]] <- c("1","2")
	a[[2]] <- c("A","B")
	alleles <- c("A","T","G","C","-")
	idx <- 3
	for (i in alleles) {
	for (j in alleles) {
		if (i==j) next;
		a[[idx]] <- c(i,j)
		idx <- idx + 1
	}
	}
	a[[idx]] <- c("2","1")
	idx <- idx + 1
	a[[idx]] <- c("B","A")
	a
}
alleleID.codes <- function() {
	a <- alleleID.alleles()
	out <- c("OPPA")
	idx <- 1
	for (i in a) {
		out[idx] <- paste(i[1],i[2],sep="")
		idx <- idx + 1
	}
	out
}
alleleID.codes.reverse <- function() {
	a <- alleleID.alleles()
	out <- c("OPPA")
	idx <- 1
	alleles <- c("A","T","G","C","-")
	reval <- c("T","A","C","G","-")
	for (i in a) {
		x <- match(alleles,i)
		if (sum(!is.na(x)) == 2) {
			names(reval) <- alleles
			xx <- reval[i]
			out[idx] <- paste(xx[1],xx[2],sep="")
		} else {
			out[idx] <- paste(i[2],i[1],sep="")
		}
		idx <- idx + 1
	}
	out
}
alleleID.raw2char.matrix <- function() {
	alleles <- alleleID.alleles()
	idx <- 1
	coding <- c("OPPA")
	for (i in alleles) {
		coding[idx] <- paste(i[1],i[1],sep="/");idx <- idx + 1
		coding[idx] <- paste(i[1],i[2],sep="/");idx <- idx + 1
		coding[idx] <- paste(i[2],i[2],sep="/");idx <- idx + 1
	}
	coding <- matrix(coding,nrow=3)
	coding <- t(coding)
	out <- coding
	colnames(out) <- c("0","1","2")
	rownames(out) <- as.character(as.raw(c(1:dim(out)[1])))
	out
}
alleleID.char2raw <- function() {
	coding <- alleleID.codes()
	out <- as.raw(c(1:length(coding)))
	names(out) <- coding
	out
}
alleleID.raw2char <- function() {
	x <- alleleID.char2raw()
	out <- names(x)
	names(out) <- x
	out
}
alleleID.revstrand <- function() {
	coding1 <- alleleID.codes()
	coding2 <- alleleID.codes.reverse()
	out1 <- out2 <- as.raw(c(1:length(coding1)))
	names(out1) <- coding1
	names(out2) <- coding2
	oo <- out1[names(out2)]
	names(oo) <- out2
	oo
}

alleleID.reference <- function() {
	codes <- alleleID.codes()
	out <- codes
	x <- strsplit(codes,"")
	for (i in c(1:length(codes))) {
		out[i] <- x[[i]][1]
	}
	names(out) <- codes
	out
}

alleleID.effective <- function() {
	codes <- alleleID.codes()
	out <- codes
	x <- strsplit(codes,"")
	for (i in c(1:length(codes))) {
		out[i] <- x[[i]][2]
	}
	names(out) <- codes
	out
}


