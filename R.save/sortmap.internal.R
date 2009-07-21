"sortmap.internal" <-
function(chrom,map,delta=1) {
	chnum <- chrom.char2num(chrom)

	ix <- order(chnum,map)

	map <- map[ix]

   off <- c(0,map[1:(length(map)-1)])
	off <- map - off
	off[which(off<=0)] <- delta
	cummap <- cumsum(off)

	out <- list()
	out$ix <- ix
	out$cummap <- cummap
	out$chnum <- chnum
	out
}

chrom.char2num <- function(chrom) {
	chrom <- as.character(chrom)
   chdesU <- unique(chrom)
   chnumU <- suppressWarnings(as.numeric(chdesU))
   chCHU <- sort(chdesU[is.na(chnumU)])
	maxch <- max(chnumU,na.rm=T)
   j <- 1
	for (i in chCHU) {
		chrom[which(chrom==i)] <- as.character(maxch+j)
		j <- j + 1
	}
	chnum <- as.numeric(chrom)
	chnum
}
