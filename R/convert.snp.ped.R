"convert.snp.ped" <- function(pedfile,mapfile,outfile,format="premakeped",traits=1,extendedmap=FALSE,strand="u",bcast=10000000) {

	alcodes <- alleleID.codes()

	poss <- c("u","+","-")
	intstrand <- which(poss == strand)-1
	if (length(intstrand)!=1) {cat("strand argument must be one of",poss,"\n");stop();}

# strand == 0 -> unknown ('u')
# strand == 1 -> plus ('+')
# strand == 2 -> minus ('-')

	if (!is.logical(extendedmap)) stop("extendedmap argument must be TRUE or FALSE");
	if (extendedmap) intstrand <- 3;

# extendedmap == TRUE -> five columns (chr,name,pos,strand,coding) expected in map-file
# passed through intstrand <- 3 !

	posf <- c("premakeped","mach")
	intfmt <- which(posf == format)-1
	if (length(intfmt)!=1) {cat("format argument must be one of",posf,"\n");stop();}

# intfmt = 0 -> premakeped, merlin
# else -> mach

	.C("convert_snp_merlin",as.character(pedfile),as.character(mapfile),as.character(outfile),as.integer(intstrand),as.integer(bcast),as.character(alcodes),as.integer(length(alcodes)),as.integer(intfmt),as.integer(traits),PACKAGE="GenABEL")
    	return(invisible(0))
}

