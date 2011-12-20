### --- Test setup ---
#
# regression tests
#

if(FALSE) {
	## Not really needed, but can be handy when writing tests
	library(RUnit)
	library(GenABEL)
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source(paste("../inst/unitTests/shared_functions.R"))
#source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.exports <- function()
{
	data(ge03d2.clean)
	dta <- ge03d2.clean[1:200,autosomal(ge03d2.clean)[1:3000]]
	gkin <- ibs(dta[,1:2000],w="freq")
	h2 <- polygenic(height,data=dta,kin=gkin)
	xOld <- mmscore(h2,dta,cppFun="old")
	xNew <- mmscore(h2,dta,cppFun="new")
	checkEquals(results(xNew),results(xOld))
}
