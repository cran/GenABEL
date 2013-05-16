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
	nTestIds <- sample(c(10:min(100,nids(ge03d2.clean))),1)
	nTestSnps <- sample(c(10:min(1000,nsnps(ge03d2.clean))),1)
	dta <- ge03d2.clean[sort(sample(1:nids(ge03d2.clean),nTestIds)),
			sort(sample(1:nsnps(ge03d2.clean),nTestSnps))]
	
	export.plink(dta,filebasename="tmpOld",dpieceFun="old")
	export.plink(dta,filebasename="tmpNew",dpieceFun="new")
	xO <- read.table(file="tmpOld.ped",head=FALSE,strings=FALSE)
	xN <- read.table(file="tmpNew.ped",head=FALSE,strings=FALSE)
	checkIdentical(xN,xO)
	xO <- read.table(file="tmpOld.map",head=FALSE,strings=FALSE)
	xN <- read.table(file="tmpNew.map",head=FALSE,strings=FALSE)
	checkIdentical(xN,xO)
	xO <- read.table(file="tmpOld.map.ext",head=FALSE,strings=FALSE)
	xN <- read.table(file="tmpNew.map.ext",head=FALSE,strings=FALSE)
	checkIdentical(xN,xO)
	
	unlink("tmpOld*")
	unlink("tmpNew*")
	
	export.merlin(dta,pedfile="tmpOld.ped",datafile="tmpOld.dat",
			mapfile="tmpOld.map",dpieceFun="old")
	export.merlin(dta,pedfile="tmpNew.ped",datafile="tmpNew.dat",
			mapfile="tmpNew.map",dpieceFun="new")
	xO <- read.table(file="tmpOld.ped",head=FALSE,strings=FALSE)
	xN <- read.table(file="tmpNew.ped",head=FALSE,strings=FALSE)
	checkIdentical(xN,xO)
	xO <- read.table(file="tmpOld.map",head=FALSE,strings=FALSE)
	xN <- read.table(file="tmpNew.map",head=FALSE,strings=FALSE)
	checkIdentical(xN,xO)
	xO <- read.table(file="tmpOld.map.ext",head=FALSE,strings=FALSE)
	xN <- read.table(file="tmpNew.map.ext",head=FALSE,strings=FALSE)
	checkIdentical(xN,xO)
	xO <- read.table(file="tmpOld.dat",head=FALSE,strings=FALSE)
	xN <- read.table(file="tmpNew.dat",head=FALSE,strings=FALSE)
	checkIdentical(xN,xO)
	
	unlink("tmpOld*")
	unlink("tmpNew*")	
	
	export.plink(dta,filebasename="tmpTrans",transpose=TRUE)
	convert.snp.tped(tpedfile="tmpTrans.tped",tfamfile="tmpTrans.tfam",out="tmpTrans.raw")
	xBack <- load.gwaa.data(gen="tmpTrans.raw",phe="tmpTrans.phe",id="IID")
	strand(xBack) <- strand(dta)
	phdata(xBack) <- phdata(xBack)[,-1]
	sameCode <- which(coding(dta) == coding(xBack))
	checkIdentical(xBack[,sameCode],dta[,sameCode])
	swappedCode <- which(refallele(dta)==effallele(xBack) & effallele(dta)==refallele(xBack))
	for (i in swappedCode) {
		gtDta <- c(2,1,0)[as.vector(as.numeric(dta[,i]))+1]
		gtBack <- as.vector(as.numeric(xBack[,i]))
		print(checkIdentical(gtBack,gtDta))
	}
	
	unlink("tmpTrans*")
	
}

test.export.merlin.bug2525 <- function()
{
	data(srdta)
	export.merlin(
			srdta[,1:2],dpieceFun="new",
			mapfile="tmpNew.map",pedfile="tmpNew.ped",datafile="tmpNew.dat"
	)
	export.merlin(
			srdta[,1:2],dpieceFun="old",
			mapfile="tmpOld.map",pedfile="tmpOld.ped",datafile="tmpOld.dat"
	)
	xO <- read.table(file="tmpOld.ped",head=FALSE,strings=FALSE)
	xN <- read.table(file="tmpNew.ped",head=FALSE,strings=FALSE)
	checkIdentical(xN,xO)
}
