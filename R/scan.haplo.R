"scan.haplo" <- 
function(y,data,snpsubset,idsubset,n.slide=2,bcast=25,simulate=FALSE,...) {
	if (class(data) != "gwaa.data") stop("wrong type of data argument, must be gwaa.data")
	if (missing(snpsubset)) snpsubset <- data@gtdata@snpnames
	if (missing(idsubset)) idsubset <- data@gtdata@idnames
	gtdata <- data@gtdata[idsubset,snpsubset]
	idsubset <- gtdata@idnames
	phdata <- data@phdata[match(idsubset,data@phdata$id),]
	attach(phdata,warn.conflicts=FALSE)
	cc <- get(y)
	detach(phdata)
        if (length(levels(as.factor(cc)))<2) stop("cc status is monomorphic!") 
	oldmap <- gtdata@map
	nsnps <- gtdata@nsnps
        Pv <- rep(1,nsnps-n.slide+1)
	name <- rep("",nsnps-n.slide+1)
	map <- rep(-1,nsnps-n.slide+1)
	for (i in 1:(nsnps-n.slide+1)) {
	cumm <- oldmap[i]
	name[i] <- snpsubset[i];
	for (j in (i+1):(i+n.slide-1)) {
		name[i]=paste(name[i],snpsubset[j],sep="-");
		cumm <- cumm + oldmap[j]
	}
	map[i] <- cumm/n.slide
	}
	formula <- match.call()
	family <- "haplo.score.slide test"
	ids <- idsubset
	for (j in 1:(nsnps-n.slide+1)) {
	  data <- as.hsgeno.snp.data(gtdata[,j:(j+n.slide-1)])
	  tmpo <- haplo.score.slide(y=cc,geno=data,n.slide=n.slide,simulate=simulate,...)
	if (simulate) 
	  	Pv[j] <- tmpo$df$global.p.sim
	else 
	 	 Pv[j] <- tmpo$df$score.global.p
	  if (j/bcast == round(j/bcast)) print(j)
	}
	out <- list(P1df=Pv,name=name,formula=formula,family=family,map=map,ids=ids)
	class(out) <- "scan.gwaa"
	out
}
