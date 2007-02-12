"scan.haplo.2D" <- 
function(y,data,snpsubset,idsubset,bcast=25,simulate=FALSE,...) {
	if (class(data) != "gwaa.data") stop("wrong type of data argument, must be gwaa.data")
	if (missing(snpsubset)) snpsubset <- data@gtdata@snpnames
	if (missing(idsubset)) idsubset <- data@gtdata@idnames
	tdta <- data[idsubset,snpsubset]
	gtdata <- tdta@gtdata
	phdata <- tdta@phdata
	snpsubset <- tdta@gtdata@snpnames
	idsubset <- tdta@gtdata@idnames
	rm(tdta)
	attach(phdata,warn.conflicts=FALSE)
	cc <- get(y)
	detach(phdata)
        if (length(levels(as.factor(cc)))<2) stop("cc status is monomorphic!") 
	oldmap <- gtdata@map
	nsnps <- gtdata@nsnps
        Pv <- matrix(rep(NA,(nsnps*nsnps)),nrow=nsnps)
	name <- snpsubset
	map <- oldmap
	formula <- match.call()
	family <- "haplo.score.slide test"
	ids <- idsubset
	donan <- 0
	for (j1 in 1:(nsnps-1)) {
	for (j2 in (j1+1):nsnps) {
	  twonam <- c(snpsubset[j1],snpsubset[j2])
	  data <- as.hsgeno.snp.data(gtdata[,twonam])
	  tmpo <- haplo.score.slide(y=cc,geno=data,n.slide=2,simulate=simulate,...)
	if (simulate) 
		Pv[j1,j2] <- tmpo$df$global.p.sim
	else 
	  	Pv[j1,j2] <- tmpo$df$score.global.p
	  donan <- donan + 1
	  if (donan/bcast == round(donan/bcast)) print(100*donan/((nsnps-1)*nsnps/2),dig=2)
	}
	}
	rownames(Pv) <- name #[length(name):1]
	colnames(Pv) <- name
  	Pint1df <- rep(NA,length(Pv))
  	Pint2df <- Pint1df
	out <- list(P1df=Pv,Pint1df=Pint1df,Pint2df=Pint2df,P2df=Pv,name=name,formula=formula,family=family,map=map,ids=ids)
	class(out) <- "scan.gwaa.2D"
	out
}
