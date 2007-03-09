"summary.snp.data" <-
	function(object,...) {
		if (class(object) != "snp.data") stop("wrong object class, should be snp.data")
		res <- .C("snp_summary_exhwe",as.raw(object@gtps),as.integer(object@nids),as.integer(object@nsnps), out = double(object@nsnps*7), PACKAGE="GenABEL")$out
		dim(res) <- c(object@nsnps,7)
		res <- as.data.frame(res)
		res <- cbind(res,as.factor(object@chromosome))
		if (any(object@chromosome == "X")) {
		  oX <- object[(object@male != 1),(object@chromosome == "X")]
		  resX <- .C("snp_summary_exhwe",as.raw(oX@gtps),as.integer(oX@nids),as.integer(oX@nsnps), out = double(oX@nsnps*7), PACKAGE="GenABEL")$out
		  vec <- (object@chromosome == "X")
		  res[vec,7] <- resX[(oX@nsnps*6+1):(oX@nsnps*7)]
		  rm(oX,vec,resX);gc(verbose=FALSE)
		}
		rownames(res) <- object@snpnames
		colnames(res) <- c("NoMeasured","CallRate","Q.2","P.11","P.12","P.22","Pexact","Chromosome")
		res
	}
