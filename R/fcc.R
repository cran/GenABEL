"fcc" <-
function(data,pheno) {
	if (class(data) == "gwaa.data") data <- data@gtdata else if (class(data) != "snp.data") stop ("Wrong argument (should be \"gwaa.data\" or \"snp.data\")")
	out <- .C("fastcc",as.raw(data@gtps),as.integer(pheno),as.integer(data@nids),as.integer(data@nsnps), chi2 = double(6*data@nsnps), PACKAGE="GenABEL")$chi2
	out
}

