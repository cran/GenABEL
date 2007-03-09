"hom" <- 
function(data,snpsubset,idsubset,weight="no") {
  if (class(data) == "gwaa.data") {
    data <- data@gtdata
  }
  if (class(data) != "snp.data") {
    stop("data argument should have gwaa.data-class or snp.data-class")
  }
  wargs <- c("no","freq")
  if (!(match(weight,wargs,nomatch=0)>0)) {
    out <- paste("weight argument should be one of",wargs,"\n")
    stop(out)
  }

  if (!missing(snpsubset)) data <- data[,snpsubset]
  if (!missing(idsubset)) data <- data[idsubset,]

  totsnps <- data@nsnps
  totids <- data@nids
  out <- .C("hom",as.raw(data@gtps),as.integer(data@nids),as.integer(data@nsnps),as.integer(match(weight,wargs)-1),sout = double(2*data@nids), PACKAGE="GenABEL")$sout
  dim(out) <- c(data@nids,2)
  rownames(out) <- data@idnames
  colnames(out) <- c("NoMeasured","Hom")
  out[,2] <- out[,2]/out[,1]
  out
}
