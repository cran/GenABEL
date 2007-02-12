"heteroz" <- 
function(x,snpsubset,idsubset, bcast=10) {
  if (class(x) == "gwaa.data") {
    gdta <- x@gtdata
  } else if (class(x) == "snp.data") {
    gdta <- x
  } else {
    stop("argument should have gwaa.data-class or snp.data-class")
  }
  if (missing(snpsubset)) snpsubset <- gdta@snpnames
  if (missing(idsubset)) idsubset <- gdta@idnames
  gdta <- gdta[idsubset,snpsubset]
  totsnps <- length(snpsubset)
  totids <- length(idsubset)
  out <- new.env()
  for (i in 1:totids) {
    y <- as.vector(as.numeric(gdta[idsubset[i],]))
    y <- replace(y,y==2,0)
    out$mean[i] <- mean(y,na.rm=T)
    out$sd[i] <- sd(y,na.rm=T)
    nmea <- sum(!is.na(y))
    out$nmea[i] <- nmea
    out$call[i] <- nmea/totsnps
    out$sem[i] <- out$sd[i]/sqrt(nmea)
    if (i/bcast == round(i/bcast)) print(100*i/totids,dig=2);
  }
  out$idsubset <- idsubset
  out$snpsubset <- snpsubset
  out <- as.list(out)
  out
}
