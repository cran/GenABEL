"show.ncbi" <-
function(region) {
  qry <- ""
  if (length(region)<=1) qry <- region else {
  	qry <- region[1]
  	for (i in 2:length(region)) qry <- paste(qry,region[i],sep="+OR+");
  }
  url <- paste("http://www.ncbi.nlm.nih.gov/mapview/map_search.cgi?taxid=9606&query=",qry,"&qchr=&strain=All&advsrch=off",sep="");
  browseURL(url);
}
