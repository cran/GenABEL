"show.ncbi" <-
function(region) {
  qry <- ""
  for (i in 1:(length(region)-1)) qry <- paste(region[i],region[i+1],sep="+OR+");
  url <- paste("http://www.ncbi.nlm.nih.gov/mapview/map_search.cgi?taxid=9606&query=",qry,"&qchr=&strain=All&advsrch=off",sep="");
  browseURL(url);
}
