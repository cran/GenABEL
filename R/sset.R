"sset" <-
function(data,nsnps,nids,list) {
  result <-.Call("sset_call",as.raw(data), 
               as.integer(nsnps), 
               as.integer(nids), 
               as.integer(list), 
               as.integer(length(list)))
  return(result)
}
