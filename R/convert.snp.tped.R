"convert.snp.tped" <-
  function(tpedfile,tfamfile,outfile,bcast=10000) {
    .C("convert_snp_tped",as.character(tpedfile),as.character(tfamfile),as.character(outfile),as.integer(bcast),PACKAGE="GenABEL")
    return(invisible(0))
  }
