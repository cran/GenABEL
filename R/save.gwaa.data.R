"save.gwaa.data" <-
function(data,phenofile = "pheno.dat", genofile = "geno.raw",human=FALSE) {
	print("writing phenotypes...")
	write.table(data@phdata,file=phenofile,sep=" ",na="NA",row.names=FALSE,col.names=TRUE)
	if (human==TRUE && missing(genofile)) genofile="geno.dat"
	ofile <- file(genofile,"w")
	print("writing id names...")
	cat(data@gtdata@idnames,"\n",file=ofile)
	print("writing snp names...")
	cat(data@gtdata@snpnames,"\n",file=ofile)
	print("writing chromosome data...")
	cat(as.character(data@gtdata@chromosome),"\n",file=ofile)
	print("writing map data...")
	cat(data@gtdata@map,"\n",file=ofile)
	print("writing genotypic data...")
	if (human) {
		idta <- as.double(data@gtdata)+1
		idta <- replace(idta,is.na(idta),0)
		write(idta,file=ofile,ncolumns=data@gtdata@nids)
	} else 
		write(data@gtdata@gtps,file=ofile,ncolumns=data@gtdata@nbytes)
#	for (i in 1:data@gtdata@nsnps) {
#		cat(data@gtdata@gtps[,i],"\n",file=ofile)
#	}
	close(ofile)
}

