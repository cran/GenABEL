"merge.snp.data" <-
function(set1, set2, whichset=1, replacena="no") {

#check whether set1 and set2 have required data class
if (class(set1)!="snp.data")
	{
	stop("Wrong data class: the first argument should be snp.data")
	}																	    

if (class(set2)!="snp.data")
	{
	stop("Wrong data class: the second argument should be snp.data")
	}																	    



if( !(whichset == 1 || whichset == 2) ) 
	{
	stop("Wrong whichset vallue. Input variable whichset must be 1 (default) or 2.")
	}

if( !(replacena == "no" || replacena == "yes") )
	{
	stop("Wrong ignorna vallue. Input variable replacena must be \"yes\" (default) or \"no\".")
	}









#What set has minimum amount of IDs and what is this amount? 
#if(set1@nids > set2@nids)	min_ids_amount <- set2@nids	else min_ids_amount <- set1@nids

cat("replacena=",replacena, "\n")

which_id_intersect_in_set1 <- which(is.element(set1@idnames,set2@idnames)*c(1:set1@nids) > 0) 
#cat("1\n")
which_id_intersect_in_set2 <- which(is.element(set2@idnames,set1@idnames)*c(1:set2@nids) > 0)
#Numbers in array "which_id_intersect_in_set1" (in set1) is intersected with numbers in array "which_id_intersect_in_set2" (in set2)
#For example: set1@idnames[which_id_intersect_in_set1[1]] and set2@idnames[which_id_intersect_in_set2[1]] is same ID.
#cat("2\n")

num_ids_intersected <- length(which_id_intersect_in_set1)	#amount of intersected subset


#cat("3\n")
ids_intersected <- c(which_id_intersect_in_set1, which_id_intersect_in_set2) #Final array which will be passed to c-function
#cat("4\n")







#Below is same procedure for SNPs.
#What set has minimum amount of SNPs and what is this amount?
#if(set1@nsnps > set2@nsnps)	min_snp_amount <- set2@nsnps	else min_snp_amount <- set1@nsnps;



which_snp_intersect_in_set1 <- which(is.element(set1@snpnames,set2@snpnames)*c(1:set1@nsnps) > 0)
which_snp_intersect_in_set2 <- which(is.element(set2@snpnames,set1@snpnames)*c(1:set2@nsnps) > 0)

#cat("which_snp_intersect_in_set1=", which_snp_intersect_in_set1, "\n")
#cat("which_snp_intersect_in_set2=", which_snp_intersect_in_set2, "\n")


num_snps_intersected <- length(which_snp_intersect_in_set1)
snp_intersected <- c(which_snp_intersect_in_set1, which_snp_intersect_in_set2) 




#cat("\n num_snps_intersected", num_snps_intersected, "\n")
#cat("\n snp_intersected", snp_intersected, "\n")



#cat("set1@nids =",set1@nids, "\n")
#cat("set2@nids =",set2@nids, "\n")
#cat("num_ids_intersected =", num_ids_intersected, "\n")
#cat("ceiling(set1@nids + set2@nids - num_ids_intersected)/4 =",ceiling(set1@nids + set2@nids - num_ids_intersected)/4, "\n")

id_amount_in_new_array <- set1@nids + set2@nids - num_ids_intersected;
byte_amount_for_one_snp_in_new_array <- ceiling((id_amount_in_new_array)/4)

snps_amount_in_new_array <- set1@nsnps + set2@nsnps - num_snps_intersected
byte_number_in_new_array <- byte_amount_for_one_snp_in_new_array * (snps_amount_in_new_array)


#byte_number_in_new_array = 10000;
cat("byte_number_in_new_array =", byte_number_in_new_array, "\n")

#new_array_raw <- raw(byte_number_in_new_array)





new_array_raw <-  .C("fast_merge_C_",
								as.raw(set1@gtps), as.integer(set1@nids), as.integer(set1@nsnps), 
								as.raw(set2@gtps), as.integer(set2@nids), as.integer(set2@nsnps),
							  as.integer(num_ids_intersected), as.integer(num_snps_intersected), as.integer(snp_intersected), as.integer(ids_intersected),
#								as.integer(set1@strand), as.integer(set2@strand),
								as.integer(whichset), as.character(replacena),
								return_val = raw(byte_number_in_new_array))$return_val




#cat("as.raw(set1@gtps)=", as.raw(set1@gtps), "\n")
#cat("as.raw(set2@gtps)=", as.raw(set2@gtps), "\n")
#cat("new_array_raw=", new_array_raw, "\n")


#cat("set1@gtps=", set1@gtps, "\n")


dim(new_array_raw) <- c(byte_amount_for_one_snp_in_new_array, snps_amount_in_new_array)


new_array_raw <- new("snp.mx",new_array_raw);gc(verbose=FALSE)




if(whichset == 1)
	{
	tmp_logic_array <- c()
	#Form array with ID names (without intersection)
	#-------------------------------------------------------------------
	tmp_logic_array[c(1:length(set2@idnames))] = TRUE
	tmp_logic_array[which_id_intersect_in_set2] = FALSE

	set2_idnames_without_intersected <- set2@idnames[tmp_logic_array]
	ids_array <- c(set1@idnames, set2_idnames_without_intersected)
	#-------------------------------------------------------------------

	#Form array with males (without intersection)
	#-------------------------------------------------------------------
	set2_male_without_intersected <- set2@male[tmp_logic_array]
	male_array <- c(set1@male, set2_male_without_intersected)
	#-------------------------------------------------------------------


	#Form array with SNP names (without intersection)
	#-------------------------------------------------------------------
	tmp_logic_array <- c()
	tmp_logic_array[c(1:length(set2@snpnames))] = TRUE
	tmp_logic_array[which_snp_intersect_in_set2] = FALSE

	set2_snps_without_intersected <- set2@snpnames[tmp_logic_array]
	snp_array <- c(set1@snpnames, set2_snps_without_intersected)
	#-------------------------------------------------------------------




	#Form array with chromosomes (without intersection)
	#-------------------------------------------------------------------
	set2_chromosome_without_intersected <- as.character(set2@chromosome[tmp_logic_array])
	chrom_array <- factor(c(as.character(set1@chromosome), set2_chromosome_without_intersected))
#	cat("class(chrom_array)=", class(chrom_array), "\n")
	#chrom_array <- new()
	#new_array_raw <- new("snp.mx",new_array_raw);gc(verbose=FALSE)
	#-------------------------------------------------------------------


	#Form array with map (without intersection)
	#-------------------------------------------------------------------
	set2_map_without_intersected <- set2@map[tmp_logic_array]
	map_array <- c(set1@map, set2_map_without_intersected)
	#-------------------------------------------------------------------






	#-------------------------------------------------------------------
	set2_coding_without_intersected <- set2@coding[tmp_logic_array]
	coding_array <- new("snp.coding", c(set1@coding, set2_coding_without_intersected))
	#-------------------------------------------------------------------


	#-------------------------------------------------------------------
	set2_strand_without_intersected <- set2@strand[tmp_logic_array]
	strand_array <- new("snp.strand", c(set1@strand, set2_strand_without_intersected))
	#-------------------------------------------------------------------

	}
else
	{
	tmp_logic_array <- c()
	#Form array with ID names (without intersection)
	#-------------------------------------------------------------------
	tmp_logic_array[c(1:length(set1@idnames))] = TRUE
	tmp_logic_array[which_id_intersect_in_set1] = FALSE

	set1_idnames_without_intersected <- set1@idnames[tmp_logic_array]
	ids_array <- c(set1_idnames_without_intersected, set2@idnames)
	#-------------------------------------------------------------------

	#Form array with males (without intersection)
	#-------------------------------------------------------------------
	set1_male_without_intersected <- set1@male[tmp_logic_array]
	male_array <- c(set1_male_without_intersected, set2@male)
	#-------------------------------------------------------------------


	#Form array with SNP names (without intersection)
	#-------------------------------------------------------------------
	tmp_logic_array <- c()
	tmp_logic_array[c(1:length(set1@snpnames))] = TRUE
	tmp_logic_array[which_snp_intersect_in_set1] = FALSE

	set1_snps_without_intersected <- set1@snpnames[tmp_logic_array]
	snp_array <- c(set1_snps_without_intersected, set2@snpnames)
	#-------------------------------------------------------------------




	#Form array with chromosomes (without intersection)
	#-------------------------------------------------------------------
	set1_chromosome_without_intersected <- set1@chromosome[tmp_logic_array]
	chrom_array <- factor(c(set1_chromosome_without_intersected, set2@chromosome))
#	cat("class(chrom_array)=", class(chrom_array), "\n")
	#chrom_array <- new()
	#new_array_raw <- new("snp.mx",new_array_raw);gc(verbose=FALSE)
	#-------------------------------------------------------------------


	#Form array with map (without intersection)
	#-------------------------------------------------------------------
	set1_map_without_intersected <- set1@map[tmp_logic_array]
	map_array <- c(set1_map_without_intersected, set2@map)
	#-------------------------------------------------------------------






	#-------------------------------------------------------------------
	set1_coding_without_intersected <- set1@coding[tmp_logic_array]
	coding_array <- new("snp.coding", c(set1_coding_without_intersected, set2@coding))
	#-------------------------------------------------------------------


	#-------------------------------------------------------------------
	set1_strand_without_intersected <- set1@strand[tmp_logic_array]
	strand_array <- new("snp.strand", c(set1_strand_without_intersected, set2@strand))
	#-------------------------------------------------------------------
	
	
	}





mearged_set_snp_data  <- snp.data(nids=id_amount_in_new_array, rawdata=new_array_raw, idnames=ids_array,
			snpnames=snp_array, chromosome=chrom_array, map=map_array, coding=coding_array,
			strand=strand_array, male=male_array)





#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#rm(rdta,ids,mnams,chrom,pos,coding,strand);gc(verbose=FALSE)
return(mearged_set_snp_data)

}











#-------------------------------------------------
#	a <- snp.data(nids=nids,rawdata=rdta,idnames=ids,snpnames=mnams,chromosome=chrom,map=pos,coding=coding,strand=strand,male=newdta$sex)
#	Exapmle is in R/load.gwaa.data.R
# Create object with class snp.data. rdta - object with class snp.mx

#		.C("qtscore",as.raw(data@gtdata@gtps),as.double(resid),as.integer(bin),as.integer(data@gtdata@nids),as.integer(data@gtdata@nsnps), as.integer(nstra), as.integer(strata), chi2 = double(6*data@gtdata@nsnps), PACKAGE="GenABEL")$chi2
#-------------------------------------------------
