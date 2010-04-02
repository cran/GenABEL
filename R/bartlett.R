#=====================================================================================
#
#       Filename:  R/bartlett.R
#
#    Description:  Bartlett's test of variance homogeneity
#
#        Version:  1.0
#        Created:  
#       Revision:  none
#       
#
#         Author:  Maksim V. Struchalin
#        Company:  ErasmusMC, Epidemiology, The Netherlands.
#          Email:  m.struchalin@erasmusmc.nl
#
#=====================================================================================


"bartlett" <-
function(trait, data, analysis_type) {#, trait.type = "gaussian") {


#if(class(data) != "gwaa.data" && class(data) != "snp.data")
#	{
#	stop("Wrong data class: data input parameter must be \"gwaa.data\" or \"snp.data\"")
#	}																	    

#class_formula <- class(formula)

#if(class_formula != "formula" || class_formula != "numeric" || class_formula != "integer" || class_formula != "double")
#	{
#	stop("Wrong formula class: formula must have one of the following types: \"formula\", \"numeric\", \"integer\", \"double\"")
#	}

#if(class(data) == "snp.data" && class(formula) == "formula" )
#	{
#	stop("Wrong types of input parameters. data can not have type \"snp.data\" and formula - \"formula\" ")
#	}

#if(class(trait.type) != "character")
#	{
#	stop("Wrong types of trait.type. It must be \"character\".")
#	}


#if(class(trait.type) != "gaussian")
#	{
#	stop("error: gaussian is only support")
#	}







if(length(trait) != data@gtdata@nids)
	{
	stop("Length of trait does not match to ids amount")
	}



is_trait_na <- is.na(trait)


trait[is_trait_na] <- 0


is_trait_na <- as.integer(is_trait_na)




print("start calculating")

retrun_val <- .C("bartlett_test_R", as.raw(data@gtdata@gtps), as.integer(data@gtdata@nids), as.integer(data@gtdata@nsnps),
			 														as.double(trait), as.integer(is_trait_na),
																	chi2 = double(data@gtdata@nsnps), analysis_type)


chi2 <- retrun_val$chi2


print("all snps done")

chi2[chi2 == -1] <- NA


return(chi2)

}
