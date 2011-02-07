#=====================================================================================
#
#       Filename:  chi2-CG.R
#
#    Description:  Functions for performing chi2 test for Detecting Rare Recessive Alleles.
#
#        Version:  1.0
#        Created:  27-May-2010
#       Revision:  none
#       
#
#         Author:  Maksim V. Struchalin
#        Company:  ErasmusMC, Epidemiology, The Netherlands.
#          Email:  m.struchalin@erasmusmc.nl
#
#=====================================================================================


"cocohet" <- function(data, trait, window, return_all_result=TRUE, manhettan_filename="manhettan_plot.jpeg", test="CHI2", min_expected_cut_off=-1)
{
#Correction is performed in case expected value one of the margin in contigency table less 5. Available values for chi2_correction:
# "CHI2"
# "YATES"
# "FISHER"


if(!is(data,"snp.data"))
	{
	stop("Wrong data class: the argument \"data\" should has type \"snp.data\".")
	}																	    



if(!is(trait, "integer") && !is(trait, "numeric"))
	{
	stop("Wrong data class: the argument \"trait\" should has type \"integer\" or \"numeric\".")
	}																	    

if(dim(table(trait)) != 2)
	{
	stop("\"trait\" must be binary")
	}


if(!is(window, "integer") && !is(window, "numeric"))
	{
	stop("Wrong data class: the argument \"window\" should has type \"integer\" or \"numeric\".")
	}																	    



snp_number <- data@nsnps


if(window >= snp_number)
	{
	print(paste("warning: \"window\" value ", window, " is bigger than total snp number. Automaticly reduced to maximum possible."))
	window <- snp_number - 1
	print(paste("New window value is ", window, ".\n"))
	}

if(window==0) {stop("Input parameter window can not be zero.\n")}



output <- .Call("interaction_rare_recesive_allele_C_", as.raw(data@gtps), as.integer(data@nids), as.integer(snp_number),
			 																								 as.integer(trait),
																											 as.integer(window), 
																											 return_all_result,
																											 test,
																											 as.integer(min_expected_cut_off))



results <- list()

results$chi2_max <- output[1:(data@nsnps-1)]
if(return_all_result)
	{
	results$chi2_all <- matrix(output[data@nsnps:length(output)], ncol=window, nrow=data@nsnps-1, byrow=T)
	}



print("plotting manhettan plot picture...")
bitmap(manhettan_filename, type="jpeg")
results_pval <- pchisq(results$chi2_max, df=1, lower.tail=F)
results_list <- list()
results_list$P1df <- results_pval
results_list$map <- data@map[1:(data@nsnps-1)]
results_list$chromosome <- data@chromosome[1:(data@nsnps-1)]
class(results_list) <- "scan.gwaa"
print("1")
try(plot.scan.gwaa(results_list))
print("2")
dev.off()
print(paste("The plot is saved into file ", manhettan_filename, sep=""))



snp_position_top <- which.max(results$chi2_max)
snp_name_top <- data@snpnames[snp_position_top]
snp_chromosome_top <- data@chromosome[snp_position_top]
snp_map_top <- data@map[snp_position_top]

print(paste("The top snp has name ", snp_name_top, ",is located in chromosome ", snp_chromosome_top, ", position ", snp_map_top, sep=""))

#if(return_all_result)
#	{
#	print(paste("Plotting the manhettan plot for interaction bewteen snp ", snp_name_top, " and snps on the left/right within a window", sep=""))
#	look_at_snp_chi2s_interactions_rare_recesive_alleles(results$chi2_all, data, snp_position_top)
#	}

results
}






"look_at_snp_chi2s_interactions_rare_recesive_alleles" <- function(chi2_all, data, snp_position, manhettan_filename="one_snp_manhettan.jpeg")
{

if(!is(data,"snp.data"))
	{
	stop("Wrong data class: the argument \"data\" should has type \"snp.data\".")
	}																	    


if(!is(chi2_all,"matrix"))
	{
	stop("Wrong data class: the argument \"chi2_all\" should has type \"matrix\".")
	} 




window <- dim(chi2_all)[2]
snp_number <- dim(chi2_all)[1] + 1

if(window == 0 | snp_number==0) stop("Wrong chi2 matrix.")

results <- list()
results$right_window_chi2s <- chi2_all[snp_position,]
results$left_window_chi2s <- c()

left_border <- snp_position - window
if(left_border <=0 ) left_border <- 1

right_border <- snp_position + window
if(right_border > snp_number-1) right_border <- snp_number-1

for(i in 1:left_border)
	{
	results$left_window_chi2s <- c(results$left_window_chi2s, chi2_all[snp_position-i, i])
	}

bitmap(manhettan_filename, type="jpeg")

print(snp_position)		
results_chi2 <- c(results$left_window_chi2s, results$right_window_chi2s)
print(length(results$left_window_chi2s))
print(length(results$right_window_chi2s))
print(length(results_chi2))
print("________________")
results_pval <- pchisq(results_chi2, df=1, lower.tail=F)

max_position_local <- which.max(results_chi2)
max_position <- max_position_local + snp_position

results_list <- list()
results_list$P1df <- results_pval
results_list$map <- c(data@map[left_border:(snp_position-1)], data@map[(snp_position+1):right_border])
results_list$chromosome <- as.factor(c(as.character(data@chromosome[left_border:(snp_position-1)]),
			 	as.character(data@chromosome[(snp_position+1):right_border])))
class(results_list) <- "scan.gwaa"


print(length(results_list$P1df))
print(length(results_list$map))
print(length(results_list$chromosome))

try(plot.scan.gwaa(results_list))

dev.off()

print(paste("The plot is saved into file ", manhettan_filename, sep=""))

print(paste("Maximum chi2 is ", results_chi2[max_position_local], " for a snp in position ", max_position, sep=""))


results_chi2
}



