#' Export GenABEL data in PLINK format
#' 
#' Export GenABEL data in PLINK format. This function is 
#' a simple wrapper to \code{\link{export.merlin}} function 
#' with specific arguments + few lines of code to 
#' export phenotypes
#' 
#' @param data GenABEL data object of 'gwaa.data'-class to 
#' be exported
#' 
#' @param filebasename base file name for exported data, 
#' extensions '.ped', '.map' and '.phe' (for phenotype file) 
#' are added for specific output files
#' 
#' @param phenotypes NULL (no phenotypes exported), "all" for 
#' all phenotypes or a vector of character with names of phneotypes 
#' to be exported 
#' 
#' @param ... arguments passed to \code{\link{export.merlin}}
#' 
#' @author Yurii Aulchenko
#' 
#' @keywords IO
#' 

"export.plink" <- function(data, filebasename, phenotypes= "all", ...) {
	
	if (!is.null(phenotypes)) {
		phef <- paste(filebasename,".phe",sep="")
		phed <- phdata(data)
		phed <- data.frame(FID=seq(1:dim(phed)[1]),IID=phed[,"id"],
				phed[,which(names(phed) != "id")])
		if (phenotypes != "all") {
			phed <- phed[,c("FID","IID",phenotypes)]
		}
		write.table(phed,file=phef,row.names=FALSE,col.names=TRUE,quote=FALSE,sep=" ")
	}
	
	pedf <- paste(filebasename,".ped",sep="")
	mapf <- paste(filebasename,".map",sep="")
	
	export.merlin(data,pedfile=pedf,datafile=NULL,
			mapfile=mapf,format="plink", extendedmap=FALSE, ... )
	
}