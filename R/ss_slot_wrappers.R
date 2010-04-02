get_slot <- function(data,slotname)
{
    if (class(data) == "gwaa.data") return(slot(data@gtdata,slotname))
    else if (class(data) == "snp.data") return(slot(data,slotname))
    else stop("data should be of class 'snp.data' or 'gwaa.data'")
}

nsnps <- function(data) {return(get_slot(data,"nsnps"))}
nids <- function(data) {return(get_slot(data,"nids"))}
idnames <- function(data) {return(get_slot(data,"idnames"))}
snpnames <- function(data) {return(get_slot(data,"snpnames"))}
map <- function(data) {return(get_slot(data,"map"))}
chromosome <- function(data) {return(as.character(get_slot(data,"chromosome")))}
strand <- function(data) {return(as.character(get_slot(data,"strand")))}
coding <- function(data) {return(as.character(get_slot(data,"coding")))}
male <- function(data) {return(get_slot(data,"male"))}


gtdata <- function(data) {
	if (class(data) == "gwaa.data") return(data@gtdata)
	else if (class(data) == "snp.data") return(data)
	else stop("data should be of class 'snp.data' or 'gwaa.data'")
}
