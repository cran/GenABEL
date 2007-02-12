"catable" <-
function(data,categories = c(quantile(data,c(0.01,0.1,0.5,0.9,0.99))), cumulative = F) {
	if (class(data)!="numeric") stop("data should be numeric vector")
	if (class(categories)!="numeric") stop("categories should be numeric vector")
	ouv <- rep(NA,length(data))

	categories <- sort(categories)
	outmat <- matrix(rep(0,2*(length(categories)+1)),nrow=2)
	tot <- sum(!is.na(data))
	outmat[1,1] <- sum(data<=categories[1])
	outmat[2,1] <- outmat[1,1]/tot
	for (i in 1:(length(categories)-1)) {
		outmat[1,i+1] <- sum(data>categories[i] & data<=categories[i+1])
		outmat[2,i+1] <- outmat[1,i+1]/tot
	}
	outmat[1,length(categories)+1] <- sum(data>categories[length(categories)])
	outmat[2,length(categories)+1] <- outmat[1,length(categories)+1]/tot

	if (cumulative) {
		for (i in 2:(length(categories)+1)) {
			outmat[1,i] <- outmat[1,i] + outmat[1,i-1]
			outmat[2,i] <- outmat[2,i] + outmat[2,i-1]
		}
	}
	cnams <- rep("",length(categories)+1)
	cnams[1] <- paste("X<=",categories[1],sep="")
	for (i in 1:(length(categories)-1)) {
		if (cumulative) 
			cnams[i+1] <- paste("X<=",categories[i+1],sep="")
		else 
			cnams[i+1] <- paste(categories[i],"<X<=",categories[i+1],sep="")
	}
	cnams[length(categories)+1] <- paste("X>",categories[length(categories)],sep="")

	colnames(outmat) <- cnams
	rownames(outmat) <- c("No","Prop")
	outmat
}
