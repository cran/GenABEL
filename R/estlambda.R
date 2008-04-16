"estlambda" <-
function(data,plot=TRUE,proportion=1.0, ...) {
	data <- data[which(!is.na(data))]
	if (proportion>1.0 || proportion<=0) stop("proportion argument should be greater then zero and less than or equal to one")
	ntp <- round(proportion*length(data))
	if (ntp<1) stop("no valid measurments")
	if (ntp==1) {
		warning(paste("One measurment, Lambda = 1 returned"))
		return(list(estimate=1.0,se=999.99))
	}
	if (ntp<10) warning(paste("number of points is too small:",ntp))
	if (min(data)<0) stop("data argument has values <0")
	if (max(data)<=1) {
		lt16 <- (data < 1.e-16)
		if (any(lt16)) {
			warning(paste("Some probabilities < 1.e-16; set to 1.e-16"))
			data[lt16] <- 1.e-16
		}
		data <- qchisq(1-data,1)
	}
	data <- sort(data)
	ppoi <- ppoints(data)
	ppoi <- sort(qchisq(1-ppoi,1))
	data <- data[1:ntp]
	ppoi <- ppoi[1:ntp]
	s <- summary(lm(data~offset(ppoi)))$coeff
	if (plot) {
		lim <- c(0,max(data,ppoi,na.rm=T))
#		plot(ppoi,data,xlim=lim,ylim=lim,xlab="Expected",ylab="Observed", ...)
		plot(ppoi,data,xlab="Expected",ylab="Observed", ...)
		abline(a=0,b=1)
		abline(a=0,b=(1+s[1,1]),col="red")
	}
	out <- list()
	out$estimate <- s[1,1]+1.0
	out$se <- s[1,2]
	out
}
