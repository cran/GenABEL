"estlambda" <-
function(data,plot=TRUE,proportion=1.0) {
	if (min(data)<0) stop("data argument has values <0")
	if (proportion>1.0 || proportion<=0) stop("proportion argument should be greater then zero and less than or equal to one")
	ppoi <- ppoints(data)
	ntp <- round(proportion*length(data))
	if (ntp<10) stop("too few P-values provided")
	if (max(data)<=1) {
		data[data<1e-16] <- 1e-16
		data <- qchisq(1-data,1)
	}
	data <- sort(data)
	ppoi <- sort(qchisq(1-ppoi,1))
	s <- summary(lm(data[1:ntp]~offset(ppoi[1:ntp])))$coeff
	if (plot) {
		lim <- c(0,max(data,ppoi))
		plot(ppoi,data,xlim=lim,ylim=lim,xlab="Expected",ylab="Observed")
		abline(a=0,b=1)
		abline(a=0,b=(1+s[1,1]),col="red")
	}
	out <- list()
	out$estimate <- s[1,1]+1.0
	out$se <- s[1,2]
	out
}
