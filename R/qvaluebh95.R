"qvaluebh95" <-
function(p,fdrate=0.1) {
	ps <- sort(p)
	rnkp <- match(p,sort(p))
	ll <- length(p)
	mul <- fdrate*c(1:length(p))/ll
	pass <- (ps <= mul)
	if (length(ps[pass])) pmin <- max(ps[pass]) else pmin <- -1
	out <- list()
	out$significant <- (p <= pmin)
# are these right Q-values???
	x <- rep(0,ll)
	for (i in 1:ll) x[i] <- p[i]*ll/sum(p<=p[i])
	for (i in 1:ll) out$qvalue[i] <- min(x[p>=p[i]])
	out
}

