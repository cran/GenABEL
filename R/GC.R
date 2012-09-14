#' Genomic control for various model of inheritance using VIF
#' 
#' This function estimates the genomic controls
#' for diffrent model (recessive,dominant,additive etc.),
#' using VIF. VIF coefficients are estimated
#' by optimizing diffrent error functions: regress,
#' median and ks.test.
#'  
#' @param data Input vector of Chi square statistic
#' @param method Function of error to be optimized. Can be
#' "regress", "median" or "ks.test"
#' @param p Input vector of allele frequencies
#' @param x Model of inheritance (0 for recessive,0.5 for additive, 1 for dominant)
#' @param index.filter Indexes for variables that will be use for analisis in data vector
#' @param n size of the sample
#' @param clust For developers only
#' @param vart0 For developers only
#' @param tmp For developers only
#' @param min.p For developers only
#' @param CA For developers only
#' @param p.table For developers only
#' 
#' @return A list with elements
#' \item{Zx}{output vector corrected Chi square statistic}
#' \item{vv}{output vector of VIF}
#' \item{exeps}{output vector of exepsons (NA)}
#' \item{calrate}{output vector of calrate}
#' \item{F}{F}
#' \item{K}{K}
#' 
#' @author Yakov Tsepilov
#' 
#' @examples
#' data(ge03d2)
#' set.seed(1)
#' ge03d2 <- ge03d2[sample(1:nids(ge03d2),200),1:1000]
#' qts=mlreg(phdata(ge03d2)$dm2~1,data=ge03d2,gtmode = "dominant")
#' chi2.1df=results(qts)$chi2.1df
#' s=summary(ge03d2)
#' freq=s$Q.2
#' result=GC(p=freq,x=1,method = "median",CA=FALSE,data=chi2.1df,n=nids(ge03d2))
#' 
#' @keywords htest
#'

GC = function(data=0,p,x,method = "regress",clust=0,vart0=0,tmp=0,min.p=0.05,CA=FALSE,index.filter=0,n,p.table=0){
#
	if (length(index.filter)<=1){
			ind.function=(1:length(data))
		} else ind.function=index.filter
	#
	if (!(method=="regress" | method=="median" | method=="ks.test")){
		print("Error. I do not know this method");
		break;
	}
	#
	if (CA){
		ololo=min.p

		p.table[1,1,]->p00
		p.table[1,2,]->p01
		p.table[1,3,]->p02
		p.table[1,4,]->n0

		p.table[2,1,]->p10
		p.table[2,2,]->p11
		p.table[2,3,]->p12
		p.table[2,4,]->n1

		p.table[3,1,]->p0
		p.table[3,2,]->p1
		p.table[3,3,]->p2
		p.table[3,4,]->n
		
		Zx=0;
		exeps=0;
		Dx=(p12-p02)+x*(p11-p01);
		i=1;
		while (i<=length(n)){
			vart=((p2[i]+p1[i]*x^2)-(p2[i]+x*p1[i])^2)*((1/n1[i])+(1/n0[i]));
			if (vart==0) {Zx[i]=NA; exeps[i]=T; } else{
			Zx[i]=((Dx[i])^2)/(vart); exeps[i]=F}
			i=i+1;
		}

		#Zx=CASforVIF(p.table,p,x)
		inf=which(!is.na(Zx));
		notinf=which(is.na(Zx));
	} else {
		Zx=data
		inf=which(!is.na(Zx));
		notinf=which(is.na(Zx));
	}
	#
	Clust <- function(k){
		#l=which(data1@phdata$dm2==k)
		#data1 <- data1[l,inf]
		#data1 <- Xfix(data1)
		#detach(df@phdata)
		#data1.gkin <- ibs(data1[,sample(which(data1@gtdata@chromosome!="X"),1000)],weight="freq")
		data.dist <- as.dist(0.5-tmp)
		data.mds <- cmdscale(data.dist)
		return(data.mds)
	}
	#
	ClustN <- function(k,data1.mds){
		km <- kmeans(data1.mds,centers=k,nstart=1000)
		i=1;
		  cl=0
		while (i<=k){
			cl[i] <- length(which(km$cluster==i))
			i=i+1;
		}
		return (cl);
	}
	#

	VIF <- function(p,N,F,K){
		q=1-p;
		Var = ((2*p*(1-F)*q*x^2+F*p+(1-F)*p^2)-(2*p*q*(1-F)*x+F*p+(1-F)*p^2)^2);
		S=(1+F)*(1+2*F);
		Cov1=((6*F^3*p+11*(1-F)*F^2*p^2+(1-F)^3*p^4+6*(1-F)^2*F*p^3)/S)-((1-F)*p^2 + F*p)^2; 
		Cov2=x*(((4*(2*(1-F)*F^2*p*q+(1-F)^3*p^3*q+3*(1-F)^2*F*p^2*q))/S)-2*((1-F)*p^2+F*p)*(2*(1 - F)*p*q)); 
		Cov3=x^2*((4*((1-F)^3*p^2*q^2+F*(1-F)*p*q))/S-(2*(1-F)*p*q)^2);
		Cov=Cov1+Cov2+Cov3;
		if (vart0==1){VarT0=N*Var-N*Cov;}
		if (vart0==0){VarT0=N*Var;}
		VarT=VarT0+Cov*K;
		#VarT=VarT0+Cov*(K-N);
		VIF=VarT/VarT0;
		return(VIF)
	}
	#
	GC_VIF_nlm = function(FK){
	###
		F=FK[1];
		K=FK[2];
		#if ((F>1) || (F<0) || (K<0) || (K>(n[1]^2))) {dMedian=(length(Zx)*100 + abs(K));} else{
		vector_vif=VIF(p,n,F,K);
		Zxl=Zx/vector_vif;
		Zxl=sort(Zxl[ind.function]);
		if (method=="ks.test"){
		dMedian=-log(ks.test(Zxl,"pchisq",df=1)$p.value);
		}
		if (method=="median"){
		dMedian=abs(qchisq(.5,1)-median(Zxl));
		}
		if (method=="regress")
		{dMedian=sum((Zxl-Chi2)^2);}
		#}
		return(1*dMedian)
	}
	###
	GC_VIF = function(F,K){
	###
		#F=FK[1];
		#K=FK[2];
		#if ((F>1) || (F<0) || (K<0) || (K>(n[1]^2))) {dMedian=(length(Zx)*100 + abs(K));} else{
		vector_vif=VIF(p,n,F,K);
		Zxl=Zx/vector_vif;
		Zxl=sort(Zxl[ind.function]);
		if (method=="ks.test"){
		dMedian=-log(ks.test(Zxl,"pchisq",df=1)$p.value);
		}
		if (method=="median"){
		dMedian=abs(qchisq(.5,1)-median(Zxl));
		}
		if (method=="regress")
		{dMedian=sum((Zxl-Chi2)^2);}
		#}
		return(1*dMedian)
	}

	###
		#Dx=(p12-p02)+x*(p11-p01);
		F=0.5; K=n[1];
	#
	#i=1;
	#j=0;
	#	Zx=0;
	#	while (i<=length(n)){
	#		vart=((p2[i]+p1[i]*x^2)-(p2[i]+x*p1[i])^2)*((1/n1[i])+(1/n0[i]));
	#		if (vart==0){j=j+1; Zx[i]=NA;} else{Zx[i]=Dx[i]^2/vart;}
	#		i=i+1;
	#	}


	#length(Zx)==length(inf)+length(notinf)
	#N=length(n);
	#if (j!=0){N=N-j;}
	N_inf=sum(!is.na(data[ind.function]))
	Chi2=(0:(N_inf-1))/N_inf;
	Chi2=qchisq(Chi2,1)
			#dta <- sort(Zx)
			#ppp <- ppoints(dta)
			#Chi2=qchisq(ppp,1);

			#delta=1/N;
			#x1=0;
			#i=1;
			#Chi2=0;
			#while (i<=N){
			#	x2=x1+(delta/2);
			#	Chi2[i]=qchisq(x2,1);
			#	x1=x1+delta;
			#	i=i+1;
			#}
	#####
	if (clust==1){
	data.mds0 <- Clust(0)
	#plot(data.mds0)
	k=2;
	cl0 <- ClustN(k,data.mds0)
	data.mds1 <- Clust(1)
	#plot(data.mds1)
	k=2;
	cl1 <- ClustN(k,data.mds1)
	if (length(cl1)>length(cl0)){cl0[(length(cl0)+1):length(cl1)]=0}
	if (length(cl0)>length(cl1)){cl1[(length(cl1)+1):length(cl0)]=0}
	S=sum(cl0);
	R=sum(cl1);
	K=sum((cl1*S-cl0*R)^2)
	K1=K/(R*S)
	#####
	#opt_1=optimize(GC_VIF,c(0,n[1]^2),tol=0.0001,F=0.2);
	opt_2=optimize(GC_VIF,c(0,1),tol=0.0001,K=K1);
	#opt_3=optimize(GC_VIF,c(0,n[1]^2),tol=0.0001,F=opt_2$minimum);
	F=opt_2$minimum;
	K=K1;
		#j=0.01;
		#min=10000;
		#while (j<1){
		#	opt_1=optimize(GC_VIF,c(0,n[1]^2),tol=0.0001,F=j);
		#	opt_2=optimize(GC_VIF,c(0,1),tol=0.0001,K=opt_1$minimum);
		#	if (opt_1$objective<min){	
		#		K=opt_1$minimum;
		#		F=j;
		#		min=opt_1$objective;
		#	}
		#	if (opt_2$objective<min){	
		#		K=opt_1$minimum;
		#		F=opt_2$minimum;
		#		min=opt_2$objective;
		#	}
		#	j=j+0.1;	
		#}	
		#FK=nlm(GC_VIF ,p=c(F=0.5, K=680))$estimate;
		#F=FK[1];
		#K=FK[2];
	}
	if (clust==0){
		#diag(tmp) <- diag(tmp)-0.5;
		#F1=abs(mean(tmp[lower.tri(tmp,diag=T)]));
		#l1=which (tmp>=0 & lower.tri(tmp,diag=T));
		#F1=mean(tmp[l1]);
		if ((median(Zx,na.rm=T)/qchisq(.5,1))>=1){
			lda=median(Zx,na.rm=T)/qchisq(.5,1)+0.05;
			c=(lda-1)*(n[1]/2)
			F1=1/(c)^0.5
			K1=c*(1+F1)/F1
		} else{
			F1=0.5
			K1=n[1]
		}
		#F1=0.1225
		#K1=245.254
		#K1=279.394
		
		#opt_1=optimize(GC_VIF,c(0,n[1]^2),tol=0.0001,F=F1);
		#K1=opt_1$minimum;
		
		####opt_1=optimize(GC_VIF,lower=10,upper=n[1]^2,tol=0.01,F=F1,maximum=T);
			
		FK=nlm(GC_VIF_nlm,c(F=F1,K=K1))$estimate;
		F=FK[1];
		K=FK[2];
	}
	vector_vif=VIF(p,n,F,K);
	Zx=Zx/vector_vif;
	#Zx[notinf]=NA;
	exeps <- is.na(Zx);
	#out <- as.data.frame(Zx)
	out=list()
	out$Zx<-Zx
	out$vv <-vector_vif;
	out$exeps <-exeps
	out$calrate <- n/max(n)
	out$F=F
	out$K=K
	out
	#Zx
	#return(Zx)
}
#ZxVIF=GC(p12,p02,p11,p01,p,x,n1,n0,n);
