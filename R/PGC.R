#' Polynomial genomic control
#' 
#' This function estimates the genomic controls
#' for diffrent model and degrees of freedom,
#' using polinom function. Polinomic coefficients are estimated
#' by optimizing diffrent error functions: regress,
#' median and ks.test.
#' 
#' @param data Input vector of Chi square statistic
#' @param method Function of error to be optimized. Can be
#' "regress", "median" or "ks.test"
#' @param p Input vector of allele frequencies
#' @param df Number of dergees of freedom
#' @param pol.d Polinom degree
#' @param plot If true, functon makes plot of lambda from allele frequencies
#' @param start.corr For regress method use it only when you want to make calculations faster
#' @param index.filter Index of variables in data vector, that will be used in analysis, 
#' if zero - all variables will be used
#' 
#' @return A list with elements
#' \item{data}{Output vector corrected Chi square statistic}
#' \item{b}{Polinom coefficients}
#' 
#' @author Yakov Tsepilov
#' 
#' @examples
#' data(ge03d2)
#' qts=qtscore(phdata(ge03d2)$dm2, ge03d2)
#' chi2.1df=results(qts)$chi2.1df
#' s=summary(ge03d2)
#' MAF=s$Q.2
#' result=PGC(data=chi2.1df,method="regress",p=MAF,df=1, pol.d=3, plot=FALSE, start.corr=FALSE)
#'
#' @keywords htest
#' 

PGC = function(data,method="regress",p,df, pol.d=3, plot=TRUE, index.filter=0,start.corr=FALSE){
	if (length(index.filter)<=1){
		ind.function=(1:length(data))
	} else ind.function=index.filter
    pol.m = function(p, b, pol.d){
        out=0;
        for (i in (1:pol.d)){
            out=(p^i)*b[i]+out;
        }     
        out=out+b[pol.d+1];
        return(out);
    }     
    
     if (!(method=="regress" | method=="median" | method=="ks.test")){
          print("Error. I do not know this method");
          break;
     }
     if (start.corr){
          data=data/(1*median(data[ind.function],na.rm=T)/qchisq(.5,df));
    }
     N=sum(!is.na(data[ind.function]))
     #dta <- sort(data)
     #ppp <- ppoints(dta)
     #Chi2=qchisq(ppp,df);
    
    #delta=1/N;
    #x1=0;
    #i=1;
    #Chi2=0;
    #while (i<=N){
    #x2=x1+(delta/2);
    #Chi2[i]=qchisq(x2,df);
    #x1=x1+delta;
    #i=i+1;
    #}
    Chi2=(0:(N-1))/N;
    Chi2=qchisq(Chi2,df)
    
     ppp=seq(from=min(p),to=max(p),length=1000)
    #inf=which(!is.na(data))
    f2 = function(b){
        Zx=data
              
        #Chi2=sort(rchisq(N,df))
        
        #for (i in (1:N)){
        #    j=inf[i];
        #    if (pol.d==3) {
        #        Zx[j]=Zx[j]/(b[1]*p[j]+b[2]*p[j]^2+b[3]*p[j]^3+b[4]);}
        #    if (pol.d==2){
          #            Zx[j]=Zx[j]/(b[1]*p[j]+b[2]*p[j]^2+b[3]);}
        #}
        #if (pol.d==4) {
        #    Zx=Zx/(b[1]*p+b[2]*p^2+b[3]*p^3+b[4]*p^4+b[5]);}
        #if (pol.d==3) {
        #    Zx=Zx/(b[1]*p+b[2]*p^2+b[3]*p^3+b[4]);}
        #if (pol.d==2){
              #Zx=Zx/(b[1]*p+b[2]*p^2+b[3]);}
        #}
        Zx=Zx/(pol.m(p,b,pol.d))     
        Zxl=sort(Zx[ind.function])
        if (method=="ks.test"){
			F=-log(ks.test(Zxl,"pchisq",df=df)$p.value);
        }
        if (method=="median"){
			F=abs(median(Zxl)-qchisq(.5,df))
        }
        if (method=="regress"){
			F=sum((Zxl-Chi2)^2);
        }
        return(F);
    }
    #b=c(-0.1807326,-0.2753122,-0.2110807,1.8401476);
    #if (pol.d==4) {
    #    b=c(0,0,0,0,1)
         #    nlm=nlm(f2,b)
         #    b=nlm$estimate;
     #    if (plot) {plot(sort(b[1]*ppp+b[2]*ppp^2+b[3]*ppp^3+b[4]*ppp^4+b[5]),typ="l");}
      #   data=data/(b[1]*p+b[2]*p^2+b[3]*p^3+b[4]*p^4+b[5]);
    #}
    #if (pol.d==3) {
    #    b=c(0,0,0,1)
         #nlm=nlm(f2,b)
         #b=nlm$estimate;
     #    if (plot) {plot(sort(b[1]*ppp+b[2]*ppp^2+b[3]*ppp^3+b[4]),typ="l");}
      #   data=data/(b[1]*p+b[2]*p^2+b[3]*p^3+b[4]);
    #}
    # if (pol.d==2) {
    #     b=c(0,0,1)
         #nlm=nlm(f2,b)
         #b=nlm$estimate;
    #     if (plot) {plot(sort(b[1]*ppp+b[2]*ppp^2+b[3]),typ="l");}
    #     data=data/(b[1]*p+b[2]*p^2+b[3]);
    #}   

    b=rep(0,pol.d+1);
    b[pol.d+1]=1;
    nlm=nlm(f2,b)
    b=nlm$estimate;
    if (plot) {plot(pol.m(ppp,b,pol.d),typ="l");}
    vv=pol.m(p,b,pol.d);
    data_check=data/vv;
    data_check=data_check[!is.na(data_check)];
    if (sum(data_check<0)==0 & sum(data_check>100)==0){ data=data/vv;}
    out=list()
    out$data=data
    out$b=b
    return(out)
	#data
}

