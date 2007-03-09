"?" <- function(...) {
	a<-readline("press ENTER to continue")
}

?"LOAD THE DATA"
data(ge03d2)
attach(ge03d2@phdata)

?"DESCRIBE TRAIT DATA"
descriptives.trait(ge03d2)

?"DESCRIBE TRAIT DATA, COMPARE CASES AND CONTROLS"
descriptives.trait(ge03d2,by=dm2)

?"DESCRIBE MARKER DATA IN ALL"
descriptives.marker(ge03d2)
?"DESCRIBE MARKER DATA IN CONTROLS"
descriptives.marker(ge03d2,ids=(dm2==0))
?"DESCRIBE MARKER DATA IN CASES"
descriptives.marker(ge03d2,ids=(dm2==1))

?"GENERATE QQ PLOT FOR HWE P-VALUES IN CONTROLS"
?"GET SUMMARY FOR SNPS IN CONTROLS"
s <- summary(ge03d2@gtdata[(dm2==0),])
?"LOOK UP FIRST 10 SNPs"
s[1:10,]
?"GENERATE QQ PLOT"
estlambda(s[,"Pexact"])

?"LET US FIRST TRY ANALYSIS OF RAW DATA"
an0 <- qtscore("dm2~CRSNP",ge03d2)
?"PLOT THE ANALYSIS RESULTS"
plot(an0)
?"CHECK IF POPULATION IS HOMOGENIOUS (Lambda=1)"
an0$lambda
?"DESCRIBE TOP 10 SNPs"
descriptives.scan(an0)
?"ACCESS EXPERIMENT-WIDE SIGNIFICANCE"
an0.e <- emp.qtscore("dm2~CRSNP",ge03d2)
?"CHECK THE RESULTS"
descriptives.scan(an0.e)
descriptives.scan(an0.e,sort="Pc1df")

?"QUALITY CONTROL, FIRST PASS"
qc1 <- check.marker(ge03d2,p.level=0)
?"DETAILED SUMMARY OF ERRORS"
summary(qc1)
?"DATA SET 1 IS RELATIVELY CLEAN"
data1 <- ge03d2[qc1$idok,qc1$snpok]

?"dDETACH PHENOTYPIC DATA"
detach(ge03d2@phdata)

?"CHECK DATA1 FOR STRONG genetic OUTLIERS; FIRST COMPUTE IBS USING ALL AUTOSOMAL MARKERS"
data1.ibs <- ibs(data1[,data1@gtdata@chromosome!="X"])
?"COMPUTE EUCLIDIAN DISTANCE AS 1-IBS"
data1.dist <- as.dist(1-data1.ibs)
?"PERFORM PRINCIPLE COMPONENTS ANALYSIS OF IBS (FIST TWO COMP-S)"
data1.mds <- cmdscale(data1.dist)
?"PLOT THE RESULTS; YOU WILL SEE TWO DISTINCT CLUSTERS"
plot(data1.mds)
?"IDENTIFY SMALLER CLUSTER"
km <- kmeans(data1.mds,centers=2,nstart=1000)
cl1 <- names(which(km$cluster==1))
cl2 <- names(which(km$cluster==2))
if (length(cl1) > length(cl2)) cl1 <- cl2;
cl1
?"PAINT THE OUTLIERS IN RED"
points(data1.mds[cl1,],pch=19,col="red")
?"EXCLUDE THE OUTLIERS"
data2 <- data1[!(data1@gtdata@idnames %in% cl1),]

?"ATTCH THE DATA"
attach(data2@phdata)
?"CHECK THE DATA AGAIN, NOW CHECK HWE (ONLY IN CONTROLS)"
qc2 <- check.marker(data2,hweids=(dm2==0))
summary(qc2)
data2 <- data2[qc2$idok,qc2$snpok]
?"DATA SET SHOULD BE CLEAN EXCEPT FOR FEW SPORADIC X-ERRORS"
qc3 <- check.marker(data2)
?"FIX SPORADIC X-ERRORS"
data2 <- Xfix(data2)
?"CHECK THAT DATA SET IS INDEED CLEAN"
qc2 <- check.marker(data2)

?"DESCRIBE THE FINAL ANALYSIS DATA"
descriptives.trait(data2,by=dm2)
descriptives.marker(data2)
?"GENERATE QQ PLOT FOR HWE P-VALUES IN CONTROLS"
s <- summary(data2@gtdata[(dm2==0),])
estlambda(s[,"Pexact"])

?"PERFORM GWA ANALYSIS"
data2.qt <- qtscore("dm2~CRSNP",data2)
data2.qt$lambda
plot(data2.qt)
descriptives.scan(data2.qt)
?"ANALYSE EMPIRICAL GWA SIGNIFICANCE"
data2.qte <- emp.qtscore("dm2~CRSNP",data2)
descriptives.scan(data2.qte,sort="Pc1df")

?"TRY FURTHER ADJUSTMENT FOR GENETIC STRATAS (PCA)"
?"FIRST, IDENTIFY MARKERS FOR PCA: HIGH CALL, MAF"
pca.qc <- check.marker(data2,maf=0.1,callrate=0.98,p.level=0,ibs.mrk=0,het.fdr=0.00001)
summary(pca.qc)
?"EXCLUDE THE MARKERS SHOWING ASSOCIATION WITH P <= 0.01"
mrks.pca <- !(pca.qc$snpok %in% data2.qt$snpnames[data2.qt$P1df<=0.01])
?"SELECT RANDOM 20%"
mrks.pca <- sample(mrks.pca,round(length(mrks.pca)*0.2),replace=FALSE)
?"HOW MANY MARKERS ARE THERE?"
length(mrks.pca)
?"COMPUTE FIRST PRINCIPAL COMPONENT OF IBS"
data2.pca <- cmdscale(as.dist(1-ibs(data2[,mrks.pca])),k=1)
?"THE DISTRIBUTION OF PC1"
hist(data2.pca)
?"BIND PC1 TO THE PHEOTYPIC DATA OF DATA2"
data2@phdata$pc1 <- data2.pca[,1]
?"PERFORM GWA ANALYSIS ADJUSTING FOR PC1"
data2.qtp <- qtscore("dm2~pc1+CRSNP",data2)
?"CHECK IF THERE IS STILL ANY INFLATION"
data2.qtp$lambda
descriptives.scan(data2.qtp)
?"ANALYSE EMPIRICAL GWA SIGNIFICANCE"
data2.qtpe <- emp.qtscore("dm2~pc1+CRSNP",data2)
descriptives.scan(data2.qtpe)

?"TRY ADJUSTMENT FOR SEX AND AGE"
data2.qta <- qtscore("dm2~sex+age+pc1+CRSNP",data2)
data2.qta$lambda
plot(data2.qta)
descriptives.scan(data2.qta)
?"EMPIRICAL GW SIGNIFICANCE WITH ADJUSTMENT"
data2.qtae <- emp.qtscore("dm2~sex+age+pc1+CRSNP",data2)
descriptives.scan(data2.qtae)

?"STRATIFIED ANLYSIS WITH OBESE CASES"
data2.qts <- qtscore("dm2~sex+age+pc1+CRSNP",data2,ids=((bmi>30 & dm2==1) | dm2==0))
data2.qts$lambda
plot(data2.qts)
descriptives.scan(data2.qts)
?"EMPIRICAL GW SIGNIFICANCE IN STRATIFIED ANLYSIS"
data2.qtse <- emp.qtscore("dm2~sex+age+pc1+CRSNP",data2,ids=((bmi>30 & dm2==1) | dm2==0))
descriptives.scan(data2.qtse,sort="Pc1df")

?"TRY REPLICATION IN OTHER DATA SET"
?"SELECT TOP 10 SNPs FROM ADJUSTED ANALYSIS"
top10 <- rownames(descriptives.scan(data2.qta))
top10
?"TRY TO REPLICATE IN THE SMALL DATA SET"
data(ge03d2c)
confirm <- qtscore("dm2~sex+age+CRSNP",ge03d2c[,top10])
descriptives.scan(confirm)
?"IS EMPIRICAL EXPERIMENT-WISE P-VALUE ALSO OK?"
confirm.e <- emp.qtscore("dm2~sex+age+CRSNP",ge03d2c[,top10],times=10000,bcast=100)
descriptives.scan(confirm.e)
?"CONFIRMED SNPs"
csnps <- rownames(descriptives.scan(confirm))[c(1,3,4)]
csnps
?"CHECK WHAT ARE THE CONFIRMED SNPs ON NCBI"
show.ncbi(csnps)

?"DO REGIONAL ANALYSIS"
?"SELECT SNPS IN 300 KB REGION AROUND SNP3"
snp <- csnps[1]
snppos <- data2@gtdata@map[which(data2@gtdata@snpnames==snp)]
reg <- snp.names(data2,begin=snppos-150000,end=snppos+150000)
reg
?"ONE SNP ANALYSES"
reg.qt <- qtscore("dm2~sex+age+pc1+CRSNP",data2[,reg],trait.type="binomial")
print(min(reg.qt$P1df))
reg.glm <- scan.glm("dm2~sex+age+pc1+CRSNP",data=data2[,reg],family=binomial())
print(min(reg.glm$P1df))
?"HAPLOTYPE ANALYSIS IN SLIDING WINDOW"
reg.h2 <- scan.haplo("dm2~sex+age+pc1+CRSNP",data2[,reg],trait.type="binomial")
print(min(reg.h2$P1df))
reg.h3 <- scan.haplo("dm2~sex+age+pc1+CRSNP",data2[,reg],n.slide=3,trait.type="binomial")
print(min(reg.h3$P1df))
?"PLOT RESULTS"
minp <- min(reg.qt$P1df,reg.glm$P1df,reg.h2$P1df,reg.h3$P1df)
print(minp)
plot(reg.qt,ylim=c(0,ceiling(-log10(minp))))
add.plot(reg.glm,cex=2)
add.plot(reg.h2,col="green",typ="l")
add.plot(reg.h3,col="blue",typ="l")
?"DO 2D HAPLOTYPE ANALYSIS IN 50 KB REGION OF THE PEAK"
bpos <- reg.h2$map[which.min(reg.h2$P1df)]
sreg <- snp.names(data2,begin=bpos-25000,end=bpos+25000)
sreg
sreg.h2D <- scan.haplo.2D("dm2~age+sex+pc1+CRSNP",data2[,sreg],trait.type="binomial")
print(min(sreg.h2D$P1df,na.rm=T))
plot(sreg.h2D)
?"DO LD ANALYSIS AND UPDATE THE PLOT"
sreg.LD <- LD(as.genotype(data2@gtdata[,sreg]))
image(sreg.h2D$map,sreg.h2D$map,t(sreg.LD$"D'"),add=T)

