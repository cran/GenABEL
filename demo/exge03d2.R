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

?"TRY ADJUSTMENT FOR SEX AND AGE"
data2.qta <- qtscore("dm2~sex+age+CRSNP",data2)
data2.qta$lambda
plot(data2.qta)
descriptives.scan(data2.qta)
?"EMPIRICAL GW SIGNIFICANCE WITH ADJUSTMENT"
data2.qtae <- emp.qtscore("dm2~sex+age+CRSNP",data2)
descriptives.scan(data2.qtae)

?"CHECK WHAT HAPPENS IF ADJUSTMENT IS MADE FOR BMI"
data2.qtabmi <- qtscore("dm2~sex+age+bmi+CRSNP",data2)
data2.qtabmi$lambda
plot(data2.qtabmi)
descriptives.scan(data2.qtabmi)
?"EMPIRICAL GW SIGNIFICANCE WITH ADJUSTMENT"
data2.qtabmie <- emp.qtscore("dm2~sex+age+bmi+CRSNP",data2)
descriptives.scan(data2.qtabmie)

?"STRATIFIED ANLYSIS WITH OBESE CASES"
data2.qts <- qtscore("dm2~sex+age+CRSNP",data2,ids=((bmi>30 & dm2==1) | dm2==0))
data2.qts$lambda
plot(data2.qts)
descriptives.scan(data2.qts)
?"EMPIRICAL GW SIGNIFICANCE IN STRATIFIED ANLYSIS"
data2.qtse <- emp.qtscore("dm2~sex+age+CRSNP",data2,ids=((bmi>30 & dm2==1) | dm2==0))
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
csnps <- rownames(descriptives.scan(confirm))[1:3]
csnps
?"CHECK WHAT ARE THE CONFIRMED SNPs ON NCBI"
show.ncbi(csnps)

