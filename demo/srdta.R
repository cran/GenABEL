"?" <- function(...) invisible(0)
#
?"loading the data"
#
?"first, you can convert data from ped-files to internal format"
adpr <- paste(.Library,"/GenABEL/exdata/",sep="")
genofile <- paste(adpr,"pedin.18",sep="")
cat("location of ped-file:",genofile,"\n")
mapfile <- paste(adpr,"map.18",sep="")
cat("location of map-file:",mapfile,"\n")
phfile <- paste(adpr,"phenos.18",sep="")
cat("location of pheno-file:",phfile,"\n")
convert.snp.ped(pedfile=genofile,mapfile=mapfile,outfile="chr18.raw",bcast=50)
a<-readline("PRESS <Enter> TO CONTINUE...")
?"These file can be read in using load.gwaa.data:"
srdta <- load.gwaa.data(geno="chr18.raw",pheno=phfile)
a<-readline("PRESS <Enter> TO CONTINUE...")
?"Now, let us convert and read data from other human-readable format:"
genofile <- paste(adpr,"srgenos.dat",sep="")
cat("location of file with human-readable genotypes:",genofile,"\n")
a<-readline("PRESS <Enter> TO CONTINUE...")
convert.snp.text(infile=genofile,outfile="srgenos.raw",bcast=200)
a<-readline("PRESS <Enter> TO CONTINUE...")
phenofile <- paste(adpr,"srphenos.dat",sep="")
cat("location of file with phenotypes:",phenofile,"\n")
srdta <- load.gwaa.data(geno="srgenos.raw",pheno=phenofile)
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# manipulating phenotypic data
#
names(srdta@phdata)
summary(srdta@phdata)
a<-readline("PRESS <Enter> TO CONTINUE...")
hist(srdta@phdata$age)
a<-readline("PRESS <Enter> TO CONTINUE...")
summary(glm(srdta@phdata$qt2 ~ srdta@phdata$age + srdta@phdata$sex))
a<-readline("PRESS <Enter> TO CONTINUE...")
attach(srdta@phdata)
summary(glm(qt2 ~ age + sex))
a<-readline("PRESS <Enter> TO CONTINUE...")
summary(glm(bt ~ age + sex, family=binomial()))
detach(srdta@phdata)
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# phenotypic QC
#
check.trait("age",srdta)
a<-readline("PRESS <Enter> TO CONTINUE...")
check.trait("qt1",srdta)
a<-readline("PRESS <Enter> TO CONTINUE...")
alltra <- names(srdta@phdata)
alltra
check.trait(alltra,srdta@phdata[srdta@phdata$bt==0,],fdrate=0.2,graph=FALSE)
a<-readline("PRESS <Enter> TO CONTINUE...")
hist(srdta@phdata$qt2)
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta@phdata$qt2 <- replace(srdta@phdata$qt2,(srdta@phdata$qt2==888),NA)
check.trait("qt2",srdta)
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(srdta@phdata$age,srdta@phdata$qt2)
a<-readline("PRESS <Enter> TO CONTINUE...")
summary(glm(srdta@phdata$qt2 ~ srdta@phdata$age + srdta@phdata$sex))
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# Manipulating genetic data
#
srdta@gtdata@nids
srdta@gtdata@nsnps
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta@gtdata@idnames[1:12]
srdta@gtdata@male[1:12]
srdta@gtdata@male[c("p1","p2","p3","p4")]
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta@gtdata@snpnames[1:4]
srdta@gtdata@chromosome[1:4]
srdta@gtdata@map[1:4]
a<-readline("PRESS <Enter> TO CONTINUE...")
n4 <- c("rs18","rs655")
n4
srdta@gtdata@map[n4]
n4 <- c("rs18","rs65")
n4
srdta@gtdata@map[n4]
srdta@gtdata@chromosome[n4]
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta[1:12,1:4]
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta@phdata[1:12,]
srdta@gtdata[1:12,1:4]
a<-readline("PRESS <Enter> TO CONTINUE...")
as.numeric(srdta@gtdata[1:12,1:4])
a<-readline("PRESS <Enter> TO CONTINUE...")
as.numeric(srdta@gtdata[c("p1","p3","p4"),c("rs18","rs65")])
a<-readline("PRESS <Enter> TO CONTINUE...")
as.character(srdta@gtdata[c("p1","p3","p4"),c("rs18","rs65")])
a<-readline("PRESS <Enter> TO CONTINUE...")
as.genotype(srdta@gtdata[c("p1","p3","p4"),c("rs18","rs65")])
a<-readline("PRESS <Enter> TO CONTINUE...")
as.hsgeno(srdta@gtdata[c("p1","p3","p4"),c("rs18","rs65")])
a<-readline("PRESS <Enter> TO CONTINUE...")
srdta@gtdata@snpnames[1:20]
srdta@gtdata@map[1:20]
snp.names(srdta,end=31000)
nams <- snp.names(srdta,end=31000)
nams
srdta@gtdata@map[nams]
srdta@gtdata@chromosome[nams]
a<-readline("PRESS <Enter> TO CONTINUE...")
snp.names(srdta,begin=5000,end=31000)
a<-readline("PRESS <Enter> TO CONTINUE...")
snp.names(srdta,end="rs143")
snp.names(srdta,begin="rs29",end="rs143")
a<-readline("PRESS <Enter> TO CONTINUE...")
snp.names(srdta,begin="rs29",end="rs143")
maprs114 <- srdta@gtdata@map["rs114"]
maprs114
snp.names(srdta,begin=maprs114-12000,end=maprs114+12000)
snp.names(srdta,begin=maprs114-12000,end="rs130")
a<-readline("PRESS <Enter> TO CONTINUE...")
snn <- snp.names(srdta,begin=maprs114-12000,end="rs130")
snn
ids <- srdta@gtdata@idnames[1:4]
ids
as.numeric(srdta@gtdata[ids,snn])
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# genetic data QC
#
summary(srdta@gtdata)
a<-readline("PRESS <Enter> TO CONTINUE...")
summary(srdta@gtdata[,1:10])
a<-readline("PRESS <Enter> TO CONTINUE...")
summary(srdta@gtdata[srdta@phdata$bt==0,1:10])
a<-readline("PRESS <Enter> TO CONTINUE...")
mc <- check.marker(srdta)
summary(mc)
HWE.show(snps=mc$nohwe,data=srdta)
plot(mc)
a<-readline("PRESS <Enter> TO CONTINUE...")
mc <- check.marker(data=srdta,call=0.94,maf=0.01,p.level=0.01,bcast=200)
summary(mc)
plot(mc)
a<-readline("PRESS <Enter> TO CONTINUE...")
HWE.show(snps=mc$nohwe,data=srdta)
a<-readline("PRESS <Enter> TO CONTINUE...")
mc <- check.marker(data=srdta,call=0.94,maf=0.01,p.level=0.01,bcast=200,hweids=(srdta@phdata$bt==0))
summary(mc)
plot(mc)
a<-readline("PRESS <Enter> TO CONTINUE...")
mc <- check.marker(data=srdta,call=0.94,maf=0.01,p.level=0.01,bcast=200,mincon=.90)
summary(mc)
plot(mc)
a<-readline("PRESS <Enter> TO CONTINUE...")
mc1 <- snp.subset(mc,snp.names(srdta,begin=300000,end=500000,chrom="1"))
plot(mc1)
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# running analysis
#
mc <- check.marker(data=srdta,call=0.94,maf=(5/srdta@gtdata@nids),fdr=0.05,hweids=(srdta@phdata$bt==0))
summary(mc)
mymrk <- mc$ok
a<-readline("PRESS <Enter> TO CONTINUE...")
res.fcc <- ccfast("bt",data=srdta,snps=mymrk)
min(res.fcc$P1df)
-log10(min(res.fcc$P1df))
which(res.fcc$P1df<=0.01)
plot(res.fcc)
a<-readline("PRESS <Enter> TO CONTINUE...")
res.boocc <- emp.ccfast("bt",data=srdta,snps=mymrk)
min(res.boocc$P1df)
plot(res.boocc)
a<-readline("PRESS <Enter> TO CONTINUE...")
res.score <- qtscore("bt~CRSNP",data=srdta,snps=mymrk)
min(res.score$P1df)
-log10(min(res.score$P1df))
which(res.score$P1df<=0.01)
plot(res.score)
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(res.fcc)
points(res.score$map,-log10(res.score$P1df),col="red")
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(res.fcc$P1df,res.score$P1df,pch=20)
a<-readline("PRESS <Enter> TO CONTINUE...")
res.score <- qtscore("bt ~ sex + age+CRSNP",data=srdta,snps=mymrk)
min(res.score$P1df)
-log10(min(res.score$P1df))
which(res.score$P1df<=0.01)
plot(res.score)
a<-readline("PRESS <Enter> TO CONTINUE...")
res.boosc <- emp.qtscore("bt~CRSNP",data=srdta,snps=mymrk)
min(res.boosc$P1df)
plot(res.boosc)
a<-readline("PRESS <Enter> TO CONTINUE...")

# analysis of a region

index <- which.min(res.score$P1df)
index
osnp <- res.score$name[index]
osnp
pososnp <- srdta@gtdata@map[osnp]
pososnp
pososnp <- res.score$map[osnp]
pososnp
# wrong -- it will get SNP names from all data, including these which did not pass QC!
reg <- snp.names(res.score,begin=pososnp-50000,end=pososnp+50000) 
reg
r1.fcc <- snp.subset(res.fcc,snps=reg)
min(r1.fcc$P1df)
r1.score <- snp.subset(res.score,snps=reg)
min(r1.score$P1df)
r1.1.glm <- scan.glm("bt ~ sex + age + CRSNP",family=binomial(),data=srdta,snps=reg,bcast=20)
min(r1.1.glm$P1df)
min(r1.1.glm$P2df)
ymax <- -log10(min(r1.1.glm$P1df))+.1
ymax
plot(r1.fcc,ylim=c(0,ymax))
add.plot.scan.gwaa(r1.score,df=1,col="red",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(r1.fcc,ylim=c(0,ymax))
points(r1.score$map,-log10(r1.score$P1df),col="red",type="l")
points(r1.1.glm$map,-log10(r1.1.glm$P1df),col="green",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(r1.fcc,ylim=c(0,ymax))
points(r1.score$map,-log10(r1.score$P1df),col="red",type="l")
add.plot.scan.gwaa(r1.1.glm,df=1,col="green",type="l")
add.plot.scan.gwaa(r1.1.glm,df=2,col="blue",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")

g2d <- scan.glm.2D("bt ~ sex + age + CRSNP",family=binomial(),data=srdta,snps=reg[10:30],bcast=20)
plot(g2d)
a<-readline("PRESS <Enter> TO CONTINUE...")

# haplotypic analysis for bt

library(haplo.stats)
x.adj <- matrix(c(srdta@phdata$sex,srdta@phdata$age),ncol=2)
srdta@phdata[1:5,]
x.adj[1:5,]
a<-readline("PRESS <Enter> TO CONTINUE...")
r1.2.hap <- scan.haplo("bt",data=srdta,snps=reg,trait="binomial",x.adj=x.adj,n.slide=2)
min(r1.2.hap$P1df)
r1.3.hap <- scan.haplo("bt",data=srdta,snps=reg,trait="binomial",x.adj=x.adj,n.slide=3)
min(r1.3.hap$P1df)
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(r1.1.glm)
points(r1.2.hap$map,-log10(r1.2.hap$P1df),col="red",type="l")
points(r1.3.hap$map,-log10(r1.3.hap$P1df),col="green",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")
sr1 <- snp.names(r1.1.glm,beg=2150000,end=2170000,chrom="1")
plot(snp.subset(r1.1.glm,snp=sr1))
points(r1.2.hap$map,-log10(r1.2.hap$P1df),col="red",type="l")
points(r1.3.hap$map,-log10(r1.3.hap$P1df),col="green",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")
sr1 <- snp.names(r1.1.glm,beg=2175000,end=2200000,chrom="1")
plot(snp.subset(r1.1.glm,snp=sr1))
points(r1.2.hap$map,-log10(r1.2.hap$P1df),col="red",type="l")
points(r1.3.hap$map,-log10(r1.3.hap$P1df),col="green",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")

h2d <- scan.haplo.2D("bt",data=srdta,snps=reg[10:30],trait="binomial",x.adj=x.adj)
plot(h2d)
a<-readline("PRESS <Enter> TO CONTINUE...")


# analysis of trait qt2

attach(srdta@phdata)
check.trait("qt2",data=srdta)
summary(glm(qt2~age+sex))
a<-readline("PRESS <Enter> TO CONTINUE...")
r.qts <- qtscore("qt2~CRSNP",data=srdta,snps=mymrk)
min(r.qts$P1df)
plot(r.qts)
a<-readline("PRESS <Enter> TO CONTINUE...")
r.qts <- qtscore("qt2~CRSNP",strata=srdta@phdata$sex,data=srdta,snps=mymrk)
min(r.qts$P1df)
plot(r.qts)
a<-readline("PRESS <Enter> TO CONTINUE...")
r.aqts <- qtscore("qt2~age+sex+CRSNP",data=srdta,snps=mymrk)
min(r.aqts$P1df)
which.min(r.aqts$P1df)
plot(r.aqts)
add.plot.scan.gwaa(r.aqts,col="red",cex=2.)
a<-readline("PRESS <Enter> TO CONTINUE...")
r.booqts <- emp.qtscore("qt2~age+sex+CRSNP",data=srdta,snps=mymrk,times=200)
min(r.booqts$P1df)
which.min(r.booqts$P1df)
plot(r.booqts)
a<-readline("PRESS <Enter> TO CONTINUE...")

#investigation of the region
osnp <- r.booqts$name[which.min(r.booqts$P1df)]
osnp
pososnp <- srdta@gtdata@map[osnp]
pososnp
reg <- snp.names(r.booqts,beg=pososnp-20000,end=pososnp+20000,chr="1")
reg
a<-readline("PRESS <Enter> TO CONTINUE...")
rr.qts <- snp.subset(r.aqts,snps=reg)
min(rr.qts$P1df)
rr.glm <- scan.glm("qt2~age+sex+CRSNP",data=srdta,snps=reg)
min(rr.glm$P1df)
plot(rr.glm)
points(rr.qts$map,-log10(rr.qts$P1df),col="red")
a<-readline("PRESS <Enter> TO CONTINUE...")
r.2.hap <- scan.haplo("qt2",data=srdta,snps=reg,x.adj=x.adj,n.slide=2)
min(r.2.hap$P1df)
r.3.hap <- scan.haplo("qt2",data=srdta,snps=reg,x.adj=x.adj,n.slide=3)
min(r.3.hap$P1df)
a<-readline("PRESS <Enter> TO CONTINUE...")
plot(rr.glm)
add.plot.scan.gwaa(rr.qts,col="cyan")
add.plot.scan.gwaa(r.2.hap,col="red",type="l")
add.plot.scan.gwaa(r.3.hap,col="green",type="l")
a<-readline("PRESS <Enter> TO CONTINUE...")


#
# Using gwaa data in other packages -- haplo.stats
#
osnp
index <- which(srdta@gtdata@snpnames==osnp)
index
r5 <- srdta@gtdata@snpnames[(index-2):(index+2)]
r5
haplo.score(qt2,as.hsgeno(srdta@gtdata[,r5]),x.adj=x.adj)
a<-readline("PRESS <Enter> TO CONTINUE...")

#
# Using gwaa data in other packages -- genetics
#
gdta <- as.genotype(srdta@gtdata[,reg])
library(genetics)
ld <- LD(gdta)
plot(ld)
a<-readline("PRESS <Enter> TO CONTINUE...")

