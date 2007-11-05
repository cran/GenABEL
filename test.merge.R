library(GenABEL)
data(ge03d2)
# merging an object with sub-object 
a <- ge03d2@gtdata[1:10,1:10]
b <- a[c("id107","id33","id25"),c("rs1646456","rs9308393","rs1292700","rs2847446")]
c <- merge(a,b)
xos <- as.numeric(a)
xts <- as.numeric(c)
xos <- xos[sort(rownames(xos)),sort(colnames(xos))]
xts <- xts[sort(rownames(xts)),sort(colnames(xts))]
xos == xts
# merging an sub-object with object 
a <- ge03d2@gtdata[1:10,1:10]
b <- a[c("id107","id33","id25"),c("rs1646456","rs9308393","rs1292700","rs2847446")]
c <- merge(b,a)
xos <- as.numeric(a)
xts <- as.numeric(c)
xos <- xos[sort(rownames(xos)),sort(colnames(xos))]
xts <- xts[sort(rownames(xts)),sort(colnames(xts))]
xos == xts

