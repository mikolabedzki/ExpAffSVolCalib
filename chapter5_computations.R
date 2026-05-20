datafolder = "C:/Users/mikolaj.labedzki/Downloads/marketdata"
scriptsfolder = "C:/Users/mikolaj.labedzki/Documents/Rscripts"
workfolder = "C:/Users/mikolaj.labedzki/Documents/dissertation_R"
#############################
## Intro                   ##
#############################
setwd(workfolder)
options("scipen"=100, "digits"=7)#to force switch of scientific notation
#load("ICMfinalExpo.RData")
#save.image("ICMfinal.RData")
library(Rcpp); Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(paste0(scriptsfolder,"/nl_expaffsvol.cpp"))
source(paste0(scriptsfolder,"/aff_pricing.R"))
#library(TTR);library(tseries);
library(xtable)#latex export

tau=c(1/12,2/12,3/12,6/12,1,2)

datatable=read.table(paste0(datafolder,"/FX/EURUSDall.txt"), header = TRUE, sep="\t");pair="EURUSD"
datatable$Date = as.Date(datatable$Date)
n=dim(datatable)[1]
tau1 = sweep(matrix(1,n,6),MARGIN=2,tau,`*`)
P1=matrix(NA,n,6);VIX=P1;P2=P1;P3=P1;P4=P1;P5=P1;P6=P1;SS=P1;CC=P1;mu=P1;Q1=P1;Q2=P1;Q3=P1;Q4=P1;l1=P1;l2=P1;skew=P1;kurt=P1;level=P1;S_Leif=P1[,1];C_Leif=P1[,1]
data.load(datatable,pair)

S = as.numeric(datatable[,32])
vol = datatable[,2:31]/100
rates = get_rates(datatable,pair)
rd = rates[,1:6]
rf = rates[,7:12]
fwd = rates[,13:18]
F = S+fwd
options = get_options(vol,F,rd,rf,tau)
K = options[,1:30]
mktPrice = options[,31:60]
vega = options[,61:90]
remove(rates);remove(options)

#S = om*rho/4
#C = om^2*(2-5*rho^2)/24/vol
#omrho = 4*slope
#om2 = 12*curv*vol + 5*8*slope^2
#omegaA = sqrt(12*C_Leif*sqrt(nu0) + 5*8*S_Leif^2)
#rhoA = 4*S_Leif/omegaA

#############################
## Tables Printing         ##
#############################
tab0 = rbind(summarypy(hestonTSpar3m[,1]),summarypy(hestonTSpar3m[,2]),summarypy(hestonTSpar3m[,3]),summarypy(sqrt(hestonTSpar3m[,4]/5)))
rownames(tab0)=c("$\\nu_0$","$\\theta$","$\\kappa$","MAE")
cfa0.table <- xtable(tab0)
digits(cfa0.table) <- 4
write.textable(cfa0.table,"r5.t1.txt")  

tab1=rbind(summarypy(durr[[1]]),summarypy(hist3m[[1]])[-7],summarypy(lab[[1]]),summarypy(perfectCal3alt[,6]),summarypy(perfectCal3falt[,6]),summarypy(perfectCal3szalt[,6]))
rownames(tab1)=c("Durrleman-Heston","Hist.-Heston","ICM-Heston","Calib.-Heston", "Calib.-Heston-Feller", "Calib.-Schöbel-Zhu")
cfa1.table <- xtable(tab1)
digits(cfa1.table) <- 4
write.textable(cfa1.table,"r5.t2.txt")

dat1=as.data.frame(cbind(perfectCal3szalt[,6],lab[[1]]))
dat1=as.data.frame(cbind(perfectCalOUOU2[,11],perfectCalBates2[,11]))
dat1=as.data.frame(cbind(perfectCalOUOU1[,11],perfectCalBates1[,11]))
dat1=as.data.frame(cbind(perfectCalBates1f[,11],perfectCalOUOU1[,11]))
dat1=as.data.frame(cbind(perfectCalBates2f[,11],perfectCalOUOU2[,11]))
dat1=as.data.frame(cbind(perfectCalBates1f[,11],perfectCalBates2f[,11]))
dat1=as.data.frame(cbind(lab2df1[[1]],lab2do1[[1]]))
dat1=as.data.frame(cbind(par5d,par5l))
dat1=as.data.frame(cbind(errMSEd,errMSEl))
dat1=as.data.frame(cbind(rhoLmed,rhoP))
dat1=as.data.frame(cbind(omegaLmed,homega3m))
ndat1=reshape(dat1, direction="long", varying=list(names(dat1)), v.names="Value")
oneway.test(Value~time, ndat1, var.equal = FALSE)
t.test(Value~time, data = ndat1, paired = TRUE)

test1=cbind(omegaDA,homega3m,omegaLmed,omegaP,rhoDA,hrho3m,rhoLmed,rhoP)
test1=test1[-(1:62),]
#dtest1=diff(test1)
ctab1=cor(test1)
rownames(ctab1)=c("$\\omega$ Durr.","$\\omega$ Hist.","$\\omega$ ICM","$\\omega$ Calib.","$\\rho$ Durr.","$\\rho$ Hist.","$\\rho$ ICM","$\\rho$ Calib.")
colnames(ctab1)=rownames(ctab1)
sumtab=rbind(summarypy(omegaDA),summarypy(homega3m),summarypy(omegaLmed),summarypy(omegaP),summarypy(rhoDA),summarypy(hrho3m),summarypy(rhoLmed),summarypy(rhoP))
rownames(sumtab)=rownames(ctab1)

cfa2.table <- xtable(ctab1)
cfa3.table <- xtable(sumtab)
digits(cfa2.table) <- 3
digits(cfa3.table) <- 3
write.textable(cfa3.table,"r5.t3.txt")
write.textable(cfa2.table,"r5.t4.txt")

tab4=rbind(summarypy(par1d),summarypy(par1l),summarypy(par4d),summarypy(par4l),summarypy(par5d),summarypy(par5l),summarypy(errMSEd),summarypy(errMSEl))
rownames(tab4)=c("$\\theta$ Durr.","$\\theta$ ICM","$\\kappa$ Durr.","$\\kappa$ ICM","$\\nu_0$ Durr.","$\\nu_0$ ICM","RMSE Durr.","RMSE ICM")
cfa4.table <- xtable(tab4)
digits(cfa4.table) <- 4
write.textable(cfa4.table,"r5.t5.txt")

tab6=rbind(summarypy(lab2db1[[1]]),summarypy(lab2db2[[1]]),summarypy(lab2do1[[1]]),summarypy(lab2do2[[1]]),summarypy(lab2df1[[1]]),summarypy(lab2df2[[1]]))
modelsmethods=c("Bates/2-stage","Bates/ICM/EVP","OUOU/2-stage","OUOU/ICM/MEVP","BatesFeller/2-stage","BatesFeller/ICM/MEVP")
rownames(tab6)=modelsmethods
cfa6.table <- xtable(tab6)
digits(cfa6.table) <- 4
write.textable(cfa6.table,"r5.t6.txt")

tab7=rbind(summarypy(perfectCalBates1[,11]),summarypy(perfectCalBates2[,11]),summarypy(perfectCalOUOU1[,11]),summarypy(perfectCalOUOU2[,11]),summarypy(perfectCalBates1f[,11]),summarypy(perfectCalBates2f[,11]))
rownames(tab7)=modelsmethods
cfa7.table <- xtable(tab7)
digits(cfa7.table) <- 5
write.textable(cfa7.table,"r5.t7.txt")

round(c(41,87,60,119,212,207,337)*252/c(30.84,31.32,33.37,43.77,57.79,65.91,76.93)/1000*100,1)

onepage:
i, iii, vii-xii, 1-2, 3, 4, 6, 8, 23, 118, 12, 16, 28, 30, 37, 89, 105-106, 121-122, 125-126, 107, 109

1-2 heston 5 parameters
JEL: C13 Estimation: General, C52 Model Evaluation, Validation, and Selection, G13 Contingent Pricing,

#############################
## Plot Printing           ##
#############################

myaxis = seq(as.Date("2010-07-22"), as.Date("2015-08-31"), by = 91)
datatable$Date = as.Date(datatable$Date)
mywidth = 15
myheight = 10

pdf("../images/rhos.pdf", width=mywidth, height=myheight)
plot(rhoDA~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", type = "b", lty=1, pch=1, ylab="rho", xlab="Date", ylim=c(-0.75,0.2), cex=0.5, cex.lab=1.4, cex.sub=1.4)
points(datatable$Date,hrho3m,col=2, pch = 2, cex=0.5);lines(datatable$Date,hrho3m, col=2, lty = 2, cex=0.5);
points(datatable$Date,rhoLmed,col=3, pch = 3, cex=0.5);lines(datatable$Date,rhoLmed, col=3, lty = 3, cex=0.5)
points(datatable$Date,rhoP, col=4, pch = 4, cex=0.5);lines(datatable$Date,rhoP, col=4, lty = 4, cex=0.5)
abline(h=0,lty=2)
legend("topleft", c("Durrleman","Hist. 3m","ICM","Calib."), lwd=1.4, lty=1:4, pch=c(1,2,3,4), col=c(1,2,3,4), bty="n", cex=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

pdf("../images/omegas.pdf", width=mywidth, height=myheight)
plot(omegaDA~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", type = "b", lty=1, pch=1, col=1, ylab="omega", ylim=c(0.0,0.85), xlab="Date", cex=0.5, cex.lab=1.4, cex.sub=1.4)
points(datatable$Date,homega3m,col=2, pch = 2, cex=0.5);lines(datatable$Date,homega3m,col=2, lty = 2, cex=0.5)
points(datatable$Date,omegaLmed, col=3, pch = 3, cex=0.5);lines(datatable$Date,omegaLmed, col=3, lty = 3, cex=0.5)
points(datatable$Date,omegaP, col=4, pch = 4, cex=0.5);lines(datatable$Date,omegaP, col=4, lty = 4, cex=0.5)
legend("topleft", c("Durrleman","Hist. 3m","ICM","Calib."), lwd=1.4, lty=1:4, pch=c(1,2,3,4), col=c(1,2,3,4), bty="n", cex=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

pdf("../images/durrvslab.pdf", width=mywidth, height=myheight)
plot(durr[[1]]~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", pch=1, lty=1, t="b", ylab=c("RMSE"), xlab="Date", ylim=c(0.0005,0.028), cex=0.5, cex.lab=1.4, cex.sub=1.4)
points(datatable$Date,hist3m[[1]], col=2, pch=2, cex=0.5);lines(datatable$Date,hist3m[[1]], col=2, lty=2, cex=0.5)
points(datatable$Date,lab[[1]], col=3, pch=3, cex=0.5);lines(datatable$Date,lab[[1]], col=3, lty=3, cex=0.5)
points(datatable$Date,perfectCal3alt[,6], col=4, pch=4, cex=0.5);lines(datatable$Date,perfectCal3alt[,6], col=4, lty=4, cex=0.5)
points(datatable$Date,perfectCal3falt[,6], col=5, pch=5, cex=0.5);lines(datatable$Date,perfectCal3falt[,6], col=5, lty=5, cex=0.5)
points(datatable$Date,perfectCal3szalt[,6], col=6, pch=6, cex=0.5);lines(datatable$Date,perfectCal3szalt[,6], col=6, lty=6, cex=0.5)
legend("topleft", c("Durrleman-Heston","Hist.-Heston","ICM-Heston","Calib.-Heston", "Calib.-Heston-Feller", "Calib.-Schöbel-Zhu"), lwd=1.4, pch=1:6, lty=1:6, col=c(1,2,3,4,5,6), bty="n", cex=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

library(latticeExtra)
library(reshape)
wide = as.data.frame(cbind(c(" 1m"," 2m"," 3m"," 6m","12m","24m"),durr[[2]]-lab[[2]]))
colnames(wide) = c("tenor","-10p","-25p","50c","25c","10c")
wide.out = melt(wide, measure.vars = c("-10p","-25p","50c","25c","10c"), variable.name = "delta", value.name = "error")
wide.out$value=as.numeric(as.character(wide.out$value))
wide.out2=wide.out
wide.out2$value=-wide.out$value

mat=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4)
mat=matrix(c(1,0.2,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4)
pdf("../images/better1.pdf", width=14, height=7)
print(cloud(value~tenor+variable, wide.out, main="max(RMSE(Durr)-RMSE(ICM),0)", panel.3d.cloud=panel.3dbars, distance = 0.3, R.mat = mat, xbase=0.5, ybase=0.5, col.facet='green', scales=list(arrows=FALSE, col=1), par.settings = list(axis.line = list(col = "transparent")), zlim=c(0,0.006), aspect = c(1,0.3,1), par.box = list(col="grey"), drape=FALSE, shade = FALSE, ylab=c("delta"), zlab=c("RMSE")), split = c(1,1,2,1), more = TRUE)
print(cloud(value~tenor+variable, wide.out2, main="max(RMSE(ICM)-RMSE(Durr),0)", panel.3d.cloud=panel.3dbars, distance = 0.3, R.mat = mat, xbase=0.5, ybase=0.5, col.facet='red', scales=list(arrows=FALSE, col=1), par.settings = list(axis.line = list(col = "transparent")), zlim=c(0,0.006), aspect = c(1,0.3,1), par.box = list(col="grey"), drape=FALSE, shade = FALSE, ylab=c("delta"), zlab=c("RMSE")), split = c(2,1,2,1))
dev.off()

pdf("../images/betterA.pdf", width=7, height=7)
print(cloud(value~tenor+variable, wide.out, panel.3d.cloud=panel.3dbars, distance = 0.3, R.mat = mat, xbase=0.5, ybase=0.5, col.facet='green', scales=list(arrows=FALSE, col=1), par.settings = list(axis.line = list(col = "transparent")), zlim=c(0,0.006), aspect = c(1,0.3,1), par.box = list(col="grey"), drape=FALSE, shade = FALSE, ylab=c("delta"), zlab=c("RMSE")), split = c(1,1,1,1), more = TRUE)
dev.off()

pdf("../images/betterB.pdf", width=7, height=7)
print(cloud(value~tenor+variable, wide.out2, panel.3d.cloud=panel.3dbars, distance = 0.3, R.mat = mat, xbase=0.5, ybase=0.5, col.facet='red', scales=list(arrows=FALSE, col=1), par.settings = list(axis.line = list(col = "transparent")), zlim=c(0,0.006), aspect = c(1,0.3,1), par.box = list(col="grey"), drape=FALSE, shade = FALSE, ylab=c("delta"), zlab=c("RMSE")), split = c(1,1,1,1))
dev.off()

#pdf("omega.pdf", width=12, height=6)
#plot((par2d), t="b", pch=1, col=1, lty=1, ylab="Risk", xlab="Date", cex=0.5, cex.lab=1.4, cex.sub=1.4)
#points(par2l, col=2, pch = 2, cex=0.5);lines(par2l, col=2, lty = 2, cex=0.5)
#legend("topleft", c("Omega Durr.","Omega ICM"), lty=c(1,2), pch=c(1,2), col=c(1,2), bty="n", cex=1.4, lwd=1.4)
#grid()
#dev.off()

#pdf("rho.pdf", width=12, height=6)
#plot((par3d), t="b", pch=1, col=1, lty=1, ylab="Risk", xlab="Date", cex=0.5, cex.lab=1.4, cex.sub=1.4)
#points(par3l, col=2, pch = 2, cex=0.5);lines(par3l, col=2, lty = 2, cex=0.5)
#legend("topleft", c("Rho Durr.","Rho ICM"), lty=c(1,2), pch=c(1,2), col=c(1,2), bty="n", cex=1.4, lwd=1.4)
#grid()
#dev.off()

pdf("../images/theta.pdf", width=mywidth, height=myheight-2.98)
plot(par1d~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", t="b", pch=1, col=1, lty=1, ylab="Risk", xlab="Date", cex=0.5, cex.lab=1.4, cex.sub=1.4)
points(datatable$Date,par1l, col=2, pch = 2, cex=0.5);lines(datatable$Date,par1l, col=2, lty = 2, cex=0.5)
legend("topleft", c("Theta Durr.","Theta ICM"), lty=c(1,2), pch=c(1,2), col=c(1,2), bty="n", cex=1.4, lwd=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

pdf("../images/kappa.pdf", width=mywidth, height=myheight-2.98)
plot(par4l~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", t="b", pch=2, col=2, lty=2, ylab="Risk", xlab="Date", cex=0.5, cex.lab=1.4, cex.sub=1.4)
points(datatable$Date,par4d, col=1, pch = 1, cex=0.5);lines(datatable$Date,par4d, col=1, lty = 1, cex=0.5)
legend("topleft", c("Kappa Durr.","Kappa ICM"), lty=c(1,2), pch=c(1,2), col=c(1,2), bty="n", cex=1.4, lwd=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

pdf("../images/nu0.pdf", width=mywidth, height=myheight-2.98)
plot(par5d~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", t="b", pch=1, col=1, lty=1, ylab="Risk", xlab="Date", cex=0.5, cex.lab=1.4, cex.sub=1.4)
points(datatable$Date,par5l, col=2, pch = 2, cex=0.5);lines(datatable$Date,par5l, col=2, lty = 2, cex=0.5)
legend("topleft", c("Nu0 Durr.","Nu0 ICM"), lty=c(1,2), pch=c(1,2), col=c(1,2), bty="n", cex=1.4, lwd=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

pdf("../images/rmse.pdf", width=mywidth, height=myheight-2.98)
plot(errMSEd~datatable$Date, xlim=c(as.Date("2010-08-30"),as.Date("2015-07-10")), xaxt="n", t="b", pch=1, col=1, lty=1, ylab="RMSE", xlab="Date", cex=0.5, cex.lab=1.4, cex.sub=1.4)
points(datatable$Date,errMSEl, col=2, pch = 2, cex=0.5);lines(datatable$Date,errMSEl, col=2, lty = 2, cex=0.5)
legend("topleft", c("RMSE Durr.","RMSE ICM"), lty=c(1,2), pch=c(1,2), col=c(1,2), bty="n", cex=1.4, lwd=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

#pdf("Feller.pdf", width=12, height=6)
#plot((cumFeller),ylab="Cumulated share of breaks of Feller condition ", xlab="Date", lty=0, cex=0, cex.lab=1.4, cex.sub=1.4)
#lines(cumFeller)
#lines(diffinv(as.numeric((test.fellerL<0))),col=2,lty=2)
#legend("topleft", c("Durrleman","ICM"), lty=1:2, col=c('black','red'), bty="n", cex=1.2)
#dev.off()

# Part 1 : "OUOU / Equal Var.","OUOU / Fitted Ratio"

pdf("../images/rmse2d_start.pdf", width=mywidth, height=myheight)
plot(lab2df2[[1]]~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", t="b", pch=6, lty=6, ylab=c("RMSE"), xlab="Date", ylim=c(0.001,0.031), cex=0.5, cex.lab=1.4, cex.sub=1.4, col=6)
points(datatable$Date,lab2db2[[1]], col=2, pch=2, cex=0.5);lines(datatable$Date,lab2db2[[1]], col=2, lty=2, cex=0.5)
points(datatable$Date,lab2do1[[1]], col=3, pch=3, cex=0.5);lines(datatable$Date,lab2do1[[1]], col=3, lty=3, cex=0.5)
points(datatable$Date,lab2do2[[1]], col=4, pch=4, cex=0.5);lines(datatable$Date,lab2do2[[1]], col=4, lty=4, cex=0.5)
points(datatable$Date,lab2df1[[1]], col=5, pch=5, cex=0.5);lines(datatable$Date,lab2df1[[1]], col=5, lty=5, cex=0.5)
points(datatable$Date,lab2db1[[1]], col=1, pch=1, cex=0.5);lines(datatable$Date,lab2db1[[1]], col=1, lty=1, cex=0.5)
legend("topright", modelsmethods, pch=(1:6), lty=(1:6), col=1:6, bty="n", cex=1.4, lwd=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

# Part 2: Comparison of better starting point and full calibration vs 2phase calibration simple starting point vs Bates old full calib

#x=5
#x1=(1+x*rnorm(n)/100)
#x2=(1+x*rnorm(n)/100)
#x3=(1+x*rnorm(n)/100)
#x4=(1+x*rnorm(n)/100)
#x5=(1+x*rnorm(n)/100)
#x6=(1+x*rnorm(n)/100)

pdf("../images/rmse2d_end.pdf", width=mywidth, height=myheight)
plot(perfectCalBates1f[,11]~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", t="b", pch=5, lty=5, ylab=c("RMSE"), xlab="Date", cex=0.5, cex.lab=1.4, cex.sub=1.4, col=5, ylim=c(0.0003,0.0142))
points(datatable$Date,perfectCalBates2[,11], col=2, pch=2, cex=0.5);lines(datatable$Date,perfectCalBates2[,11], col=2, lty=2, cex=0.5)
points(datatable$Date,perfectCalOUOU1[,11], col=3, pch=3, cex=0.5);lines(datatable$Date,perfectCalOUOU1[,11], col=3, lty=3, cex=0.5)
points(datatable$Date,perfectCalOUOU2[,11], col=4, pch=4, cex=0.5);lines(datatable$Date,perfectCalOUOU2[,11], col=4, lty=4, cex=0.5)
points(datatable$Date,perfectCalBates1[,11], col=1, pch=1, cex=0.5);lines(datatable$Date,perfectCalBates1[,11], col=1, lty=1, cex=0.5)
points(datatable$Date,perfectCalBates2f[,11], col=6, pch=6, cex=0.5);lines(datatable$Date,perfectCalBates2f[,11], col=6, lty=6, cex=0.5)
legend("topright", modelsmethods, pch=(1:6), lty=(1:6), col=1:6, bty="n", cex=1.4, lwd=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

convergence1=read.table("proper/nl_out1b.txt",sep=";",fill=TRUE)
convergence2=read.table("proper/nl_out2b.txt",sep=";",fill=TRUE)
convergence1o=read.table("proper/nl_out1o.txt",sep=";",fill=TRUE)
convergence2o=read.table("proper/nl_out2o.txt",sep=";",fill=TRUE)
convergence1f=read.table("proper/nl_out1f.txt",sep=";",fill=TRUE)
convergence2f=read.table("proper/nl_out2f.txt",sep=";",fill=TRUE)

convergence1o=t(matrix(rep(as.numeric(unlist(t(convergence1o))),6),161,n+5))[-((n+1):(n+5)),]
convergence2o=t(matrix(rep(as.numeric(unlist(t(convergence2o))),6),161,n+5))[-((n+1):(n+5)),]
convergence1f=t(matrix(rep(as.numeric(unlist(t(convergence1f))),6),161,n+5))[-((n+1):(n+5)),]
convergence2f=t(matrix(rep(as.numeric(unlist(t(convergence2f))),6),161,n+5))[-((n+1):(n+5)),]
convergence1=t(matrix(rep(as.numeric(unlist(t(convergence1))),6),161,n+5))[-((n+1):(n+5)),]
convergence2=t(matrix(rep(as.numeric(unlist(t(convergence2))),6),161,n+5))[-((n+1):(n+5)),]

#fix(convergence1)
convergence1[,161] = as.numeric(rep(NA,n))
convergence2[,161] = as.numeric(rep(NA,n))
convergence1o[,161] = as.numeric(rep(NA,n))
convergence2o[,161] = as.numeric(rep(NA,n))
convergence1f[,161] = as.numeric(rep(NA,n))
convergence2f[,161] = as.numeric(rep(NA,n))
convergence1=addatend(convergence1,perfectCalBates1[,11])
convergence2=addatend(convergence2,perfectCalBates2[,11])
convergence1o=addatend(convergence1o,perfectCalOUOU1[,11])
convergence2o=addatend(convergence2o,perfectCalOUOU2[,11])
convergence1f=addatend(convergence1f,perfectCalBates1f[,11])
convergence2f=addatend(convergence2f,perfectCalBates2f[,11])
convergence1=cbind(lab2db1[[1]],convergence1)
convergence2=cbind(lab2db2[[1]],convergence2)
convergence1o=cbind(lab2do1[[1]],convergence1o)
convergence2o=cbind(lab2do2[[1]],convergence2o)
convergence1f=cbind(lab2df1[[1]],convergence1f)
convergence2f=cbind(lab2df2[[1]],convergence2f)

pdf("../images/rmse2d_conv_log.pdf", width=mywidth, height=myheight)
x=(0:161)*10
x1=1:162
plot(x,log(colMeans(convergence2f[,x1])), t="b", lty=6, pch=6,col=6, ylab=c("Natural logarithm of mean RMSE"), xlab="Number of iterations", cex=0.5, cex.lab=1.4, cex.sub=1.4, ylim=c(-6.7,-4.5)) # ylim=c(0.0015,0.0111))ylim=c(-6.7,-4.5)
points(x,log(colMeans(convergence2[,x1])), col=2, pch=2, cex=0.5);lines(x,log(colMeans(convergence2[,x1])),col=2,lty=2)
points(x,log(colMeans(convergence1o[,x1])), col=3, pch=3, cex=0.5);lines(x,log(colMeans(convergence1o[,x1])),col=3,lty=3)
points(x,log(colMeans(convergence2o[,x1])), col=4, pch=4, cex=0.5);lines(x,log(colMeans(convergence2o[,x1])),col=4,lty=4)
points(x,log(colMeans(convergence1[,x1])), col=1, pch=1, cex=0.5);lines(x,log(colMeans(convergence1[,x1])),col=1,lty=1)
points(x,log(colMeans(convergence1f[,x1])), col=5, pch=5, cex=0.5);lines(x,log(colMeans(convergence1f[,x1])),col=5,lty=5)
legend("topright", modelsmethods, pch=(1:6), lty=(1:6), col=1:6, bty="n", cex=1.4, lwd=1.4)
grid()
dev.off()

pdf("../images/rmse2d_conv.pdf", width=mywidth, height=myheight)
x=(0:161)*10
x1=1:162
plot(x,(colMeans(convergence2f[,x1])), t="b", lty=6, pch=6,col=6, ylab=c("Natural logarithm of mean RMSE"), xlab="Number of iterations", cex=0.5, cex.lab=1.4, cex.sub=1.4, ylim=c(0.0015,0.0111)) #ylim=c(-6.7,-4.5)
points(x,(colMeans(convergence2[,x1])), col=2, pch=2, cex=0.5);lines(x,(colMeans(convergence2[,x1])),col=2,lty=2)
points(x,(colMeans(convergence1o[,x1])), col=3, pch=3, cex=0.5);lines(x,(colMeans(convergence1o[,x1])),col=3,lty=3)
points(x,(colMeans(convergence2o[,x1])), col=4, pch=4, cex=0.5);lines(x,(colMeans(convergence2o[,x1])),col=4,lty=4)
points(x,(colMeans(convergence1[,x1])), col=1, pch=1, cex=0.5);lines(x,(colMeans(convergence1[,x1])),col=1,lty=1)
points(x,(colMeans(convergence1f[,x1])), col=5, pch=5, cex=0.5);lines(x,(colMeans(convergence1f[,x1])),col=5,lty=5)
legend("topright", modelsmethods, pch=(1:6), lty=(1:6), col=1:6, bty="n", cex=1.4, lwd=1.4)
grid()
dev.off()

#2d bar plot
#library(mapplots)
#barplot2D()
#library(ggplot2)
#ggplot(as.data.frame(durr0[[2]]),aes(strike,tenor))+geom_hex(bins=30)

#grid(nx = NULL, ny = nx, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
#dates <- as.character(datatable[,1])
#dates = as.Date(dates, "%Y-%m-%d")
#plot(dates, rhoLmed, xaxt = "n", type = "l")
#axis(1, dates, format(dates, "%b %d"), cex.axis = .7)

#require(ggplot2)
#test = as.data.frame(cbind(dates,rhoLmed))
#test$dates = as.Date(test$dates, "%Y-%m-%d")
#test$rhoLmed = as.numeric(as.character(test$rhoLmed))
#ggplot( data = test, aes( dates, rhoLmed )) + geom_line() 

#############################
## Variance Term Structure ##
#############################
#cbind(volofvol(VIX[,1],63),volofvol(VIX[,2],63),volofvol(VIX[,3],63),volofvol(VIX[,4],63),volofvol(VIX[,5],63),volofvol(VIX[,6],63))
#hrho1mm=rowMeans(hrho1m)
#homega1mm=rowMeans(homega1m)

#kappa.init = 2
#nu0 = VIX[,1]^2
#theta = VIX[,6]^2
#kt = kappa.init*tau1
#ey2int2=((exp(-kt)*(((2*kt-3)*exp(kt)+4*kt)*theta-4*kt*nu0))-((exp(-2*kt))*((2*exp(2*kt)-4*exp(kt)-1)*theta-2*nu0*exp(2*kt)+2*nu0)))/(2*kappa.init^3)
#exyint2 = kappa.init^(-2)*(exp(-kt)*(2*theta-nu0-kt*(nu0-theta))+(nu0+theta*(kt-2)))
#VIXadded2 = sqrt(VIX^2 - omega*rho*exyint2/tau1 + omega^2*ey2int2/4/tau1)

#omega = c(rep(homega1m[21],20),homega1m[-(1:20)])
#rho = c(rep(hrho1m[21],20),hrho1m[-(1:20)])
#omega1 = cbind(omega,omega,omega,omega,omega,omega)
#rho1 = cbind(rho,rho,rho,rho,rho,rho)
#VIXadded1m = sqrt(VIX^2*(1 - omega1*rho1*tau1/2 + omega1^2*tau1^2/3/4))
#hrho1ma=hcorr(log(S),VIXadded1m[,1]^2,21)
#homega1ma=volofvol(VIXadded1m[,1],21)

hrho1m=hcorr(log(S),VIX[,1]^2,21)
hrho3m=hcorr(log(S),VIX[,1]^2,63)
#cbind(hcorr(logspot,VIX[,1]^2,63),hcorr(logspot,VIX[,2]^2,63),hcorr(logspot,VIX[,3]^2,63),hcorr(logspot,VIX[,4]^2,63),hcorr(logspot,VIX[,5]^2,63),hcorr(logspot,VIX[,6]^2,63))
homega1m=volofvol(VIX[,1],21)
homega3m=volofvol(VIX[,1],63)
  
omega = c(VIX[1:62,6],homega3m[-(1:62)])
rho = c(rep(-0.1,62),hrho3m[-(1:62)])
omega1 = cbind(omega,omega,omega,omega,omega,omega)
rho1 = cbind(rho,rho,rho,rho,rho,rho)
VIXadded3m = VIX*sqrt(1 - omega1*rho1*tau1/2 + omega1^2*tau1^2/3/4)

#VIXadded3m2 = sqrt(RC2/tau1 - omega1*rho1*RC2*tau1/2 + omega1^2/4*RC2*tau1^2/3)
#hrho3ma=hcorr(log(S),VIXadded3m[,1]^2,63)
#homega3ma=volofvol(VIXadded3m[,1],63)
#hestonTSpar1m=varTSfitAltC(VIXadded1m^2,1,tau,2,6,2,TRUE,barrier=1,levelup=1.5,leveldn=2/3)
#hestonTSpar3m0=varTSfitC(VIXadded3m^2,1,tau,2,6,2,TRUE)
#summary(cumax(log(hestonTSpar3m[,2])))
hestonTSpar3m=varTSfitAltC(VIXadded3m^2,1,tau,2,6,2,TRUE,barrier=0.4,levelup=2,leveldn=0.6);

#summary(sqrt(hestonTSpar3m))
#summary(hestonTSpar3m2)
#plot(hestonTSpar3m[,1])
#plot(hestonTSpar3m[,2])

nu0=hestonTSpar3m[,1]
theta=hestonTSpar3m[,2]
kappa=hestonTSpar3m[,3]
kappa1=cbind(kappa,kappa,kappa,kappa,kappa,kappa)
kt=sweep(kappa1,MARGIN=2,tau,`*`)
#k1 = kappaReg(VIX[,1]^2)
#M=(as.numeric(datatable[,19])/100-as.numeric(datatable[,4])/100)/(5/12)
#kappaD = kappaDurr(theta,nu0,homega3mmed,hrho3mmed,M)
#summary(cbind(kappa,kappa2))

#library(emdbook) #lambert W function
#a = (VIX[,3]^2-VIX[,6]^2)/(VIX[,1]^2-VIX[,6]^2)
#k = (a*lambertW_base(-exp(-1/a)/a)+1)/a/0.25
#plot(ts(EMA(pmin(k,rep(10,n))),10))
#lines(kappa,col=2)
#summary(pmin(k,rep(5,n)))

##Schoebel-Zhu

hrho3msz=hcorr(log(S),VIX[,1],63)
homega3msz=sqrt(EMA(c(0,VIX[-1,1]-VIX[-n,1])^2,63)*252)
omegasz = c(VIX[1:62,6]/2,homega3msz[-(1:62)])
rhosz = c(rep(-0.1,62),hrho3msz[-(1:62)])
omegasz1 = cbind(omegasz,omegasz,omegasz,omegasz,omegasz,omegasz)
rhosz1 = cbind(rhosz,rhosz,rhosz,rhosz,rhosz,rhosz)

szTSpar = varTSfitAltC(VIXadded3m^2,2,tau,2,6,1.1,TRUE,max_it=8000,barrier=0.4,levelup=2,leveldn=0.6,omega=omegasz)

#summary(szTSpar)
#szTSpar[1:5,]
#summary(sqrt(hestonTSpar3m))
#summary(shzTSpar)
#plot(szTSpar[,2])

## from I. Clark
## E[vol] = sqrt(E[var])*(1-1/8*varofvar/E[var]^2)
## W = vol^2*T
## from gatheral
## var[W_T] = theta*T * omega^2/kappa^2 + O(T)

ekt=exp(-kt)
timefactor = 2*(nu0-theta)*(1-2*kt*ekt-ekt^2) + theta*(4*ekt-3+2*kt-ekt^2)
varofvar = omega1^2/2/tau1^2/kappa1^3*timefactor
convadj = (1-1/8*varofvar/VIXadded3m^4)
convadj[,1]=pmax(convadj[,1],rep(0.00001,n))
VolIndex = VIXadded3m*convadj

#############################
## Perfect calibrated p&o  ##
#############################
#test.feller=feller(theta,homega3m,hrho3m,kappa,nu0);summary(test.feller>0)
#w=which(test.feller<0)
#homega3mf = homega3m
#homega3mf[w] = sqrt(0.9*2*kappa[w]*theta[w])

myparams = cbind(hestonTSpar3m[,2],omega,rho,hestonTSpar3m[,3],hestonTSpar3m[,1])
myparamsf = myparams
myparamsf[,2] = pmin(omega,sqrt(1.99*myparamsf[,1]*myparamsf[,4]))
myparamsz=cbind(szTSpar[,2],omegasz,rhosz,szTSpar[,3],szTSpar[,1])

t0=Sys.time(); perfectCal3 = opterCv(myparams,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600); Sys.time()-t0;
t0=Sys.time(); perfectCal3f = opterCv(myparamsf,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,fellerMust=1); Sys.time()-t0;
t0=Sys.time(); perfectCal3sz = opterCv(myparamsz,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=2,err=TRUE,max_it=1600); Sys.time()-t0;
t0=Sys.time(); perfectCalOUOUsym = opterCv(as.matrix(myparamsz),mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=5,err=TRUE,max_it=1600); Sys.time()-t0;

summary(perfectCalOUOUsym)
perfectCalOUOUsym[170,]
which(perfectCalOUOUsym[,6]>1)
summary(myparamsz)
i=170
nlOpter(myparamsz[i,],mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,5,1,5,TRUE,1600)
opterC(as.numeric(myparamsz[170,]),1,6,170,pair,opType=5,errType=1,modType=5,max_it=1600,err=TRUE)

perfectCal3alt = calibClean(perfectCal3,myparams,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,barrier=0.4,levelup=2)
perfectCal3falt = calibCleanF(perfectCal3f,myparams,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,fellerMust=1,barrier=0.4,levelup=2,leveldn=100,n)
perfectCal3falt[,2]=perfectCal3falt[,2]*0.999
test.fellerP=feller(perfectCal3falt);summary(test.fellerP>0)

plot(log(perfectCal3falt[,1]))
plot(log(perfectCal3falt[,2]))
plot(log(perfectCal3falt[,5]))
which((perfectCal3falt[,1])>1)

omegaP = perfectCal3alt[,2]
rhoP = perfectCal3alt[,3]

#perfectCal3faltf = calibNA(perfectCal3falt,myparams,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,c(2,1/2,1,2,1))

calibOutliers2(perfectCal3falt,0.4)

#Schoebel-Zhu
#perfect cal from historic estimates

perfectCal3szalt = calibClean(perfectCal3sz,myparamsz,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=2,err=TRUE,max_it=1600,barrier=0.4,levelup=2)
calibOutliers2(perfectCal3sz,0.4)
plot((perfectCal3szalt[,2]))
plot((perfectCal3sz[,2]))

plot((perfectCalOUOUsym[,5]))

which((perfectCal3sz[,2])>0.7)

summary(perfectCal3sz)
szTSperf = szTSvar(cbind(perfectCal3sz[,5],perfectCal3sz[,1],perfectCal3sz[,4],perfectCal3sz[,2]))

summary(szTSperf)

#perfectCal3=read.table("perfectCal3.txt")
#perfectCal3f=read.table("perfectCal3f.txt")
#perfectCal3sz=read.table("perfectCal3sz.txt")
#labMSE = read.table("labMSE.txt")
#labMAE = read.table("labMAE.txt")
#labMAPE = read.table("labMAPE.txt")
#durrMSE = read.table("durrMSE.txt")
#durrMAE = read.table("durrMAE.txt")
#durrMAPE = read.table("durrMAPE.txt")

#############################
## Starting Points Making  ##
#############################
rhoD=cbind(smile2rho(sqrt(nu0),SS[,1],CC[,1]),smile2rho(sqrt(nu0),SS[,2],CC[,2]),smile2rho(sqrt(nu0),SS[,3],CC[,3]),smile2rho(sqrt(nu0),SS[,4],CC[,4]),smile2rho(sqrt(nu0),SS[,5],CC[,5]),smile2rho(sqrt(nu0),SS[,6],CC[,6]))
omegaD=cbind(smile2omega(sqrt(nu0),SS[,1],CC[,1]),smile2omega(sqrt(nu0),SS[,2],CC[,2]),smile2omega(sqrt(nu0),SS[,3],CC[,3]),smile2omega(sqrt(nu0),SS[,4],CC[,4]),smile2omega(sqrt(nu0),SS[,5],CC[,5]),smile2omega(sqrt(nu0),SS[6],CC[6]))
omegaDmed = sqrt(rowMedians(omegaD^2))
rhoDmed = rowMedians(rhoD)
omegaDA = smile2omega(VIXadded3m[,1],SS[,1],CC[,1])
rhoDA = smile2rho(VIXadded3m[,1],SS[,1],CC[,1])

#R1 = log(P1 +1 +1/2*P2+1/6*P3+1/24*P4)
RC2 = (P2-P1^2)
RC3 = (P3-3*P1*P2+2*P1^3)
RC4 = (P4-4*P1*P3+6*P1^2*P2-3*P1^4)
RC5 = (P5-5*P4*P1+10*P3*P1^2-10*P2*P1^3+5*P1^5-P1^5)
RC6 = (P6 - 6*P5*P1 + 15*P4*P1^2 - 20*P3*P1^3 + 15*P2*P1^4 - 5*P1^6)

#Ev = theta + (nu0-theta)*exp(-kt)
#A = theta*tau1 + (nu0-theta)*(1-exp(-kt))/kappa1 #=Eintv
#A=nu0*tau1
#A = VIX^2*tau1
A = RC2
B = sqrt(1+A)

ey2e0 = B^(-7)*(4*B^9-8*B^8+4*B^7+(4*B^7-8*B^6+4*B^4)*RC2+(4*B^2-4*B^4)*RC3+(5-4*B^2)*RC4)
#-5*RC5)
ey2e2 = RC6/B^10 + (2*RC5)/B^8 + (5*RC4)/B^6 - (4*RC4)/B^5 + (4*RC3)/B^5 - (4*RC3)/B^3 + (4*RC2)/B^3 + 4*B^2 - (8*RC2)/B - 8*B + 4*RC2 + 4
exye2 = 2*RC2*(1/B-1)+RC3/B^3+RC4/B^5 + (RC6/2/B^10 + (RC5)/B^8 + (5/2*RC4)/B^6 - (2*RC4)/B^5 + (2*RC3)/B^5 - (2*RC3)/B^3 + (2*RC2)/B^3 + 2*B^2 - (4*RC2)/B - 4*B + 2*RC2 + 2)
exye0 = B^(-7)*(2*B^9 - 4*B^8 + 2*B^7 + (2*B^4 - 2*B^6)*RC2 + (2*B^2 - B^4)*RC3 + (10/4 - B^2)*RC4)
exye0 = (8*B^9-16*B^8+8*B^7-8*B^6*RC2-4*B^4*RC3+8*B^4*RC2-4*B^2*RC4+8*B^2*RC3+10*RC4)/(4*B^7)

exye = (2*B^8-4*B^7+2*B^6+ (2*B^3-2*B^5)*RC2+(2*B^2-B^3)*RC3+1/2*RC4)/B^6
ey2e = (4*B^8-8*B^7+4*B^6+ (4*B^6-8*B^5+4*B^3)*RC2+(4*B^2-4*B^3)*RC3+RC4)/B^6

#skosnosc=((4*B^3-4*B^2)*Eintv+(2*omega*rho*exyint-omega^2*ey2int)*B^3)/2/Eintv^(3/2)
#kurtoza=(4*B^6-8*B^5+8*B^4-4*B^3)*Eintv-4*B^8+8*B^7+(4*omega*rho*exyint-omega^2*ey2int-4)*B^6+(2*omega^2*ey2int-4*omega*rho*exyint)*B^5/Eintv^2

k1=kappa
nu0proxy = RC2[,1]/tau1[,1]
exyint = k1^(-2)*(exp(-k1*tau1)*(2*theta-nu0-k1*tau1*(nu0-theta))+(nu0+theta*(k1*tau1-2)))
ey2int = ((theta*(2*k1*tau1 - 5) + 2*nu0) + 4*exp(-k1*tau1)*(theta*(k1*tau1 + 1) - k1*tau1*nu0) + (theta - 2*nu0)*exp(-2*k1*tau1))/(2*k1^3)
ey2intL0 = nu0proxy*tau1^3/3
ey2intL = RC2*tau1^2/3
exyintL0 = nu0proxy*tau1^2/2
exyintL = RC2*tau1/2

omegaL2 = (ey2e/ey2intL)
omegaLmed = sqrt(rowMedians(ey2e/ey2intL))
omegaLmid = sqrt(apply(omegaL2,1,max)/2+apply(omegaL2,1,min)/2)
omegaLmed2 = sqrt(rowMedians(ey2e)/rowMedians(ey2intL))
omegaLm = sqrt(rowMeans(ey2e/ey2intL))
omrhoL = exye/exyintL
rhoLmed = rowMedians(exye/exyintL)/omegaLmed
rhoLmid = (apply(omrhoL,1,max)/2+apply(omrhoL,1,min)/2)/omegaLmid
rhoLm = rowMeans(exye/exyintL)/omegaLm
rhoLmed2 = rowMedians(exye)/rowMedians(exyintL)/omegaLmed2

sigma=datatable[,c(4,9,14,19,24,29)]/100; #sigma=VIX;
k2=RC2;k3=l1*(k2)^(3/2);k4=l2*k2^2
IRC2=k2;IRC3=k3;IRC4=k4+3*k2^2
mean(IRC3/RC3)

#############################
## Specific point
#############################
plot(ts(omegaLmed))
m=350; a=1; b=6
omegaP[m]*rhoP[m]
omegaL0[m,];rhoL0[m,]
omegaL[m,]*rhoL[m,]
median(omegaLO[m,]);median(rhoLO[m,])
median(omegaL[m,]);median(rhoL[m,])
omegaD[m,];rhoD[m,];
opt=opter(c(theta[m],omegaD[m,1],rhoD[m,1],kappa[m],nu0[m]),"Heston","E5",a,b,m,pair="EURUSD");opt
opt2=opter(c(theta[m],omegaD[m,1],rhoD[m,1],kappa[m],nu0[m]),"Heston","E2",a,b,m,pair="EURUSD");opt2
opt2=opter(c(theta[m],homega3mmed[m],hrho3mmed[m],kappa[m],nu0[m]),"Heston","E2",a,b,m,pair="EURUSD");opt2
data.load1(datatable,m)

j=3
plotLines(y=vol,x=K,j)
yhat1=vol[j,3]+(K[j,]/S-1)*SS[m,j]+(K[j,]/S-1)^2*CC[m,j]/2
points(y=yhat1,x=K[j,],col=3);lines(y=yhat1,x=K[j,],col=3);

yy=fitsmile(c(theta[m],omegaD[m,1],rhoD[m,1],kappa[m],nu0[m]),S,vol,rd,rf,0,tau)
points(y=yy[j,],x=K[j,],col=2);lines(y=yy[j,],x=K[j,],col=2);

yyL=fitsmile(c(theta[m],omegaL[m,3],rhoL[m,3],kappa[m],nu0[m]),S,vol,rd,rf,0,tau)
yyL=fitsmile(c(theta[m],omegaLmed2[m],rhoLmed2[m],kappa[m],nu0[m]),S,vol,rd,rf,0,tau)
points(y=yyL[j,],x=K[j,],col=4);lines(y=yyL[j,],x=K[j,],col=4);

#############################
## Fit testing             ##
#############################
#hist1m=tester(cbind(hestonTSpar1m[,2],homega1m,hrho1m,hestonTSpar1m[,3],hestonTSpar1m[,1]),F,K,rd,mktPrice,vega,tau,model=1);summary(ts(hist1m[[1]]))
hist3m=tester(cbind(hestonTSpar3m[,2],homega3m,hrho3m,hestonTSpar3m[,3],hestonTSpar3m[,1]),F,K,rd,mktPrice,vega,tau,model=1);summary(ts(hist3m[[1]]))
#paramsDurr = cbind(perfectCal3alt2[,1],omegaDA,rhoDA,perfectCal3alt2[,4],perfectCal3alt2[,5])
#paramsDurr = cbind(VIXadded3m[,6]^2,omegaDA,rhoDA,2,VIXadded3m[,1]^2)
paramsDurr = cbind(hestonTSpar3m[,2],omegaDA,rhoDA,hestonTSpar3m[,3],hestonTSpar3m[,1])
durr=tester(paramsDurr,F,K,rd,mktPrice,vega,tau,model=1);summary(ts(durr[[1]]))
#paramsLmed = cbind(VIXadded3m[,6]^2,omegaLmed,rhoLmed,2,VIXadded3m[,1]^2)
#paramsLmed = cbind(perfectCal3alt2[,1],omegaLmed,rhoLmed,perfectCal3alt2[,4],perfectCal3alt2[,5])
paramsLmed = cbind(hestonTSpar3m[,2],omegaLmed,rhoLmed,hestonTSpar3m[,3],hestonTSpar3m[,1])
lab=tester(paramsLmed,F,K,rd,mktPrice,vega,tau,model=1);summary(ts(lab[[1]]))

#t0=Sys.time(); perfectCalDurr = opterCv(paramsDurr,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600); Sys.time()-t0;
#t0=Sys.time(); perfectCalLmed = opterCv(paramsLmed,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600); Sys.time()-t0;

lab_err <- data.frame(err=lab[[1]])
durr_err <- data.frame(err=durr[[1]])
lab_err$met <- 'lab'
durr_err$met <- 'durr'
metLengths <- rbind(lab_err, durr_err)
library(ggplot2)
ggplot(metLengths, aes(err, fill = met)) + geom_density(alpha = 0.2)

plot(ts(lab1-durr0))
lab1mix=tester(datatable,tau,hestonTSpar[,1],hestonTSpar[,2],hestonTSpar[,3],(rho1w+rho1)/2,omega1)
summary(lab1mix)
plot(ts(lab1mix-lab1w))
plot(ts(log(lab1mix*100)))
sd(lab1mix)


#############################
## Fit testing SZ          ##
#############################
#from VolIndex
# 1. policzyć Vol używając Var i convexity adj
# 2. skalibrować vol term struct standardowo
# 3. policzyć uproszczone omegaSZ i rhoSZ używając nu0 z powyższej
szTSpar1 = varTSfitAltC(VolIndex,1,tau,2,6,0.95,TRUE,max_it=8000,barrier=1,levelup=1.6,leveldn=0.4)
nu0SZ=szTSpar1[,1]
nu0SZ1 = cbind(nu0SZ,nu0SZ,nu0SZ,nu0SZ,nu0SZ,nu0SZ)
omegaSZmat = sqrt(( -nu0SZ1^2/tau1 + sqrt(nu0SZ1^4 + ey2e/2)/tau1))
omegaSZ = (sqrt(rowMeans( -nu0SZ1^2/tau1 + sqrt(nu0SZ1^4 + ey2e/2)/tau1)))
omegaSZ1 = matrix(rep(omegaSZ,6),n,6)
rhoSZ = rowMedians(exye / (2*nu0SZ1^2*tau1 + omegaSZ^2*tau1^2))/omegaSZ
sz1=tester(cbind(szTSpar1[,2],omegaSZ,rhoSZ,szTSpar1[,3],szTSpar1[,1]),F,K,rd,mktPrice,vega,tau,model=2);summary(ts(sz1[[1]]))
sz1=tester(cbind(szTSpar1[,2],omegaLmed/2,rhoLmed,szTSpar1[,3],szTSpar1[,1]),F,K,rd,mktPrice,vega,tau,model=2);summary(ts(sz1[[1]]))

plot2(szTSpar[,1],sqrt(hestonTSpar3m[,1]))

omegaSZ = sqrt(rowMeans(ey2e/4/RC2))
#rhoSZ = rowMedians(exye/2/RC2)/omegaSZ

### druga możliwość
# 1. Policzyć omegaLmed/2
# 2. skalibrować var-sz term struct używając danych z omegaLmed/2
# 3. policzyć uproszczone omegaSZ i rhoSZ używając nu0 z powyższej
nu0SZ=VIX[,1]
nu0SZ1 = cbind(nu0SZ,nu0SZ,nu0SZ,nu0SZ,nu0SZ,nu0SZ)
#Emv2 = VIX^2*tau1
Emv2 = (nu0SZ1^2*tau1 + omegasz1^2*tau1^2/2)
ey2SZ = 4*omegasz1^2*Emv2
exySZ = 2*omegasz1*rhosz1*Emv2 
VIXaddedSZ = sqrt(VIX^2 - exySZ/tau1 + ey2SZ/4/tau1)
VIXaddedSZ = sqrt(1 - 2*omegasz1*rhosz1 + 4*omegasz1^2/4)

szTSpar2 = varTSfitAltC(VIXaddedSZ^2,2,tau,2,6,0.95,TRUE,max_it=8000,barrier=1,levelup=1.6,leveldn=0.4,omega=omegaLmed/2)
nu0SZ=szTSpar2[,1]
nu0SZ1 = cbind(nu0SZ,nu0SZ,nu0SZ,nu0SZ,nu0SZ,nu0SZ)
omegaSZ2 = (sqrt(rowMeans( -nu0SZ1^2/tau1 + sqrt(nu0SZ1^4 + ey2e/2)/tau1)))
omegaSZ21 = matrix(rep(omegaSZ2,6),n,6)
rhoSZ2 = rowMedians(exye / (2*nu0SZ1^2*tau1 + omegaSZ2^2*tau1^2))/omegaSZ2
sz2=tester(cbind(szTSpar2[,2],omegaSZ2,rhoSZ2,szTSpar2[,3],szTSpar2[,1]),F,K,rd,mktPrice,vega,tau,model=2);summary(ts(sz2[[1]]))
summary(cbind(omegaSZ2,rhoSZ2))
summary(cbind(omegaSZ2,omegaSZ))

#############################
## Test Feller Condition   ##
#############################
test.fellerP=feller(perfectCal3alt);summary(test.fellerP>0)
test.fellerL=feller(paramsLmed);summary(test.fellerL>0)
test.fellerD=feller(paramsDurr);summary(test.fellerD>0)

plot2(omegaLmed,omegaD[,1])

cumFeller = as.numeric((test.fellerP<0))
cumFeller[1:62] = rep(0,62)
cumFeller = cumsum(cumFeller)/(1:n)

outcome = as.numeric(summary(test.fellerP>0)[2:3])
outcome[2]/sum(outcome) #=0.3131393

#lab
test.feller=feller(cbind(hestonTSpar[,2],omegaL2/3,rhoL,hestonTSpar[,3],hestonTSpar[,1]));summary(test.feller>0)
#zatem ekspercko wybieram 3
w=which(test.fellerL<0)
omegaL2[w] = sqrt(0.9*2*hestonTSpar[w,3]*hestonTSpar[w,2])

#durr
test.feller=feller(cbind(hestonTSpar[,2],0.59*omegaD2,rhoD,hestonTSpar[,3],hestonTSpar[,1]));summary(test.feller>0)
#zatem ekspercko wybieram 4
w=which(test.fellerD<0)
omegaD2[w] = sqrt(0.9*2*hestonTSpar[w,3]*hestonTSpar[w,2])
#############################
## Test calibration risk   ##
#############################
t0=Sys.time()
labMSE = opterCv(paramsLmed,mktPrice,vega,F,K,rd,tau,opType=3,errType=1,modType=1,err=TRUE,max_it=1600)
labMAE = opterCv(paramsLmed,mktPrice,vega,F,K,rd,tau,opType=3,errType=2,modType=1,err=TRUE,max_it=1600)
labMAPE = opterCv(paramsLmed,mktPrice,vega,F,K,rd,tau,opType=3,errType=3,modType=1,err=TRUE,max_it=1600)
Sys.time()-t0

t0=Sys.time()
durrMSE = opterCv(paramsDurr,mktPrice,vega,F,K,rd,tau,opType=3,errType=1,modType=1,err=TRUE,max_it=1600)
durrMAE = opterCv(paramsDurr,mktPrice,vega,F,K,rd,tau,opType=3,errType=2,modType=1,err=TRUE,max_it=1600)
durrMAPE = opterCv(paramsDurr,mktPrice,vega,F,K,rd,tau,opType=3,errType=3,modType=1,err=TRUE,max_it=1600)
Sys.time()-t0

par1d=calibRisk(1,durrMSE,durrMAE,durrMAPE)
par1l=calibRisk(1,labMSE,labMAE,labMAPE)
par4d=calibRisk(4,durrMSE,durrMAE,durrMAPE)
par4l=calibRisk(4,labMSE,labMAE,labMAPE)
par5d=calibRisk(5,durrMSE,durrMAE,durrMAPE)
par5l=calibRisk(5,labMSE,labMAE,labMAPE)
errMSEd = durrMSE[,6]
errMSEl = labMSE[,6]

summary(durrMAE)

summary(cbind(par1d,par1l,par4d,par4l,par5d,par5l))
summary(cbind(errMSEd,errMSEl))

#############################
## 2-factor models tester  ##
#############################
## 1. OUOU ##
ouousymTS = varTSfitAltC(VIXadded3m^2/2,2,tau,2,6,1.1,TRUE,max_it=8000,barrier=0.4,levelup=2,leveldn=0.6,omega=omegaLmed/2/sqrt(2))
omega1=omegaLmed/2/sqrt(2)*(sqrt(1-rhoLmed^2)+rhoLmed)
omega2=omegaLmed/2/sqrt(2)*(sqrt(1-rhoLmed^2)-rhoLmed)

myparams1o = factor2from1(perfectCalOUOUsym[,-6],c(1,1,1,1,1))
myparams2o = cbind(ouousymTS[,2],omega1,0.99,ouousymTS[,3],ouousymTS[,1],ouousymTS[,2],omega2,-0.99,ouousymTS[,3],ouousymTS[,1])

lab2do1=tester(myparams1o,F,K,rd,mktPrice,vega,tau,model=4);summary(ts(lab2do1[[1]]))
lab2do2=tester(myparams2o,F,K,rd,mktPrice,vega,tau,model=4);summary(ts(lab2do2[[1]]))

## 2. Bates no Feller ##
myparams1 = factor2from1(perfectCal3alt[,-6],c(0.5,1,1,1,0.5))
myparams2 = factor2from1(paramsLmed,c(0.5,1,1,1,0.5))

lab2db1=tester(myparams1,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2db1[[1]]))
lab2db2=tester(myparams2,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2db2[[1]]))

## 3. Bates Feller ##
myparams1f = factor2from1(perfectCal3falt[,-6],c(0.5,1,1,1,0.5))
myparams1f[,2] = pmin(myparams1f[,2],sqrt(1.99*myparams1f[,1]*myparams1f[,4]))
myparams1f[,7] = myparams1f[,2] 

myparams2f = factor2from1(paramsLmed,c(0.5,1,1,1,0.5))
myparams2f[,3] = rep(0.99,n)
myparams2f[,8] = rep(-0.99,n)
myparams2f[,2] = pmin(omegaLmed*(sqrt(1-rhoLmed^2)+rhoLmed),sqrt(1.99*myparams2f[,1]*myparams2f[,4]))
myparams2f[,7] = pmin(omegaLmed*(sqrt(1-rhoLmed^2)-rhoLmed),sqrt(1.99*myparams2f[,1]*myparams2f[,4]))

lab2df1=tester(myparams1f,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2df1[[1]]))
lab2df2=tester(myparams2f,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2df2[[1]]))

## 3b. Bates Feller - other ##
myparams3f = factor2from1(paramsLmed,c(0.5,1,1,1,0.5))
myparams3f[,2] = pmin(myparams3f[,2],sqrt(1.99*myparams3f[,1]*myparams3f[,4]))
myparams3f[,7] = myparams3f[,2] 

myparams3f = factor2from1(paramsLmed,c(0.5,1,1,1,0.5))
myparams3f[,2] = pmin(myparams3f[,2],sqrt(1.99*myparams3f[,1]*myparams3f[,4]))
myparams3f[,7] = myparams3f[,2] 

myparams4f = factor2from1(perfectCal3alt[,-6],c(0.5,1,1,1,0.5))
myparams4f[,2] = pmin(myparams4f[,2],sqrt(1.99*myparams4f[,1]*myparams4f[,4]))
myparams4f[,7] = myparams4f[,2] 

myparamsff = myparams
myparamsff[,2] = pmin(omega,sqrt(0.99*myparamsff[,1]*myparamsff[,4]))
t0=Sys.time(); perfectCalFellerx2 = opterCv(as.matrix(myparamsff),mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,fellerMust=2); Sys.time()-t0;
myparams5f = factor2from1(perfectCalFellerx2[,-6],c(1,1,1,1,1))

myparams6f = factor2from1(myparams,c(0.5,1,1,1,0.5))
myparams6f[,2] = pmin(myparams6f[,2],sqrt(1.99*myparams6f[,1]*myparams6f[,4]))
myparams6f[,7] = myparams6f[,2] 

lab2df3=tester(myparams3f,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2df3[[1]]))
lab2df4=tester(myparams4f,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2df4[[1]]))
lab2df5=tester(myparams5f,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2df5[[1]]))
lab2df6=tester(myparams6f,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2df6[[1]]))

## final testing
t0=Sys.time(); perfectCalOUOU1 = opterCv(myparams1o,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=4,err=TRUE,max_it=1600,logfile=paste0(workfolder,"/nl_out1o.txt")); Sys.time()-t0;
t0=Sys.time(); perfectCalOUOU2 = opterCv(myparams2o,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=4,err=TRUE,max_it=1600,logfile=paste0(workfolder,"/nl_out2o.txt")); Sys.time()-t0;
t0=Sys.time(); perfectCalBates1 = opterCv(myparams1,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,logfile=paste0(workfolder,"/nl_out1b.txt")); Sys.time()-t0;
t0=Sys.time(); perfectCalBates2 = opterCv(myparams2,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,logfile=paste0(workfolder,"/nl_out2b.txt")); Sys.time()-t0;
t0=Sys.time(); perfectCalBates1f = opterCv(myparams1f,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,fellerMust=1,logfile=paste0(workfolder,"/nl_out1f.txt")); Sys.time()-t0;
t0=Sys.time(); perfectCalBates2f = opterCv(myparams2f,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,fellerMust=1,logfile=paste0(workfolder,"/nl_out2f.txt")); Sys.time()-t0;

## other testing
t0=Sys.time(); perfectCalBates3f = opterCv(myparams3f,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,fellerMust=1,logfile=paste0(workfolder,"/nl_out3f.txt")); Sys.time()-t0;
t0=Sys.time(); perfectCalBates4f = opterCv(myparams4f,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,fellerMust=1,logfile=paste0(workfolder,"/nl_out4f.txt")); Sys.time()-t0;
t0=Sys.time(); perfectCalBates5f = opterCv(myparams5f,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,fellerMust=1,logfile=paste0(workfolder,"/nl_out5f.txt")); Sys.time()-t0;
t0=Sys.time(); perfectCalBates6f = opterCv(myparams6f,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,fellerMust=1,logfile=paste0(workfolder,"/nl_out6f.txt")); Sys.time()-t0;

#W metodzie z przeksztalceniem ucina się srednio mniej
2*summary(omegaLmed-pmin(omegaLmed,sqrt(1.99*myparams3f[,1]*myparams3f[,4])))
summary(omegaLmed*(sqrt(1-rhoLmed^2)+rhoLmed)-pmin(omegaLmed*(sqrt(1-rhoLmed^2)+rhoLmed),sqrt(1.99*myparams2f[,1]*myparams2f[,4])))+summary(omegaLmed*(sqrt(1-rhoLmed^2)-rhoLmed)-pmin(omegaLmed*(sqrt(1-rhoLmed^2)-rhoLmed),sqrt(1.99*myparams2f[,1]*myparams2f[,4])))

perfectCalBates3falt = calibNA(perfectCalBates3f,myparams3f,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600)

## 4. Simple starting points

simple = cbind(VIXadded3m[,6]^2/2,sqrt(2)*VIXadded3m[,6],0.99,2.1,VIXadded3m[,1]^2/2)
simplepars1 = cbind(simple,sweep(simple,2,c(1,1,-1,1,1),`*`))
lab2ds1=tester(simplepars1,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2ds1[[1]]))

simple = cbind(VIXadded3m[,6]/sqrt(2),sqrt(2)*VIXadded3m[,6]/2,0.99,1.1,VIXadded3m[,1]/sqrt(2))
simplepars2 = cbind(simple,sweep(simple,2,c(1,1,-1,1,1),`*`))
lab2ds2=tester(simplepars2,F,K,rd,mktPrice,vega,tau,model=4);summary(ts(lab2ds2[[1]]))

#perfectCalOUOUsym=read.table("perfectCalOUOUsym.txt")
#perfectCalOUOU1=read.table("perfectCalOUOU1.txt")
#perfectCalOUOU2=read.table("perfectCalOUOU2.txt")
perfectCalBates1=read.table("perfectCalBates1.txt")
#perfectCalBates2=read.table("perfectCalBates2.txt")
perfectCalBates3=read.table("perfectCalBates3.txt")
perfectCalBates1f=read.table("perfectCalBates1f.txt")
#perfectCalBates2f=read.table("perfectCalBates2f.txt")
perfectCalBates3f=read.table("perfectCalBates3f.txt")
perfectCalBates4f=read.table("perfectCalBates4f.txt")
perfectCalBates5f=read.table("perfectCalBates5f.txt")
perfectCalFellerx2=read.table("perfectCalFellerx2.txt")

#perfectCalBates1 = rbind(perfectCalBates1_1[1:667,],perfectCalBates1_2[668:n,])
#perfectCalBates1f = rbind(perfectCalBates1f_1[1:667,],perfectCalBates1f_2[668:n,])

i=11
perfectCalOUOU2[i,]
t0=Sys.time();test=nlOpter(myparams2o[i,],mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,10,1,4,TRUE,1600,0,paste0(workfolder,"/pop1.txt"));Sys.time()-t0; test
perfectCalOUOU2[i,]=test

summary(cbind(perfectCalBates1[,11],perfectCalBates2[,11],perfectCalBates1f[,11],perfectCalBates2f[,11],perfectCalBates3f[,11],perfectCalOUOU1[,11],perfectCalOUOU2[,11]))

cbind(summary(perfectCalBates2f[,11]),summary(perfectCalBates6f[,11]),summary(perfectCalBates5f[,11]),summary(perfectCalBates3f[,11]),summary(perfectCalBates1f[,11]),summary(perfectCalFellerx2[,6]))

plot2(perfectCalBates1f[,11],perfectCalBates5f[,11])

plot(ts(perfectCal3alt2[,6]))
lines(ts(perfectCal3sz[,6]),col=2)

plot(convergence2o[,1],ylim=c(0,0.03))
points(convergence2o[,100],col=2)
points(convergence2[,100],col=3)
i=310

summary(perfectCalBates1[,11])

t0=Sys.time()
nlOpter(myparams3o[i,],mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,10,1,4,TRUE,1600)
Sys.time()-t0

plot(colMeans(convergence1[,91:111])-colMeans(convergence2[,91:111]))
mean(convergence1[,95])-mean(convergence2[,95])
median(convergence1[,160])-median(convergence2[,160])

#############################
## Model error regression  ##
#############################
datatable2=datatable/100
level=datatable2$EURATM1M
slope=datatable2$EURATM2Y-datatable2$EURATM1M
curv=datatable2$EURATM2Y+datatable2$EURATM1M-2*datatable2$EURATM6M
RR1M=datatable2$EUR10P1M-datatable2$EUR10C1M
RR3M=datatable2$EUR10P3M-datatable2$EUR10C3M
RR1Y=datatable2$EUR10P1Y-datatable2$EUR10C1Y
RR2Y=datatable2$EUR10P2Y-datatable2$EUR10C2Y
BF1M=datatable2$EUR10P1M+datatable2$EUR10C1M-2*datatable2$EURATM1M
BF3M=datatable2$EUR10P3M+datatable2$EUR10C3M-2*datatable2$EURATM3M
BF1Y=datatable2$EUR10P1Y+datatable2$EUR10C1Y-2*datatable2$EURATM1Y
BF2Y=datatable2$EUR10P2Y+datatable2$EUR10C2Y-2*datatable2$EURATM2Y

lm1a = lm(lab[[1]]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm1b = lm(durr[[1]]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm1 = lm(perfectCal3alt[,6]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm1 = lm(perfectCal3sz[,6]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm2a = lm(perfectCalBates1[,11]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm2b = lm(perfectCalBates2[,11]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm2c = lm(perfectCalBates1f[,11]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm2d = lm(perfectCalBates2f[,11]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm2e = lm(perfectCalOUOU1[,11]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm2f = lm(perfectCalOUOU2[,11]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)

library(texreg)
screenreg(list(lm1a, lm1b, lm2a, lm2b, lm2c, lm2d, lm2e, lm2f))
tab10 = texreg(list(m1, m2, m3, m4, m5, m6a, m6b), digits = 1, dcolumn = FALSE, booktabs = FALSE,use.packages = FALSE, label = "tab:10", caption = "Models explaining the probability of increases and decreases of implied volatility for the EURUSD 1M option in 1 quarter horizon",float.pos = "h!")
tab10 = gsub("hestonTSpar\\[, 1]", "$\\\\nu_{0,H}$", tab10)
tab10 = gsub("shzTSpar2\\[, 1]", "$\\\\nu_{0,S}^2$", tab10)
tab10 = gsub("hestonTSpar\\[, 2]", "$\\\\theta_H$", tab10)
tab10 = gsub("shzTSpar2\\[, 2]", "$\\\\theta_S^2$", tab10)
tab10 = gsub("kappa2shz", "$\\\\xi_S$", tab10)
tab10 = gsub("kappa2", "$\\\\xi_H$", tab10)
write.table(c(tab10,AUROC_stats),paste0(workfolder,"/r3.t10.txt"))

#############################
## 2-factor models other   ##
#############################
# 1. Comparison of my 2 different starting points for Bates no Feller
omegas = ey2e/ey2intL
omrho= exye/exyintL
rho1=apply(omrho/sqrt(omegas),1,absmax)
rho2=apply(omrho/sqrt(omegas),1,absmin)
wr=(rhoLmed-rho2)/(rho1-rho2)
myparams3 = params2factor(omegas,omrho,wr,VIXadded3m^2,mod=1,a=0)
myparams3[,2] = omegaLmed
myparams3[,7] = omegaLmed
myparams3[,3] = rho1
myparams3[,8] = rho2

omrho1=apply(omrho,1,absmax)
omrho2=apply(omrho,1,absmin)
omega1=apply(sqrt(omegas),1,max)
omega2=apply(sqrt(omegas),1,min)
wo=(omegaLmed^2-omega2^2)/(omega1^2-omega2^2)
wro=(rhoLmed*omegaLmed-omrho2)/(omrho1-omrho2)
b=wo/wro
omega1b = omega1*sqrt(b)
omega2b = omega2*sqrt((1-b*wro)/(1-wro))
rho1b = omrho1/omega1b
rho2b = omrho2/omega2b
myparams4 = params2factor(omegas,omrho,wro,VIXadded3m^2,mod=1,a=0)
myparams4[,2] = omega1b
myparams4[,7] = omega2b
myparams4[,3] = rho1b
myparams4[,8] = rho2b

myparams3=myparams1
w=1/2;m=0.9
myparams3[,2] = omegaLmed*(sqrt((rhoLmed^2-m^2)*w^2+(m^2-rhoLmed^2)*w)+w*rhoLmed)/m/w
myparams3[,7] = omegaLmed*(sqrt((rhoLmed^2-m^2)*w^2+(m^2-rhoLmed^2)*w)-w*rhoLmed)/m/w
myparams3[,3] = rep(0.9,n)
myparams3[,8] = rep(-0.9,n)

lab2db1=tester(myparams1,F,K,rd,mktPrice,vega,tau,model=3)
lab2db2=tester(myparams2,F,K,rd,mktPrice,vega,tau,model=3)
lab2db3=tester(myparams3,F,K,rd,mktPrice,vega,tau,model=3)
summary(ts(lab2db1[[1]]))
summary(ts(lab2db2[[1]]))
summary(ts(lab2db3[[1]]))

#omegaSZ = sqrt(rowMedians(ey2e/4/RC2))
#rhoSZ = rowMedians(exye/2/RC2)/omegaSZ
myparams1o = params2factor(omegas,omrho,rep(0.5,n),VIXadded3m^2,mod=2,a=0,omega)
myparams1o[,2] = omegaSZ
myparams1o[,7] = omegaSZ
myparams1o[,3] = rhoSZ
myparams1o[,8] = rhoSZ
lab2do1=tester(myparams1o,F,K,rd,mktPrice,vega,tau,model=4);summary(ts(lab2do1[[1]]))


nu0SZ=szTSpar[,1]
nu0SZ1 = cbind(nu0SZ,nu0SZ,nu0SZ,nu0SZ,nu0SZ,nu0SZ)
omegas = -nu0SZ1^2/tau1 + sqrt(nu0SZ1^4 + ey2e/2)/tau1
omrho = exye / (2*nu0SZ1^2*tau1 + omegas*tau1^2)
#omegas = ey2e/4/RC2
#omrho = exye/2/RC2
rho1=apply(omrho/sqrt(omegas),1,absmax)
rho2=apply(omrho/sqrt(omegas),1,absmin)
omrho1=apply(omrho,1,absmax)
omrho2=apply(omrho,1,absmin)
omega1=apply(sqrt(omegas),1,max)
omega2=apply(sqrt(omegas),1,min)
wo=(omegaSZ^2-omega2^2)/(omega1^2-omega2^2)
wro=(rhoSZ*omegaSZ-omega2*rho2)/(omega1*rho1-omega2*rho2)
wro2=(rhoSZ*omegaSZ-omrho2)/(omrho1-omrho2)
b=wo/wro
omega1b = omega1*sqrt(b)
omega2b = omega2*sqrt((1-b*wro)/(1-wro))
rho1b = omega1*rho1/omega1b
rho2b = omega2*rho2/omega2b
myparams2o = params2factor(omegas,omrho,wro,VIXadded3m^2,mod=2,a=0,omega)
#myparams2o = params2factor(omegas,omrho,wro,VIXaddedSZ^2,mod=3,a=0,omega)
myparams2o[,2] = omega1b
myparams2o[,7] = omega2b
myparams2o[,3] = rho1b
myparams2o[,8] = rho2b

omega1=sqrt(apply(omegas,1,max)*(1+a))
omega2=sqrt(apply(omegas,1,min)*(1-a))
omega1=matrix(rep(omega1,6),n,6)
omega2=matrix(rep(omega2,6),n,6)
w = (omegaLmed^2 - omega2^2)/(omega1^2-omega2^2)
w = (omegas - omega2^2)/(omega1^2-omega2^2)
#omegaLmed*rhoLmed = w*rho1*omega1 + (1-w)*rho2*omega2
rho2 = apply(omrho/sqrt(omegas),1,signabsmax)*apply(abs(omrho/sqrt(omegas)),1,max)
w = apply(w,1,median)
rho1 = (omegaLmed*rhoLmed - (1-w)*rho2*omega2 ) / w / omega1  
#rho2 = (omegaLmed*rhoLmed - w*rho1*omega1)/(1-w)/omega2 
myparams6 = params2factor(omegas,omrho,w,VIXadded3m^2,mod=1,a=0)
myparams6[,2] = omega1
myparams6[,7] = omega2
myparams6[,3] = rho1
myparams6[,8] = rho2
plot(ts(w))
lines(w,col=1)

rho1=apply(omrho/sqrt(omegas),1,absmax)
rho2=apply(omrho/sqrt(omegas),1,absmin)
wr=(rhoLmed-rho2)/(rho1-rho2)
myparams7 = params2factor(omegas,omrho,wr,VIXadded3m^2,mod=1,a=0)
myparams7[,2] = omegaLmed
myparams7[,7] = omegaLmed
myparams7[,3] = rho1
myparams7[,8] = rho2

omega1=apply(sqrt(omegas),1,max)
omega2=apply(sqrt(omegas),1,min)
wo=(omegaLmed^2-omega2^2)/(omega1^2-omega2^2)
myparams8 = params2factor(omegas,omrho,wo,VIXadded3m^2,mod=1,a=0)
myparams8[,2] = omega1
myparams8[,7] = omega2
myparams8[,3] = rhoLmed
myparams8[,8] = rhoLmed

omrho1=apply(omrho,1,absmax)
omrho2=apply(omrho,1,absmin)
wro=(rhoLmed*omegaLmed-rho2*omega2)/(rho1*omega1-rho2*omega2)
myparams8 = params2factor(omegas,omrho,wro,VIXadded3m^2,mod=1,a=0)
myparams8[,2] = omega1
myparams8[,7] = omega2
myparams8[,3] = omrho1/omega1
myparams8[,8] = omrho2/omega2

omrho1=apply(omrho,1,absmax)
omrho2=apply(omrho,1,absmin)
wro=(rhoLmed*omegaLmed-omrho2)/(omrho1-omrho2)
b=wo/wro
omega1b = omega1*sqrt(b)
omega2b = omega2*sqrt((1-b*wro)/(1-wro))
rho1b = omrho1/omega1b
rho2b = omrho2/omega2b
myparams11 = params2factor(omegas,omrho,wro,VIXadded3m^2,mod=1,a=0)
myparams11[,2] = omega1b
myparams11[,7] = omega2b
myparams11[,3] = rho1b
myparams11[,8] = rho2b
summary(rho2b)


summary(1/sqrt((1-b*wro)/(1-wro)))
plot2(rho1b,1/sqrt(b))

wo=(omegaLmed^2-omega2^2)/(omega1^2-omega2^2)
wro=(rhoLmed*omegaLmed-omega2*rho2)/(omega1*rho1-omega2*rho2)
b=wro/wo
rho1b = pmax(rho1*b,rep(-1,n))
rho2b = rho2*(1-b*wo)/(1-wo)
myparams10 = params2factor(omegas,omrho,wo,VIXadded3m^2,mod=1,a=0)
myparams10[,2] = omega1
myparams10[,7] = omega2
myparams10[,3] = rho1b
myparams10[,8] = rho2b


weights = weights1(omegas,omrho,a=0)
weights3 = weights1(omegas,omrho,a=0)
myparams5 = params2factor(omegas,omrho,weights3,VIXadded3m^2,mod=1,a=0)
myparams1 = params2factor(omegas,omrho,rep(0.5,n),VIXadded3m^2,mod=1,a=0)
myparams2 = params2factor(omegas,omrho,weights,VIXadded3m^2,mod=1,a=0)
myparams3 = myparams1
myparams1[,2] = omega1
myparams1[,7] = omega2
myparams1[,3] = rho1
myparams1[,8] = rho2
myparams4 = myparams1

#onlyomrho = opterCv26(cbind(hestonTSpar3m[,2],homega3m,hrho3m,hestonTSpar3m[,3],hestonTSpar3m[,1]),mktPrice,vega,F,K,rd,tau,opType=2,errType=1,modType=1,err=FALSE,max_it=1600)
omegasP = onlyomrho[,1:6]^2
omrhoP = onlyomrho[,1:6]*onlyomrho[,7:12]
w6 = weights2(omegasP,omrhoP,a=0)
w6 = weights2(omegas,omrho,a=0)
myparams5 = params2factor(omegas,omrho,w6,VIXadded3m^2,mod=1,a=0)

t0=Sys.time()
perfectCalBatesSM = opterCv(cbind(hestonTSparA2[,2]/2,omegaLmed,rhoLmed,hestonTSparA2[,3],hestonTSparA2[,1]/2,hestonTSparA2[,2]/2,omegaLmed,-rhoLmed,hestonTSparA2[,3],hestonTSparA2[,1]/2),mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600)
Sys.time()-t0;
write.table(perfectCalBatesICM, paste0(workfolder,"/BatesICM.txt"))

t0=Sys.time()
perfectCalBatesICM = opterCv(cbind(hestonTSparA2[,2]/2,omega1,0.99,hestonTSparA2[,3],hestonTSparA2[,1]/2,hestonTSparA2[,2]/2,omega2,-0.99,hestonTSparA2[,3],hestonTSparA2[,1]/2),mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600)
Sys.time()-t0;

#############################
## SZ                      ##
#############################
szTSpar=varTSfit(sqrt(VIX),tau,FALSE,2,2,5)
nu0SZ=szTSpar[,1]
thetaSZ=szTSpar[,2]
kappaSZ=szTSpar[,3]
hrho3msz=hrho(logspot,VIX[,1],63)
homega3msz=sqrt(EMA(c(0,VIX[-1,1]-VIX[-n,1])^2,63)*252)
omega = c(rep(homega3msz[63],62),homega3msz[-(1:62)])
rho = c(rep(hrho3msz[63],62),hrho3msz[-(1:62)])
exyintSZ = 2*(omega^2*tau1^3+3*nu0SZ^2*tau1^2)/6
ey2intSZ = 4*(omega^2*tau1^4+4*nu0SZ^2*tau1^3)/12
addon = - omega*rho*exyintSZ/tau1 + omega^2*ey2intSZ/4/tau1
VIXadded = sqrt(VIX^2 - omega*rho*exyintSZ/tau1 + omega^2*ey2intSZ/4/tau1)
plot(ts(rowMedians(VIX)-rowMedians(VIXadded)))
plot2(homega3msz,homega3m)
szTSpar0=varTSfit(sqrt(VIXadded),tau,FALSE,2,2,5)
summary(cbind(szTSpar[,3],szTSpar0[,3]))
nu0SZ=szTSpar0[,1]
thetaSZ=szTSpar0[,2]
kappaSZ=szTSpar0[,3]

# o2 (a o2 + b) = c
# a x2 + b x - c = 0
x = -b/2/a +- sqrt(b^2+4*a*c)/2/a

a = 1/3*tau1^4
b = 4/3*VIX[,1]^2*tau1^3
c = ey2e
omegaSZ = sqrt(rowMedians(-2*nu0SZ/tau1 + sqrt(16/9*nu0SZ^2*tau1^6+4/3*tau1^4*(ey2e))/(2/3*tau1^4)))
rhoSZ = rowMedians(exye/(2*(omegaSZ^2*tau1^3+3*nu0SZ*tau1^2)/6))/omegaSZ
labSZ=tester(datatable,pair,tau,nu0SZ,thetaSZ,kappaSZ,rhoSZ,omegaSZ,model="SZ");summary(ts(labSZ[[1]]))

omegaSZ2 = sqrt(rowMedians(ey2e/(1/3*tau1^2*(2*RC2+2*VIX[,1]^2*tau1))))
rhoSZ2 = rowMedians(exye/(1/3*tau1*(2*RC2+VIX[,1]^2*tau1)))/omegaSZ2

plot(ts(rhoSZ))
plot2(rhoSZ,rhoSZ2)
plot2(exye)
plot2(omegaSZ,omegaLmed2)
plot2(omegaSZ,omegaSZ2)

labSZ=tester(datatable,pair,tau,sqrt(nu0),sqrt(theta),kappa,rhoLmed2,omegaLmed2/2,model="SZ");summary(ts(labSZ[[1]]))
labSZ=tester(datatable,pair,tau,sqrt(nu0),sqrt(theta),kappa,rhoSZ2,omegaSZ2,model="SZ");summary(ts(labSZ[[1]]))
labSZ=tester(datatable,pair,tau,sqrt(nu0),sqrt(theta),kappa,rhoSZ,omegaSZ,model="SZ");summary(ts(labSZ[[1]]))

labSZ=tester(datatable,pair,tau,nu0SZ,thetaSZ,kappaSZ,rhoLmed2,omegaLmed2/2,model="SZ");summary(ts(labSZ[[1]]))
summary(perfectCalSZ[,6])

perfectCalSZ = matrix(0,1,7)
for (i in 1:62) {
	perfectCalSZ=rbind(perfectCalSZ,rep(NA,7))
}
for (i in 63:n) {
	opt=opter(c(thetaSZ[i],homega3msz[i],hrho3msz[i],kappaSZ[i],nu0SZ[i]),"SZ","ESZ2",1,6,i,pair,"MSE")
	perfectCalSZ=rbind(perfectCalSZ,opt)
}
perfectCalSZ=perfectCalSZ[-1,]
omegaPsz = perfectCalSZ[,2]
rhoPsz = perfectCalSZ[,3]

plot2(omegaPsz,omegaLmed2/2)
plot2(rhoPsz,rhoSZ)

opt5sz = matrix(0,1,7)
for (i in 1:(k/1)) {
	j=i#*20
	m=days[j]
	opt=opter(c(sqrt(theta[m]),omegaLmed2[m]/2,rhoLmed2[m],kappa[m],sqrt(nu0[m])),"SZ","ESZ5",1,6,m,pair,"MSE")
	opt5sz = rbind(opt5sz,opt)
}
opt5sz=opt5sz[-1,]

summary(perfectCalSZ)
summary(opt5)

#############################
## All other               ##
#############################
theta = 0.016443845
sigma = 0.207525426
rho = -0.398326265
kappa = 1.309675845
nu0 = 0.006120888

params = log(c(0.016443845,0.207525426,exp(atanh(-0.398326265)),1.309675845,0.006120888))
data.load1(datatable,602);
a=1;b=6

t0=Sys.time()

E5(params,data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun="MSE")

Sys.time()-t0

opter(c(0.012,0.1,-0.3,2,0.007),"Heston","E5",1,6,602,pair,"MSE")
sourceCpp("nlheston.cpp")
data.load1(datatable,602)
K1=as.numeric(t(K));vega1=as.numeric(t(vega));mktPrice1=as.numeric(t(mktPrice))
params=c(0.016443845,0.207525426,-0.398326265,1.309675845,0.006120888)
t0=Sys.time()
nlHeston(params,mktPrice1,vega1,S,K1,rd,rf,tau)
Sys.time()-t0


sqrt(E5(log(c(0.010423325104222933, 0.23081429767001443,exp(atanh(-0.38386799582959491)),5.4652722019197775, 0.0053903675937874837)),data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun="MSE"))
VIX[602,6]^2
hrho3m[602]
homega3m[602]

library(Rcpp)
Sys.getenv("PATH") 

plugin = Rcpp.plugin.maker(
include.before = "",
include.after = "",
#LinkingTo = unique(c(package, "Rcpp")),
#Depends = unique(c(package, "Rcpp")),
libs = paste0(workfolder,"/VectorTester.lib"),
Makevars = NULL,
Makevars.win = NULL,
package = "Rcpp"
)


registerPlugin("plugin1", 
                function() {
                    list(env = list(PKG_CXXFLAGS=paste0("-I",path.expand(scriptsfolder))))
                }
              )

registerPlugin("plugin1", plugin)
sourceCpp("MyFunc.cpp", verbose=TRUE, rebuild=TRUE)


Sys.setenv("PKG_CXXFLAGS"="-I/usr/include")
Sys.setenv("PKG_LIBS"="-L/usr/lib/x86_64-linux-gnu/ -lm -lmpc -lgmp -lmpfr")

system.time (output <- myFunc(df)) # see Rcpp function below
test(1.5,2.5)

# wrapper function to invoke helloA1
dyn.load("helloA1.so")
helloA1 <- function() {
  result <- .Call("helloA1")
}

greeting <- helloA1()
class(greeting)

cppFunction("int g(int n) { if (n < 2) return(n); return(g(n-1) + g(n-2)); }")
## Using it on first 11 arguments
sapply(0:10, g)
#set R_HOME=C:\Program Files\R\R-3.1.1\bin


alpha=5
tau=1
rd=0.01
K=1.3
F=1.25
S=1.25/exp(0.01)
hestongg(1.25/exp(0.01),1.3,2,0.01,0.2,-0.3,0,0.01,0.01,0,0,1)
batesa(1.25/exp(0.01),1.3,2,0.01,0.2,-0.3,0,0.01,0.01,0,0,1)
hestonaa2(1.25/exp(0.01),1.3,2,0.01,0.2,-0.3,0,0.01,0.01,0,0,1)

target=integrate(function(u) Re(exp(-1i*u*log(1.3/1.25))*((1-1i/u)/(1+u^2))*CFaa(u,kappa=2,theta=0.01,sigma=0.2,rho=-0.3,lambda=0,nu0=0.01,t=0,T=1)),1e-16,100)$value

u = 20+(0:120)*2/3
fx=Re(exp(-1i*u*log(1.3/1.25))*((1-1i/u)/(1+u^2))*CFaa(u,kappa=2,theta=0.01,sigma=0.2,rho=-0.3,lambda=0,nu0=0.01,t=0,T=1))
(trapz(u,fx)-target)*10^5

target2=integrate(function(u) Re(exp(-1i*u*log(1.3/1.25))*((1-1i/u)/(1+u^2))*CFaa(u,kappa=2,theta=0.01,sigma=0.2,rho=-0.3,lambda=0,nu0=0.01,t=0,T=1)),10,20)$value

u = 10+(0:40)/4
fx=Re(exp(-1i*u*log(1.3/1.25))*((1-1i/u)/(1+u^2))*CFaa(u,kappa=2,theta=0.01,sigma=0.2,rho=-0.3,lambda=0,nu0=0.01,t=0,T=1))
(trapz(u,fx)-target2)*10^5

target1=integrate(function(u) Re(exp(-1i*u*log(1.3/1.25))*((1-1i/u)/(1+u^2))*CFaa(u,kappa=2,theta=0.01,sigma=0.2,rho=-0.3,lambda=0,nu0=0.01,t=0,T=1)),1e-16,10+1e-16)$value

u = 1e-16+(0:60)/6
fx=Re(exp(-1i*u*log(1.3/1.25))*((1-1i/u)/(1+u^2))*CFaa(u,kappa=2,theta=0.01,sigma=0.2,rho=-0.3,lambda=0,nu0=0.01,t=0,T=1))
(trapz(u,fx)-target1)*10^5

u=c(1e-16,(1:75)/5,15+(0:45)/3,30+(0:70)*1)
fx=Re(exp(-1i*u*log(1.4/1.25))*((1-1i/u)/(1+u^2))*CFaa(u,kappa=2,theta=0.09,sigma=1.2,rho=-0.3,lambda=0,nu0=0.09,t=0,T=1/12))
(trapz(u,fx)-target)*10^5	

#0;inf fx dx
#-inf;inf f(logx) dx/x

target = integrate(function(u) Re(exp(-1i*exp(u)*log(1.4/1.25))*((1-1i/exp(u))/(1+exp(u)^2))*exp(u)*CFaa(exp(u),kappa=2,theta=0.09,sigma=1.2,rho=-0.3,lambda=0,nu0=0.09,t=0,T=1/12)),-36,6)$value

length(u)
u=linspacing(-17,5,0.4)
u=linspacing(-18.4206,5.8,0.4)
u=linspacing(-17,4.9,0.4)
fx=Re(exp(-1i*exp(u)*log(1.4/1.25))*((1-1i/exp(u))/(1+exp(u)^2))*exp(u)*CFaa(exp(u),kappa=2,theta=0.09,sigma=1.2,rho=-0.3,lambda=0,nu0=0.09,t=0,T=1/12))
(trapz(u,fx)-target)*10^5
trapz(u,fx)

u=exp(linspacing(-12,5,2/3))
fx=Re(exp(-1i*u*log(1.4/1.25))*((1-1i/u)/(1+u^2))*CFaa(u,kappa=2,theta=0.09,sigma=1.2,rho=-0.3,lambda=0,nu0=0.09,t=0,T=1/12))

exp(-0.01/12)*(1.25 - 1.4 * (1/2 + 1/pi * target ))
hestonAttari(c(0.09,1.2,-0.3,2,0.09), 1.25, 1.4, 0.01, 1/12)
hestonAttari2(c(0.09,1.2,-0.3,2,0.09), 1.25, 1.4, 0.01, 1/12)


1.25*exp(0.3^2*1/12/2 - 0.3*sqrt(1/12)*qnorm(0.10))
summary(vol)

system.time(for (k in 1:3000) hestonaa(1.25/exp(0.01),1.3,2,0.01,0.2,-0.3,0,0.01,0.01,0,0,1))
system.time(for (k in 1:3000) hestonaa2(1.25/exp(0.01),1.3,2,0.01,0.2,-0.3,0,0.01,0.01,0,0,1))

library(rbenchmark)
benchmark( hestonAttari(c(0.01,0.2,-0.3,2,0.01), 1.25, 1.3, 0.01, 1), hestonAttari2(c(0.01,0.2,-0.3,2,0.01), 1.25, 1.3, 0.01, 1), replications = 5000)
system.time(for (k in 1:4000) hestonAttari2(c(0.01,0.2,-0.3,2,0.01), 1.25, 1.3, 0.01, 1))
system.time(for (k in 1:4000) hestonAttari(c(0.01,0.2,-0.3,2,0.01), 1.25, 1.3, 0.01, 1))

szHeston(c(0.01,0.2,-0.3,2,0.01), F, K, rd, tau)-szCarr(c(0.01,0.2,-0.3,2,0.01), F, K, rd, tau, 3)
sz2d(S,K,0.01,0.2,-0.8,4,0.01,0.01,0.2,0.8,2,0.01,rd,0,0,tau)
sz2dHeston(c(0.01,0.2,-0.8,4,0.01,0.01,0.2,0.8,2,0.01), 1.25, 1.3, 0.01, 1)
sz2dCarr(c(0.01,0.2,-0.8,4,0.01,0.01,0.2,0.8,2,0.01), 1.25, 1.3, 0.01, 1,5)

CFaSZ2ddouble(10, F, c(0.01,0.2,-0.8,4,0.01,0.01,0.2,0.8,2,0.01), 0, tau)

cfSZd(10,S,rd,0,0,tau,0.01,0.2,-0.8,4,0.01,0.01,0.2,0.8,2,0.01)


### lets investigate integral convergence

u = (0:1000)/10+1e-3
u = (0:300)/3+1e-3
fx = Re(exp(-1i*u*log(1.5/1.25))*((1-1i/u)/(1+u^2))*exp(-1i*u*log(1.25))*CFa(u,S=1.25/exp(0.01),kappa=2,theta=0.01,sigma=0.2,rho=-0.3,lambda=0,nu0=0.01,rd=0.01,rf=0,t=0,T=1,type=2))
u1 = exp((-18:9)*2/3)
fx1 = Re(exp(-1i*u*log(1.5/1.25))*((1-1i/u)/(1+u^2))*exp(-1i*u*log(1.25))*CFa(u,S=1.25/exp(0.01),kappa=2,theta=0.01,sigma=0.2,rho=-0.3,lambda=0,nu0=0.01,rd=0.01,rf=0,t=0,T=1,type=2))

plot(u,(fx),t="l")
plot(u1,fx1,col=1)
lines(u1,(fx1),col=2)

plot(diff(-fx))
analyze = diff(-fx)/0.001
mean(analyze[1:50])
mean(analyze[50:100])
mean(analyze[100:200])

plot(ts(abs(analyze[0:100])))

#0-2:0,1
#2-4:0,2
#4-8:0,4
#8-100:1

plot(u,(fx),xlim=c(0,100),t="p")
plot(log(fx))
plot(1/(1+exp(-fx)))
lines(u,fx,col=2)

x = (0:1000)/10+1e-3
plot((1/x[-1]))





alpha < sqrt(theta^2/omega^4+2/omega^2/nu0) - theta/omega^2 - 1
alpha < 
sqrt(0.01^2/0.4^4+2/0.4^2/0.01) - 0.01/0.4^2 - 1

i=1300

sourceCpp("nlheston.cpp")
system.time(for (i in 1:10) nlOpter(c(0.01,0.2,-0.3,2,0.01),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,5,1,1,TRUE))

CFa(-(3+1)*1i,S,2,0.01,0.2,-0.3,0,0.01,rd,0,0,1,type=2)
S^(3+1)

hestonaa(F=1.25,K=1.3,kappa=1,theta=VIX[i,1]^2,sigma=1.2*VIX[i,1],rho=-0.6,lambda=0,nu0=VIX[i,1]^2,rd=rd[4],rf=rf[4],t=0,T=tau[4])

##########################
## Controled experiment ##
##########################

#normal conditions
#round(colMedians(perfectCal3alt),3)
nc_th = 0.015
nc_om = 0.30
nc_rh = -0.4
nc_ka = 1.6
nc_nu = 0.002
nc_rd = 0.0045
nc_rf = 0.0015
nc_S = 1.32
nc_F = nc_S*exp((nc_rd-nc_rf)*tau)
nc_vart = tau*(nc_th + (nc_nu-nc_th)*(1-exp(-nc_ka*tau))/nc_ka/tau) #* (1 - 0.5*nc_rh*nc_om*tau + nc_om^2/12*tau^2)
deltas = c(0.06,0.23,0.5,0.73,0.86)
#start proposition
#kappa=1,theta=VIX[i,1]^2,sigma=1.2*VIX[i,1],rho=-0.6,lambda=0,nu0=VIX[i,1]^2

nc_prices = matrix(NA,6,5)
nc_vol = matrix(NA,6,5)
nc_vega = matrix(NA,6,5)
nc_delta = matrix(NA,6,5)
nc_K = matrix(NA,6,5)
for(i in 1:6){#tenors
  nc_K[i,] = nc_F[i]*exp(nc_vart[i]/2 + sqrt(nc_vart[i])*qnorm(deltas/exp(-as.numeric(nc_rf)*tau[i])))
  for(j in 1:5){#smile
    nc_prices[i,j] = hestonAttari(c(nc_th,nc_om,nc_rh,nc_ka,nc_nu), nc_F[i], nc_K[i,j], nc_rd, tau[i])
    nc_vol[i,j] = impVol(nc_prices[i,j],0,tau[i],nc_S,nc_K[i,j],nc_rd,nc_rf,type=1)
    #hestonaa1(S=nc_S,K=nc_K[j],kappa=nc_ka,theta=nc_th,sigma=nc_om,rho=nc_rh,lambda=0,nu0=nc_nu,rd=nc_rd,rf=nc_rf,t=0,T=tau[i])
    d1 = (log(nc_F[i]/nc_K[i,j]) + (nc_vol[i,j]^2*tau[i]/2))/(sqrt(tau[i])*nc_vol[i,j])
    d2 = d1 - sqrt(tau[i])*nc_vol[i,j]
    nc_vega[i,j] = nc_K[i,j]*exp(-nc_rd*tau[i])*sqrt(tau[i])*dnorm(d2)
    nc_delta[i,j] = pnorm(d1)
  }
}

plot(ts(nc_vol[6,]))
nc_opt = nlOpter(c(0.095^2,2*0.05,-0.25,2,0.05^2),as.numeric(t(nc_prices/nc_vega)),as.numeric(t(nc_vega)),nc_F,as.numeric(t(nc_K)),rep(nc_rd,6),tau,5,1,1,TRUE)
round(nc_opt,3)

out = test_methods(nc_S,nc_K,nc_vol,nc_rd,nc_rf,nc_prices/nc_vega,nc_vega)
out2 = test_methods(nc_S,nc_K,nc_vol,nc_rd,nc_rf,nc_prices/nc_vega,nc_vega,2)
out[1:5]
out[6:10]
controltab = rbind(c(nc_th,nc_om,nc_rh,nc_ka,nc_nu,0),c(out[6:10],out2[[3]]),c(out[1:5],out2[[1]]))
rownames(controltab) = c("Benchmark","ICM","Durrleman")
colnames(controltab) = c("theta","omega","rho","kappa","nu0","RMSE")

cfaX.table <- xtable(controltab)
digits(cfaX.table) <- 4
write.textable(cfaX.table,"appendix.t1.txt")

######

t0=Sys.time(); perfectDurr = opterCv(paramsDurr,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=51,logfile="nl_t1.txt"); Sys.time()-t0;
t0=Sys.time(); perfectLmed = opterCv(paramsLmed,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,logfile="nl_t2.txt"); Sys.time()-t0;

perfectDurr = read.table("perfectDurr.txt",sep=",",fill=TRUE,header=T)[,-1]
perfectLmed = read.table("perfectLmed.txt",sep=",",fill=TRUE,header=T)[,-1]

convDurr=read.table("nl_outDurr.txt",sep="\t",fill=TRUE)
convLmed=read.table("nl_outLmed.txt",sep="\t",fill=TRUE)

#convDurr=t(matrix(rep(as.numeric(unlist(t(convDurr))),6),161,n+5))[-((n+1):(n+5)),]
#convLmed=t(matrix(rep(as.numeric(unlist(t(convLmed))),6),161,n+5))[-((n+1):(n+5)),]
naDurr = is.na(convDurr)*1
naLmed = is.na(convLmed)*1
stepsDurr = as.numeric(colSums(naDurr))
stepsLmed = c(as.numeric(colSums(naLmed)),rep(n,23))

#convDurr[,161] = as.numeric(rep(NA,n))
#convLmed[,161] = as.numeric(rep(NA,n))
convDurr=addatend(convDurr,perfectDurr[,6])
convLmed=addatend(convLmed,perfectLmed[,6])
convDurr=cbind(durr[[1]],convDurr)
convLmed=cbind(lab[[1]],convLmed)

pdf("../images/rmse1d_conv_log.pdf", width=mywidth, height=myheight)
x=(0:61)*10
x1=1:62
plot(x,log(colMeans(convDurr[,x1])), t="b", lty=1, pch=1,col=1, ylab=c("Natural logarithm of mean RMSE"), xlab="Number of iterations", cex=0.5, cex.lab=1.8, cex.sub=1.4, cex.axis=2) # ylim=c(0.0015,0.0111))ylim=c(-6.7,-4.5) #ylim=c(-6.7,-4.5)
points(x,log(colMeans(convLmed[,x1])), col=2, pch=2, cex=0.5);lines(x,log(colMeans(convLmed[,x1])),col=2,lty=2)
legend("topright", c("Durrleman","ICM"), pch=(1:2), lty=(1:2), col=1:2, bty="n", cex=2, lwd=1.4)
grid()
dev.off()

pdf("../images/rmse1d_conv_steps.pdf", width=mywidth, height=myheight)
x=(0:84)*10
x1=1:85
plot(x,stepsDurr[x1]/1333, t="b", lty=1, pch=1,col=1, ylab=c("Share of days with calibration finished below given number of iterations"), xlab="Number of iterations", cex=0.5, cex.lab=1.8, cex.sub=1.4, cex.axis=2) # ylim=c(0.0015,0.0111))ylim=c(-6.7,-4.5) #ylim=c(-6.7,-4.5)
points(x,stepsLmed[x1]/1333, col=2, pch=2, cex=0.5);lines(x,stepsLmed[x1]/1333,col=2,lty=2)
legend("topleft", c("Durrleman","ICM"), pch=(1:2), lty=(1:2), col=1:2, bty="n", cex=2, lwd=1.4)
grid()
dev.off()

pdf("../images/rmse1d_conv_steps_rel.pdf", width=mywidth, height=myheight)
x=(0:84)*10
x1=1:85
plot(x,(stepsLmed[x1]-stepsDurr[x1])/1333, t="b", lty=1, pch=1,col=1, ylab=c("Share of days with calibration finished below given number of iterations"), xlab="Number of iterations", cex=0.5, cex.lab=1.8, cex.sub=1.4, cex.axis=2) # ylim=c(0.0015,0.0111))ylim=c(-6.7,-4.5) #ylim=c(-6.7,-4.5)
legend("topleft", c("ICM-Durrleman"), pch=(1:2), lty=(1:2), col=1:2, bty="n", cex=2, lwd=1.4)
grid()
dev.off()



