datafolder = "C:/Users/mikolaj.labedzki/Downloads/marketdata"
scriptsfolder = "C:/Users/mikolaj.labedzki/Documents/Rscripts"
workfolder = "C:/Users/mikolaj.labedzki/Documents/dissertation_R"
library(xtable)#latex export
library(moments)
datatable=read.table(paste0(datafolder,"/FX/EURUSDall.txt"), header = TRUE, sep="\t")
n=dim(datatable)[1]
m=dim(datatable)[2]-1
d.datatable=datatable[-1,-1]-datatable[-n,-1]
#summary(datatable)
#summary(d.datatable)

summ.table=xtable(cbind(colMeans(d.datatable),apply(d.datatable,2,sd),apply(d.datatable,2,skewness),apply(d.datatable,2,kurtosis)))
colnames(summ.table)=c("Mean","Std. dev.","Skewness","Kurtosis")
digits(summ.table) <- 3
digits(summ.table)[2] <- 4
write.table(print(summ.table,floating=FALSE),paste0(workfolder,"/r3.t1.txt"))
#d.datatable=d.datatable[,-(6:10)]
smile1m=d.datatable[,1:5]
smile3m=d.datatable[,11:15]
smile6m=d.datatable[,16:20]
smile1y=d.datatable[,21:25]
smile2y=d.datatable[,26:30]
atmtstructure=d.datatable[,c(3,13,18,23,28)]
#remove 2M tenor
d.datatable = d.datatable[,-c(6:10)]

##################
## PCA
##################

#Test whole surf
pca1 = prcomp(d.datatable)
pca1.var=summary(pca1)
pca1.table <- xtable(pca1.var$importance[,1:8])
digits(pca1.table)[(1:9)] <- 3
write.table(print(pca1.table,floating=FALSE),paste0(workfolder,"/r3.t2.txt"))

#                     PC1     PC2     PC3     PC4     PC5     PC6    PC7    PC8     PC9
#Standard deviation     1.8845 0.42385 0.32593 0.24249 0.23079 0.17635 0.1674 0.1209 0.08953
#Proportion of Variance 0.8748 0.04426 0.02617 0.01448 0.01312 0.00766 0.0069 0.0036 0.00197
#Cumulative Proportion  0.8748 0.91909 0.94526 0.95975 0.97287 0.98053 0.9874 0.9910 0.99300
#                     PC1     PC2     PC3     PC4
#Standard deviation     0.9578 0.15603 0.07808 0.03112
#Proportion of Variance 0.9664 0.02565 0.00642 0.00102
#Cumulative Proportion  0.9664 0.99206 0.99848 0.99950
#                          PC1    PC2     PC3     PC4
#Standard deviation     0.6072 0.1650 0.09897 0.05004
#Proportion of Variance 0.9032 0.0667 0.02400 0.00613
#Cumulative Proportion  0.9032 0.9699 0.99387 1.00000

#biplot(pca1)
#explained=pca1$sdev/sum(pca1$sdev)
#surf4=sum(explained[1:4])
#loss=round((explained[-m]-explained[-1])*100,2)
#plot(loss)
#plot(explained[1:10])

#Test atmtstructure
pca3 = prcomp(atmtstructure)
pca3.var=summary(pca3)
pca3.table <- xtable(pca3.var$importance[,1:5])
digits(pca3.table)[(1:6)] <- 3
write.table(print(pca3.table,floating=FALSE),paste0(workfolder,"/r3.t3.txt"))

#explained=pca3$sdev/sum(pca3$sdev)
#term4=sum(explained[1:4])

#Test smiles
pca2 = prcomp(smile1m)
pca2a = prcomp(smile3m)
pca2b = prcomp(smile6m)
pca2c = prcomp(smile1y)
pca2d = prcomp(smile2y)
avgPCA=(summary(pca2)$importance+summary(pca2a)$importance+summary(pca2b)$importance+summary(pca2c)$importance+summary(pca2d)$importance)/5
pca2.table = xtable(summary(pca2)$importance);digits(pca2.table)[(1:6)] <- 3
pca2a.table = xtable(summary(pca2a)$importance);digits(pca2a.table)[(1:6)] <- 3
pca2b.table = xtable(summary(pca2b)$importance);digits(pca2b.table)[(1:6)] <- 3
pca2c.table = xtable(summary(pca2c)$importance);digits(pca2c.table)[(1:6)] <- 3
pca2d.table = xtable(summary(pca2d)$importance);digits(pca2d.table)[(1:6)] <- 3
avgPCA.table = xtable(avgPCA);digits(avgPCA.table)[(1:6)] <- 3
write.table(rbind(print(pca2.table),print(pca2a.table),print(pca2b.table),print(pca2c.table),print(pca2d.table),print(avgPCA.table)),paste0(workfolder,"/r3.t4.txt")

#explained=pca2$sdev/sum(pca2$sdev)
#smile4=sum(explained[1:3])

#plot(pca3.var$rotation[,1], ylim=c(-1,1), t="l")
#lines(pca3.var$rotation[,2], col="red")
#lines(pca3.var$rotation[,3], col="green")
#plot(pca2.var$rotation[,1], ylim=c(-1,1), t="l")
#lines(pca2.var$rotation[,2], col="red")
#lines(pca2.var$rotation[,3], col="green")

#pcized = t( t(pca1$rot) %*% ( t(d.datatable) - apply(d.datatable, 2, mean) ) )
#plot(ts(diffinv(pcized[,4])))
#plot(ts(datatable[,3]))

pca1rot=t(matrix(data=summary(pca1)$rotation[,1],5,5))
rownames(pca1rot)=c("1M","3M","6M","1Y","2Y")
colnames(pca1rot)=c("10P","25P","ATM","25C","10C")
pca1rot.table <- xtable(pca1rot);digits(pca1rot.table) <- 3

pca2rot=t(matrix(data=summary(pca1)$rotation[,2],5,5))
rownames(pca2rot)=c("1M","3M","6M","1Y","2Y")
colnames(pca2rot)=c("10P","25P","ATM","25C","10C")
pca2rot.table <- xtable(pca2rot);digits(pca2rot.table) <- 3

pca3rot=t(matrix(data=summary(pca1)$rotation[,3],5,5))
rownames(pca3rot)=c("1M","3M","6M","1Y","2Y")
colnames(pca3rot)=c("10P","25P","ATM","25C","10C")
pca3rot.table <- xtable(pca3rot);digits(pca3rot.table) <- 3

pca4rot=t(matrix(data=summary(pca1)$rotation[,4],5,5))
rownames(pca4rot)=c("1M","3M","6M","1Y","2Y")
colnames(pca4rot)=c("10P","25P","ATM","25C","10C")
pca4rot.table <- xtable(pca4rot);digits(pca4rot.table) <- 3

write.table(rbind(print(pca1rot.table),print(pca2rot.table),print(pca3rot.table),print(pca4rot.table)),paste0(workfolder,"/r3.t5.txt"))

##################
## CFA
##################
fit <- princomp(d.datatable, cor=TRUE)
fit <- princomp(d.datatable, cor=FALSE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

# Maximum Likelihood Factor Analysis
fit <- factanal(d.datatable,4,rotation="varimax", scores = "regression")
nx = colSums(fit$loadings^2)
p = nrow(fit$loadings)
cfa2 = rbind(nx,nx/p,cumsum(nx/p))
#rownames(cfa2)=c("Ładunek", "Udział w wariancji","Udział skumulowany")
rownames(cfa2)=c("SS loadings", "Proportion Var","Cumulative Var")
colnames(cfa2)=c("F1","F2","F3","F4")
cfa2.table <- xtable(cfa2);digits(cfa2.table) <- 3
write.table(print(cfa2.table,floating=FALSE),paste0(workfolder,"/r3.t6.txt"))

#print(fit, digits=2, cutoff=.3, sort=TRUE)
spectable=matrix(data=fit$uniquenesses,5,5)
rownames(spectable)=c("10P","25P","ATM","25C","10C")
colnames(spectable)=c("1M","3M","6M","1Y","2Y")
cfa1.table <- xtable(spectable)
digits(cfa1.table) <- 3
write.table(print(cfa1.table,floating=FALSE),paste0(workfolder,"/r3.t7.txt"))

## 3 factors, not used
fit <- factanal(d.datatable, 3, rotation="varimax", scores = "regression")
print(fit, digits=2, cutoff=.3, sort=TRUE)
cfa3=t(matrix(data=c(8.55,6.90,6.65,0.34,0.28,0.27,0.34,0.62,0.88),3,3))
rownames(cfa3)=c("Ładunek", "Udział w wariancji","Udział skumulowany")
colnames(cfa3)=c("F1","F2","F3")
cfa3.table <- xtable(cfa3)
digits(cfa3.table) <- 3
print(cfa3.table,floating=FALSE)

fit <- factanal(EUR2MO~EUR1MO,2,data=atmtstructure, rotation="varimax")

# plot factor 1 by factor 2 
load <- fit$loadings[,1:2] 
plot(load,type="n") # set up plot 
text(load,labels=names(atmtstructure),cex=.7) # add variable names

# Determine Number of Factors to Extract
library(nFactors)
ev <- eigen(cor(d.datatable)) # get eigenvalues
ap <- parallel(subject=nrow(d.datatable),var=ncol(d.datatable),rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

# PCA Variable Factor Map 
library(FactoMineR)
result <- PCA(d.datatable) # graphs generated automatically


##################
## Autoregression
##################
datatable=read.table(paste0(datafolder,"/FX/EURUSDall.txt"), header = TRUE, sep="\t");pair="EURUSD"
source(paste0(scriptsfolder,"aff_pricing.R"))
n=dim(datatable)[1]
P1=matrix(NA,n,6);VIX=P1;P2=P1;P3=P1;P4=P1;P5=P1;P6=P1;SS=P1;CC=P1;mu=P1;Q1=P1;Q2=P1;Q3=P1;Q4=P1;slope=P1;curv=P1;skew=P1;kurt=P1
data.load(datatable,pair)

v=VIX[,1]
dv=c(NA,diff(v));lv=lagger(v);invlv=1/lv
h=v^2;dh=c(NA,diff(h));lh=lagger(h)

lm1=lm(dh~lh);lm1.out=summary(lm1)
lm2=lm(dh~lh+lv);lm2.out=summary(lm2) 
lm3=lm(dv~lv+invlv-1);lm3.out=summary(lm3) 
lm4=lm(dv~lv);lm4.out=summary(lm4)#== m4=lm(v~lv)

library(texreg)
#screenreg(list(lm1, lm2, lm3, lm4), digits=4)
tab9 = texreg(list(lm1, lm2, lm3, lm4), digits=4, dcolumn = FALSE, booktabs = FALSE,use.packages = FALSE, label = "tab:9", caption = "Estimation results for models with alternative functional forms.",float.pos = "h!")
fstat = round(c(lm1.out$fstatistic[1],lm2.out$fstatistic[1],lm3.out$fstatistic[1],lm4.out$fstatistic[1]),4)
F_stats = paste("F Statistics",fstat[1],fstat[2],fstat[3],fstat[4],sep=" & ")
write.table(c(tab9,F_stats),paste0(workfolder,"/r3.t9.txt"))

library(stargazer)
tab9 = stargazer(lm1, lm2, lm3, lm4, title="Results", align=TRUE)

err1=dh-c(NA,predict(lm1));
err2=dh-c(NA,predict(lm2));
err3=dv-c(NA,predict(lm3));
err4=dv-c(NA,predict(lm4));
lm1e=lm(log(err1)~log(lh));summary(lm1e)
lm2e=lm(log(err2)~log(lh));summary(lm2e)
lm4e=lm(log(err4)~log(lv));summary(lm4e)
shapiro.test(err1/v)
shapiro.test(err2/v)
shapiro.test(err3)
shapiro.test(err4)

##################
## Logit
##################
tau=c(1/12,2/12,3/12,6/12,1,2)
change1M=datatable$EURATM3M[-(1:21)]-datatable$EURATM3M[-((n-20):n)]
change3M=datatable$EURATM1M[-(1:63)]-datatable$EURATM1M[-((n-62):n)]
change6M=datatable$EURATM1M[-(1:126)]-datatable$EURATM1M[-((n-125):n)]
change1Y=datatable$EURATM3M[-(1:252)]-datatable$EURATM3M[-((n-251):n)]
#plot(change3M,t="l")
#plot(change1Y,t="l")
change1M.s=as.numeric(change1M>0)
change3M.s=as.numeric(change3M>0)
change6M.s=as.numeric(change6M>0)
change3M.s=as.numeric(change3M>0)
datatable2=cbind(datatable[-((n-62):n),])/100
n2 = dim(datatable2)[1]
level=datatable2$EURATM1M
slope=datatable2$EURATM2Y-datatable2$EURATM1M
curv=datatable2$EURATM2Y+datatable2$EURATM1M-2*datatable2$EURATM6M
skew1M=datatable2$EUR25P1M-datatable2$EUR25C1M
skew6M=datatable2$EUR25P6M-datatable2$EUR25C6M
RR1M=datatable2$EUR25P1M-datatable2$EUR25C1M
RR3M=datatable2$EUR25P3M-datatable2$EUR25C3M
RR1Y=datatable2$EUR25P1Y-datatable2$EUR25C1Y
RR2Y=datatable2$EUR25P2Y-datatable2$EUR25C2Y
BF1M=datatable2$EUR25P1M+datatable2$EUR25C1M-2*datatable2$EURATM1M
BF3M=datatable2$EUR25P3M+datatable2$EUR25C3M-2*datatable2$EURATM3M
BF1Y=datatable2$EUR25P1Y+datatable2$EUR25C1Y-2*datatable2$EURATM1Y
BF2Y=datatable2$EUR25P2Y+datatable2$EUR25C2Y-2*datatable2$EURATM2Y

termstruct=cbind(datatable2$EURATM1M,datatable2$EURATM2M,datatable2$EURATM3M,datatable2$EURATM6M,datatable2$EURATM1Y,datatable2$EURATM2Y)
#P1=matrix(NA,n,6);VIX=P1;P2=P1;P3=P1;P4=P1;P5=P1;P6=P1;SS=P1;CC=P1;mu=P1;Q1=P1;Q2=P1;Q3=P1;Q4=P1;slope=P1;curv=P1;skew=P1;kurt=P1
#data.load(datatable2*100,pair)

library(Rcpp)
sourceCpp(paste0(scriptsfolder,"/nlheston.cpp"))

hestonTSpar=hestonTSfitC(termstruct^2,tau,2,6,2)
hestonTSpar=hestonTSfitC(VIXadded^2,tau,1,6,2)
hestonTSpar0=hestonTSfitC(VIX^2,tau,2,6,2)
shzTSpar=hestonTSfitC(termstruct,tau,2,6,1.01)

hestonTSpar=varTSfitAltC(termstruct^2,1,tau,a=2,b=6,kappa=2,barrier=1000)
shzTSpar=varTSfitAltC(termstruct,1,tau,2,6,1.01,barrier=1000)

summary(hestonTSpar)
summary(perfectCal[,1])
i=1321
which(perfectCal[,1]>400)
plot2(hestonTSpar[,2],hestonTSpar0[,2])

kappa2=(1-exp(-hestonTSpar[,3]))/hestonTSpar[,3]
kappa2shz=(1-exp(-shzTSpar[,3]))/shzTSpar[,3]

shzTSpar2=shzTSpar^2

m0=glm(change3M.s~level+slope+curv+BF1M+skew1M)
#summary(m0)

m1=glm(change3M.s~level+slope+curv,family=binomial());m1.out=summary(m1)
m2=glm(change3M.s~level+slope+curv+BF1M,family=binomial());m2.out=summary(m2)
m3=glm(change3M.s~(hestonTSpar[,1])+(hestonTSpar[,2])+kappa2,family=binomial());m3.out=summary(m3)
m4=glm(change3M.s~shzTSpar2[,1]+shzTSpar2[,2]+kappa2shz,family=binomial());m4.out=summary(m4)
m5=glm(change3M.s~(hestonTSpar[,1])+(hestonTSpar[,2])+kappa2+shzTSpar2[,1]+shzTSpar2[,2]+kappa2shz,family=binomial());m5.out=summary(m5)
#m6=glm(change3M.s~(hestonTSpar[,1])+(hestonTSpar[,2])+kappa2+shzTSpar2[,1]+shzTSpar2[,2]+kappa2shz+level+slope+curv+BF1M,family=binomial());m6.out=summary(m6)
m6a=glm(change3M.s~(hestonTSpar[,1])+(hestonTSpar[,2])+kappa2+level+slope+curv+BF1M,family=binomial());m6.out=summary(m6a)
m6b=glm(change3M.s~shzTSpar2[,1]+shzTSpar2[,2]+kappa2shz+level+slope+curv+BF1M,family=binomial());m6.out=summary(m6b)

m6s=glm(change3M.s~sqrt(hestonTSpar[,1])+sqrt(hestonTSpar[,2])+kappa2+shzTSpar[,1]+shzTSpar[,2]+kappa2shz+level+slope+curv+BF1M,family=binomial());m6s.out=summary(m6s)

#library(ROCR)
test=as.data.frame(cbind(change3M.s,(hestonTSpar[,1]),(hestonTSpar[,2]),kappa2,shzTSpar2[,1]+shzTSpar2[,2]+kappa2shz,level,slope,curv,BF1M,RR3M,RR1Y,BF3M,BF1Y))
#tests=as.data.frame(cbind(change3M.s,sqrt(hestonTSpar[,1]),sqrt(hestonTSpar[,2]),kappa2,shzTSpar[,1]+shzTSpar[,2]+kappa2shz,level,slope,curv,BF1M,RR3M,RR1Y,BF3M,BF1Y))
score1<-predict(m1,type='response',test)
score2<-predict(m2,type='response',test)
score3<-predict(m3,type='response',test)
score4<-predict(m4,type='response',test)
score5<-predict(m5,type='response',test)
#score6<-predict(m6,type='response',test)
score6a<-predict(m6a,type='response',test)
score6b<-predict(m6b,type='response',test)
pred1<-prediction(score1,change3M.s)
pred2<-prediction(score2,change3M.s)
pred3<-prediction(score3,change3M.s)
pred4<-prediction(score4,change3M.s)
pred5<-prediction(score5,change3M.s)
#pred6<-prediction(score6,change3M.s)
pred6a<-prediction(score6a,change3M.s)
pred6b<-prediction(score6b,change3M.s)
perf1 <- performance(pred1,"tpr","fpr")
perf2 <- performance(pred2,"tpr","fpr")
perf3 <- performance(pred3,"tpr","fpr")
perf4 <- performance(pred4,"tpr","fpr")
perf5 <- performance(pred5,"tpr","fpr")
#perf6 <- performance(pred6,"tpr","fpr")
perf6a <- performance(pred6a,"tpr","fpr")
perf6b <- performance(pred6b,"tpr","fpr")

pdf("../images/logit.pdf", width=12, height=9)
plot(perf1, lty=1, cex.lab=1.2, cex.sub=1.2, lwd=2)
lines((1:n2)/n2,(1:n2)/n2, col="gray", lty=2, lwd=2) 
lines(perf2@x.values[[1]],perf2@y.values[[1]], col=2, lty=2, lwd=2)
lines(perf3@x.values[[1]],perf3@y.values[[1]], col=3, lty=3, lwd=2)
lines(perf4@x.values[[1]],perf4@y.values[[1]], col=4, lty=4, lwd=2)
lines(perf5@x.values[[1]],perf5@y.values[[1]], col=5, lty=5, lwd=2)
lines(perf6a@x.values[[1]],perf6a@y.values[[1]], col=6, lty=6, lwd=2)
lines(perf6b@x.values[[1]],perf6b@y.values[[1]], col="darkorange", lty=7, lwd=2)
legend("bottomright", c("Model 1","Model 2","Model 3","Model 4","Model 5","Model 6","Model 7"), lty=1:7, col=c(1:6,"darkorange"), bty="n", cex=1.5, lwd=1.4)
grid()
dev.off()

png(paste("logit.png"), width = 640, height = 480)
plot(perf1, lty=1, cex.lab=1.2, cex.sub=1.2, lwd=2)
lines((1:n2)/n2,(1:n2)/n2, col="gray", lty=2, lwd=2) 
lines(perf2@x.values[[1]],perf2@y.values[[1]], col=2, lty=2, lwd=2)
lines(perf3@x.values[[1]],perf3@y.values[[1]], col=3, lty=3, lwd=2)
lines(perf4@x.values[[1]],perf4@y.values[[1]], col=4, lty=4, lwd=2)
lines(perf5@x.values[[1]],perf5@y.values[[1]], col=5, lty=5, lwd=2)
lines(perf6a@x.values[[1]],perf6a@y.values[[1]], col=6, lty=6, lwd=2)
lines(perf6b@x.values[[1]],perf6b@y.values[[1]], col="darkorange", lty=7, lwd=2)
legend("bottomright", c("Model 1","Model 2","Model 3","Model 4","Model 5","Model 6","Model 7"), lty=1:7, col=c(1:6,"darkorange"), bty="n", cex=1.5, lwd=2)
dev.off()

#library(texreg)
screenreg(list(m1, m2, m3, m4, m5, m6a, m6b))
tab10 = texreg(list(m1, m2, m3, m4, m5, m6a, m6b), digits = 1, dcolumn = FALSE, booktabs = FALSE,use.packages = FALSE, label = "tab:10", caption = "Models explaining the probability of increases and decreases of implied volatility for the EURUSD 1M option in 1 quarter horizon",float.pos = "h!")
tab10 = gsub("hestonTSpar\\[, 1]", "$\\\\nu_{0,H}$", tab10)
tab10 = gsub("shzTSpar2\\[, 1]", "$\\\\nu_{0,S}^2$", tab10)
tab10 = gsub("hestonTSpar\\[, 2]", "$\\\\theta_H$", tab10)
tab10 = gsub("shzTSpar2\\[, 2]", "$\\\\theta_S^2$", tab10)
tab10 = gsub("kappa2shz", "$\\\\xi_S$", tab10)
tab10 = gsub("kappa2", "$\\\\xi_H$", tab10)
auroc = round(unlist(c(performance(pred1,"auc")@y.values,performance(pred2,"auc")@y.values,performance(pred3,"auc")@y.values,performance(pred4,"auc")@y.values,performance(pred5,"auc")@y.values,performance(pred6a,"auc")@y.values,performance(pred6b,"auc")@y.values)),2)
AUROC_stats = paste("AUROC",auroc[1],auroc[2],auroc[3],auroc[4],auroc[5],auroc[6],auroc[7],sep=" & ")
write.table(c(tab10,AUROC_stats),"r3.t10.txt")

performance(pred6a,"auc")@y.values
performance(pred6b,"auc")@y.values

library(stargazer)
stargazer(m1, m2, m3, m4, m5, m6a, m6b, title="Results", align=TRUE)
