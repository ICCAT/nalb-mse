## Script to produce the figures of the paper (SCRS-OEM)
## Adapted by Gorka Merino (AZTI)  gmerino@azti.es

library(FLCore)
library(mpb)
library(plyr)
library(reshape)
library(FLBRP)

# Directories, scenarios and filenames:

theme_set(theme_bw(10))

dirMy ="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT"
dirDat=paste(dirMy,"/RObjects",sep="")
Dir="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/MFCL scenarios 240/"
scenarios=basename(list.dirs(Dir, full.names=FALSE, recursive=FALSE))

load(paste(dirDat,"om240.RData",sep="/"))

setwd(dirDat)
rels=as.numeric(stock(om[[201]])/refpts(eql[[201]])["msy", "biomass"])
write.csv(rels, file="OM_stock.csv")

biol=ldply(eql, function(x) as.data.frame(x[["catch.sel","m","mat","stock.wt"]]))
biol=transform(biol,
               Quantity=factor(qname,labels=c("Selectivity","M","Maturity","Mass")),
               .id     =names(om)[X1])

# Figure 1: Operating Model
source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/RScripts/Multiplot.R")


x11()    #### "SCRS2017_091_Figure1.tiff", res=300, unit="in", 6,7)
par(mfrow=c(3,1), mar=c(4,4,1,1))
a=ggplot()+
  geom_line(aes(1930:2014, as.numeric(rec(om[[201]]))/1000000), lwd=2)+
  ylab("Recruitment (Millions)")+xlab("")
b=ggplot()+
  geom_line(aes(1930:2014, as.numeric(stock(om[[201]]))/1000000), lwd=2)+
  ylab("Stock biomass (k tonnes)")+xlab("")
c=ggplot()+
  geom_line(aes(1930:2014, as.numeric(catch(om[[201]]))/1000000), lwd=2)+
  ylab("Catch (k tonnes)")+xlab("")
d=ggplot()+
  geom_line(aes(1930:2014, as.numeric(fbar(om[[201]]))), lwd=2)+
  ylab("fbar")+xlab("Years")

multiplot(a,b,c,d,cols=1)

## 


## Figure 2) Overall selectivity pattern
names(eql)=names(om)
biol=ldply(eql, function(x) as.data.frame(x[["catch.sel","m","mat","stock.wt"]]))
biol=transform(biol,
               Quantity=factor(qname,labels=c("Selectivity","M","Maturity","Mass")),
               .id     =names(om)[X1])

plotbiol=ggplot(biol)+
  geom_line(aes(age,data,group=.id))+ 
  facet_wrap(~Quantity,scale="free_y")+ 
  theme(legend.position="none")+
  scale_colour_discrete(name="Scenario")+
  xlab("Age") +ylab("")

bioldat=plotbiol$data

x=bioldat$.id
bioldat$M=substr(x, unlist(gregexpr("M", x)), unlist(gregexpr("_h", x))-1)
bioldat$h=substr(x, unlist(gregexpr("h", x)), unlist(gregexpr("_q", x))-1)
bioldat$dynq=substr(x, unlist(gregexpr("q", x)), length(x))
bioldat$origin=substr(x, 1, 4)
bioldat$Mage="fixed"
bioldat$Mage[5761:7200]="age"

x11()    #### "SCRS2017_091_Figure2.tiff", 4, 3, unit="in", res=300)

ggplot(subset(bioldat,  Quantity=="Selectivity" & h=="h6" & origin=="Base" & M=="M03"))+
  geom_line(aes(age, data/max(data), col=M, group=.id, lty=Mage))+
  ylab("")+xlab("Age")+facet_wrap(~Quantity, scale="free")+
  theme(legend.position="none")+
  guides(col=guide_legend(title=""), lty=guide_legend(title=""))+
  scale_linetype_manual(values=c("solid")) # Change linetypes

## 
####

# print(names(om)[i], i)
## 
## Generate CPUEs and catch for the MP
## saa.RData contains the outputs from MFCL for fisheries' selectivity
load("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/RObjects/saa.RData")
# Source the OEM
source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/RScripts/OEM.R")
# 1) cpue1=f(stock)
# 2) cpue2=f(overall cpue)
# 3) cpue3=f(overall age corrected cpue)
# 4) cpue4=f(5 cpue)
# 5) cpue5=f(5 cpue with CT 10%)

# 1) cpue1=f(stock)
i=201
cpue1=oem(window(om[["Base_M03_h6_q0"]],start=1930), 
          fish.dependent=FALSE)

# 2) cpue2=f(overall cpue)
i=201
cpue2=oem(window(om[["Base_M03_h6_q0"]],start=1930), 
          fish.dependent=TRUE,
          sel=FLQuant(FLQuant(1,dimnames=dimnames(harvest(window(om[["Base_M03_h6_q0"]],start=1930))))))

# 3) cpue3=f(overall age corrected cpue)
i=201
cpue3=oem(window(om[["Base_M03_h6_q0"]],start=1930), 
          fish.dependent=TRUE,
          sel=FLQuant(subset(biol, .id==names(om)[i] & Quantity=="Selectivity")$data*rnorm(dim(stock.n(window(om[["Base_M03_h6_q0"]],start=1930)))[1], 1, 0.15), 
                      dimnames=dimnames(harvest(window(om[["Base_M03_h6_q0"]],start=1930)))))
# 4) cpue4=f(5 cpue)
i=201

cpue4=FLQuants(llply(saa[c(1,7,10,11)], function(x) 
  mpb:::oem(window(om[["Base_M03_h6_q0"]],start=1981),sel=FLQuant(c(x),dimnames=dimnames(stock.n(window(om[["Base_M03_h6_q0"]],start=1981)))))))

# 5) cpue5=f(5 cpue with CT 10%)
i=201
cpue5=FLQuants(llply(saa[c(1,7,10,11)], function(x) 
  mpb:::oem(window(om[["Base_M03_h6_q0"]],start=1981),bias=FLPar(q=0.10),sel=FLQuant(c(x),dimnames=dimnames(stock.n(window(om[["Base_M03_h6_q0"]],start=1981)))))))

cpue5[c(1,2,4)]=cpue4[c(1,2,4)]



x11()    #### "SCRS2017_091_Figure3.tiff", 6,4, unit="in", res=300)

ggplot()+
  geom_line(aes(1930:2014, as.numeric(cpue1)/max(as.numeric(cpue1))), lwd=1, col="black")+
  geom_line(aes(1930:2014, as.numeric(cpue2)/max(as.numeric(cpue2))), lwd=1, col="red")+
  geom_line(aes(1930:2014, as.numeric(cpue3)/max(as.numeric(cpue3))), lwd=1, col="blue")+
  xlab("Years")+ylab("Norm Index")+ylim(0,1)+
  annotate("text", label = "Index 1~Stock", x = 1994, y = 0.96, size = 4, colour = "black")+
  annotate("text", label = "Index 2~Catch/fbar", x = 1996, y = 0.88, size = 4, colour = "blue")+
  annotate("text", label = "Index 3~Sel x Catch/fbar", x = 1998, y = 0.80, size = 4, colour = "red")

## 

## Generate CPUE with all cpues
cpues=data.frame(year=1930:2014, model=rep("Stock", length(1930:2014)), data=as.numeric(cpue1)/as.numeric(max(cpue1)), fleet=rep("1 Overall Stock", length(1930:2014)))
cpues=rbind(cpues, data.frame(year=1930:2014, model=rep("CPUE", length(1930:2014)), data=as.numeric(cpue2)/as.numeric(max(cpue2)), fleet=rep("2 Overall CPUE", length(1930:2014))))
cpues=rbind(cpues, data.frame(year=1930:2014, model=rep("CPUE-Sel", length(1930:2014)), data=as.numeric(cpue3)/as.numeric(max(cpue3)), fleet=rep("3 Overall Age CPUE", length(1930:2014))))
cpues=rbind(cpues, data.frame(year=1981:2014, model=rep("ind CPUE-Sel", length(1981:2014)), data=as.numeric(cpue4[[1]])/max(as.numeric(cpue4[[1]])), fleet=rep("4A Spain BB", length(1981:2014))))
cpues=rbind(cpues, data.frame(year=1981:2014, model=rep("ind CPUE-Sel", length(1981:2014)), data=as.numeric(cpue4[[2]])/max(as.numeric(cpue4[[2]])), fleet=rep("4B Japan LL byc", length(1981:2014))))
cpues=rbind(cpues, data.frame(year=1981:2014, model=rep("ind CPUE-Sel", length(1981:2014)), data=as.numeric(cpue4[[4]])/max(as.numeric(cpue4[[4]])), fleet=rep("4C KoreaPaCu LL4", length(1981:2014))))
cpues=rbind(cpues, data.frame(year=1981:2014, model=rep("ind CPUE-Sel", length(1981:2014)), data=as.numeric(cpue4[[3]])/max(as.numeric(cpue4[[3]])), fleet=rep("4D China Taipei LL late", length(1981:2014))))
cpues=rbind(cpues, data.frame(year=1981:2014, model=rep("ind CPUE-Sel q10%", length(1981:2014)), data=as.numeric(cpue5[[3]])/max(as.numeric(cpue5[[3]])), fleet=rep("5D China Taipei 10%", length(1981:2014))))


## Save CPUEs:
setwd(dirDat)
save(cpues, file="SimulatedCPUEs.RData")
write.csv(cpues, file="SimulatedCPUEs.csv")

## Now restart and run mpb to do stock assessments


#### --------------------------

library(ggplot2)
library(reshape)
library(dplyr)
library(plyr)
library(FLCore)
library(FLBRP)
library(ggplotFL)
library(diags)
library(mpb)

# Source additional functions
source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/RScripts/Multiplot.R")
source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/RScripts/OEM.R")
source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/RScripts/diags-MFCL.R")

# 1st. Build the operating models from MFCL scenarios.

# Directories, scenarios and filenames:

theme_set(theme_bw(10))

dirMy ="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/"
dirDat=paste(dirMy,"/RObjects",sep="")
### paste(dirMy,"/Figures/check091",sep="")
Dir="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/MFCL scenarios 240/"
scenarios=basename(list.dirs(Dir, full.names=FALSE, recursive=FALSE))

rodFn=FLife:::rodFn



flt=c("ESP_BBrec","EsFr_TR","EsFr_BBear","PRT_BB",
      "JPN_LLtrg","JPN_LLtra","JPN_LLbyc",
      "TAI_LL1","TAI_LL2","TAI_LL3",
      "KrPaCu_LL","Other_SU")

nms=c("Spain BB","Spain France\nTroll","Spain France\nEarly BB",
      "Portugal\nBB","Japan\nLL Target","Japan\nLL tra","Japan\nLL Bycatch",
      "Chinese-Taipei\nLL Early","Chinese-Taipei\nLL Mid","Chinese-Taipei\nLL Late",
      "KoreaPaCu LL","Other surface")

nm2=paste("flt",1:12,sep="")

u=diags.mfcl("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/MFCL scenarios 240/Base_M03_h6_q0/plot-07.par.rep")

u=transform(u,fleet=factor(name,labels=nms,levels=c("flt1","flt2", "flt3", "flt4",
                                                    "flt5","flt6", "flt7", "flt8",
                                                    "flt9","flt10","flt11","flt12")))

u=transform(subset(u,!is.na(residual)&residual!=0&name%in%paste("flt",c(1,7,10,11),sep="")),
            Month=factor(month))


x11()    #### "SCRS2017_091_Figure4.tiff", 7,6, unit="in", res=300)

ggplot(ddply(u,.(fleet,month),transform,obs=obs/mean(obs)))+
  geom_point(aes(year,obs, col=Month, fill=fleet),  pch=21)+
  xlim(1980,2013)+
  geom_smooth(aes(year,obs),se=FALSE)+
  facet_wrap(~fleet, ncol=2, scale="free")+
  xlab("Year")+ylab("CPUEs used to fit MFCL runs in 2013")+
  theme_bw()+theme(legend.position="none")

## 

x11()    #### "SCRS2017_091_Figure6.tiff", 7,6, unit="in", res=300)
dgs=ddply(u,.(Month,fleet), with, data.frame(dev=diags:::stdz(effDev.value),year=year))
ggplot(dgs)+
  geom_abline(aes(slope=0,intercept=0))+
  geom_point(aes(year,dev,fill=fleet,col=Month),pch=21)+
  # geom_smooth(aes(year,dev),se=FALSE)+
  stat_smooth(aes(year,dev, col=fleet),method="lm", fill="salmon",pch=21)+
  facet_wrap(~fleet, ncol=2, scale="free")+
  xlab("Year")+ylab("Residuals of fit from MFCL runs in 2013")+
  theme_bw()+theme(legend.position="none")
## 

# Figure 

x11()    #### "SCRS2017_091_Figure8.tiff", 6,6, unit="in", res=300)
par(mar=c(4,4,1,2))
plot(1:15, as.numeric(saa[[1]]), type="l", xlab="Age class", ylab="Selectivity", ylim=c(0,0.8))
grid()
lines(1:15, as.numeric(saa[[7]]), type="l", col=2)
lines(1:12, as.numeric(saa[[10]]), type="l", col="darkgreen")
lines(1:15, as.numeric(saa[[11]]), type="l", col=4)
legend(8,0.8, c("Spain BB", "Japan LL byc", "China Taipei LL", "KoreaPaCu LL"), col=c(1,2,"darkgreen",4), lty=1)

## 


## Now use the input data of 2016 Madeira
# Source an additional function from 'kobe'

source('c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/Rscripts/biodyn-kobe.R')

options(digits=4)
theme_set(theme_bw(10))

# -----------------------------------------------------------------------
# 1) Upload and check input CPUE series
# -----------------------------------------------------------------------

v=read.csv(paste(dirDat,"data.csv", sep="/" ), sep=";")
names(v)=c("name", "year", "index", "cv", "ref")

v=subset(v, name!="sp_troll")
v=subset(v, name!="catch")
v=subset(v, name!="jpll_byctcL")

v  =v[!duplicated(v[,c("name","year")]),c("name","year","index")]

## Prettify names
nms=c("China Taipei late","Japan bycatch","Esp Baitboat",
      "US cont", 
      "Ven LL")
names(nms)=sort(unique(v$name))

nm2=c("Chinese Taipei late Longline","Japan bycatch Longline",
      "Spanish Baitboat",
      "US continuity Longline", 
      "Venezuela Longline")
names(nm2)=sort(unique(v$name))

v  =transform(v,name=factor(name,levels=sort(unique(v$name)),labels=nm2))

cpue=FLQuants(dlply(v,.(name),with, 
                    as.FLQuant(data.frame(year=year,data=index))))
names(cpue)=nm2

# Look at the CPUE's -------------------------------------



x11()    #### "SCRS2017_091_Figure5.tiff", 7,6, unit="in", res=300)

mpb:::plotIndex(cpue)+ylab("CPUE index")+facet_wrap(~qname, ncol=3, scale="free")

## 

load(paste(dirDat, "TFGO.RData", sep="/"))

res=ldply(bds2, function(x) x@diags)



x11()    #### "SCRS2017_091_Figure7.tiff", 6,6, unit="in", res=300)
hey=subset(res, .id==names(bds2)[1])
levels(hey$name)=c("China Taipei LL", "Spain BB", "Japan bycatch LL", "US cont", "Venezuela LL")
hey=na.omit(hey)
ggplot(hey)+
  geom_hline(aes(yintercept=0),col="black")+
  geom_point(aes(year,residual),fill="salmon",col="grey",pch=21)+
  stat_smooth(aes(year,residual, col=name),method="lm", fill="salmon",pch=21)+
  facet_wrap(~name,ncol=2, scale="free")+theme_bw()+
  xlab("Year")+ylab("Residuals of fit from biodyn base case in 2016")+
  theme(legend.position="none")
## 

##### Ahora realizar ajustes sobre el OM con el MP usando diferentes CPUEs.

catch=as.FLQuant(read.csv(paste(dirDat, "/catch.csv",sep=""), sep=";"))

cpue =FLQuants(dlply(read.csv(paste(dirDat, "/SimulatedCPUEs.csv", sep="")),
                     .(fleet),with,as.FLQuant(data.frame(year=year,data=data))))


cpue=window(cpue, start=1981)

cpue[4]=window(cpue[4], start=1981)
cpue[5]=window(cpue[5], start=1988, end=2012)
cpue[6]=window(cpue[6], start=1987)
cpue[7]=window(cpue[7], start=1999)
cpue[8]=window(cpue[8], start=1999)



x11()    #### "SCRS2017_091_Figure9.tiff", 7,6, unit="in", res=300)
mpb:::plotIndex(cpue[4:8])+ylab("Simulated CPUE index")+facet_wrap(~qname, ncol=3, scale="free")+xlim(1981, 2014)
## 

## Now do the stock assessment with the Biomass Dynamic model

bd=biodyn("pellat", params=FLPar(r=0.3,k=5.5e5,b0=1,p=0.001),
          catch = catch)

# Create the list of objects that inlcudes the Base CAse and Sensitivities

# Specify controls and do fits for each scenario 
set.seed(500)
bds=mpb:::biodyns(list("Stock"=bd,
                       "CPUE"=bd,
                       "CPUE~Sel"=bd,
                       "Fleet"=bd,
                       "Fleet China-Tai 10%"=bd,
                       "China-Tai 10% only"=bd))

params(bds[[1]])=FLPar(r=.1,k=3.2e6,p=0.001,b0=1) 
setParams(bds[[1]])=cpue[1]
setControl(bds[[1]])=params(bds[[1]])
bds[[1]]=fit(bds[[1]],cpue[1])

params(bds[[2]])=FLPar(r=.1,k=3.2e6,p=0.001,b0=1) 
setParams(bds[[2]])=cpue[2]
setControl(bds[[2]])=params(bds[[2]])
bds[[2]]=fit(bds[[2]],cpue[2])

params(bds[[3]])=FLPar(r=.1,k=3.2e6,p=0.001,b0=1) 
setParams(bds[[3]])=cpue[3]
setControl(bds[[3]])=params(bds[[3]])
bds[[3]]=fit(bds[[3]],cpue[3])

params(bds[[4]])=FLPar(r=0.1,k=3.2e6, p=0.001,b0=1) 
setParams(bds[[4]])=cpue[4:7]
setControl(bds[[4]])=params(bds[[4]])
bds[[4]]=fit(bds[[4]],cpue[4:7])

params(bds[[5]])=FLPar(r=.1,k=3.2e6,p=0.001,b0=1) 
setParams(bds[[5]])=cpue[c(4,5,6,8)]
setControl(bds[[5]])=params(bds[[4]])
bds[[5]]=fit(bds[[5]],cpue[c(4,5,6,8)])

params(bds[[6]])=FLPar(r=.1,k=3.2e6,p=0.001,b0=1) 
setParams(bds[[6]])=cpue[8]
setControl(bds[[6]])=params(bds[[4]])
bds[[6]]=fit(bds[[6]],cpue[8])

## Figures of Stock Status

omstock=read.csv(paste(dirDat, "OM_stock.csv", sep="/"))

theme_set(theme_bw(10))


x11()    #### "SCRS2017_091_Figure10.tiff", 5,5, unit="in", res=300)

MadeStock=read.csv("c:/use/Gorka/1_AZTI_tuna/ICCAT/Species Groups/2016/ALBACORE/Haritz deterministic Madeira OK.csv")

par(mar=c(4,4,1,1))
plot(1930:2014, omstock[,2], 
     ylim=c(0, 3), type="l", ylab="MP estimated B/Bmsy", lwd=2,xlab="Years")
grid()
lines(1930:2015, as.numeric(stock(bds[[1]])/bmsy(bds[[1]])),col=2, lwd=1)
lines(1930:2015, as.numeric(stock(bds[[2]])/bmsy(bds[[2]])), col="darkgreen", lwd=1)
lines(1930:2015, as.numeric(stock(bds[[3]])/bmsy(bds[[3]])), col=4, lwd=1)
lines(1930:2015, as.numeric(stock(bds[[4]])/bmsy(bds[[4]])), col="violet", lwd=1)
lines(1930:2015, as.numeric(stock(bds[[5]])/bmsy(bds[[5]])), col="green", lwd=1)
#lines(1930:2015, as.numeric(stock(bds[[6]])/bmsy(bds[[6]])), col="orange", lwd=1)

lines(1930:2014, omstock[,2], lwd=4)
lines(1930:2015, MadeStock$stock, col="gray", lwd=4)

abline(h=1, lty=2)

legend(1968, 3, c("OM", "2016 SA", "Index 1~Stock", "Index 2~Catch/fbar",
                  "Index 3~Sel x Catch/fbar", "Index 4", "Index 5"), #, "China Taipei 10% only"), 
       col=c(1, "gray", 2,"darkgreen",4, "violet", "green", "orange"), cex=0.8,  border="white", lty=1, 
       lwd=c(4,4,2,2,2,2,2))  #,2))

## 

## Do residuals of fits
res=ldply(bds, function(x) x@diags)

dats=subset(res, .id==names(bds)[4])
levels(dats$name)=c("Spain BB", "Japan LL byc", "KoreaPaCu LL", "China Taipei 0%")
dats=na.omit(dats)

a=ggplot(dats)+
  geom_hline(aes(yintercept=0),col="orange")+
  # geom_smooth(aes(year,residual,col=name),fill="salmon",pch=21)+
  stat_smooth(aes(year,residual,col=name),method="lm", fill="salmon",pch=21)+
  geom_point(aes(year,residual,col=name),fill="salmon",pch=21)+
  facet_wrap(~name,ncol=1, scale="free")+theme_bw()+theme(legend.position="none")+
  ggtitle("Index 4 (no bias in China Taipei LL)")+ylab("Residuals of simulated Indices to MP")

dats=subset(res, .id==names(bds)[5])
levels(dats$name)=c("Spain BB", "Japan LL byc", "KoreaPaCu LL", "China Taipei 10%")
dats=na.omit(dats)

b=ggplot(dats)+
  geom_hline(aes(yintercept=0),col="orange")+
  # geom_smooth(aes(year,residual,col=name),fill="salmon",pch=21)+
  stat_smooth(aes(year,residual,col=name),method="lm", fill="salmon",pch=21)+
  geom_point(aes(year,residual,col=name),fill="salmon",pch=21)+
  facet_wrap(~name,ncol=1, scale="free")+theme_bw()+theme(legend.position="none")+
  ggtitle("Index 5 (10% bias in China Taipei LL)")+ylab("")

x11()    #### "SCRS2017_091_Figure11.tiff",7,7,unit="in", res=300)
multiplot(a,b, cols=2)
## 

## Now projections: 

library(FLCore) 
library(FLash) 
library(FLBRP) 
library(biodyn) ## ------------------------------------------------------------------------ 

# 1st. Build the operating models from MFCL scenarios.

# Directories, scenarios and filenames:

theme_set(theme_bw(10))

## Stock-Recruitment relationships

# Catch projection to 2040
iYr=2015
fwdYr=2035
nit=1
CVc=0

rDev2=function(runnames, iYr, fwdYr, CV){
  for (scen in 1:length(runnames)){
    srDev[[paste(runnames[scen])]]=exp(rnorm(1, FLQuant(0, dimnames=list(year=iYr:fwdYr)), CV))
  }
  res=srDev
}

srDev=FLQuants()
# srDev=rDev(scenarios, 2010, 2030, 0.2) #This generates 100 iterations.
srDev=rDev2(scenarios, iYr, fwdYr, 0.25) #This generates 1 iteration.#

omp1=om

i=201
iYr=2015
Fmsy=FLQuant(c(FLBRP:::refpts(eql[[i]])['msy','harvest']),
             dimnames=list(year=iYr:fwdYr,
                           iter=1))
omp1[[i]]=fwdWindow(omp1[[i]], end=fwdYr, refp[[i]])
omp1[[i]]=FLash:::fwd(omp1[[i]], f=Fmsy, sr=refp[[i]], sr.residuals=srDev[[i]])


setwd(dirDat)

save(omp1, file="omforMP.RData")
bmsy=as.numeric(refpts(eql[[201]])["msy", "biomass"])
stock=as.numeric(stock(omp1[[201]])/bmsy)
write.table(stock, file="stock_FUTURE.csv", sep=";")
catch=data.frame(year=1930:2035, data=as.numeric(catch(omp1[[201]])))
write.table(catch[,1:2], file="catch4proj.csv", sep=";")

# ---------------------------------------

# Directories, scenarios and filenames:

theme_set(theme_bw(10))

dirMy ="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/"
dirDat=paste(dirMy,"/RObjects",sep="")
### paste(dirMy,"/Figures/check091",sep="")
Dir="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/MFCL scenarios 240/"
scenarios=basename(list.dirs(Dir, full.names=FALSE, recursive=FALSE))



# If not loaded:

load("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/RObjects/om240.RData")
load("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/RObjects/omforMP.RData")
load("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/RObjects/saa.RData")

# Source the OEM

biol=ldply(eql, function(x) as.data.frame(x[["catch.sel","m","mat","stock.wt"]]))
biol=transform(biol,
               Quantity=factor(qname,labels=c("Selectivity","M","Maturity","Mass")),
               .id     =names(om)[X1])

i=201
cpue1=oem(window(omp1[["Base_M03_h6_q0"]],start=1981), 
          fish.dependent=FALSE)
# 2) cpue2=f(overall cpue)
cpue2=oem(window(omp1[["Base_M03_h6_q0"]],start=1981), 
          fish.dependent=TRUE,
          sel=FLQuant(FLQuant(1,dimnames=dimnames(harvest(window(omp1[["Base_M03_h6_q0"]],start=1981))))))
# 3) cpue3=f(overall age corrected cpue)

cpue3=oem(window(omp1[["Base_M03_h6_q0"]],start=1981), 
          fish.dependent=TRUE,
          sel=FLQuant(subset(biol, .id==names(om)[i] & Quantity=="Selectivity")$data*rnorm(dim(stock.n(window(omp1[["Base_M03_h6_q0"]],start=1981)))[1], 1, 0.15), 
                      dimnames=dimnames(harvest(window(omp1[["Base_M03_h6_q0"]],start=1981)))))
# 4) cpue4=f(5 cpue)

cpue4=FLQuants(llply(saa[c(1,7,10,11)], function(x) 
  mpb:::oem(window(omp1[["Base_M03_h6_q0"]],start=1981),sel=FLQuant(c(x),dimnames=dimnames(stock.n(window(omp1[["Base_M03_h6_q0"]],start=1981)))))))

# 5) cpue5=f(5 cpue with CT 10%)

cpue5=FLQuants(llply(saa[c(1,7,10,11)], function(x) 
  mpb:::oem(window(omp1[["Base_M03_h6_q0"]],start=1981),bias=FLPar(q=0.03),sel=FLQuant(c(x),dimnames=dimnames(stock.n(window(omp1[["Base_M03_h6_q0"]],start=1981)))))))

cpue5[c(1,2,4)]=cpue4[c(1,2,4)]

## Generate CPUE with all cpues

cpues=data.frame(year=1981:2035, model=rep("Stock", length(1981:2035)), Rec=rep("srDev=0.2", length(1981:2035)), data=as.numeric(cpue1)/max(as.numeric(cpue1)), fleet=rep("1 Overall Stock", length(1981:2035)))
cpues=rbind(cpues, data.frame(year=1981:2035, model=rep("CPUE", length(1981:2035)), Rec=rep("srDev=0.2", length(1981:2035)), data=as.numeric(cpue2)/max(as.numeric(cpue2)), fleet=rep("2 Overall CPUE", length(1981:2035))))
cpues=rbind(cpues, data.frame(year=1981:2035, model=rep("CPUE-Sel", length(1981:2035)), Rec=rep("srDev=0.2", length(1981:2035)),data=as.numeric(cpue3)/max(as.numeric(cpue3)), fleet=rep("3 Overall Age CPUE", length(1981:2035))))
cpues=rbind(cpues, data.frame(year=1981:2035, model=rep("ind CPUE-Sel", length(1981:2035)), Rec=rep("srDev=0.2", length(1981:2035)),data=as.numeric(cpue4[[1]])/max(as.numeric(cpue4[[1]])), fleet=rep("4A Spain BB", length(1981:2035))))
cpues=rbind(cpues, data.frame(year=1981:2035, model=rep("ind CPUE-Sel", length(1981:2035)), Rec=rep("srDev=0.2", length(1981:2035)),data=as.numeric(cpue4[[2]])/max(as.numeric(cpue4[[2]])), fleet=rep("4B Japan LL byc", length(1981:2035))))
cpues=rbind(cpues, data.frame(year=1981:2035, model=rep("ind CPUE-Sel", length(1981:2035)), Rec=rep("srDev=0.2", length(1981:2035)),data=as.numeric(cpue4[[4]])/max(as.numeric(cpue4[[4]])), fleet=rep("4C KoreaPaCu LL4", length(1981:2035))))
cpues=rbind(cpues, data.frame(year=1981:2035, model=rep("ind CPUE-Sel", length(1981:2035)), Rec=rep("srDev=0.2", length(1981:2035)),data=as.numeric(cpue4[[3]])/max(as.numeric(cpue4[[3]])), fleet=rep("4D China Taipei LL late", length(1981:2035))))
cpues=rbind(cpues, data.frame(year=1981:2035, model=rep("ind CPUE-Sel q10%", length(1981:2035)), Rec=rep("srDev=0.2", length(1981:2035)),data=as.numeric(cpue5[[3]])/max(as.numeric(cpue5[[3]])), fleet=rep("5D China Taipei 3%", length(1981:2035))))


a=ggplot(cpues)+
  geom_point(aes(year, data/max(data), col=model))+
  geom_smooth(aes(year, data/max(data), col=model))+
  theme(legend.position="none")+
  facet_wrap(~fleet, scale="free", ncol=3)+ggtitle("srDev=0.1")


## Save CPUEs:
setwd(dirDat)
save(cpues, file="SimulatedCPUEs_forProjections.RData")
write.csv(cpues, file="SimulatedCPUEs_projections.csv")

### Plot OMp1
load("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/RObjects/omforMP.RData")



a1=ggplot()+
  geom_line(aes(1930:2035, as.numeric(rec(omp1[[201]]))/1000000), lwd=1)+
  ylab("Recruitment (Millions)")+xlab("")+ggtitle("")
b1=ggplot()+
  geom_line(aes(1930:2035, as.numeric(stock(omp1[[201]]))/1000000), lwd=1)+
  ylab("Stock biomass (k tonnes)")+xlab("")
c1=ggplot()+
  geom_line(aes(1930:2035, as.numeric(catch(omp1[[201]]))/1000000), lwd=1)+
  ylab("Catch (k tonnes)")+xlab("")
d1=ggplot()+
  geom_line(aes(1930:2035, as.numeric(fbar(omp1[[201]]))), lwd=1)+
  ylab("fbar")+xlab("Years")


x11()    #### "SCRS2017_091_Figure12.tiff", res=300, unit="in", 6,7)

multiplot(a1,b1,c1,d1,cols=1)

## 

##################################################################################################
## Restart here

library(mpb)
library(plyr)
library(reshape)
#library(FLBRP)

dirMy ="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/"
dirDat=paste(dirMy,"/RObjects",sep="")
### paste(dirMy,"/Figures/check091",sep="")
Dir="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/MFCL scenarios 240/"
scenarios=basename(list.dirs(Dir, full.names=FALSE, recursive=FALSE))


catch=as.FLQuant(read.csv("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/RObjects/catch4proj.csv", sep=";"))

catch=catch/1000

cpue =FLQuants(dlply(read.csv("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/RObjects/SimulatedCPUEs_projections.csv", sep=","),
                     .(fleet),with,as.FLQuant(data.frame(year=year,data=data))))

cpue=window(cpue, start=1981)

cpue[4]=window(cpue[4], start=1981)
cpue[5]=window(cpue[5], start=1988)
cpue[6]=window(cpue[6], start=1987)
cpue[7]=window(cpue[7], start=1999)
cpue[8]=window(cpue[8], start=1999)

# Create a first bd object with approx values:


x11()    #### "SCRS2017_091_Figure13.tiff", 7,6, unit="in", res=300)
mpb:::plotIndex(cpue)+ylab("Simulated CPUE index")+facet_wrap(~qname, ncol=3, scale="free")+xlim(1981, 2035)
## 

# Now simulate stock assessment with Biodyn

bd=biodyn("pellat", params=FLPar(r=0.3,k=5.5e5,b0=1,p=0.001),
          catch = catch)

# Create the list of objects that inlcudes the Base CAse and Sensitivities

# Specify controls and do fits for each scenario 

bdsf=mpb:::biodyns(list("Stock"=bd,
                        "CPUE"=bd,
                        "CPUE~Sel"=bd,
                        "Fleet"=bd,
                        "Fleet China-Tai 5%"=bd))

params(bdsf[[1]])=FLPar(r=.1,k=3.2e6,p=0.001,b0=1) 
setParams(bdsf[[1]])=cpue[1]
setControl(bdsf[[1]])=params(bdsf[[1]])
bdsf[[1]]=fit(bdsf[[1]],cpue[1])

params(bdsf[[2]])=FLPar(r=.1,k=3.2e6,p=0.001,b0=1) 
setParams(bdsf[[2]])=cpue[2]
setControl(bdsf[[2]])=params(bdsf[[2]])
bdsf[[2]]=fit(bdsf[[2]],cpue[2])

params(bdsf[[3]])=FLPar(r=.1,k=3.2e6,p=0.001,b0=1) 
setParams(bdsf[[3]])=cpue[3]
setControl(bdsf[[3]])=params(bdsf[[3]])
bdsf[[3]]=fit(bdsf[[3]],cpue[3])

params(bdsf[[4]])=FLPar(r=0.1,k=3.2e6, p=0.001,b0=1) 
setParams(bdsf[[4]])=cpue[4:7]
setControl(bdsf[[4]])=params(bdsf[[4]])
bdsf[[4]]=fit(bdsf[[4]],cpue[4:7])

params(bdsf[[5]])=FLPar(r=.1,k=3.2e6,p=0.001,b0=1) 
setParams(bdsf[[5]])=cpue[c(4,5,6,8)]
setControl(bdsf[[5]])=params(bdsf[[4]])
bdsf[[5]]=fit(bdsf[[5]],cpue[c(4,5,6,8)])


## Figures of Stock Status



theme_set(theme_bw(10))

omstock=read.csv(paste(dirDat, "/stock_FUTURE.csv", sep=""), sep=";")

x11()    #### "SCRS2017_091_Figure14.tiff", 5,5, unit="in", res=300)

par(mar=c(4,4,1,1))
plot(1930:2035, omstock[,1], 
     ylim=c(0, 3), type="l", ylab="MP estimated B/Bmsy", lwd=2,xlab="Years")


grid()
lines(1930:2036, as.numeric(stock(bdsf[[1]])/bmsy(bdsf[[1]])),col=2, lwd=1)
lines(1930:2036, as.numeric(stock(bdsf[[2]])/bmsy(bdsf[[2]])), col="darkgreen", lwd=1)
lines(1930:2036, as.numeric(stock(bdsf[[3]])/bmsy(bdsf[[3]])), col=4, lwd=1)
lines(1930:2036, as.numeric(stock(bdsf[[4]])/bmsy(bdsf[[4]])), col="violet", lwd=1)
lines(1930:2036, as.numeric(stock(bdsf[[5]])/bmsy(bdsf[[5]])), col="green", lwd=1)

lines(1930:2035, omstock[,1], lwd=4)

legend(1980, 3, c("OM", "Index 1~Stock", "Index 2~Catch/fbar",
                  "Index 3~Sel x Catch/fbar", "Index 4", "Index 5"), 
       col=c(1, 2,"darkgreen",4, "violet", "green"), cex=0.8,  border="white", lty=1, 
       lwd=c(4,4,2,2,2,2))

## 



## Do residuals of fits
res=ldply(bdsf, function(x) x@diags)


dats=subset(res, .id==names(bdsf)[4])
levels(dats$name)=c("Spain BB", "Japan LL byc", "KoreaPaCu LL", "China Taipei")
dats=na.omit(dats)
a=ggplot(dats)+
  geom_hline(aes(yintercept=0),col="orange")+
  # geom_smooth(aes(year,residual,col=name),fill="salmon",pch=21)+
  stat_smooth(aes(year,residual,col=name),method="lm", fill="salmon",pch=21)+
  geom_point(aes(year,residual,col=name),fill="salmon",pch=21)+
  facet_wrap(~name,ncol=1, scale="free")+theme_bw()+theme(legend.position="none")+
  ggtitle("Index 4 (no bias in China Taipei LL)")+ylab("Residuals of simulated Indices to MP")

dats=subset(res, .id==names(bdsf)[5])
levels(dats$name)=c("Spain BB", "Japan LL byc", "KoreaPaCu LL", "China Taipei 3%")
dats=na.omit(dats)
b=ggplot(dats)+
  geom_hline(aes(yintercept=0),col="orange")+
  # geom_smooth(aes(year,residual,col=name),fill="salmon",pch=21)+
  stat_smooth(aes(year,residual,col=name),method="lm", fill="salmon",pch=21)+
  geom_point(aes(year,residual,col=name),fill="salmon",pch=21)+
  facet_wrap(~name,ncol=1, scale="free")+theme_bw()+theme(legend.position="none")+
  ggtitle("Index 5 (10% bias in China Taipei LL)")+ylab("")

x11()    #### "SCRS2017_091_Figure15.tiff",7,7,unit="in", res=300)
multiplot(a,b, cols=2)
## 


## Final Figures: MP fitting models with error:
# Convenient to restart:

library(mpb)
library(plyr)
library(reshape)


# Source additional functions
source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/RScripts/Multiplot.R")
source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/RScripts/OEM.R")
source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/RScripts/diags-MFCL.R")

dirMy ="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/"
dirDat=paste(dirMy,"/RObjects",sep="")
### paste(dirMy,"/Figures/check091",sep="")

catch=as.FLQuant(read.csv(paste(dirDat, "/catch.csv", sep=""), sep=";"))

cpue =FLQuants(dlply(read.csv(paste(dirDat, "/SimulatedCPUEs.csv",sep=""), sep=","),
                     .(fleet),with,as.FLQuant(data.frame(year=year,data=data))))

cpue=window(cpue, start=1981)

cpue[4]=window(cpue[4], start=1981)
cpue[5]=window(cpue[5], start=1988, end=2012)
cpue[6]=window(cpue[6], start=1987)
cpue[7]=window(cpue[7], start=1999)
cpue[8]=window(cpue[8], start=1999)



## Now fit iteratively: --- I'm sure there is an easier way to do this.

options(digits=2)

cv=as.numeric(c(0.0, 0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4))
d=0

cpuev=cpue

for (j in 1:length(cv)){   ## alternatives for cv
  
  for (i in 1:20){
    
bd=biodyn("pellat", params=FLPar(r=0.3,k=5.5e5,b0=1,p=0.001),
          catch = catch)

bdsc=mpb:::biodyns(list("Index 4 outCV"=bd))

params(bdsc[[1]])=FLPar(r=0.1,k=2.2e6, p=0.001,b0=1) 

cpuev[[4]]=cpue[[4]]*rlnorm(length(cpue[[4]]), 0, cv[j])
cpuev[[5]]=cpue[[5]]*rlnorm(length(cpue[[5]]), 0, cv[j])
cpuev[[6]]=cpue[[6]]*rlnorm(length(cpue[[6]]), 0, cv[j])
cpuev[[7]]=cpue[[7]]*rlnorm(length(cpue[[7]]), 0, cv[j])

setParams(bdsc[[1]])=cpuev[4:7]
setControl(bdsc[[1]])=params(bdsc[[1]])
bdsc[[1]]=fit(bdsc[[1]],cpuev[4:7])

if (i==1 & j==1) {results=data.frame(CV=rep(paste("CV=", cv[j]), length(1930:2015)),
                   iter=rep(i, length(1930:2015)), year=1930:2015, 
                   relB=as.numeric(stock(bdsc[[1]])/refpts(bdsc[[1]])$bmsy),
                   relF=c(NA,as.numeric(harvest(bdsc[[1]])/refpts(bdsc[[1]])$fmsy)),
                   r=rep(as.numeric(params(bdsc[[1]])$r), length(1930:2015)), 
                   K=rep(as.numeric(params(bdsc[[1]])$k), length(1930:2015)),
                   MSY=rep(as.numeric(refpts(bdsc[[1]])$msy), length(1930:2015))) 
            }else{
                results=rbind(results, data.frame(CV=rep(paste("CV=", cv[j]), length(1930:2015)), iter=rep(i, length(1930:2015)), 
                                    year=1930:2015, relB=as.numeric(stock(bdsc[[1]])/refpts(bdsc[[1]])$bmsy),
                                    relF=c(NA,as.numeric(harvest(bdsc[[1]])/refpts(bdsc[[1]])$fmsy)),
                                    r=rep(as.numeric(params(bdsc[[1]])$r), length(1930:2015)), 
                                    K=rep(as.numeric(params(bdsc[[1]])$k), length(1930:2015)),
                                    MSY=rep(as.numeric(refpts(bdsc[[1]])$msy), length(1930:2015))))
  
                }

d=d+1
print(paste("iteration ", i, ", CV ", cv[j], ". Run ", d, " of ", 20*length(cv), sep=""))
  }
}

yref=subset(results, CV==paste("CV=", cv[1]))$year
sref=subset(results, CV==paste("CV=", cv[1]))$relB
fref=subset(results, CV==paste("CV=", cv[1]))$relF
refs=data.frame(yref=yref, sref=sref, fref=fref)



x11()    #### "SCRS2017_091_Figure16.tiff",6,8, unit="in", res=300)

a=ggplot(results)+
  geom_line(aes(year, relB, group=iter))+facet_wrap(~CV)+ylab("MP: B/Bmsy")+xlab("Years")+
  geom_line(aes(yref, sref), data=refs, col="green", lwd=1)+
  stat_smooth(aes(year, relB), col="blue")+
  theme(legend.position="none")+
  geom_hline(yintercept = 1, col=2, lty=2)

b=ggplot(results)+
  geom_line(aes(year, relF, group=iter))+facet_wrap(~CV)+ylab("MP: F/Fmsy")+xlab("Years")+
  stat_smooth(aes(yref, fref), data=refs, col="green", lwd=1)+
  stat_smooth(aes(year, relF), col="blue")+
  theme(legend.position="none")+
  geom_hline(yintercept = 1, col=2, lty=2)

multiplot(a,b)

## 


## Compare a random fit from Figure 16 to Madeira SA. Fits and residuals. 
