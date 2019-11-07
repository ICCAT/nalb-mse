#------------------------------------------------------------------------##
# Script to read the MFCL runs and create a grid of 132 OMs:
# 10 from 2013 x 3 M multipliers (0.9, 1.0 and 1.1) x 3 steepness (0.7, 0.8, 0.9) x 2 q trends (0,1%)
# Names are: Scenario 2013 x M x h x q: 
# Eg. "Alt1M10h75q0" is "Alt1" from 2013, M=Mx1.0, h=0.75 and qtrend =0%.

# libraries:
# install.packages("FLCore", repos="http://flr-project.org/R",dependencies=TRUE) ## 
# install.packages("FLash", repos="http://flr-project.org/R") ## 
# install.packages("FLBRP", repos="http://flr-project.org/R") ## 
# install.packages("ggplotFL",repos="http://flr-project.org/R") ## 
# install.packages("diags", repos="http://flr-project.org/R") ## ## 
# library(devtools) ## 
# install_github("biodyn","laurieKell",dependencies=TRUE) #, echo=TRUE)

library(FLCore) 
library(FLash) 
library(FLBRP) 
library(biodyn) ## ------------------------------------------------------------------------ 

# 1st. Build the operating models from MFCL scenarios.

# Directories, scenarios and filenames:

Dir="c:/use/Gorka//1_AZTI_tuna/ICCAT/1.ALB/MSE/MFCL runs/2017 grid/"
scenarios=basename(list.dirs(Dir, full.names=FALSE, recursive=FALSE))
files=c("plot-07.par.rep","07.par")

# Read output from Assessment model:

# functions
library(R4MFCL)
source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/RScripts/Functions_readMFCL.R")

# Operating models as a form of FLStock  
om=FLStocks()
om=condOM("mfcl",Dir, scenarios, files) # read OM

x11()
plot(om)+theme(legend.position="none")

## Stock-Recruitment relationships

sr=FLSRs(llply(om,as.FLSR,model="bevholt"))
sr=FLSRs(llply(sr, fmle,control=list(silent=TRUE)))

for (i in 1:length(scenarios)){
  
  setwd(paste(Dir, scenarios[i], sep=""))
  
  x=scenarios[i]
  origin=substr(x, 1, unlist(gregexpr("_M", x))-1)
  
  if (origin!="Alt5") sr2=getBHSR("plot-07.par.rep")
  if (origin=="Alt5") sr2=getBHSR("plot-08.par.rep")
  
  params(sr[[i]])["a"]=sr2$alpha
  params(sr[[i]])["b"]=sr2$beta*1000
  
}

names(sr)=names(om)

plot(sr[[1]])

## Reference Points

refp=FLBRPs()
refp=refP(scenarios)

# Rec Deviates

rDev2=function(runnames, iYr, fwdYr, CV){
  for (scen in 1:length(runnames)){
    srDev[[paste(runnames[scen])]]=exp(rnorm(1, FLQuant(0, dimnames=list(year=iYr:fwdYr)), CV))
  }
  res=srDev
}

srDev=FLQuants()
# srDev=rDev(scenarios, 2010, 2030, 0.0) #This generates 100 iterations.
srDev=rDev2(scenarios, 2012, 2014, 0.0) #This generates 1 iteration.

## eql

eql=FLBRPs(mlply(names(sr),function(i) 
  brp(FLBRP(om[[i]],sr=sr[[i]],nyears=dims(om[[i]])$year))))
names(eql)=names(om)

plot(eql[[1]],refpts=TRUE)+theme_bw(10)+ylab("")+xlab("")+theme(legend.position="none")

# Catch projection to 2015

## Ct catch projections

fwdYr=2014
nit=1
CVc=100000
CVc=0

om1=om

for (scen in 1:length(scenarios)) { 
  
  iYr=2012
  if (scen>=121 & scen<=144) iYr=2008
  
  TAC=c(25679972.3, 24633579.2, 26660423.9) #runif(1, min(catch(om[[scen]]), na.rm=TRUE), max(catch(om[[scen]]), na.rm=TRUE))  ## Here random
  TAC=FLQuant(TAC, dimnames=list(year=iYr:fwdYr)) 
  
  if (scen>=121 & scen<=144) {
    TAC=c(20038926.6, 15374980.61, 19509094.54, 20038926.64, 25679972.3, 24633579.2, 26660423.9) #runif(1, min(catch(om[[scen]]), na.rm=TRUE), max(catch(om[[scen]]), na.rm=TRUE))  ## Here random
    TAC=FLQuant(TAC, dimnames=list(year=iYr:fwdYr)) 
                              }
  
  om[[paste(scenarios[scen])]]=fwdWindow(om[[paste(scenarios[scen])]], 
                                         end=fwdYr, refp[[paste(scenarios[scen])]])
  catch = TAC
  
  om[[paste(scenarios[scen])]]=FLash:::fwd(om[[paste(scenarios[scen])]], catch=catch, 
                                           sr=refp[[paste(scenarios[scen])]], sr.residuals=srDev[[paste(scenarios[scen])]])
  
print(names(om)[scen])
  }


# save OMs an object:
setwd("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/RObjects/")
save(om, eql, refp, sr, rDev2, file="om240.RData")



x11()
plot(om)

# Explore the range of uncertainties on OMs.
refp=eql

RPdat=data.frame()

RPdat=data.frame(OM=names(om)[1], 
                 origin=substr(scenarios[1], 1, unlist(gregexpr("_", scenarios[1]))[1]-1),
                 M=substr(scenarios[1], unlist(gregexpr("_", scenarios[1]))[1]+1, unlist(gregexpr("_", scenarios[1]))[2]-1),
                 h=substr(scenarios[1], unlist(gregexpr("_", scenarios[1]))[2]+1, unlist(gregexpr("_", scenarios[1]))[3]-1), 
                 q=substr(scenarios[1], unlist(gregexpr("_", scenarios[1]))[3]+1, unlist(gregexpr("_", scenarios[1]))[3]+2), 
                 B0=as.numeric(refpts(refp[[1]])["virgin", "biomass"]),
                 SSB0=as.numeric(refpts(refp[[1]])["virgin", "ssb"]),
                 MSY=as.numeric(refpts(refp[[1]])["msy", "yield"]),
                 BMSY=as.numeric(refpts(refp[[1]])["msy", "biomass"]),
                 SMSY=as.numeric(refpts(refp[[1]])["msy", "ssb"]),
                 FMSY=as.numeric(refpts(refp[[1]])["msy", "harvest"]),
                 B15=as.numeric(tail(stock(om[[1]]), 1)/
                                  refpts(refp[[1]])["msy", "biomass"]),
                 S15=as.numeric(tail(ssb(om[[1]]), 1)/
                                  refpts(refp[[1]])["msy", "ssb"]),
                 F15=as.numeric((tail(catch(om[[1]]), 1))/(tail(stock(om[[1]]), 1))/
                                  refpts(refp[[1]])["msy", "harvest"]))

for (i in 2:length(scenarios)){
  
RPdat=rbind(RPdat, data.frame(OM=names(om)[i], 
                                origin=substr(scenarios[i], 1, unlist(gregexpr("_", scenarios[i]))[1]-1),
                                M=substr(scenarios[i], unlist(gregexpr("_", scenarios[i]))[1]+1, unlist(gregexpr("_", scenarios[i]))[2]-1),
                                h=substr(scenarios[i], unlist(gregexpr("_", scenarios[i]))[2]+1, unlist(gregexpr("_", scenarios[i]))[3]-1), 
                                q=substr(scenarios[i], unlist(gregexpr("_", scenarios[i]))[3]+1, unlist(gregexpr("_", scenarios[i]))[3]+2), 
                              B0=as.numeric(refpts(refp[[i]])["virgin", "biomass"]),
                              SSB0=as.numeric(refpts(refp[[i]])["virgin", "ssb"]),
                              MSY=as.numeric(refpts(refp[[i]])["msy", "yield"]),
                              BMSY=as.numeric(refpts(refp[[i]])["msy", "biomass"]),
                              SMSY=as.numeric(refpts(refp[[i]])["msy", "ssb"]),
                              FMSY=as.numeric(refpts(refp[[i]])["msy", "harvest"]),
                              B15=as.numeric(tail(stock(om[[i]]), 1)/
                                               refpts(refp[[i]])["msy", "biomass"]),
                              S15=as.numeric(tail(ssb(om[[i]]), 1)/
                                               refpts(refp[[i]])["msy", "ssb"]),
                              F15=as.numeric((tail(catch(om[[i]]), 1))/(tail(stock(om[[i]]), 1))/
                                               refpts(refp[[i]])["msy", "harvest"])))
            
}


## Exploratory figures to reflect the range of uncertainty covered with the grid:

# 1) Expansion of uncertainty for each of the original scenarios
B0=ggplot(RPdat)+geom_line(aes(as.factor(origin), B0/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), B0/1000000, col=M, pch=h, fg=q))+
  #geom_point(aes(as.factor(origin), B0/1000000), size=2, col="gray30", data=subset(RPdat, M=="M03"))+
  geom_hline(yintercept=1105.477, col="blue", lty=2)+
  ylim(0, 1200)+
  xlab("")+ylab("B0 (th tonnes)")

MSY=ggplot(RPdat)+geom_line(aes(as.factor(origin), MSY/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), MSY/1000000, col=M, pch=h, fg=q), size=2)+
  geom_hline(yintercept=37.0817632, col="blue", lty=2)+
  ylim(0, max(RPdat$MSY)/1000000)+
  xlab("")+ylab("MSY (th tonnes)")

BMSY=ggplot(RPdat)+geom_line(aes(as.factor(origin), BMSY/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), BMSY/1000000, col=M, pch=h, fg=q), size=2)+
  #geom_point(aes(as.factor(origin), BMSY/1000000), size=2, col="gray30", data=subset(RPdat, M=="M03"))+
  geom_hline(yintercept=406.89, col="blue", lty=2)+
  ylim(0, 500)+
  xlab("")+ylab("BMSY (th tonnes)")

FMSY=ggplot(RPdat)+geom_line(aes(as.factor(origin), FMSY), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), FMSY), size=2, col="gray30")+
  geom_point(aes(as.factor(origin), FMSY), size=2, col="gray30", data=subset(RPdat, M=="M03"))+
  geom_hline(yintercept=0.097, col="blue", lty=2)+
  ylim(0, max(RPdat$FMSY))+
  xlab("")+ylab("FMSY")

b15=ggplot(RPdat)+geom_line(aes(as.factor(origin), B15), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), B15), size=2, col="gray30")+
  geom_point(aes(as.factor(origin), B15), size=2, col="gray30", data=subset(RPdat, M=="M03"))+
  ylim(0, max(RPdat$B15))+geom_hline(yintercept=1, col="red")+
  geom_hline(yintercept=1.4032, col="blue", lty=2)+
  xlab("")+ylab("B(2015)/BMSY")

f15=ggplot(RPdat)+geom_line(aes(as.factor(origin), F15), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), F15), size=2, col="gray30")+
  geom_point(aes(as.factor(origin), F15), size=2, col="gray30", data=subset(RPdat, M=="M03"))+
  ylim(0, max(RPdat$F15))+geom_hline(yintercept=1, col="red")+
  geom_hline(yintercept=0.518, col="blue", lty=2)+
  xlab("")+ylab("F(2015)/FMSY")

source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/RScripts/Multiplot.R")

## 

# In blue results from 2016 SA with biodyn
setwd("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/RObjects/OMfigs//")

tiff("RPs covered by OM grid.tiff",7, 7, unit="in", res=300)
multiplot(B0, MSY, b15, BMSY, FMSY, f15, cols=2)
dev.off()



## 
