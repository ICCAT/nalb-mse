#------------------------------------------------------------------------##
# Script to read the MFCL runs and create a grid of 240 OMs:
# 10 from 2013 x 3 M multipliers (0.9, 1.0 and 1.1) of 0.3 + 4 steepness (0.75w, 0.7, 0.8, 0.9) x 2 q trends (0,1%)
# Scripts adapted by Gorka Merino (AZTI)  gmerino@azti.es
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
library(biodyn)
library(plyr)

# 1st. Build the operating models from MFCL scenarios.

## ----------------------------------------------------------------------------------------------------------------
# Directories, scenarios and filenames: ---------------------------------------------------------------------------
dirMy ="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Files ICCAT/"  # Change this to where the user installs the files
Dir=paste(dirMy, "MFCL scenarios 240", sep="/")
dirDat=paste(dirMy,"/RObjects",sep="")

## ----------------------------------------------


scenarios=basename(list.dirs(Dir, full.names=FALSE, recursive=FALSE))
files=c("plot-07.par.rep","07.par")
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

# Read output from Assessment model:

# Additional functions ------------------------------------------------------------

source(paste(dirMy, "RScripts/Functions_readMFCL.R", sep="/"))
source(paste(dirMy, "RScripts/Multiplot.R", sep="/"))

## *************************************************************************************
## *************************************************************************************
# 1) Condition Operating models as a form of FLStock from MFCL runs  

om=FLStocks()
om=condOM("mfcl",Dir, scenarios, files) # read OM

x11()
plot(om)+theme(legend.position="none")

## Stock-Recruitment relationships

sr=FLSRs(llply(om,as.FLSR,model="bevholt"))
sr=FLSRs(llply(sr, fmle,control=list(silent=TRUE)))

for (i in 1:length(scenarios)){
  
  setwd(paste(Dir, scenarios[i], sep="/"))
  
  x=scenarios[i]
  origin=substr(x, 1, unlist(gregexpr("_M", x))-1)
  
  if (origin!="Alt5") sr2=getBHSR("plot-07.par.rep")
  if (origin=="Alt5") sr2=getBHSR("plot-08.par.rep")
  
  params(sr[[i]])["a"]=sr2$alpha
  params(sr[[i]])["b"]=sr2$beta*1000
  
}

names(sr)=names(om)

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

# Catch projection to 2015

## Ct catch projections

fwdYr=2014
nit=1
CVc=100000
CVc=0

# om1=om

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
                                           sr=refp[[paste(scenarios[scen])]], 
                                           sr.residuals=srDev[[paste(scenarios[scen])]])

print(names(om)[scen])

}



# save OMs an object:
setwd(paste(dirMy, "/RObjects/", sep="/"))

save(om, eql, refp, sr, rDev2, file="om240.RData")

####################################################################

## Figures for SCRS/2017/092


gridir=paste(dirMy, "MFCL scenarios 240/", sep="/")
scenarios=list.dirs(gridir, recursive=FALSE)
scenario=list.dirs(gridir, recursive=FALSE, full.names = FALSE)

for (i in 1:length(scenarios)){
  
  setwd(scenarios[i])
  
  print(i)
  
  if (i>=97 & i<=120){
    for (j in 1:8){  
      parfile=readLines(paste("0", j, ".par", sep=""))
      lineOb=which(grepl("# Objective function value", parfile))
      linePar=which(grepl("# The number of parameters", parfile))
      ObV=as.numeric(parfile[lineOb+1])
      ParV=as.numeric(parfile[linePar+1])
      
      plotpar=readLines(paste("plot-0", 8, ".par.rep", sep=""))
      lineSt=which(grepl("# Beverton-Holt stock-recruitment relationship report", plotpar))
      stpns=as.numeric(substr(plotpar[lineSt+1], 55, 100))
      
      runsDiags=rbind(runsDiags, data.frame(scenario=scenario[i], Phase=j, ObFun=ObV, Npars=ParV, stpns=stpns))
      
    }}else{
      
      for (j in 1:7){  
        parfile=readLines(paste("0", j, ".par", sep=""))
        lineOb=which(grepl("# Objective function value", parfile))
        linePar=which(grepl("# The number of parameters", parfile))
        ObV=as.numeric(parfile[lineOb+1])
        ParV=as.numeric(parfile[linePar+1])
        
        plotpar=readLines(paste("plot-0", 7, ".par.rep", sep=""))
        lineSt=which(grepl("# Beverton-Holt stock-recruitment relationship report", plotpar))
        stpns=as.numeric(substr(plotpar[lineSt+1], 55, 100))
        
        
        if (i==1 & j==1) {runsDiags=data.frame(scenario=scenario[i], Phase=j, ObFun=ObV, Npars=ParV, stpns=stpns) } 
        else {
          runsDiags=rbind(runsDiags, data.frame(scenario=scenario[i], Phase=j, ObFun=ObV, Npars=ParV, stpns=stpns))
        }
      }
    }
}

x=runsDiags$scenario
runsDiags$M=substr(x, unlist(gregexpr("M", x)), unlist(gregexpr("_h", x))-1)
runsDiags$dynq=substr(x, unlist(gregexpr("q", x)), length(x))
runsDiags$origin=substr(x, 1, unlist(gregexpr("_M", x))-1)


setwd(dirFigs)

x11()     ## ("SCRS2017_092_Figure1.x11()     ## ", 8, 4, unit="in", res=300)

ggplot(runsDiags)+
  geom_line(aes(Phase, ObFun, group=scenario, col=origin))+
  guides(col=guide_legend(title=""))+
  ylab("Objective function")+
  theme_bw()

#

x11()     ## ("SCRS2017_092_Figure2.x11()     ## ", 8, 4, unit="in", res=300)

ggplot(runsDiags)+
  geom_line(aes(Phase, Npars, group=scenario, col=origin))+
  ylab("Number of parameters")+
  guides(col=guide_legend(title=""))+
  theme_bw()

#


x=runsDiags$scenario
runsDiags$h=substr(x, unlist(gregexpr("h", x)), unlist(gregexpr("_q", x))-1)

datsh6=subset(runsDiags, h=="h6")
datsh7=subset(runsDiags, h=="h7")
datsh8=subset(runsDiags, h=="h8")
datsh9=subset(runsDiags, h=="h9")

x11()     ## ("SCRS2017_092_Figure12.x11()     ## ", 7, 2, unit="in", res=300)

par(mfrow=c(1,4), mar=c(4,3,3,0.5))

hist(rbeta(10000, 3.5, 2.1), xlim=c(0,1), main="mode 0.75, sd = 0.15", xlab="", ylab="", col="gray", breaks=30, font.main=1, freq=FALSE)
for (i in 1:length(datsh6$scenario)){
  abline(v=datsh6$stpns, col="red", lwd=0.5)
}

hist(rbeta(10000, 36.8, 22.5), xlim=c(0,1), main="mode 0.7, sd = 0.05", xlab="", ylab="", col="gray", breaks=30, font.main=1, freq=FALSE)
for (i in 1:length(datsh7$scenario)){
  abline(v=datsh7$stpns, col="darkgreen", lwd=0.5)
}

hist(rbeta(10000, 35.7, 12.6), xlim=c(0,1), main="mode 0.8, sd = 0.05", xlab="", ylab="", col="gray", breaks=30, font.main=1, freq=FALSE)
for (i in 1:length(datsh8$scenario)){
  abline(v=datsh8$stpns, col="lightblue", lwd=0.5)
}

hist(rbeta(10000, 26.9, 4.7), xlim=c(0,1), main="mode 0.9, sd = 0.05", xlab="", ylab="", col="gray", breaks=30, font.main=1, freq=FALSE)
for (i in 1:length(datsh9$scenario)){
  abline(v=datsh9$stpns, col="violet", lwd=0.5)
}

#

## Figure 3**  Equilibrium curves
eqdat=plot(eql,refpts=FALSE)$data
x=eqdat$.id
eqdat$M=substr(x, unlist(gregexpr("M", x)), unlist(gregexpr("_h", x))-1)
eqdat$h=substr(x, unlist(gregexpr("h", x)), unlist(gregexpr("_q", x))-1)
eqdat$dynq=substr(x, unlist(gregexpr("q", x)), length(x))
eqdat$origin=substr(x, 1, unlist(gregexpr("_M", x))-1)

# For base case
eqplot0bc=ggplot(subset(eqdat, origin=="Base" & M=="M03" & h=="h6" & dynq=="q0"))+
  geom_line(aes(x, y, col=origin, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("2013 base case only")

eqplot1bc=ggplot(subset(eqdat, origin=="Base" & h=="h6" & dynq=="q0"))+
  geom_line(aes(x, y, col=M, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("2013 base case ~ M")

eqplot2bc=ggplot(subset(eqdat, origin=="Base" &  M=="M03" & dynq=="q0"))+
  geom_line(aes(x, y, col=h, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("2013 base case ~ Steepness")

eqplot3bc=ggplot(subset(eqdat,origin=="Base" & M=="M03" & h=="h6"))+
  geom_line(aes(x, y, col=dynq, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("2013 base case ~ Dynamic q")

# For 2013
eqplot0o=ggplot(subset(eqdat, M=="M03" & h=="h6" & dynq=="q0"))+
  geom_line(aes(x, y, col=origin, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("2013 runs only")

eqplot1o=ggplot(subset(eqdat, h=="h6" & dynq=="q0"))+
  geom_line(aes(x, y, col=M, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("2013 runs ~ Natural Mortality")

eqplot2o=ggplot(subset(eqdat, M=="M03" & dynq=="q0"))+
  geom_line(aes(x, y, col=h, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("2013 runs ~ Steepness")

eqplot3o=ggplot(subset(eqdat, M=="M03" & h=="h6"))+
  geom_line(aes(x, y, col=dynq, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("2013 runs ~ Dynamic catchability")

# For 2017

eqplot0=ggplot(eqdat)+
  geom_line(aes(x, y, col=origin, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("grid ~ 2013 runs")

eqplot1=ggplot(eqdat)+
  geom_line(aes(x, y, col=M, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("grid ~ Natural Mortality")

eqplot2=ggplot(eqdat)+
  geom_line(aes(x, y, col=h, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("grid ~ Steepness")

eqplot3=ggplot(eqdat)+
  geom_line(aes(x, y, col=dynq, group=.id))+
  facet_wrap(~pnl, scale="free", ncol=1)+xlab("")+ylab("")+
  theme(legend.position="top")+
  guides(col=guide_legend(title=""))+
  ggtitle("grid ~ Dynamic catchability")


x11()     ## ("SCRS2017_092_Figure3.x11()     ## ", 12.8, 7, unit="in", res=300)
multiplot(eqplot0bc, eqplot1bc, eqplot2bc, eqplot3bc, cols=4)
#

x11()     ## ("SCRS2017_092_Figure4.x11()     ## ", 12.8, 7, unit="in", res=300)
multiplot(eqplot0o, eqplot1o, eqplot2o, eqplot3o, cols=4)
#

x11()     ## ("SCRS2017_092_Figure5.x11()     ## ", 12.8, 7, unit="in", res=300)
multiplot(eqplot0, eqplot1, eqplot2, eqplot3, cols=4)
#


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
# 2013
M13=ggplot(subset(bioldat, Quantity=="Selectivity" & h=="h6" & dynq=="q0"))+
  geom_line(aes(age, data, col=M, group=.id))+
  theme(legend.position="top")+ggtitle("2013 ~ Natural Mortality")+
  guides(col=guide_legend(title=""))+
  ylab("Selectivity")+xlab("Age")+facet_wrap(~origin, ncol=1)

h13=ggplot(subset(bioldat, Quantity=="Selectivity" & dynq=="q0" & M=="M03"))+
  geom_line(aes(age, data, col=h, group=.id))+ylab("")+
  theme(legend.position="top")+ggtitle("2013 ~ Steepness")+
  guides(col=guide_legend(title=""))+
  xlab("Age")+facet_wrap(~origin, ncol=1)

q13=ggplot(subset(bioldat, Quantity=="Selectivity" & h=="h6" & M=="M03"))+
  geom_line(aes(age, data, col=dynq, group=.id))+ylab("")+
  theme(legend.position="top")+ggtitle("2013 ~ Catchability")+
  guides(col=guide_legend(title=""))+
  xlab("Age")+facet_wrap(~origin, ncol=1)

x11()     ## ("SCRS2017_092_Figure8.x11()     ## ", 7, 8, unit="in", res=300)
multiplot(M13, h13, q13, cols=3)
#

plotbiol=ggplot(biol)+
  geom_line(aes(age,data,group=.id))+ 
  facet_wrap(~Quantity,scale="free_y")+ 
  theme(legend.position="none")+
  scale_colour_discrete(name="Scenario")+
  xlab("Age") +ylab("")


## There is an error in the SCRS paper, there are two Figure 4.

x11()     ## ("SCRS2017_092_Figure4B.x11()     ## ", 4, 3, unit="in", res=300)

ggplot(subset(bioldat,  Quantity=="M" & h!="h6"))+
  geom_line(aes(age, data, col=M, group=.id, lty=Mage))+
  ylab("")+xlab("Age")+facet_wrap(~Quantity, scale="free")+
  guides(col=guide_legend(title=""), lty=guide_legend(title=""))+
  scale_linetype_manual(values=c("dashed", "solid")) # Change linetypes

#
####

# Play wi


x11()     ## ("SCRS2017_092_Figure9.x11()     ## ", 6, 5, unit="in", res=300)
plot(om)+theme(legend.position="none")+ scale_color_manual(values=c(rep("gray", 240)))
#

## Reference Points
srs=sr

x11()     ## ("SCRS2017_092_Figure10.x11()     ## ", 6, 5, unit="in", res=300)
plot(srs[[193]])   # "Base_h6_M03_q0"
#

# Explore the range of uncertainties on OMs.
# Figure 2 Time series of recruits, SSB, biomass, mean F, F apex and Catch

scenarios=basename(list.dirs(Dir, full.names=FALSE, recursive=FALSE))

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


## oms: Contains the numbers of the OMs finally used in SCRS/2017/093

## Exploratory figures to reflect the range of uncertainty covered with the grid:
oms=seq(1,240,2)  ## Select only OMs with q0 (non dynamic q).
x=seq(194, 216, 2)   ## Add Robustness CASE: Base case + dynamic q.
oms=c(oms, x) 


# 1) Expansion of uncertainty for each of the original scenarios
B0=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), B0/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), B0/1000000, color=M, shape=h), size=2.5)+
  # geom_point(aes(as.factor(origin), B0/1000000, shape=h, col=M), fill="yellow", size=1.8, data=subset(RPdat, q=="q1"))+
  #geom_point(aes(as.factor(origin), B0/1000000), size=2, col="gray30", data=subset(RPdat, M=="M03"))+
  #geom_point(aes(as.factor(origin), B0/1000000, shape=h), color="white", size=2, data=subset(RPdat, q=="q1"))+
  #  geom_hline(yintercept=1105.477, col="blue", lty=2)+
  ylim(0, 1200)+guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  xlab("")+ylab("B0 (th tonnes)")

MSY=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), MSY/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), MSY/1000000, col=M, shape=h), size=2.5)+
  #  geom_hline(yintercept=37.0817632, col="blue", lty=2)+
  ylim(0, max(RPdat$MSY)/1000000)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  xlab("")+ylab("MSY (th tonnes)")

BMSY=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), BMSY/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), BMSY/1000000, col=M, shape=h), size=2.5)+
  #geom_hline(yintercept=406.89, col="blue", lty=2)+
  ylim(0, 500)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  xlab("")+ylab("BMSY (th tonnes)")

FMSY=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), FMSY), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), FMSY, col=M, shape=h), size=2.5)+
  #  geom_hline(yintercept=0.097, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  ylim(0, max(RPdat$FMSY))+
  xlab("")+ylab("FMSY")

b15=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), B15), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), B15, col=M, shape=h), size=2.5)+
  ylim(0,2)+geom_hline(yintercept=1, col="red")+
  #  geom_hline(yintercept=1.4032, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  xlab("")+ylab("B(2015)/BMSY")

f15=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), F15), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), F15, col=M, shape=h), size=2.5)+
  ylim(0, 3.5)+geom_hline(yintercept=1, col="red")+
  #  geom_hline(yintercept=0.518, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  xlab("")+ylab("F(2015)/FMSY")


setwd(dirFigs)

x11()     ## ("SCRS2017_092_Figure6.x11()     ## ", 8, 8, unit="in", res=300)

multiplot(B0, MSY, b15, BMSY, FMSY, f15, cols=2)

#
## 


BMSYh=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), BMSY/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), BMSY/1000000, col=h, shape=q), size=2)+
  #geom_hline(yintercept=406.89, col="blue", lty=2)+
  ylim(0, 500)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("red", "blue", "darkgreen", "violet"))+ # Change linetypes
  xlab("")+ylab("BMSY (th tonnes)")

FMSYh=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), FMSY), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), FMSY, col=h, shape=q), size=2)+
  #  geom_hline(yintercept=0.097, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("red", "blue", "darkgreen", "violet"))+ # Change linetypes
  ylim(0, max(RPdat$FMSY))+
  xlab("")+ylab("FMSY")

BMSYq=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), BMSY/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), BMSY/1000000, col=q, shape=M), size=2)+
  #geom_hline(yintercept=406.89, col="blue", lty=2)+
  ylim(0, 500)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("green3", "brown"))+ # Change linetypes
  xlab("")+ylab("BMSY (th tonnes)")

FMSYq=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), FMSY), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), FMSY, col=q, shape=M), size=2)+
  #  geom_hline(yintercept=0.097, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("green3", "brown"))+ # Change linetypes
  ylim(0, max(RPdat$FMSY))+
  xlab("")+ylab("FMSY")

b15q=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), B15), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), B15, col=q, shape=M), size=2)+
  ylim(0,2)+geom_hline(yintercept=1, col="red")+
  #  geom_hline(yintercept=1.4032, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("green3", "brown"))+ # Change linetypes
  xlab("")+ylab("B(2015)/BMSY")

f15q=ggplot(RPdat[oms,])+geom_line(aes(as.factor(origin), F15), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), F15, col=q, shape=M), size=2)+
  ylim(0, 3.5)+geom_hline(yintercept=1, col="red")+
  #  geom_hline(yintercept=0.518, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("green3", "brown"))+ # Change linetypes
  xlab("")+ylab("F(2015)/FMSY")

x11()     ## ("SCRS2017_092_Figure7.x11()     ## ", 8, 6, unit="in", res=300)

multiplot(BMSYq,  b15q, FMSYq, f15q, cols=2)

#

## This figures originally were different in SCRS092 as they would also include 120 scenarios with dynamic q

# 1) Expansion of uncertainty for each of the original scenarios
B0=ggplot(RPdat)+
  geom_line(aes(as.factor(origin), B0/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), B0/1000000, color=M, shape=h), size=2.5)+
  # geom_point(aes(as.factor(origin), B0/1000000, shape=h, col=M), fill="yellow", size=1.8, data=subset(RPdat, q=="q1"))+
  #geom_point(aes(as.factor(origin), B0/1000000), size=2, col="gray30", data=subset(RPdat, M=="M03"))+
  #geom_point(aes(as.factor(origin), B0/1000000, shape=h), color="white", size=2, data=subset(RPdat, q=="q1"))+
  #  geom_hline(yintercept=1105.477, col="blue", lty=2)+
  ylim(0, 1200)+guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  xlab("")+ylab("B0 (th tonnes)")

MSY=ggplot(RPdat)+geom_line(aes(as.factor(origin), MSY/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), MSY/1000000, col=M, shape=h), size=2.5)+
  #  geom_hline(yintercept=37.0817632, col="blue", lty=2)+
  ylim(0, max(RPdat$MSY)/1000000)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  xlab("")+ylab("MSY (th tonnes)")

BMSY=ggplot(RPdat)+geom_line(aes(as.factor(origin), BMSY/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), BMSY/1000000, col=M, shape=h), size=2.5)+
  #geom_hline(yintercept=406.89, col="blue", lty=2)+
  ylim(0, 500)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  xlab("")+ylab("BMSY (th tonnes)")

FMSY=ggplot(RPdat)+geom_line(aes(as.factor(origin), FMSY), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), FMSY, col=M, shape=h), size=2.5)+
  #  geom_hline(yintercept=0.097, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  ylim(0, max(RPdat$FMSY))+
  xlab("")+ylab("FMSY")

b15=ggplot(RPdat)+geom_line(aes(as.factor(origin), B15), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), B15, col=M, shape=h), size=2.5)+
  ylim(0,2)+geom_hline(yintercept=1, col="red")+
  #  geom_hline(yintercept=1.4032, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  xlab("")+ylab("B(2015)/BMSY")

f15=ggplot(RPdat)+geom_line(aes(as.factor(origin), F15), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), F15, col=M, shape=h), size=2.5)+
  ylim(0, 3.5)+geom_hline(yintercept=1, col="red")+
  #  geom_hline(yintercept=0.518, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  xlab("")+ylab("F(2015)/FMSY")


setwd(dirFigs)

x11()     ## ("SCRS2017_092_Figure6_original.x11()     ## ", 8, 8, unit="in", res=300)

multiplot(B0, MSY, b15, BMSY, FMSY, f15, cols=2)

#
## 


BMSYh=ggplot(RPdat)+geom_line(aes(as.factor(origin), BMSY/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), BMSY/1000000, col=h, shape=q), size=2)+
  #geom_hline(yintercept=406.89, col="blue", lty=2)+
  ylim(0, 500)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("red", "blue", "darkgreen", "violet"))+ # Change linetypes
  xlab("")+ylab("BMSY (th tonnes)")

FMSYh=ggplot(RPdat)+geom_line(aes(as.factor(origin), FMSY), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), FMSY, col=h, shape=q), size=2)+
  #  geom_hline(yintercept=0.097, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("red", "blue", "darkgreen", "violet"))+ # Change linetypes
  ylim(0, max(RPdat$FMSY))+
  xlab("")+ylab("FMSY")

BMSYq=ggplot(RPdat)+geom_line(aes(as.factor(origin), BMSY/1000000), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), BMSY/1000000, col=q, shape=M), size=2)+
  #geom_hline(yintercept=406.89, col="blue", lty=2)+
  ylim(0, 500)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("green3", "brown"))+ # Change linetypes
  xlab("")+ylab("BMSY (th tonnes)")

FMSYq=ggplot(RPdat)+geom_line(aes(as.factor(origin), FMSY), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), FMSY, col=q, shape=M), size=2)+
  #  geom_hline(yintercept=0.097, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("green3", "brown"))+ # Change linetypes
  ylim(0, max(RPdat$FMSY))+
  xlab("")+ylab("FMSY")

b15q=ggplot(RPdat)+geom_line(aes(as.factor(origin), B15), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), B15, col=q, shape=M), size=2)+
  ylim(0,2)+geom_hline(yintercept=1, col="red")+
  #  geom_hline(yintercept=1.4032, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("green3", "brown"))+ # Change linetypes
  xlab("")+ylab("B(2015)/BMSY")

f15q=ggplot(RPdat)+geom_line(aes(as.factor(origin), F15), lwd=2, col="gray")+
  geom_point(aes(as.factor(origin), F15, col=q, shape=M), size=2)+
  ylim(0, 3.5)+geom_hline(yintercept=1, col="red")+
  #  geom_hline(yintercept=0.518, col="blue", lty=2)+
  guides(col=guide_legend(title=""), pch=guide_legend(title=""))+
  scale_color_manual(values=c("green3", "brown"))+ # Change linetypes
  xlab("")+ylab("F(2015)/FMSY")

x11()     ## ("SCRS2017_092_Figure7_original.x11()     ## ", 8, 6, unit="in", res=300)

multiplot(BMSYq,  b15q, FMSYq, f15q, cols=2)

#

x11()     ## ("Histogram_B2015.x11()     ## ", 8, 6, unit="in", res=300)

hist(RPdat[oms,]$B15, col="gray", xlab="Operating Models B2015/Bmsy", main="OMs Stock in 2015, 59% < Bmsy")

rect(1,-2,6,40, col="lightgreen", border="limegreen")#
rect(-1,-2,1,40, col="pink", border="red3")

hist(RPdat[oms,]$B15, col="gray", xlab="Operating Models B2015/Bmsy", main="OMs Stock in 2015, 59% < Bmsy", add=TRUE)


abline(v=1, col="red", lty=2)
abline(v=RPdat[201,]$B15, col="blue", lty=1)

text(1.25,24, "Base Case MFCL", col="blue")

#
### ---------------------------------------------------------------------------------

# Figure 11 requires a projection with R regimes:

# Catch projection to 2040

#
iYr=2015
fwdYr=2035
nit=1
CVc=1000
CVc=0

# Add increase and reduction of rec in SrDev function

rDev2=function(runnames, iYr, fwdYr, CV){
  for (scen in 1:length(runnames)){
    srDev[[paste(runnames[scen])]]=exp(rnorm(1, FLQuant(0, dimnames=list(year=iYr:fwdYr)), CV))
  }
  res=srDev
}

srDev=FLQuants()
# srDev=rDev(scenarios, 2010, 2030, 0.2) #This generates 100 iterations.
srDev=rDev2(scenarios, iYr, fwdYr, 0.0) #This generates 1 iteration.#


omp1=om
omp2=om
omp3=om

for (i in 1:length(names(om))) { 
  
  iYr=2015
  
  Fmsy=FLQuant(c(FLBRP:::refpts(eql[[i]])['msy','harvest']),
               dimnames=list(year=iYr:fwdYr,
                             iter=1))
  
  omp1[[i]]=fwdWindow(omp1[[i]], end=fwdYr, refp[[i]])
  omp2[[i]]=fwdWindow(omp2[[i]], end=fwdYr, refp[[i]])
  omp3[[i]]=fwdWindow(omp3[[i]], end=fwdYr, refp[[i]])
  
  omp1[[i]]=FLash:::fwd(omp1[[i]], f=Fmsy, sr=refp[[i]], sr.residuals=srDev[[i]])
  omp2[[i]]=FLash:::fwd(omp2[[i]], f=Fmsy, sr=refp[[i]], sr.residuals=1.2*srDev[[i]])
  omp3[[i]]=FLash:::fwd(omp3[[i]], f=Fmsy, sr=refp[[i]], sr.residuals=0.8*srDev[[i]])
  
  print(paste(names(om)[i], i))
  
}


p1=plot(omp1)+theme(legend.position="none")+ scale_color_manual(values=c(rep("gray", 240)))+ggtitle("Current regime")
p2=plot(omp2)+theme(legend.position="none")+ scale_color_manual(values=c(rep("red", 240)))+ggtitle("Higher recruitment (+20%)")
p3=plot(omp3)+theme(legend.position="none")+ scale_color_manual(values=c(rep("blue", 240)))+ggtitle("Lower recruitment (-20%)")

setwd(dirFigs)

x11()     ## ("SCRS2017_092_Figure11.x11()     ## ", 9, 6, unit="in", res=300)
multiplot(p1, p2, p3, cols=3)
#


