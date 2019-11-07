## Script to evaluate HCRs on a series of OMs with MSE 

# Load libraries
library(FLCore)  ## install.packages("FLCore", repos="http://flr-project.org/R")
library(ggplotFL)  ## install.packages("ggplotFL", repos="http://flr-project.org/R")
library(reshape)   ## install.packages("reshape")  # instalara tambien 'plyr'
library(plyr)   ## install.packages("plyr")
# library(kobe)   # install.packages("bbmle")
# install_github("lauriekell/kobe")
library(FLBRP)  ## install.packages("FLBRP", repos="http://flr-project.org/R")

library(biodyn) ## install.packages("biodyn", repos="http://flr-project.org/R")
library(popbio)  ## install.packages("popbio")
library(mpb)  ## install_github("lauriekell/mpb")

# source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/RScripts/new_oem.R")

## Define folders
dirMFCL ="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/files for the alb wg/MFCL scenarios 240//"  # Change this to where the user installs the files
scenarios=basename(list.dirs(dirMFCL, full.names=FALSE, recursive=FALSE))

dirMy ="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Rerun HCRs/"
dirDat="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/Rerun HCRs/RObjects/"
dir2018="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/2018 improvements/2018 runs/output/500 iterations/"

load(paste(dirDat,"/om240.RData",    sep=""))
load(paste(dirDat,"/saa.RData",    sep=""))
# load(paste(dirDat,"/om17.RData",    sep=""))
#load(paste(dirDat,"/om.RData",    sep=""))
#load(paste(dirDat,"/sa.RData",sep=""))

## Additional functions

# source(paste(dirMy, "/RScripts/SCRS093_FUNCTIONS.R", sep=""))
source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/2018 improvements/2018 runs/RScripts/SCRS093_functions_AdoptedHCR_exact.R")
source(paste(dirMy, "/RScripts/fitPella.R", sep=""))

##

biol=ldply(eql, function(x) as.data.frame(x[["catch.sel","m","mat","stock.wt"]]))
biol=transform(biol,
               Quantity=factor(qname,labels=c("Selectivity","M","Maturity","Mass")),
               .id     =names(om)[X1])

nits  = 500

options(digits=2)

# Create a data frame with the HCRs

HCRs=array(dim=c(4, 4))
HCRs[,1]=0.4
HCRs[,2]=0.1
HCRs[1:2,3]=0.8
HCRs[3:4,3]=1
HCRs[,4]=rep(c(0.8, 1), 2)

HCRs=as.data.frame(HCRs)

names(HCRs)=c("Blim", "Fmin", "Btrig", "Ftar")

# scenarios_MSE will be a list of all the HCRs and OMs run.

d=0

# Run MSE through the runMSE function for all om and HCRs.

# for (i in 1: length(names(om))){ # Empezar por i=2

oms=seq(1,240,2)  ## Select only OMs with q0 (non dynamic q).
x=seq(194, 216, 2)   ## Add Robustness CASE: Base case + dynamic q.
oms=c(oms, x)   ## Subscripts of the desired OMs.   

# for (z in 132:1){ 
  
#  for (j in 1:length(HCRs[,1])){
for (z in 1:132){ 

  j=3      
    i=oms[z]
    
    setwd(dir2018)
    
    namefile=paste("MSE_hcr_om", i, "_adoptedHCR_exact.RData", sep="")[]    
    
    if(file.exists(namefile)){print(namefile)}
    if(!file.exists(namefile)){
      
      
      MSE_hcr=runMSE(window(om[[i]], start=1930),  # change the start to 1950 would be more realistic but doesn't work. Check sa.
                     eql[[i]],
                     range=c("start"=2015,"end"=2015+30, interval=3),
                     hcr=c(ftar=HCRs[j, 4],
                           fmin=HCRs[j, 2],
                           blim=HCRs[j, 1],
                           btrig=HCRs[j, 3]),
                     srDev=0.4, maxF=2, 
                     uDev=0.2) 
      
      setwd(dir2018)
      save(MSE_hcr, file=paste(namefile))
      rm(MSE_hcr)
      
    }
    
  }


