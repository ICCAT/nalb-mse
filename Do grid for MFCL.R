## Create files for MFCL

# Starting from the initial scenarios, make them equal with dynamic q

library(R4MFCL)  # Requires an oldish version of R, run OK with 2.15.3

source("c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/RScripts/Functions_readMFCL.R")

mfcldir="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/MFCL runs"
dir13="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/2013_assessment"
gridir="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/files for Ai Kimoto Review/MFCL runs/"
runs0=list.dirs(dir13, full.names=TRUE, recursive=FALSE)[-11]
MAlt5dir="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/MFCL runs/Alt5M"

# ---------------------------------------------------------
# 1) Create grid of M
# ---------------------------------------------------------

for (d in 1:length(runs0)){

ll=length(strsplit(runs0[d], "/")[[1]])
nameI=strsplit(runs0[d], "/")[[1]][ll]

setwd(runs0[d])
iniOb=readLines("albN.ini")

Mnames=c("M02", "M03", "M04")

for (m in 1:3){

if (m==1) iniOb[11]="\t0.2\t\t\t\t\t\t\t\t\t\t\t\t\t\t"
if (m==2) iniOb[11]="\t0.3\t\t\t\t\t\t\t\t\t\t\t\t\t\t"
if (m==3) iniOb[11]="\t0.4\t\t\t\t\t\t\t\t\t\t\t\t\t\t"

griddirM=paste(gridir, "/", nameI, "_", Mnames[m], "_h6_q0", sep="")

dir.create(griddirM)
setwd(griddirM)

if (nameI!="Alt5") write(iniOb, "albN.ini")

if (m==1 & nameI=="Alt5") { 
  dirM=list.dirs(MAlt5dir, recursive=FALSE)[1]
  file.copy(paste(dirM, "albN.ini", sep="/"), griddirM) 
}

if (m==2 & nameI=="Alt5") { 
  dirM=list.dirs(MAlt5dir, recursive=FALSE)[2]
  file.copy(paste(dirM, "albN.ini", sep="/"), griddirM) 
}

if (m==3 & nameI=="Alt5") { 
  dirM=list.dirs(MAlt5dir, recursive=FALSE)[3]
  file.copy(paste(dirM, "albN.ini", sep="/"), griddirM) 
}

file.copy(paste(runs0[d], "doitallN.alb", sep="/"), griddirM)
file.copy(paste(runs0[d], "albN.frq", sep="/"), griddirM)
file.copy(paste(mfcldir, "mfcl.exe", sep="/"), griddirM)

if (d==10) file.copy(paste(runs0[d], "albN.tag", sep="/"), griddirM)

}
}

# ---------------------------------------------------------
# 2) x 2 dynamic q
# ---------------------------------------------------------

## Create files for MFCL

# Starting from the initial scenarios, make them equal with dynamic q

dir13="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/2013_assessment"
gridir="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/MFCL runs/maspruebas2/test5"
runs=list.dirs(gridir, full.names=FALSE)[-1]

for(d in 1:length(runs)){
  
  nameI=substr(runs[d], 1, unlist(gregexpr("q0", runs[d]))[1]-1)
  namefolderqdir=paste(nameI, "q1", sep="")
  
  frq=read.frq(paste(runs[d], "albN.frq", sep="/"))
  
  frq$mat=frq$mat[-1,]
  
  for (i in 2:length(frq$mat[,1])){
    f=as.numeric(frq$mat[i,4])
    fleet=frq$mat[-1,4]
    rr=range(frq$mat[fleet==f,1])
    iyr=min(rr)
    if (iyr<1975) iyr=1980
    fyr=max(rr)
    if (fyr<iyr) fyr=1980
    old=as.numeric(frq$mat[i,6])
    
    if (frq$mat[i,6]>=0.1 & frq$mat[i,1]>1980) frq$mat[i, 6]=as.numeric(frq$mat[i, 6])*1.01^(as.numeric(frq$mat[i, 1])-iyr) 
    
    new=as.numeric(frq$mat[i,6])
 #   print(c(as.numeric(frq$mat[i, 1])-iyr, as.numeric(frq$mat[i, 1]), iyr, f, old, new))
    
  }
  
  dir.create(namefolderqdir)
  setwd(namefolderqdir)
  
  write.frq("albN.frq", frq.obj = frq)
  
  file.copy(paste(runs[d], "doitallN.alb", sep="/"), namefolderqdir)
  file.copy(paste(runs[d], "albN.ini", sep="/"), namefolderqdir)
  file.copy(paste(runs[d], "mfcl.exe", sep="/"), namefolderqdir)
  
  if (d>=28) file.copy(paste(runs[d], "albN.tag", sep="/"), namefolderqdir)
  
}

# ---------------------------------------------------------
# 3) x 4 steepness: weak prior (h) and strict prior (sd=0.05, h7, h8, h9)
# ---------------------------------------------------------

# First, generate the 240 folders with the scenarios.

gridir="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/MFCL runs/maspruebas2/test5"
runs=list.dirs(gridir, full.names=FALSE, recursive=FALSE)

for(d in 1:length(runs)){

  for (h in 1:3){
  
  if (h==1) nameh="7"
  if (h==2) nameh="8"
  if (h==3) nameh="9"
  
  nameI=substr(runs[d], 1, unlist(gregexpr("q", runs[d]))[1]-3)
  nameI2=substr(runs[d], unlist(gregexpr("_", runs[d]))[5], nchar(runs[d]))
  namehfile=paste(nameI, nameh, nameI2, sep="")
 
  dir.create(namehfile)
  setwd(namehfile)
  
  file.copy(paste(runs[d], "albN.frq", sep="/"), namehfile)
  file.copy(paste(runs[d], "albN.ini", sep="/"), namehfile)
  file.copy(paste(runs[d], "mfcl.exe", sep="/"), namehfile)
  file.copy(paste(runs[d], "doitallN.alb", sep="/"), namehfile)
  
  if (d>=55) file.copy(paste(runs[d], "albN.tag", sep="/"), namehfile)

   }
  }

# Now modify priors for steepness: Works only OK with Linux (cluster).

gridir="c:/use/Gorka/1_AZTI_tuna/ICCAT/1.ALB/MSE/MFCL runs/maspruebas2/tes4"
runs=list.dirs(gridir, full.names=FALSE, recursive=FALSE)

for (d in 1:length(runs)){
  
  setwd(runs[d])
  
  ll=length(strsplit(runs[d], "/")[[1]])
  nameI=strsplit(runs[d], "/")[[1]][ll]
  steepness=substr(runs[d], unlist(gregexpr("h", runs[d]))[1], unlist(gregexpr("_q", runs[d]))[1]-1)
  
  doi=readLines("doitallN.alb", n=-1L)
  print(doi)
  
  l1=which(grepl(153, doi))
  l2=which(grepl(154, doi))
  
  if (steepness=="h7") {
    doi[l1] = "    2 153 368 #130        # parameters of beta distribution defining prior for"
    doi[l2] = "    2 154  225 # 50        # steepness - mode = 0.76, sd = 0.08" 
  }
  
  if (steepness=="h8") {
    doi[l1] = "    2 153 357 #130        # parameters of beta distribution defining prior for"
    doi[l2] = "    2 154  126 # 50        # steepness - mode = 0.76, sd = 0.08" 
  }
  
  if (steepness=="h9") {
    doi[l1] = "    2 153 269 #130        # parameters of beta distribution defining prior for"
    doi[l2] = "    2 154  47 # 50        # steepness - mode = 0.76, sd = 0.08" 
  }
  
  write(doi, "doitallN_mod.alb")
  
}

# The batch file to run for each scenario is "doitallN_mod.alb"