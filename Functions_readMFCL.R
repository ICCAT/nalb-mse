# Functions to read from MFCL

## Read ouput
condOM=function(model, DirDat, runnames, filenames){
  if (model=="VPA2Box") om=readVPA2Box(paste(DirDat, runname, filename,sep=""))
  if (model=="mfcl") for (scen in 1:length(runnames)){ om[[paste(runnames[scen])]]=readMFCL(paste(DirDat, runnames[scen], filenames,sep="/")) }
  res=om
}

## Stock Rec relationships

getBHSR <- function (plotrepfile = "plot.rep") {
  
  scanText<-function(string, what = character(0), ...){
    ## Like scan() but reading from a vector of character strings
    tc <- textConnection(string)
    result <- scan(tc, what = what, quiet = TRUE, ...)
    close(tc)
    return(result)}
  
  ##==============================================================
  ## Beverton-Holt stock-recruitment relationship report
  ##==============================================================
  con <- file(plotrepfile,open="rt")
  on.exit(close(con))
  
  while((tt <- readLines(con,n=1)) != "# Beverton-Holt stock-recruitment relationship report"){}
  
  a <- as.numeric(strsplit(readLines(con,n=1)," +")[[1]][c(4,7,10)])
  dat <- list(alpha=a[1],beta=a[2],steepness=a[3])
  a <- readLines(con,n=8)
  if(length(grep("# Observed spawning Biomass",a[1]))>0)
    dat$spawnB.o <- scanText(a[2],what=0)
  if(length(grep("# Observed recruitment",a[3]))>0)
    dat$recruit.o <- scanText(a[4],what=0)
  if(length(grep("# Spawning Biomass",a[5]))>0)
    dat$spawnB <- scanText(a[6],what=0)
  if(length(grep("# Predicted recruitment",a[7]))>0)
    dat$recruit.p <- scanText(a[8],what=0)
  dat}

srOM=function(model, runnames){
  
  for (scen in 1:length(runnames)){ 
    sr[[scen]]=as.FLSR(om[[scen]], model=model)
    sr[[scen]] = fmle(sr[[scen]], control = list(silent = TRUE), fixed = list(b = 0))
  } 
  res=sr
}

## Reference points

refP=function(runnames){
  for (scen in 1:length(runnames)){
    refp[[paste(runnames[scen])]]=FLBRP(om[[scen]], sr=sr[[scen]], fbar=seq(1,3,length.out=101))
    refp[[paste(runnames[scen])]]=brp(refp[[paste(runnames[scen])]])
  }
  res=refp
}

## Recruitment Deviates

rDev=function(runnames, iYr, fwdYr, CV){
  for (scen in 1:length(runnames)){
    srDev[[paste(runnames[scen])]]=exp(rnorm(100, FLQuant(0, dimnames=list(year=iYr:fwdYr)), CV))
  }
  res=srDev
}

scanText <- function(string, what=character(0),...){
  tc <- textConnection(string)
  result <- scan(tc, what=what, quiet=TRUE,...)
  close(tc)
  return(result)
}
