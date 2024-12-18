setGeneric('hcr', function(object,...) standardGeneric('hcr'))

setMethod('hcr', signature(object='biodyn'),
          function(object, 
                   params=hcrParam(ftar =0.70*refpts(object)['fmsy'],
                                   btrig=0.80*refpts(object)['bmsy'],
                                   fmin =0.01*refpts(object)['fmsy'],
                                   blim =0.40*refpts(object)['bmsy']),
                   stkYrs=max(as.numeric(dimnames(stock(object))$year)),
                   refYrs=max(as.numeric(dimnames(catch(object))$year)),
                   hcrYrs=max(as.numeric(dimnames(stock(object))$year)),
                   tac   =TRUE,
                   tacMn =TRUE,
                   bndF  =NULL, #c(1,Inf),
                   bndTac=c(0.9,1.1),
                   iaF   =TRUE, 
                   iaTac =TRUE, 
                   maxF  =2,
                   ...) {
            ## HCR
            dimnames(params)$params=tolower(dimnames(params)$params)
            params=as(params,'FLQuant')  
            #if (blim>=btrig) stop('btrig must be greater than blim')
            a=(params['ftar']-params['fmin'])/(params['btrig']-params['blim'])
            b=params['ftar']-a*params['btrig']
            
            ## Calc F
            # bug
            #val=(SSB%*%a) %+% b
            stk=FLCore::apply(stock(object)[,ac(stkYrs)],6,mean)
            
            rtn=(stk%*%a)  
            rtn=FLCore::sweep(rtn,2:6,b,'+')
            
            fmin=as(params['fmin'],'FLQuant')
            ftar=as(params['ftar'],'FLQuant')
            for (i in seq(dims(object)$iter)){
              FLCore::iter(rtn,i)[]=max(FLCore::iter(rtn,i),FLCore::iter(fmin,i))
              FLCore::iter(rtn,i)[]=min(FLCore::iter(rtn,i),FLCore::iter(ftar,i))} 
            
            rtn=window(rtn,end=max(hcrYrs))
            #dimnames(rtn)$year=min(hcrYrs)  
            if (length(hcrYrs)>1){
              rtn=window(rtn,end=max(hcrYrs))
              rtn[,ac(hcrYrs)]=rtn[,dimnames(rtn)$year[1]]}
            
            ### Bounds ##################################################################################
            ## F
            if (!is.null(bndF)){  
              
              ref=FLCore::apply(harvest(object)[,ac(refYrs-1)],6,mean)
              
              rtn[,ac(min(hcrYrs))]=qmax(rtn[,ac(min(hcrYrs))],ref*bndF[1])
              rtn[,ac(min(hcrYrs))]=qmin(rtn[,ac(min(hcrYrs))],ref*bndF[2])
              
              if (length(hcrYrs)>1)        
                for (i in hcrYrs[-1]){
                  if (iaF){
                    rtn[,ac(i)]=qmax(rtn[,ac(i)],rtn[,ac(i-1)]*bndF[1])
                    rtn[,ac(i)]=qmin(rtn[,ac(i)],rtn[,ac(i-1)]*bndF[2])
                  }else{
                    rtn[,ac(i)]=rtn[,ac(i-1)]}
                  
                  if (!is.null(maxF)) rtn=qmin(rtn,maxF)}}
            hvt=rtn
            
            
            ## TAC
            if (tac){
              
              ref=FLCore::apply(catch(object)[,ac(refYrs)],6,mean)
              
              object=window(object, end=max(as.numeric(hcrYrs)))
              # # print("check before error 74")
              object=fwd(object,harvest=harvest(object)[,ac(min(as.numeric(hcrYrs)-1))])
              # # print("check before error 76")
              rtn   = catch(fwd(object, harvest=rtn))[,ac(hcrYrs)]
              print(rtn)
              
              ref=apply(rtn[,ac(refYrs)],6,mean)  #change TAC constaint years
      #        print(ref)
              
              print("check")
              if (!is.null(bndTac)){  
                
                rtn[,ac(hcrYrs)]=qmax(rtn[,ac(hcrYrs)],ref*bndTac[1])
                rtn[,ac(hcrYrs)]=qmin(rtn[,ac(hcrYrs)],ref*bndTac[2])
                print("TACs")
                print( rtn)
                
              if (length(hcrYrs)>1)        
                  for (i in hcrYrs[-1]){
                    if (iaTac){
                      rtn[,ac(i)]=qmax(rtn[,ac(i)],rtn[,ac(i-1)]*1000000*bndTac[1])
                      rtn[,ac(i)]=qmin(rtn[,ac(i)],rtn[,ac(i-1)]*1000000*bndTac[2])
                    }else{
                      rtn[,ac(i)]=rtn[,ac(i-1)]}}
                
                if (tacMn) rtn[]=c(apply(rtn,3:6,mean))}}
            
            if (tac) rtn=list(hvt=hvt,tac=rtn,stock=stk) else rtn=list(hvt=hvt,stock=stk)
            
            return(rtn)})



runMSE<-function(om,eql,
                 hcr,
                 srDev,uCV,
                 range=c(start   =dims(om)$maxyear,
                         end     =dims(om)$maxyear+30,
                         interval=3),
                 maxF=2,
                 maxTAC=35000000,  # 50 th tons
                 uDev,
                 iaF=.1,
                 cmdOps='-maxfn 500 -iprint 0 -est'){
  
  ## OM projections
  # # print("check before error 108")
  om  =fwdWindow(om,end=range["end"],eql)
  
  ## F in longterm
  prj=FLQuant(c(FLBRP:::refpts(eql)['msy','harvest']),
              dimnames=list(year=range["start"]:range["end"],
                            iter=seq(nits)))
  
#  prj=FLQuant(c(FLBRP:::refpts(eql)['msy','harvest']*hcr['ftar']),
  #              dimnames=list(year=range["start"]:range["end"],
  #                          iter=seq(nits)))
  
  ## Add stochastcity
  if (!(srDev%in%"FLQuant"))
    set.seed(11)
    srDev=rlnorm(nits,FLQuant(0,
                              dimnames=list(year=range["start"]:range["end"])),srDev)
  # print("check before error 121")
  om =FLash:::fwd(om,f=prj, sr=eql,sr.residuals=srDev) 
  prj=om
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eql=dims(params(eql))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) 
    ("Stop, iters not '1 or n' in OM")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  ## Cut in capacity
  maxF=mean(apply(fbar(window(om,end=range["start"])),6,mean)*maxF)
  
  #### Observation Error (OEM) setup #######################
  ## Random variation for CPUE  
  #cpue=oem(window(om,end=range["start"]), uDev)
  #cpue=FLQuants(llply(saa[c(1,7,10,11)], function(x) 
  

  cpue=FLQuants(llply(saa[c(1,7,10,11)], function(x) 
    mpb:::oem(window(om,start=1981,end=range["start"]),sel=FLQuant(c(x),dimnames=dimnames(stock.n(window(om,start=1981, end=range["start"])))))))

#  cpue=oem(window(om,end=range["start"]))

  #cpue=oem(window(om,end=range["start"]))
# print("OK here1") 
  ## use data from last year
  cpue=window(cpue,end=range["start"]-1)
  #cpue[,ac(range["start"]-(range["interval"]:1))]=oem(om[,ac(range["start"]-(range["interval"]:1))], uDev)
# print("OK here2")   
#  cpue[,ac(range["start"]-(range["interval"]:1))]=oem(om[,ac(range["start"]-(range["interval"]:1))])
  cpue=window(cpue, start=1981) 
  # As done in Madeira
  cpue[[1]]=window(cpue[[1]], start=1981)
  cpue[[2]]=window(cpue[[2]], start=1988)
  cpue[[3]]=window(cpue[[3]], start=1999)
  cpue[[4]]=window(cpue[[4]], start=1987) 
  
  set.seed(11)
  cpue[[1]]=cpue[[1]]*rlnorm(length(cpue[[1]]), 0, uDev)
  set.seed(12)
  cpue[[2]]=cpue[[2]]*rlnorm(length(cpue[[2]]), 0, uDev)
  set.seed(13)
  cpue[[3]]=cpue[[3]]*rlnorm(length(cpue[[3]]), 0, uDev)
  set.seed(14)
  cpue[[4]]=cpue[[4]]*rlnorm(length(cpue[[4]]), 0, uDev)
      
  #### Management Procedure
  ## Set up assessment parameter options
  #sa@control["p","val"]=1
    
  mp=biodyn:::FLStock2biodyn(window(om,end=range["start"]-1))
 # pnms=dimnames(control(sa))$param[dimnames(control(sa))$param%in%
  #                                   dimnames(params(mp))$params]
  pnms=as.character(c("r", "k", "p", "b0"))
  # params(mp)[pnms]=apply(sa@params[pnms],1,median)
  params(mp)[pnms]=FLPar(r=.1,k=3.2e9,p=0.001,b0=1) 
  
# Add constraints to r and K. r>=0.1 and K<=2^6.
  #control(mp)["r","min"]=0.1
  #control(mp)["r","max"]=0.4
  #control(mp)["k","max"]=3500000000
  #control(mp)["k","val"]=2000000000
  #control(mp)["k","min"]=1000000000
  
  setParams( mp)=cpue
  
  setControl(mp)=params(mp)
  
# mp@control["p","val"]=0.001
  
  # print(mp@params)
  # print(mp@control)
  
  ## Loop round years
  TACsave=FLQuant(NA, 
                  dimnames=list(year=range["start"]:range["end"],
                                iter=1:nits))
  
  ## Loop round years
  for (iYr in seq(range["start"],range["end"]-range["interval"],range["interval"])){
    # print("check before error 1")
    cat('\n===================', iYr, '===================\n')
    # print("check before error 2")
    ## use data from last year
    mp=window(mp,end=iYr)
        cpue=window(cpue,end=iYr-1)
     
    #cpue[,     ac(iYr-(range["interval"]:1))]=oem(  om[,ac(iYr-(range["interval"]:1))], uDev)
    # cpue[,     ac(iYr-(range["interval"]:1))]=oem(  om[,ac(iYr-(range["interval"]:1))])
     cpue=FLQuants(llply(saa[c(1,7,10,11)], function(x) 
       mpb:::oem(window(om,start=1981,end=iYr),sel=FLQuant(c(x),dimnames=dimnames(stock.n(window(om,start=1981, end=iYr)))))))
    
     
     #cpue[[1]]=window(cpue[[1]], start=1981)
     #cpue[[2]]=window(cpue[[2]], start=1988)
     #cpue[[3]]=window(cpue[[3]], start=1999)
     #cpue[[4]]=window(cpue[[4]], start=1987) 
     cpue[[1]]=window(cpue[[1]], start=1981)
     cpue[[2]]=window(cpue[[2]], start=1988)
     cpue[[3]]=window(cpue[[3]], start=1999)
     cpue[[4]]=window(cpue[[4]], start=1987) 
     #
     
    catch(mp)[,ac(iYr-(range["interval"]:1))]=catch(om[,ac(iYr-(range["interval"]:1))])        
    # print(c(iYr, range["interval"]:1))
    ## fit
    # print("check before error 6")
    mmm<<-list(om=om,mp=mp,prj=prj,cpue=cpue)
    #mp =biodyn::fwd(mp,catch=catch(mp))
    # print("check before error 7")
    catch(mp)[is.na(catch(mp))]=0.01
    # print("check before error 8")
    #mp<<-mp
    ## print("check before error 9")
    #cpue<<-cpue
    # print("check before error 10")
    #save(mp,cpue, file="test.RData")
  #  plot(mp)
  #  mp =biodyn:::fit(mp,FLQuants(A=cpue))
  #  browser()
   # print("here")
#    cpue=cpue*rlnorm(length(cpue), 0, uDev)  ## Add observation error to the CPUE index. Error: logNormal(0, uDev)
    set.seed(11)
    cpue[[1]]=cpue[[1]]*rlnorm(length(cpue[[1]]), 0, uDev)
    set.seed(12)
    cpue[[2]]=cpue[[2]]*rlnorm(length(cpue[[2]]), 0, uDev)
    set.seed(13)
    cpue[[3]]=cpue[[3]]*rlnorm(length(cpue[[3]]), 0, uDev)
    set.seed(14)
    cpue[[4]]=cpue[[4]]*rlnorm(length(cpue[[4]]), 0, uDev)
    
    pnms=as.character(c("r", "k", "p", "b0"))
    #params(mp)[pnms]=FLPar(r=.1,k=3.2e9,p=0.001,b0=1)     # What was done in Madeira
    params(mp)[pnms]=FLPar(r=.1,k=3.2e9,p=0.001,b0=1) 

    setParams( mp)=cpue[1:4]
    setControl(mp)=params(mp)
    
     # print("here2")  
#
    #params(mp)=FLPar(r=.1,k=3.2e9,b0=1,p=0.001)
    #setParams(mp)=cpue[1:4]
    #setControl(mp)=params(mp)
  
    mp =fitPella(mp,cpue[1:4])
    # print("check error")
    #print(mp@params)
    #print(mp@control)
    
    ## HCR
    hcrPar=hcrParam(ftar =hcr["ftar" ]*fmsy(mp),
                    btrig=hcr["btrig"]*bmsy(mp),
                    fmin =hcr["fmin" ]*fmsy(mp), 
                    blim =hcr["blim" ]*bmsy(mp))
    
    #hcrOutcome=biodyn:::hcr(mp,hcrPar,
    #                       hcrYrs=iYr+seq(range["interval"]),
    #                       refYrs=iYr, #change TAC constaint years
    #                       tac =TRUE, #) #, 
    #bndTac=c(0.9,1.1))
    mp=fwd(mp,catch=catch(om)[,ac(iYr)])
    
    hcrOutcome=hcrAlb(mp, hcrPar,
                      hcrYrs=iYr+seq(range["interval"])-1,
                      #hcrYrs=iYr+seq(range["interval"]),
                      refYrs=iYr, #change TAC constaint years
                      stkYrs=iYr, #change TAC constaint years
                      ref1=ifelse(iYr==range["start"], 28000000, NA),
                      tac =TRUE,
                      bndTac=c(0.80,1.20),
                      maxTAC=50000000,
                      minTAC=1) 

    ## TACs for next year (iYtr+1) for n=interval years
    TAC  = hcrOutcome$tac
    
    TAC[]=rep(apply(TAC,6,mean)[drop=T],each=range["interval"])

    TACsave[,dimnames(TAC)$year]=TAC
    
    # Operating Model Projectionfor TAC
    # print("check before error 205")
    
    om =FLash:::fwd(om,catch=TAC, maxF=maxF,sr=eql,sr.residuals=srDev)  
    
    mmm<<-list(om=om,mp=mp,prj=prj)
    #print(omMp(om,mp))
  }
  
  return(list(om=om,mp=mp,prj=prj, tac=TACsave, cpue=cpue))}

hcrAlb=function(object, 
                params=hcrParam(ftar =0.70*refpts(object)['fmsy'],
                                btrig=0.80*refpts(object)['bmsy'],
                                fmin =0.01*refpts(object)['fmsy'],
                                blim =0.40*refpts(object)['bmsy']),
                stkYrs=max(as.numeric(dimnames(stock(object))$year)),
                refYrs=max(as.numeric(dimnames(catch(object))$year)),
                hcrYrs=max(as.numeric(dimnames(stock(object))$year)),
                ref1=NA,
                tac   =TRUE,
                tacMn =TRUE,
                bndF  =NULL, #c(1,Inf),
                bndTac=c(0.8,1.2),
                iaF   =TRUE, 
                iaTac =TRUE, 
                maxF  =2,
                maxTAC=50000000,
                minTAC=1,
                ...) {
  
  print(stkYrs)
  print(hcrYrs)
  print(refYrs)
  
  ## HCR
  dimnames(params)$params=tolower(dimnames(params)$params)
  params=as(params,'FLQuant')  
  #if (blim>=btrig) stop('btrig must be greater than blim')
  #a=(params['ftar']-params['fmin'])*refpts(object)['fmsy']/((refpts(object)['bmsy'])*(params['btrig']-params['blim']))
  #b=params['ftar']*refpts(object)['fmsy']-a*params['btrig']*refpts(object)['bmsy']

# Originally 
  a=(params['ftar']-params['fmin'])/(params['btrig']-params['blim'])
  b=params['ftar']-a*params['btrig']
  
    
  ## Calc F
  # bug
  #val=(SSB%*%a) %+% b
  stk=FLCore::apply(stock(object)[,ac(stkYrs)],6,mean)
  
  rtn=(stk%*%a)  
  rtn=FLCore::sweep(rtn,2:6,b,'+')
  
  fmin=as(params['fmin'],'FLQuant')
  ftar=as(params['ftar'],'FLQuant')
  
  for (i in seq(dims(object)$iter)){
    FLCore::iter(rtn,i)[]=max(FLCore::iter(rtn,i),FLCore::iter(fmin,i))
    FLCore::iter(rtn,i)[]=min(FLCore::iter(rtn,i),FLCore::iter(ftar,i))} 
  
  rtn=window(rtn,end=max(hcrYrs))
  #dimnames(rtn)$year=min(hcrYrs)  
  if (length(hcrYrs)>1){
    rtn=window(rtn,end=max(hcrYrs))
    rtn[,ac(hcrYrs)]=rtn[,dimnames(rtn)$year[1]]}
  
  ### Bounds ##################################################################################
  ## F
  if (!is.null(bndF)){  
    
    ref=FLCore::apply(harvest(object)[,ac(refYrs-1)],6,mean)
    
    rtn[,ac(min(hcrYrs))]=qmax(rtn[,ac(min(hcrYrs))],ref*bndF[1])
    rtn[,ac(min(hcrYrs))]=qmin(rtn[,ac(min(hcrYrs))],ref*bndF[2])
    
    if (length(hcrYrs)>1)        
      for (i in hcrYrs[-1]){
        if (iaF){
          rtn[,ac(i)]=qmax(rtn[,ac(i)],rtn[,ac(i-1)]*bndF[1])
          rtn[,ac(i)]=qmin(rtn[,ac(i)],rtn[,ac(i-1)]*bndF[2])
        }else{
          rtn[,ac(i)]=rtn[,ac(i-1)]}
        
        if (!is.null(maxF)) rtn=qmin(rtn,maxF)}}
  hvt=rtn
  
  ## TAC
  # if (tac){
  #   print("TAC")
  #   
  #   ref=FLCore::apply(catch(object)[,ac(refYrs)],6,mean)
  #   
  #   object=window(object, end=max(as.numeric(hcrYrs)))
  #   object=fwd(object,harvest=harvest(object)[,ac(min(as.numeric(hcrYrs)-1))])
  #   
  #   rtn   = catch(fwd(object, harvest=rtn))[,ac(hcrYrs)]
  #   
  #   if (!is.null(bndTac)){  
  #     rtn[,ac(min(hcrYrs))]=qmax(rtn[,ac(min(hcrYrs))],ref*bndTac[1])
  #     rtn[,ac(min(hcrYrs))]=qmin(rtn[,ac(min(hcrYrs))],ref*bndTac[2])
  #     
  #     print(rtn)
  #     if (length(hcrYrs)>1)        
  #       for (i in hcrYrs[-1]){
  #         if (iaTac){
  #           rtn[,ac(i)]=qmax(rtn[,ac(i)],rtn[,ac(i-1)]*bndTac[1])
  #           rtn[,ac(i)]=qmin(rtn[,ac(i)],rtn[,ac(i-1)]*bndTac[2])
  #         }else{
  #           rtn[,ac(i)]=rtn[,ac(i-1)]}}
  #     print(rtn)          
  #     
  #   if (tacMn) rtn[]=c(apply(rtn,3:6,mean))}
  #   }
  
  if (tac){
    print("TAC")
    
    ## TACs for target F
    object=window(object, end=max(as.numeric(hcrYrs)))
    object=fwd(object,harvest=harvest(object)[,ac(min(as.numeric(hcrYrs)))])
    rtn   = catch(fwd(object, harvest=rtn))[,ac(hcrYrs)]
    
    ## Bounds
   Bcurr=stock(object[,ac(refYrs)])/params['btrig']    ###*refpts(object)$bmsy   ## Replace with btrig
   print(paste("Bcurr is:", Bcurr, "RefYrs:", refYrs) )
   print(paste("params['btrig']:",params['btrig']))
   
   if (Bcurr<1) bndTac = c(0.00001, 3.0)
   if (Bcurr<1) minTac = 0
   
   if (refYrs==2015) bndTac = c(0.9999, 1.00001)
   if (refYrs==2018) bndTac = c(1.19999, 1.20001)
  
   print(paste("Ref years are:",  refYrs))
   
    if (!is.null(bndTac)){
      ## Reference TAC for bounds
      refYrs=refYrs
      ref=FLCore::apply(catch(object)[,ac(refYrs-1)],6,mean)  ## Rep
       #ref=FLCore::apply(catch(object)[,ac(refYrs-1)],6,mean)
      if (!is.na(ref1)) ref[]=ref1   
      print(bndTac)
      rtn=FLQuant(rep(c(apply(rtn,3:6,mean)),each=dim(rtn)[2]),dimnames=dimnames(rtn))
      rtn=qmin(rtn, maxTAC, ref*bndTac[2])
      rtn=qmax(rtn, minTAC, ref*bndTac[1])
    #  rtn=qmin(rtn, refpts(object)$msy*1.2)
    #  rtn=qmax(rtn,ref*bndTac[1])
    #  rtn=qmin(rtn,ref*bndTac[2])
      
    }
  }
  
  if (tac) rtn=list(hvt=hvt,tac=rtn,stock=stk) else rtn=list(hvt=hvt,stock=stk)
  
  print(rtn) 
  
  return(rtn)}


# return(list(om=om,mou=mou,bd=bd,mp=mp,oem=mcf(FLQuants(cpue=cpue,catch=catch(om)))))}

# runMSE_OEM<-function(om,eql,
#                     sa,hcr,
#                    srDev,
#                    range=c(start   =dims(om)$maxyear,
#                            end     =dims(om)$maxyear+30,
#                            interval=3),
#                    maxF=2,
#                    iaF=.1,
#                    cmdOps='-maxfn 500 -iprint 0 -est'){
  
  ## OM projections
  # print("check before error 226")
#  om  =fwdWindow(om,end=range["end"],eql)
  
  ## F in longterm
#  prj=FLQuant(c(FLBRP:::refpts(eql)['msy','harvest']*hcr['ftar']),
#              dimnames=list(year=range["start"]:range["end"],
#                            iter=seq(nits)))
  
  ## Add stochastcity
#  if (!(srDev%in%"FLQuant"))
    #  srDev=rlnorm(nits,FLQuant(0,
    #                            dimnames=list(year=range["start"]:range["end"])),srDev)
    
    # Modified by Gorka to account for recruitment increase 
    
#   rlnorm(nits,FLQuant(log(1.3),
#                       dimnames=list(year=range["start"]:range["end"])),srDev)
  
  # print("check before error 244")
# om =fwd(om,f=prj, sr=eql,sr.residuals=srDev)
# prj=om
  
  ## Get number of iterations in OM
# nits=c(om=dims(om)$iter, eql=dims(params(eql))$iter, rsdl=dims(srDev)$iter)
# if (length(unique(nits))>=2 & !(1 %in% nits)) 
    ("Stop, iters not '1 or n' in OM")
# if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  ## Cut in capacity
# maxF=mean(apply(fbar(window(om,end=range["start"])),6,mean)*maxF)
  
  #### Observation Error (OEM) setup #######################
  ## Random variation for CPUE  
# cpue=oem(window(om,end=range["start"]))
  
  ## use data from last year
#  cpue=window(cpue,end=range["start"]-1)
#  cpue[,ac(range["start"]-(range["interval"]:1))]=oem(om[,ac(range["start"]-(range["interval"]:1))])
  
  #### Management Procedure
  ## Set up assessment parameter options
#  mp@control["p","val"]=1
  
#  mp=biodyn:::FLStock2biodyn(window(om,end=range["start"]-1))
#  pnms=dimnames(control(sa))$param[dimnames(control(sa))$param%in%
#                                     dimnames(params(mp))$params]
#  params(mp)[pnms]=apply(sa@params[pnms],1,median)
  
#  setParams( mp)=cpue 
#  setControl(mp)=params(mp)
  
#  print(mp@params)
#  print(mp@control)
  
  ## Loop round years
#  for (iYr in seq(range["start"],range["end"]-range["interval"],range["interval"])){
#    cat('\n===================', iYr, '===================\n')
    
    ## use data from last year
#    mp=window(mp,end=iYr)
#    cpue=window(cpue,end=iYr-1)
#    cpue[,     ac(iYr-(range["interval"]:1))]=oem(  om[,ac(iYr-(range["interval"]:1))])
#   catch(mp)[,ac(iYr-(range["interval"]:1))]=catch(om[,ac(iYr-(range["interval"]:1))])        
    
    ## fit
#   mmm<<-list(om=om,mp=mp,prj=prj,cpue=cpue)
    #mp =biodyn::fwd(mp,catch=catch(mp))
#   catch(mp)[is.na(catch(mp))]=0.01
#   print("hey")
#   mp =biodyn::fit(mp,cpue,cmdOps=cmdOps)
    # print("check error")
    ## HCR
#   hcrPar=hcrParam(ftar =hcr["ftar" ]*fmsy(mp),
#                   btrig=hcr["btrig"]*bmsy(mp),
#                   fmin =hcr["fmin" ]*fmsy(mp), 
#                   blim =hcr["blim" ]*bmsy(mp))
#   hcrOutcome=biodyn::hcr(mp,hcrPar,
#                          hcrYrs=iYr+seq(range["interval"]),
#                          tac =TRUE)
    
    ## TACs for next year (iYtr+1) for n=interval years
#   TAC  =hcrOutcome$tac
#   TAC[]=rep(apply(TAC,6,mean)[drop=T],each=range["interval"])
    
    #### Operating Model Projectionfor TAC
    # print("check before error 311")
#   om =FLash:::fwd(om,catch=TAC,maxF=maxF,sr=eql,sr.residuals=srDev)  
    
#   mmm<<-list(om=om,mp=mp,prj=prj)
    #print(omMp(om,mp))
# }
  
#  return(list(om=om,mp=mp,prj=prj))}



pMeas=function(stk,brp,proxy="msy"){
  
  res=FLQuants(stock  =stock(stk)%/%FLBRP:::refpts(brp)[proxy,"biomass"],
               ssb    =ssb(  stk)%/%FLBRP:::refpts(brp)[proxy,"ssb"],
               rec    =rec(  stk)%/%FLBRP:::refpts(brp)[proxy,"rec"],
               catch  =catch(stk)%/%FLBRP:::refpts(brp)[proxy,"yield"],
               fbar   =fbar( stk)%/%FLBRP:::refpts(brp)[proxy,"harvest"],
               harvest=(catch(stk)/stock(stk))%/%(FLBRP:::refpts(brp)[proxy,"yield"]/FLBRP:::refpts(brp)[proxy,"biomass"]))
  
  model.frame(res,drop=T)}

pM<-function(om,eql,LRP){
  
  dat=pMeas(om,eql)
  
  res=with(dat, cbind(kobeSmry(stock,harvest),
                      dRate(catch, 0),
                      dRate(catch,0.05),
                      dRate(catch,0.10),
                      diags:::av(catch),
                      diags:::av(harvest),
                      count(subset(dat, year>=2010)$stock>LRP)$freq[2]/length(subset(dat, year>=2010)$stock)
  ))
  
  names(res)[10:15]=c("catch","catch5","catch10","aavCatch","aavF", "pLRP")
  
  res}

utils::globalVariables(c('eql','srDev'))

## MSE function
#' mseBiodyn
#' @description Runs a full MSE using an \code{FLStock} object as the Operating Model and \code{biodyn} as the Mangement Procedure
#'           
#' @aliases mse
#' 
#' @param om an \code{FLStock} object 
#' @param eql an \code{FLBRP} object that holds the biological parameters for use in the projections
#' @param mp an \code{biodyn} object that holds the options for the biomass dynamic assessment model
#' @param range a \code{vector} the starting and end years for the projections, and the interval for running the MP
#' @param srDev  a \code{FLQuant} with recruitment deviates
#' @param uDev an \code{FLQuant} or \code{FLQuants} with CPUE residuals
#' @param ftar a \code{numeric} with target F in HCR
#' @param fmin a \code{numeric} with minimum F in HCR
#' @param blim a \code{numeric} with biomass limit for HCR
#' @param btrig a \code{numeric} with biomass trigger (i.e. when to reduce F) in HCR 
#' @param what a \code{character} that specifies what is to be used for the reference point in the HCR, recycled as required
#' @param mult a \code{logical} that specifies whether quantity in HCR options is a multiplier or probability, recycled as required
#'
#' @return  a list of \code{data.frame}s with performance measures from OM and summaries from MP, if \code{con!=NULL} will
#' also write to a MYSQL database
#'  
#' @export
#' @rdname runMSE
#' 
#' @seealso \code{\link{biodyn}}
#' 
#' @examples
#' \dontrun{
#' sim()
#'    }
mseBiodyn<-function(om,eql,srDev,
                    control,priors,
                    start   =range(om)["maxyear"],          
                    end     =start+30,         
                    interval=3,
                    uDev  =0.2,
                    u     =oem,
                    ftar  =0.70,    btrig =0.60,
                    fmin  =0.01,    blim  =0.01,
                    what  ="msy",
                    mult  =TRUE,
                    bndF  =NULL,
                    maxF    =1.0,     
                    phaseQ=1,
                    cmdOps=paste('-maxfn 500 -iprint 0 -est'),
                    omega =1,
                    refB  =1,
                    trendQ=1+FLQuant(cumprod(rep(0,end-dims(trendQ)$minyear+1)),
                                     dimnames=list(year=dims(trendQ)$minyear:end))){
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eql=dims(params(eql))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in OM")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  ## Cut in capacity
  maxF=mean(apply(fbar(window(om,end=start)),6,max)*maxF)
  
  ## open loop feed forward
  mou=om
  
  #### Observation Error (OEM) setup #######################
  ## Random variation for CPUE  
  cpue=oem(window(om,end=start),uDev)
  #cpue=oem(window(om,end=start))
  cpue=cpue*trendQ[,dimnames(cpue)$year]
  
  ## Loop round years
  mp =NULL
  hcr=NULL
  for (iYr in seq(start,range(om,'maxyear')-interval,interval)){
    #iYr = seq(start+rcvPeriod,range(om,'maxyear')-interval,interval)[1]
    cat('\n===================', iYr, '===================\n')
    
    ## use data from last year
    cpue=window(cpue,end=iYr-1)
    cpue[,ac(iYr-(interval:1))]=oem(om[,ac(iYr-(interval:1))], uDev)
    #cpue[,ac(iYr-(interval:1))]=oem(om[,ac(iYr-(interval:1))])
    cpue[,ac(iYr-(interval:1))]=cpue[,ac(iYr-(interval:1))]*trendQ[,ac(iYr-(interval:1))]
    
    #### Management Procedure
    ## Set up assessment parameter options
    bd=FLStock2biodyn(window(om,end=iYr-1))
    pnms=dimnames(control)$param[dimnames(control)$param%in%dimnames(params(bd))$params]
    params(bd)[pnms]=control[pnms,'val']
    
    bd@priors=priors
    setParams( bd)=cpue 
    setControl(bd)=params(bd)
    bd@control[dimnames(control)$params,'phase'][]=control[dimnames(control)$params,'phase']
    bd@control['q1','phase']=phaseQ
    bd@control['q1','val']  =1
    
    ## fit
    print("check before error 1")
    bd =biodyn::fit(bd,cpue,cmdOps=cmdOps)
    print("check before error 2")
    bd =biodyn::fwd(bd,catch=catch(om)[,ac(iYr)])
    
    
    print(bd@params)
    print(bd@control)
    
    ## HCR
    hcrPar=hcrParam(ftar =ftar *fmsy(bd),
                    btrig=btrig*bmsy(bd),
                    fmin =fmin *fmsy(bd), 
                    blim =blim *bmsy(bd))
    #hcrOutcome=biodyn::hcr(bd,hcrPar,
    #                       hcrYrs=iYr+seq(interval),
    #                       bndF=bndF,
    #                       tac =TRUE)
    
    hcrOutcome=hcrAlb(bd,hcrPar,
                      hcrYrs=iYr+seq(interval),
                      bndF=bndF,
                      tac =TRUE)
    
    ## TACs for next year (iYtr+1) for n=interval years
    TAC  =hcrOutcome$tac
    TAC[]=rep(apply(TAC,6,mean)[drop=T],each=interval)
    
    #### Operating Model Projectionfor TAC
    # print("check before error 468")
    om =FLash::fwd(om,catch=TAC,maxF=maxF,sr=eql,sr.residuals=srDev)  
    
    #### Summary Statistics
    ## HCR actions, i.e. is biomass<Btrig?, what is F?, ..
    hcr =rbind(hcr,data.frame(yearHcr=min(as.numeric(dimnames(hcrOutcome$hvt)$year)),
                              #yearAss=rep(range(bd)[2],dims(bd)$iter),
                              model.frame(           hcrPar,drop=T)[,-5],
                              tac    =as.data.frame(apply(hcrOutcome$tac,6,mean),drop=T)[,'data'],
                              harvest=as.data.frame(apply(hcrOutcome$hvt,6,mean),drop=T)[,'data'],
                              stock  =as.data.frame(hcrOutcome$stock,drop=T)[,2]))
    
    ## Assessment parameters and reference points
    mp =rbind(mp,cbind(cbind(year=iYr,model.frame(params(bd))),
                       model.frame(refpts(bd))[,-4],
                       hcr))
  }
  
  ## save OM, projection without feedback, last assessment and MP summary
  return(list(om=om,mou=mou,bd=bd,mp=mp,oem=mcf(FLQuants(cpue=cpue,catch=catch(om)))))}


hcrFn=function(om,btrig,blim,ftar,fmin,start,end,interval,lag=seq(interval)){    
  
  a=(ftar-fmin)/(btrig-blim)
  b=ftar-a*btrig
  
  for (iYr in seq(start+rcvPeriod,range(om,'maxyear')-interval,interval)){
    stk=FLCore::apply(stock(om)[,ac(iYr-1)],6,mean)
    
    trgt=(stk%*%a)+b  
    
    for (i in seq(dims(om)$iter)){
      FLCore::iter(trgt,i)[]=max(FLCore::iter(trgt,i),FLCore::iter(fmin,i))
      FLCore::iter(trgt,i)[]=min(FLCore::iter(trgt,i),FLCore::iter(ftar,i))} 
    
    dmns     =dimnames(trgt)
    dmns$year=as.character(iYr+lag)
    
    trgt=FLQuant(rep(c(trgt),each=length(lag)),dimnames=dmns)
    
    # print("check before error 509")
    om=fwd(om,f=trgt,sr=eql,sr.residuals=srDev)
  }
  
  return(om)}



demoBiodyn<-function(om,mp,
                     eDev  =0.3,
                     uDev  =0.2,
                     u     =oem,
                     ftar  =0.70,    btrig =0.60,
                     fmin  =0.01,    blim  =0.01,
                     bndF  =NULL,
                     start   =range(mp)["maxyear"],          
                     end     =start+30,         
                     interval=3,
                     maxF    =2.0,     
                     cmdOps=paste('-maxfn 500 -iprint 0 -est')){
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, eql=dims(params(eql))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in OM")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))
  
  ## Cut in capacity
  maxF=mean(apply(fbar(window(om,end=start)),6,max)*maxF)
  
  
  #### Observation Error (OEM) setup #######################
  ## Random variation for CPUE  
  cpue=oem(window(om,end=start), uDev)
  #cpue=oem(window(om,end=start))
  
  ## Loop round years
  mp =NULL
  hcr=NULL
  for (iYr in seq(start,range(om,'maxyear')-interval,interval)){
    #iYr = seq(start+rcvPeriod,range(om,'maxyear')-interval,interval)[1]
    cat('\n===================', iYr, '===================\n')
    
    ## use data from last year
    cpue=window(cpue,end=iYr-1)
    cpue[,ac(iYr-(interval:1))]=oem(om[,ac(iYr-(interval:1))], uDev)
    #cpue[,ac(iYr-(interval:1))]=oem(om[,ac(iYr-(interval:1))])
    
    #### Management Procedure
    ## Set up assessment parameter options
    bd=FLStock2biodyn(window(om,end=iYr-1))
    pnms=dimnames(control)$param[dimnames(control)$param%in%dimnames(params(bd))$params]
    params(bd)[pnms]=control[pnms,'val']
    
    bd@priors=priors
    setParams( bd)=cpue 
    setControl(bd)=params(bd)
    bd@control[dimnames(control)$params,'phase'][]=control[dimnames(control)$params,'phase']
    bd@control['q1','phase']=phaseQ
    bd@control['q1','val']  =1
    
    ## fit
    bd =biodyn::fit(bd,cpue,cmdOps=cmdOps)
    # print("check before error 569")
    bd =biodyn::fwd(bd,catch=catch(om)[,ac(iYr)])
    
    ## HCR
    hcrPar=hcrParam(ftar =ftar *fmsy(bd),
                    btrig=btrig*bmsy(bd),
                    fmin =fmin *fmsy(bd), 
                    blim =blim *bmsy(bd))
    hcrOutcome=biodyn::hcr(bd,hcrPar,
                           hcrYrs=iYr+seq(interval),
                           bndF=bndF,
                           tac =TRUE)
    
    ## TACs for next year (iYtr+1) for n=interval years
    TAC  =hcrOutcome$tac
    TAC[]=rep(apply(TAC,6,mean)[drop=T],each=interval)
    
    #### Operating Model Projectionfor TAC
    # print("check before error 587")
    om =FLash::fwd(om,catch=TAC,maxF=maxF,sr=eql,sr.residuals=srDev)  
    
    #### Summary Statistics
    ## HCR actions, i.e. is biomass<Btrig?, what is F?, ..
    hcr =rbind(hcr,data.frame(yearHcr=min(as.numeric(dimnames(hcrOutcome$hvt)$year)),
                              #yearAss=rep(range(bd)[2],dims(bd)$iter),
                              model.frame(           hcrPar,drop=T)[,-5],
                              tac    =as.data.frame(apply(hcrOutcome$tac,6,mean),drop=T)[,'data'],
                              harvest=as.data.frame(apply(hcrOutcome$hvt,6,mean),drop=T)[,'data'],
                              stock  =as.data.frame(hcrOutcome$stock,drop=T)[,2]))
    
    ## Assessment parameters and reference points
    mp =rbind(mp,cbind(cbind(year=iYr,model.frame(params(bd))),
                       model.frame(refpts(bd))[,-4],
                       hcr))
  }
  
  ## save OM, projection without feedback, last assessment and MP summary
  return(list(om=om,mou=mou,bd=bd,mp=mp,oem=mcf(FLQuants(cpue=cpue,catch=catch(om)))))}


demo<-function(om,mp,pe,
               uDev  =0.3,
               oem   =function(x,uDev){
                 dmns=list(year=dimnames(stock(x))$year[-dim(stock(x))[2]])
                 res =rlnorm(dim(stock(x))[6],
                             FLQuant(0,dimnames=dmns),uDev)
                 
                 res*(stock(x)[,-1]+stock(x)[,-dim(stock(x))[2]])/2},
               ftar    =0.70,    
               btrig   =0.60,
               fmin    =0.01,   
               blim    =0.01,
               bndF    =NULL,
               start   =range(mp)["maxyear"],          
               end     =start+30,
               interval=3,
               maxF    =.75, 
               qKnown  =FALSE,
               cmdOps=paste('-maxfn 500 -iprint 0')){
  
  om=window(om,end=start+1)
  
  ## Get number of iterations in OM
  nits=dims(om)$iter
  if (dims(mp)$iter==1&nits>1) mp=propagate(mp,nits)
  
  ## Cut in capacity
  #maxF=apply(harvest(window(om,end=start)),6,max,na.rm=T)*maxF
  
  #### Observation Error (OEM) setup #######################
  ## Random variation for CPUE  
  cpue=oem(window(om,end=start), uDev)
  #cpue=oem(window(om,end=start))
  setParams( mp)=cpue  
  setControl(mp)=params(mp)  
  
  if (qKnown){
    control(mp)["q1",c("phase","val")]=c(-1,1)
    print(control(mp))
  }   
  
  ## Loop round years
  for (iYr in seq(start,end-interval,interval)){
    #iYr = seq(start+rcvPeriod,range(om,'maxyear')-interval,interval)[1]
    cat('\n===================', iYr, '===================\n')
    
    ## use data from last year
    cpue=window(cpue,end=iYr-1)
    cpue[,ac(iYr-(interval:1))]     =oem(om[,ac(iYr-(interval:0))],uDev)
    #cpue[,ac(iYr-(interval:1))]     =oem(om[,ac(iYr-(interval:0))])
    mp=window(mp,end=iYr-1)
    catch(mp)=catch(om)
    #### Management Procedure
    
    ## fit
    mp =fit(window(mp,end=iYr-1),cpue,
            cmdOps=ifelse(nits>1,cmdOps,'-maxfn 500 -iprint 0'))
    # print("check before error 664")
    mp =biodyn::fwd(mp,catch=catch(om)[,ac(iYr)],maxF=maxF)
    
    ## HCR
    hcrPar=hcrParam(ftar =ftar *fmsy(mp),
                    btrig=btrig*bmsy(mp),
                    fmin =fmin *fmsy(mp), 
                    blim =blim *bmsy(mp))
    hcrOutcome=biodyn::hcr(mp,hcrPar,
                           hcrYrs=iYr+seq(interval),
                           bndF=bndF,
                           tac =TRUE,
                           maxF=maxF)
    TAC  =hcrOutcome$tac
    TAC[]=rep(apply(TAC,6,mean)[drop=T],each=interval)
    TAC[is.na(TAC)]=10
    TAC=qmax(TAC,10)
    TAC=qmin(TAC,bmsy(mp)*fmsy(mp))
    
    #### Operating Model Projectionfor TAC
    # print("check before error 684")
    om =fwd(om,catch=TAC,pe=pe,maxF=maxF)
    
    #print((1:100)[is.na(stock(om)[,ac(iYr)])])
  }
  
  biodyns(om=om,mp=mp)}

omMp=function(om,sa,end=range(sa)["maxyear"]){
  ggplot(rbind.fill(cbind(What="MP",plot(FLQuants(sa,"stock"))$data),
                    cbind(What="OM",plot(FLQuants(om,"stock"))$data)))+
    geom_ribbon(aes(year,min=`10%`,max=`90%`,fill=What),alpha=.5)+
    geom_line(aes(year,`50%`,col=What))+
    facet_wrap(~qname,ncol=2)+
    scale_x_continuous(limits=c(1975,end))}



