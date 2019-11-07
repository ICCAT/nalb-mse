
setGeneric('oem',    function(object,...) standardGeneric('oem'))

setMethod( 'oem',   signature(object='FLStock'),
           function(object,
                    # sel       =FLQuant(FLQuant(1,dimnames=dimnames(harvest(object)))),
                    sel=FLQuant(subset(biol, .id==names(om) & Quantity=="Selectivity")$data, 
                                dimnames=dimnames(harvest(object))),
                    bias      =FLPar(q=0),
                    timing    =0.5,
                    fish.dependent=TRUE,
                    effort    =c("f","h"),
                    mass      =TRUE){
             
             timing=pmax(pmin(timing,1.0),0.0)
             stock=(stock(object)[,-dim(stock(object))[2]]*timing+stock(object)[,-1]*(1.0-timing))
             
             object=window(object,start=range(object)["minyear"],end=range(object)["maxyear"])
             
             if (fish.dependent) {
               if (effort[1]=="h")
                 E=catch(object)%/%stock
               else  
                 E=fbar(object)
               
               if (mass){
                 cpue=apply((catch.n(object)%*%catch.wt(object)%*%sel)%/%E,2:6,sum)
                                }
               else
                 cpue=apply((catch.n(object)%*%sel)%/%E,2:6,sum)
             }
             
             
             
             
             else 
               if (mass){
                 cpue=apply((stock.n(object)%*%stock.wt(object)%*%sel),2:6,sum)
                 
                                              }
             else  
               cpue=apply((stock.n(object)%*%sel),2:6,sum)
             
             #HYPERSTABILITY
             if (all(c("ref","omega")%in%dimnames(bias)$params))
               cpue=cpue%*%exp(log(stock(object)%/%bias["ref"])%*%(bias["omega"]))
             
             if ("q"%in%dimnames(bias)$params){
               cpue=cpue%*%FLQuant(cumprod(1+rep(c(bias["q"]),dim(fbar(object))[2])), dimnames=dimnames(fbar(object)))
             }
             
             cpue})