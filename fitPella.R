fitPella <- function (object, index = index, exeNm = "pella", package = "biodyn", 
          dir = tempdir(), set = biodyn:::setPella, get = getPella, cmdOps = paste("-maxfn 500 -iprint 0")) 
{
  first = TRUE
  
  slts = getSlots("biodyn")
  slts = slts[slts %in% c("FLPar", "FLQuant")]
  oldwd = getwd()
  setwd(dir)
  biodyn:::exe()
  object = list(object, index)
  bd = object[[1]]
  its = max(maply(names(slts), function(x) dims(slot(bd, x))$iter))
  its = max(its, dims(bd@control)$iter)
  nms = dimnames(params(bd))$params
  bd@vcov = FLPar(array(as.numeric(NA), dim = c(length(nms), 
                                                length(nms), its), dimnames = list(params = nms, params = nms, 
                                                                                   iter = seq(its))))
  bd@hessian = bd@vcov
  us = paste("u", seq(length(dimnames(params(bd))$params[grep("q", 
                                                              dimnames(params(bd))$params)])), sep = "")
  bd@ll = FLPar(NA, dimnames = list(params = us, iter = seq(1)))
  if (its > 1) {
    bd@stock = FLCore::iter(bd@stock, 1)
    bd@params = FLCore::iter(bd@params, 1)
    bd@objFn = FLCore::iter(bd@objFn, 1)
    bd@vcov = FLCore::iter(bd@vcov, 1)
    bd@ll = FLCore::iter(bd@ll, 1)
    bd@hessian = FLCore::iter(bd@hessian, 1)
    bd@mng = FLPar(a = 1)
    bd@mngVcov = FLPar(a = 1, a = 1)
    if (dim(bd@stock)[6] == 1) 
      bd@stock = propagate(bd@stock, iter = its, fill.iter = TRUE)
    if (dim(bd@catch)[6] == 1) 
      bd@stock = propagate(bd@catch, iter = its, fill.iter = TRUE)
    pnms <- getSlots(class(bd))
    pnames <- names(pnms)[pnms == "FLPar"]
    for (i in pnames) {
      slot(bd, i) = FLCore::iter(slot(bd, i), 1)
      slot(bd, i) <- propagate(slot(bd, i), its)
    }
  }
  cpue = object[[2]]
  bd2 = object[[1]]
  for (i in seq(its)) {
    object[[2]] = FLCore::iter(cpue, i)
    for (s in names(slts)[-(7:8)]) {
      slot(object[[1]], s) = FLCore::iter(slot(bd2, s), 
                                          i)
    }
    object[[1]] = set(object, exeNm, dir)
    system(paste(exeNm, " ", cmdOps, sep = ""))
    object[[1]] = biodyn:::getPella(object[[1]], exeNm)
    for (s in names(slts)[slts == "FLQuant"]) {
      FLCore::iter(slot(bd, s), i) = slot(object[[1]], 
                                          s)
    }
    if (its <= 1) {
      x <- file(paste(dir, "admodel.hes", sep = "/"), 
                "rb")
      nopar <- readBin(x, "integer", 1)
      H <- matrix(readBin(x, "numeric", nopar * nopar), 
                  nopar)
      try(bd@hessian@.Data[activeParams(object[[1]]), 
                           activeParams(object[[1]]), i] <- H, silent = TRUE)
      close(x)
      if (file.exists(paste(dir, "admodel.cov", sep = "/"))) 
        try(bd@vcov@.Data[activeParams(object[[1]]), 
                          activeParams(object[[1]]), i] <- cv(paste(dir, 
                                                                    "admodel.hes", sep = "/")), silent = TRUE)
      if (file.exists(paste(dir, "pella.hst", sep = "/"))) 
        bd@hst = admbProfile(paste(dir, "pella.hst", 
                                   sep = "/"))$profile
      if (file.exists(paste(dir, "lpr.plt", sep = "/"))) 
        bd@profile = mdply(data.frame(var = c("r", "k", 
                                              "bnow", "fnow", "msy", "bmsy", "fmsy", "cmsy", 
                                              "bmsy", "ffmsy", "bk", "fr", "bratio", "fratio", 
                                              "slopeb", "slopef")), function(var) {
                                                fl = paste("lp", var, ".plt", sep = "")
                                                if (file.exists(fl)) 
                                                  admbPlt(fl)
                                              })
    }
    bd@params@.Data[, i] = object[[1]]@params
    bd@control@.Data[, , i] = object[[1]]@control
    bd@ll@.Data[, i][] = unlist(c(object[[1]]@ll))
    if (file.exists("pella.std")) {
      err1 = try(mng. <- read.table("pella.std", header = T)[, 
                                                             -1])
      err2 = try(mngVcov. <- fitFn(paste(dir, "pella", 
                                         sep = "/"))$vcov)
      if (first) {
        if (any(is(err1) != "try-error")) 
          bd@mng = FLPar(array(unlist(c(mng.[, -1])), 
                               dim = c(dim(mng.)[1], 2, its), dimnames = list(param = mng.[, 
                                                                                           1], var = c("hat", "sd"), iter = seq(its))))
        if (any(is(err2) != "try-error")) 
          bd@mngVcov <- FLPar(array(unlist(c(mngVcov.)), 
                                    dim = c(dim(mng.)[1], dim(mng.)[1], its), 
                                    dimnames = list(param = dimnames(mngVcov.)[[1]], 
                                                    param = dimnames(mngVcov.)[[1]], iter = seq(its))))
        first = !first
      }
      else {
        try(if (is(err1) != "try-error") 
          bd@mng@.Data[, , i][] = unlist(c(mng.[, -1])))
        try(if (is(err2) != "try-error") 
          bd@mngVcov@.Data[, , i][] = unlist(c(mngVcov.)))
      }
    }
  }
  units(bd@mng) = "NA"
  bd = fwd(bd, catch = catch(bd)[, rev(dimnames(catch(bd))$year)[1]])
  if (length(grep("-mcmc", cmdOps)) > 0 & length(grep("-mcsave", 
                                                      cmdOps)) > 0) {
    setMCMC = function(obj, dir) {
      ps = read.psv(paste(dir, "pella.psv", sep = "/"))
      dmns = list(params = activeParams(obj), iter = seq(dim(ps)[1]))
      ps = array(t(ps), dim = unlist(llply(dmns, length)), 
                 dimnames = dmns)
      ps = FLPar(ps)
      units(ps) = "NA"
      ps
    }
    par = setMCMC(bd, dir)
    cmd = strsplit(cmdOps, ",")
    grp = unlist(gregexpr("-mcmc", cmd[[1]]))
    mcmc = sub(" +", "", cmd[[1]][grp > 0])
    mcmc = as.numeric(substr(mcmc, 6, nchar(mcmc)))
    grp = unlist(gregexpr("-mcsave", cmd[[1]]))
    mcsave = sub(" +", "", cmd[[1]][grp > 0])
    mcsave = sub(" +", "", mcsave)
    mcsave = as.numeric(substr(mcsave, 8, nchar(mcsave)))
    bd@params = propagate(bd@params[, 1], dims(par)$iter)
    bd@params[dims(par)$params, ] = par
    bd@stock = propagate(bd@stock, dim(params(bd))[2])
    bd = fwd(bd, catch = catch(bd))
    attributes(bd@params)[["mcmc"]] = mcmc
    attributes(bd@params)[["mcsave"]] = mcsave
  }
  if (its <= 1) 
    bd@diags = getDiags()
  setwd(oldwd)
  return(bd)
}