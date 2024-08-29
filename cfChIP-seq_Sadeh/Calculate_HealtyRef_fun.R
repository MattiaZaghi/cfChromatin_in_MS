EstimateScale = function(X, c, b, m, v, w = NULL)  {
  if( is.null(w) )
    w = rep(1, length(X))
  
  z = (b/c)
  z2 = z**2
  
  ll = function(x) {
    m = x*m
    v = (x**2) * v
    mu = c*m+b 
    r = m**2/v+2*z*m/v+z2/v 
    -sum(w*dnbinom(X,mu=mu, size=r, log = TRUE ))
  }
  
  Y.guess = (X-b)/c
  fit = lm(Y.guess ~ m - 1)
  
  #  print(Y.guess) 
  #  print(fit)
  opt = optimize(ll,interval = c(1E-10, 1E10))
  #  print(opt)
  c(scale = opt$minimum, ll = -opt$objective)
}

RobustEstimateScale = function(X, c, b, m, v, w = NULL) {
  if( is.null(w) )
    w = rep(1, length(X))
  xx = sapply(1:100, function(i) {
    I = sample(length(X), floor(0.75*length(X)))
    EstimateScale(X[I], c[I], b[I], m[I], v[I], w[I])[1]
  })
  #  print(sort(xx))
  median(xx)
}


#
# quad fails in examples such as ExpectedPosteriorGammaPoisson(10,.05, 1, 100, 10000, quad = TRUE, K = 100)
#
ExpectedPosteriorGammaPoisson = function(Y, c, b, m, v, K = 100, quad = FALSE) {
  alpha = m**2/v
  beta = m/v
  alpha1 = alpha+Y
  beta1 = beta+c
  
  if( qgamma(1-1/(2*K), shape = alpha1, rate = beta1) == 0 ) {
    xs = c(0,1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1)
  } else {
    if(!quad) {
      qs = c( 
        seq(0, 0.1/K, 0.01/K)+0.005/K,
        seq(0, 1/K, 0.1/K)+0.05/K,
        seq(1/K, 0.5, 1/K)+0.5/K)
      qs = c(qs, 1-qs)
      qs = qs[qs >= 0 & qs <= 1]
      qs = unique(sort(qs))
      xs = qgamma(qs, shape = alpha1, rate = beta1)
      
    }
    if(quad) {
      gq = gauss.quad(K, "laguerre", alpha1/beta1)
      xs = gq$nodes/beta1
    }
  }
  # to avoid overshoot
  while( min(xs) > 1 )
    xs = c(min(xs)/sqrt(2), xs, max(xs)*sqrt(2))
  
  xs.mid = (xs[-1]+xs[-length(xs)])/2
  ps = pgamma(xs.mid, shape = alpha1, rate = beta1, log = TRUE)
  v1 = c(ps,0)
  v0 = c(-Inf, ps)
  #  ys.old = exp(v1) - exp(v0)
  ys.log = v1+log(1 - exp(v0-v1))
  if( Y == 0 ){
    ws = exp(ys.log)
  } else {
    ws = log(1+b/(c*xs))*Y+ys.log
    ws = exp(ws-max(ws))
  }
  ws = ws/sum(ws)
  c(sum(xs*ws), sum((xs**2) * ws))
}

EstimateMeanVar = function(X, c, b, w = NULL, n0 = .1, p0 = 1, prior.w = 1) {
  
  #  p0 = 1  # product
  q0 = n0  # sum
  w0 = 0.1
  if( is.null(w) )
    w = rep(1, length(X))
  
  z = (b/c)
  z2 = z**2
  
  ll = function(x) {
    m = x[1]
    v = x[2]
    alpha = m**2/v
    beta = m/v
    mu = c*m+b  # c*alpha/beta+b
    r = m**2/v+2*z*m/v+z2/v #alpha+2*b*beta/c+(b*beta/c)**2/alpha
    -sum(w*dnbinom(X,mu=mu, size=r, log = TRUE )) -
      prior.w*((alpha-1)*log(p0)-beta*q0-n0*lgamma(alpha)+alpha*n0*log(beta)) #+ w0*(m + 1/m)
    #0
  }
  
  test.prior = function(m,v) {
    alpha = m**2/v
    beta = m/v
    -((alpha-1)*log(p0)-beta*q0-n0*lgamma(alpha)+alpha*n0*log(beta))
  }
  
  test.prior2 = function(m) {
    opt = optimize(function(v) test.prior(m,v), c(1E-4,1E10))
    opt$objective
  }
  
  ll.m = function(m) {
    opt = optimize(function(v) ll(c(m,v)), c(1E-4,1E10))
    opt$objective
  }
  
  c.mean = mean(c)
  b.mean = mean(b)
  Y.guess = pmax(0,(X-b)/c)
  mu.emp = max(mean(Y.guess),0.01)
  var.emp = max(var(Y.guess), 0.01)
  opt = optim(c(mu.emp,var.emp), ll)
  
  if( opt$par[2] < 1E-3 || opt$par[1] < 1E-3 ) {
    opt = optimize(ll.m,  c(1E-5,1E10))
    m = opt$minimum
    opt = optimize(function(v) ll(c(m,v)), c(1E-4,1E10))
    v = opt$minimum
  } else {
    m = opt$par[1]
    v = opt$par[2]
  }
  
  if(0) {
    Y = sapply(1:length(X), function(k) ExpectedPosteriorGammaPoisson(X[k], c[k], b[k], m, v))
    meanP = mean(Y[1,])
    varP = mean(Y[2,]) - meanP**2
    #need to be fixed
    if( is.na(meanP) ) {
      meanP = m
      varP = v
    }
  } else {
    meanP = m
    varP = v
  }
  
  c(mean = m, var = v, ll = -ll(c(m,v)), meanP = meanP, varP = varP)
}






cfChIP.EstimateMeanBasis = function(SigCounts, SigBackground, QQnorm) {
  SigDiff = SigCounts - SigBackground
  SigDiff[SigDiff < 0] = 0
  
  Sig.avg = rowMeans(SigDiff)    
  I = Sig.avg > 0
  J = QQnorm > 0
  Sig.avg[I] = rowMeans(t(t(SigDiff[I,J])*QQnorm[J]))
  
  Sig.avg
}

cfChIP.EstimateMeanVarianceBasis = function(SigCounts, SigBackground, QQnorm) {
  SigDiff = SigCounts - SigBackground
  SigDiff[SigDiff < 0] = 0
  
  Sig.avg = rowMeans(SigDiff)    
  Sig.var = rowMeans((SigDiff - Sig.avg)**2)
  names(Sig.avg) = rownames(SigCounts)
  names(Sig.var) = rownames(SigCounts)
  I = which(Sig.avg > 0)
  J = QQnorm > 0
  if( length(I) > 0) {
    Sig.est = sapply(I, function(w) EstimateMeanVar(SigCounts[w,J], 1/QQnorm[J], SigBackground[w,J]))
    Sig.avg[I] = Sig.est["meanP",]
    Sig.var[I] = Sig.est["varP",]
  }
  list( avg = Sig.avg, var = Sig.var )
}

cfChIP.EstimateMeanVarianceWinSig = function(LL, WinSig) {
  
  WinCounts = do.call("cbind", lapply(LL, function(l) sapply(WinSig, function(w) sum(l$Counts[w]))))
  WinBackground = do.call("cbind", lapply(LL, function(l) sapply(WinSig, function(w) sum(l$WinBackground[w]))))
  QQnorm = sapply(LL, function(l) l$QQNorm)
  names(QQnorm) = colnames(WinCounts)
  
  cfChIP.EstimateMeanVarianceBasis(WinCounts, WinBackground, QQnorm)
}

cfChIP.EstimateMeanVarianceGeneSig = function(LL, GeneSig ) {
  GeneCounts = do.call("cbind", lapply(LL, function(l) sapply(GeneSig, function(w) sum(l$GeneCounts[w]))))
  GeneBackground = do.call("cbind", lapply(LL, function(l) sapply(GeneSig, function(w) sum(l$GeneBackground[w]))))
  QQnorm = sapply(LL, function(l) l$QQNorm)
  names(QQnorm) = colnames(GeneCounts)
  
  cfChIP.EstimateMeanVarianceBasis(GeneCounts, GeneBackground, QQnorm)
}



fitNoise = function(X, thresh = 0.95, MinNumber = 50, Prior = NA, PriorStrength = 100) {
  X = X[!is.na(X)]
  if(length(X) < MinNumber) 
    return(NA)
  t = quantile(X,thresh)
  X = X[X <= t]
  T = table(X)
  N = as.integer(names(T))
  M = sum(T)
  
  alpha = 1
  beta = 0
  if( !is.na(Prior) ) {
    beta = PriorStrength
    alpha = Prior*PriorStrength
  }
  
  LL = function(l)  {
    sum( T*dpois(N,l,log = TRUE)) - M*ppois(max(N),l,log.p = TRUE) + (alpha-1)*log(l) - beta*l
  }
  m = mean(X)
  up = max(0.1,m*10)
  down = max(0.0001,m/10)
  op = optimize(LL, c(up,down),maximum = TRUE)
  op$maximum
}

fitNoise.nbin = function(X,thresh = 0.95) {
  t = quantile(X,thresh,na.rm=TRUE)
  X = X[X <= t]
  fit = fitdist(X,"nbinom")
  coef(fit)
}

fitNoise.mean = function(X,thresh = 0.95) {
  t = quantile(X,thresh)
  X = X[X <= t]
  mean(X)
}


if( !exists("Background.Regions")) {
  print("Loading Background model")
  BackgroundModel.filename = paste0(SetupDIR,"BackgroundModel.rds")
  L = readRDS(BackgroundModel.filename)
  Background.Regions.width = L$Background.Regions.width 
  Background.Regions = L$Background.Regions
  Background.Regions.num = L$Background.Regions.num
  Background.Uniq.Regions = L$Background.Uniq.Regions
  Background.Inds.width = L$Background.Inds.width
  Background.Inds = L$Background.Inds 
  Background.Inds.num = L$Background.Inds.num 
  Background.Uniq.Inds = L$Background.Uniq.Inds
  Background.RegionChr = L$Background.RegionChr
  Background.IndRegion = L$Background.IndRegion
  Background.WindowInd = L$Background.WindowInd
  Background.WindowRegion = L$Background.WindowRegion
  WinChr = L$WinChr
  Background.RegionWindows = L$Background.RegionWindows
  Background.IndWindows = L$Background.IndWindows
  rm(L)
}


Background.Ind.Strech = 
  Background.windows = TSS.windows$type == "background" & width(TSS.windows) > 4000

Background.windows.width = 5

if (TargetMod == "H3K4me3-scer") { # TODO: set a cutoff for human/yeast in params
  Background.windows = TSS.windows$type == "background" & width(TSS.windows) > 400
  Background.windows.width = 2
}

Edge.Penalty = 0




getBackground = function( mu, w ) {
  mu.genome = mu$genome
  W = TSS.windows[w]
  chr = as.character(chrom(W))
  mu.chr = mu$chr[chr]
  
  mu.region = mu$region.uniq[as.numeric(Background.WindowRegion[w])]
  
  mu.ind = mu$ind.uniq[as.numeric(Background.WindowInd[w])]
  
  len = (width(W)+Edge.Penalty)/1000
  c( ind = mu.ind*len, region = mu.region*len, chr = mu.chr*len, genome = mu.genome*len)
} 



getMultiBackground = function( mu, w ) {
  mu.genome = mu$genome
  W = TSS.windows[w]
  chr = as.character(chrom(W))
  mu.chr = mu$chr[chr]
  
  mu.region = mu$region.uniq[as.numeric(Background.WindowRegion[w])]
  
  mu.ind = mu$ind.uniq[as.numeric(Background.WindowInd[w])]
  
  len = (width(W)+Edge.Penalty)/1000
  rbind( ind = mu.ind*len, region = mu.region*len, chr = mu.chr*len, genome = mu.genome*len)
} 

getMultiBackgroundEstimate = function( mu, w ) {
  W = TSS.windows[w]
  mu.ind = mu$ind.uniq[as.numeric(Background.WindowInd[w])]
  
  len = (width(W)+Edge.Penalty)/1000
  mu.ind*len
} 

getBackgroundParallel = function( MUs, w, terse=TRUE ) {
  W = TSS.windows[w]
  len = (width(W)+Edge.Penalty)/1000
  
  i = as.numeric(Background.WindowInd[w])
  if(!terse) {
    r = as.numeric(Background.WindowRegion[w])
    chr = as.character(chrom(W))
    sapply( MUs, function(mu) {
      mu.genome = mu$genome
      mu.chr = mu$chr[chr]
      mu.region = mu$region.uniq[r]
      
      mu.ind = mu$ind.uniq[i]
      
      c( ind = mu.ind*len, region = mu.region*len, chr = mu.chr*len, genome = mu.genome*len)
    }) 
  } else
    sapply( MUs, function(mu) mu$ind.uniq[i])*len
} 

sigBackground = function(s,sig, B = Background) {
  mu = B[[s]]
  sapply(sig, function(w) getBackground(mu,w)[1])
}

sigBackgroundParallel.Jump = 1000
sigBackgroundParallel = function(sig, B = Background) {
  doWindow = function(w) getBackgroundParallel(B,w, terse = TRUE)
  if( length(sig) <= sigBackgroundParallel.Jump ) { 
    sapply(sig, doWindow)
  } else {
    breaks = c(seq(0,length(sig)-1, by = sigBackgroundParallel.Jump),length(sig))
    do.call(cbind,
            lapply(1:(length(breaks)-1), function(i) {
              print(breaks[i])
              sapply(sig[(breaks[i]+1):breaks[i+1]],doWindow)
            }))
  }
}
