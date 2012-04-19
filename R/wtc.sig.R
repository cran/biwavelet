wtc.sig <-
function (nrands=300, lag1, dt, pad=TRUE, dj=1/12, s0, J1, 
          mother=c("morlet", "paul", "dog")) {
  
  mothers=c("morlet", "paul", "dog")
  mother=match.arg(tolower(mother), mothers)
  
  ## choose a n so that largest scale have atleast some part outside the coi
  ms = s0 * (2^(J1 * dj)) / dt ## maxscale in units of samples
  n = ceiling(ms * 6)
  d1 = arima.sim(model = list(ar = lag1[1], ma = 0), n = n)
  d1 = d1 - mean(d1)
  wt1 = wt(d = d1, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1, mother = mother)
  outside.coi = matrix(0, nrow = NROW(wt1$wave), ncol = NCOL(wt1$wave))
  
  for (s in 1:length(wt1$scale)) {
    outside.coi[s, ] = wt1$period[s] <= wt1$coi
  }
  s.inv=1/t(wt1$scale)
  s.inv=matrix(rep(s.inv, n), nrow=NROW(wt1$wave))

  if (nrands < 1) {
    wtcsig=cbind(scale, 0.71)
  }
  else {
    warn=FALSE
    sig95=rep(0, length(wt1$scale))
    max.scale=1
    
    for (s in 1:length(wt1$scale)) {
      if (any(outside.coi[s, ] > 0))
        max.scale=s
      else {
        sig95[s]=NA
        if (!warn) {
          warning('Long wavelengths completely influenced by COI. 
                  Set n higher or J1 lower') 
          warn=TRUE
        }
      }
    }
    nbins=1000
    wlc=matrix(nrow=length(wt1$scale), ncol=nbins, 0)
    prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    for (r in 1:nrands) {
      ## Generate time series
      d1 = arima.sim(model = list(ar = lag1[1], ma = 0), n = n)
      d1 = d1 - mean(d1)
      d2 = arima.sim(model = list(ar = lag1[2], ma = 0), n = n)
      d2 = d2 - mean(d2)
      ## Wavelet transforms
      wt1 = wt(d = d1, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1, 
               mother = mother, do.sig=FALSE)
      wt2 = wt(d = d2, pad = pad, dj = dj, dt = dt, s0 = s0, J1 = J1, 
               mother = mother, do.sig=FALSE)
      ## Smoothed cross wavelet transform
      smooth.CW=smooth.wavelet(s.inv * (wt1$wave * Conj(wt2$wave)), dt, dj, wt1$scale)
      rsq=abs(smooth.CW)^2/
        (smooth.wavelet(s.inv * (abs(wt1$wave)^2), dt , dj, wt1$scale) *
        smooth.wavelet(s.inv * (abs(wt2$wave)^2), dt, dj, wt1$scale))     
      for (s in 1:max.scale) {
        cd=rsq[s, which(outside.coi[s, ] == 1)]
        cd=pmax(pmin(cd, 1), 0)
        cd=floor(cd * (nbins - 1)) + 1
        for (jj in 1:length(cd)) {
          wlc[s, cd[jj]] = wlc[s, cd[jj]] + 1
        }
      }
      setTxtProgressBar(prog.bar, r)
    }
    for (s in (1:max.scale)) {
      rsqy=(seq(1, nbins, 1) - 0.5)/nbins
      ptile=wlc[s, ]
      ind=which(ptile != 0)
      ptile=ptile[ind]
      rsqy=rsqy[ind]
      ptile=cumsum(ptile)
      ptile=(ptile - 0.5) / ptile[length(ptile)]
      sig95[s]=approx(x=ptile, y=rsqy, 0.95)$y
    }
    wtcsig=cbind(wt1$scale, sig95)
    
    if (any(is.na(sig95)) & !warn)
      warning('Monte Carlo randomization failed')
    
    return (wtcsig)
  }
}
