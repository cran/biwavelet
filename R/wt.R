wt <-
function (d, pad=TRUE, dt=NULL, dj=1/12, s0=2*dt, J1=NULL, max.scale=NULL, 
          mother=c("morlet", "paul", "dog"), param=-1, lag1=NULL, 
          sig.level=0.95, sig.test=0, do.sig=TRUE) {

  mothers=c("morlet", "paul", "dog")
  mother=match.arg(tolower(mother), mothers)

  if (is.null(dt) & NCOL(d) > 1) {
    if (class(d[,1])=="Date") {
      t = 1:NROW(d)
      dt = 1
    }
    else {
      dt = diff(d[, 1])[1]
      t = d[, 1]
    }
    xaxis = d[, 1]
    x = d[, 2] - mean(d[, 2])
  }
  else {
    x = d - mean(d)
    t = 1:length(x)
    xaxis = t
  }  
  n.obs = length(x)
  if (is.null(J1)) {
    if (is.null(max.scale)) {
      max.scale=(n.obs * 0.17) * 2 * dt ## automaxscale
    }
    J1=round(log2(max.scale / s0) / dj)
  }
  if (pad) {
    base2 = floor(log(n.obs) / log(2) + 0.5)
    x = c(x, rep(0, times=2^(base2 + 1) - n.obs))
  }
  n=length(x)
  k = seq(1, floor(n / 2), 1)
  k = k * ((2 * pi)/(n * dt))
  k = c(0, k, -k[seq(floor(( n - 1) / 2), 1, -1)])
  f = fft(x)
  scale = s0 * 2^(seq(0, J1, 1)*dj)
  period = scale
  wave = matrix(0, nrow=J1 + 1, ncol=n)
  wave = wave + 1i*wave
  
  for (a1 in seq(1, J1 + 1, 1)) {
    wb=wt.bases(mother, k, scale[a1], param)
    wave[a1, ]=fft(f * wb$daughter, inverse=TRUE)/length(f)
  }
  period = wb$fourier.factor * scale
  coi = wb$coi * dt * c(1e-5, seq(1, (n.obs + 1) / 2 - 1), 
                        seq(floor(n.obs / 2 - 1), 1, -1),
                        1e-5)
  wave = wave[, 1:n.obs] ## Get rid of padding before returning
  power = abs(wave)^2
  phase=atan2(Im(wave), Re(wave))
  sigma2 = var(x)
  if (do.sig) {
    signif = wt.sig (d=x, dt=dt, scale=scale, sig.test=sig.test, 
                     sig.level=sig.level, lag1=lag1, dof=-1, 
                     mother=mother, sigma2=1)$signif
    signif=matrix(signif, nrow=length(signif), ncol=1) %*% rep(1, n.obs)
    signif=power / (sigma2*signif)
  }
  else
    signif=NA
  
  results=list(coi=coi, 
               wave=wave,
               power=power,
               phase=phase,
               period=period,
               scale=scale, 
               dt=dt, 
               t=t,
               xaxis=xaxis,
               s0=s0, 
               dj=dj, 
               sigma2=sigma2,
               mother=mother,
               type = "wt", 
               signif=signif)
  class(results)="biwavelet"
  return (results)
}
