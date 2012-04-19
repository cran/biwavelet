plot.biwavelet <-
  function (x, ncol=64, xlab="Time", ylab="Period", main=NULL, 
            sig.level=0.95, plot.cb=FALSE, plot.phase=FALSE,
            type=c("power.norm", "power", "wavelet", "phase"), 
            plot.coi=TRUE, plot.sig=TRUE, bw=FALSE,
            legend.loc=c(0.885, 0.92, 0.2, 0.8), 
            legend.horiz=FALSE,
            arrow.size=0.08,
            arrow.lwd=2, arrow.cutoff=0.9, 
            form='%Y', ...) {
    if (bw) {
      bw.colors <- colorRampPalette(c("black", "white"))
      fill.colors=bw.colors(ncol)
    }
    else {
      jet.colors <-
        colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                           "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      fill.colors=jet.colors(ncol)
    }
    y.ticks = 2^(floor(log2(min(x$period))):floor(log2(max(x$period))))
    types=c("power.norm", "power", "wavelet", "phase")
    type=match.arg(tolower(type), types)
    
    if (type=="power.norm") {
      if (x$type == "xwt") {
        zvals=log2(abs(x$wave / (x$d1.sigma * x$d2.sigma)))
        zlims=c(-1, 1) * max(zvals)
        zvals [zvals < zlims[1]]=zlims[1]
        locs=pretty(range(zvals), n=5)
        leg.lab=2^locs
      }
      else if (x$type == "wtc") {
        zvals=x$rsq
        zlims=range(zvals)
        zvals [zvals < zlims[1]]=zlims[1]
        locs=pretty(range(zvals), n=5)
        leg.lab=locs
      }
      else {
        zvals=log2(abs(x$power / x$sigma2))
        zlims=c(-1, 1) * max(zvals)
        zvals [zvals < zlims[1]]=zlims[1]  
        locs=pretty(range(zvals), n=5)
        leg.lab=2^locs
      }
    }
    else if (type=="power") {
      zvals=log2(x$power)
      zlims=c(-1, 1)*max(zvals)
      zvals [zvals < zlims[1]]=zlims[1]
      locs=pretty(range(zvals), n=5)
      leg.lab=2^locs
    }
    else if (type=="wavelet") {
      zvals=(Re(x$wave))
      zlims=range(zvals)
      locs=pretty(range(zvals), n=5)
      leg.lab=locs
    }
    else if (type=="phase") {
      zvals=x$phase
      zlims=c(-pi, pi)
      locs=pretty(range(zvals), n=5)
      leg.lab=locs        
    }
    else {
      stop("type must be power, power.norm, wavelet or phase")
    }
    yvals=log2(x$period)
    image(x$t,
          yvals, 
          t(zvals), 
          zlim=zlims,
          ylim=rev(range(yvals)),
          xlab=xlab, 
          ylab=ylab,
          main=main,
          yaxt="n",
          xaxt="n",
          col=fill.colors)
    box()
    
    if (class(x$xaxis) == "Date") {
      xlocs=pretty(x$t)+1
      axis(side=1, at=xlocs, labels=format(x$xaxis[xlocs], form))
    }
    else {
      xlocs=pretty(x$t)
      axis(side=1, at=xlocs)
    }
    axis.locs=axTicks(2)
    axis(2, at=axis.locs, labels=2^axis.locs)
    
    ## Add color bar
    if (plot.cb) {
      image.plot(x$t, 
                 yvals, 
                 t(zvals), 
                 zlim=zlims,
                 ylim=rev(range(yvals)),
                 xlab=xlab, 
                 ylab=ylab,
                 col=fill.colors,
                 smallplot=legend.loc,
                 horizontal=legend.horiz,
                 legend.only=TRUE, 
                 axis.args=
                   list(at=locs, labels=format(leg.lab, dig=2)), 
                 xpd=NA)
      box()
    }
    ## COI
    if (plot.coi)
      lines(x$t, log2(x$coi), lty=2, lwd=4, col="white")
    ## sig.level contour (default is 95%)
    if (plot.sig & length (x$signif) > 1) {
      if (x$type %in% c("wt", "xwt")) {
        contour(x$t, yvals, t(x$signif), level=sig.level, col="black", lwd=4, 
                add=TRUE, drawlabels=FALSE)
      }
      else {
        contour(x$t, yvals, t(x$signif), nlevel=1, col="black", lwd=4, 
                add=TRUE, drawlabels=FALSE)
      }
    }
    ## Plot phases
    if (plot.phase) {
      a=x$phase
      ## Remove phases where power is weak
      locs=which (zvals < quantile(zvals, arrow.cutoff))
      a[locs]=NA
      x.ind=seq(max(floor(x$dt/2),1), length(x$t), length.out=40)
      y.ind=seq(max(floor(1/2),1), length(x$period), length.out=50)
      phase.plot(x$t[x.ind], log2(x$period[y.ind]), a[y.ind, x.ind], 
                 arrow.size=arrow.size, arrow.lwd=arrow.lwd)    
    }
  }
