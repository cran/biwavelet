plot.biwavelet <-
  function (x, ncol=64, xlab="Time", ylab="Period", 
            tol=0.95, plot.cb=FALSE, plot.phase=FALSE,
            type=c("power.corr.norm", "power.corr", "power.norm", "power", "wavelet", "phase"), 
            plot.coi=TRUE, plot.sig=TRUE, bw=FALSE,
            legend.loc=NULL, 
            legend.horiz=FALSE,
            arrow.size=0.08,
            arrow.lwd=2, arrow.cutoff=0.9, xlim = NULL, ylim = NULL,
            xaxt = "s", yaxt = "s", form='%Y', ...) {
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
    yrange=ylim
    y.ticks = 2^(floor(log2(min(x$period, yrange))):floor(log2(max(x$period, yrange))))
    types=c("power.corr.norm", "power.corr", "power.norm", "power", "wavelet", "phase")
    type=match.arg(tolower(type), types)
    
    if (type=="power.corr" | type=="power.corr.norm") {
      if (x$type=="wtc") {
        x$wave=x$wave.corr
        x$power=x$power.corr
      }
      else if (x$type=="xwt") {
        x$power=x$power.corr
        x$wave=x$wave.corr
      }
      else
        x$power=x$power.corr
    }
    
    if (type=="power.norm" | type=="power.corr.norm") {
      if (x$type == "xwt") {
        zvals=log2(x$power)/(x$d1.sigma*x$d2.sigma)
        zlims=range(c(-1, 1) * max(zvals))
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
        zlims=range(c(-1, 1) * max(zvals))
        zvals [zvals < zlims[1]]=zlims[1]  
        locs=pretty(range(zvals), n=5)
        leg.lab=2^locs
      }
    }
    else if (type=="power"| type=="power.corr") {
      zvals=log2(x$power)
      zlims=range(c(-1, 1)*max(zvals))
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
      stop("type must be power, power.norm, power.corr, power.corr.norm, wavelet or phase")
    }
    if (is.null(xlim))
      xlim=range(x$t)
    yvals=log2(x$period)
    if (is.null(ylim))
      ylim=range(yvals)
    else
      ylim=log2(ylim)
    image(x$t,
          yvals, 
          t(zvals), 
          zlim=zlims,
          xlim=xlim,
          ylim=rev(ylim),
          xlab=xlab, 
          ylab=ylab,
          yaxt="n",
          xaxt="n",
          col=fill.colors, ...)
    box()
    if (class(x$xaxis) == "Date") {
      if (xaxt!="n") {
        xlocs=pretty(x$t)+1
        axis(side=1, at=xlocs, labels=format(x$xaxis[xlocs], form))
      }
    }
    else {
      if (xaxt!="n") {
        xlocs=axTicks(1)
        axis(side=1, at=xlocs)        
      }
    }
    if (yaxt!="n") {
      axis.locs=axTicks(2)
      yticklab=format(2^axis.locs, dig=1)
      axis(2, at=axis.locs, labels=yticklab)
    }
    
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
        contour(x$t, yvals, t(x$signif), level=tol, col="black", lwd=4, 
                add=TRUE, drawlabels=FALSE)
      }
      else {
        tmp=x$rsq/x$signif
         contour(x$t, yvals, t(tmp), level=tol, col="black", lwd=4, 
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
