\name{biwavelet-package}
\alias{biwavelet-package}
\alias{biwavelet}
\docType{package}
\title{
Conduct univariate and bivariate wavelet analyses
}
\description{
This is a port of the WTC MATLAB package written by Aslak Grinsted and the wavelet 
program written by Christopher Torrence and Gibert P. Compo. This package can be 
used to perform univariate and bivariate (cross-wavelet, wavelet coherence) wavelet 
analyses.
}
\details{
\tabular{ll}{
Package: \tab biwavelet\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2012-04-19\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
Univariate and bivariate wavelet analyses
}
\author{
Tarik C. Gouhier

Maintainer: Tarik C. Gouhier <tarik.gouhier@gmail.com>

Code based on WTC MATLAB package written by Aslak Grinsted and the wavelet MATLAB 
program written by Christopher Torrence and Gibert P. Compo.
}
\references{
Cazelles, B., M. Chavez, D. Berteaux, F. Menard, J. O. Vik, S. Jenouvrier, and 
N. C. Stenseth. 2008. Wavelet analysis of ecological time series. 
\emph{Oecologia} 156:287-304.

Grinsted, A., J. C. Moore, and S. Jevrejeva. 2004. Application of the cross 
wavelet transform and wavelet coherence to geophysical time series. 
\emph{Nonlinear Processes in Geophysics} 11:561-566.

Torrence, C., and G. P. Compo. 1998. A Practical Guide to Wavelet Analysis. 
\emph{Bulletin of the American Meteorological Society} 79:61-78.

Torrence, C., and P. J. Webster. 1998. The annual cycle of persistence in the 
El Nino/Southern Oscillation. 
\emph{Quarterly Journal of the Royal Meteorological Society} 124:1985-2004.
}
\keyword{ wavelet }
\keyword{ wavelet coherence }
\keyword{ cross-wavelet }

\examples{
t1=cbind(1:100, rnorm(100))
t2=cbind(1:100, rnorm(100))
## Continuous wavelet transform
wt.t1=wt(t1)
## Plot power
## Make room to the right for the color bar
par(oma=c(0, 0, 0, 1), mar=c(5, 4, 4, 5) + 0.1)
plot(wt.t1, plot.cb=TRUE, plot.phase=FALSE)
## Cross-wavelet
xwt.t1t2=xwt(t1, t2)
## Plot cross wavelet power and phase difference (arrows)
plot(xwt.t1t2, plot.cb=TRUE)
## Wavelet coherence; nrands should be large (>= 1000)
wtc.t1t2=wtc(t1, t2, nrands=10)
## Plot wavelet coherence and phase difference (arrows)
## Make room to the right for the color bar
par(oma=c(0, 0, 0, 1), mar=c(5, 4, 4, 5) + 0.1)
plot(wtc.t1t2, plot.cb=TRUE)
}
