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
used to perform univariate and bivariate (cross-wavelet, wavelet coherence, wavelet clustering) 
wavelet analyses.
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

Rouyer, T., J. M. Fromentin, F. Menard, B. Cazelles, K. Briand, R. Pianet, 
B. Planque, and N. C. Stenseth. 2008. Complex interplays among population 
dynamics, environmental forcing, and exploitation in fisheries. 
\emph{Proceedings of the National Academy of Sciences} 105:5420-5425.

Rouyer, T., J. M. Fromentin, N. C. Stenseth, and B. Cazelles. 2008. 
Analysing multiple time series and extending significance testing in 
wavelet analysis. \emph{Marine Ecology Progress Series} 359:11-23. 

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

## Perform wavelet clustering of three time series
t1=cbind(1:100, sin(seq(from=0, to=10*2*pi, length.out=100)))
t2=cbind(1:100, sin(seq(from=0, to=10*2*pi, length.out=100)+0.1*pi))
t3=cbind(1:100, rnorm(100))
## Compute wavelet spectra
wt.t1=wt(t1)
wt.t2=wt(t2)
wt.t3=wt(t3)
## Store all wavelet spectra into array
w.arr=array(NA, dim=c(3, NROW(wt.t1$wave), NCOL(wt.t1$wave)))
w.arr[1, , ]=wt.t1$wave
w.arr[2, , ]=wt.t2$wave
w.arr[3, , ]=wt.t3$wave

## Compute dissimilarity and distance matrices
w.arr.dis=wclust(w.arr)
plot(hclust(w.arr.dis$dist.mat, method="ward"), sub="", main="", 
        ylab="Dissimilarity", hang=-1)
}
