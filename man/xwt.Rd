\name{xwt}
\alias{xwt}
\title{
Compute cross-wavelet
}
\description{
Compute cross-wavelet
}
\usage{
xwt (d1, d2, pad = TRUE, dj = 1/12, s0 = 2 * dt, J1 = NULL, 
     max.scale = NULL, mother = c("morlet", "paul", "dog"), 
     param = -1, lag1 = NULL, sig.level = 0.95, sig.test = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d1}{
time series 1 in matrix format (\code{n} rows x 2 columns). The first column 
  should contain the time steps and the second column should contain the values.
}
  \item{d2}{
time series 2 in matrix format (\code{n} rows x 2 columns). The first column 
  should contain the time steps and the second column should contain the values.
}
  \item{pad}{
pad the values will with zeros to increase the speed of the transform. Default is TRUE.
}
  \item{dj}{
spacing between successive scales. Default is 1/12.
}
  \item{s0}{
smallest scale of the wavelet. Default is 2*\code{dt}
}
  \item{J1}{
number of scales - 1.
}
  \item{max.scale}{
maximum scale. Computed automatically if left unspecified.
}
  \item{mother}{
type of mother wavelet function to use. Can be set to \code{morlet}, \code{dog}, 
or \code{paul}. Default is \code{mother}. Significance testing is only available
for \code{morlet} wavelet.
}
  \item{param}{
nondimensional parameter specific to the wavelet function.
}
  \item{lag1}{
vector containing the AR(1) coefficient of each time series.
}
  \item{sig.level}{
significance level. Default is 0.95.
}
  \item{sig.test}{
type of significance test. If set to 0, use a regular \eqn{\chi^2} test. 
If set to 1, then perform a time-average test. If set to 2, then do a scale-average
test.
}
}
\value{
Returns a \code{biwavelet} object containing:

\item{coi }{matrix containg cone of influence}
\item{wave }{matrix containing the cross-wavelet transform}
\item{power }{matrix of power}
\item{phase }{matrix of phases}
\item{period }{vector of periods}
\item{scale }{vector of scales}
\item{dt }{length of a time step}
\item{t }{vector of times}
\item{xaxis }{vector of values used to plot xaxis}
\item{s0 }{smallest scale of the wavelet }
\item{dj }{spacing between successive scales}
\item{d1.sigma }{standard deviation of time series 1}
\item{d2.sigma }{standard deviation of time series 2}
\item{mother }{mother wavelet used}
\item{type }{type of \code{biwavelet} object created (\code{xwt})}
\item{signif }{matrix containg significance levels}
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
\author{
Tarik C. Gouhier (tarik.gouhier@gmail.com)

Code based on WTC MATLAB package written by Aslak Grinsted.
}

\examples{
t1=cbind(1:100, rnorm(100))
t2=cbind(1:100, rnorm(100))
## Cross-wavelet
xwt.t1t2=xwt(t1, t2)
## Plot cross-wavelet and phase difference (arrows)
plot(xwt.t1t2, plot.cb=TRUE, plot.phase=TRUE)
## Real data
data(enviro.data)
## Cross-wavelet of MEI and NPGO
xwt.mei.npgo=xwt(subset(enviro.data, select=c("date", "mei")), 
        subset(enviro.data, select=c("date", "npgo")))
## Make room to the right for the color bar
par(oma=c(0, 0, 0, 1), mar=c(5, 4, 4, 5) + 0.1)    
plot(xwt.mei.npgo, plot.cb=TRUE, plot.phase=TRUE)
}