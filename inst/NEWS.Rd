\name{NEWS}
\title{News for Package 'biwavelet'}
\encoding{UTF-8}
\section{Changes in biwavelet version 0.15 (2013-04-08)}{
  \subsection{added}{
    \itemize{
      \item Users can now specify the color, line width and line type for
      the COI, significance contours and phase arrows using the \code{plot} function.
    }
  }
    \subsection{fixed}{
\itemize{
      \item Removed misleading examples showing how to compute the 'bias-corrected' wavelet coherence.
      There is no bias for the wavelet coherence function, so using the default \code{type} argument in the 
      \code{plot} function is recommended.
      \item Fixed typos in the documentation of plot.biwavelet and xwt (thanks Lei Cheng).
    }
  }
}
\section{Changes in biwavelet version 0.14 (2013-03-06)}{
  \subsection{added}{
    \itemize{
      \item As of biwavelet version 0.14, the bias-corrected wavelet and
        cross-wavelet spectra are automatically computed and plotted by default
        using the methods described by Liu et al. (2007) and Veleda et al. (2012). 
        This correction is needed because the traditional approach for computing the power
        spectrum (e.g., Torrence and Compo 1998) leads to an artificial and systematic
        reduction in power at lower periods
    }
  }
  \subsection{fixed}{
    \itemize{
      \item Plotting function now accepts traditional plotting flags such as xaxt and yaxt
        to control x and y tickmarks
    }
  }
}
