#' Compute dissimilarity between multiple wavelet spectra
#'
#' @param w.arr \code{N x p x t} array of wavelet spectra where \code{N} is the
#'   number of wavelet spectra to be compared, \code{p} is the number of periods
#'   in each wavelet spectrum and \code{t} is the number of time steps in each
#'   wavelet spectrum.
#'
#' @param quiet Do not display progress bar.
#'
#' @return Returns a list containing:
#'   \item{diss.mat}{square dissimilarity matrix}
#'   \item{dist.mat}{(lower triangular) distance matrix}
#'
#' @references
#' Rouyer, T., J. M. Fromentin, F. Menard, B. Cazelles, K. Briand, R. Pianet,
#' B. Planque, and N. C. Stenseth. 2008. Complex interplays among population
#' dynamics, environmental forcing, and exploitation in fisheries.
#' \emph{Proceedings of the National Academy of Sciences} 105:5420-5425.
#'
#' Rouyer, T., J. M. Fromentin, N. C. Stenseth, and B. Cazelles. 2008.
#' Analysing multiple time series and extending significance testing in
#' wavelet analysis. \emph{Marine Ecology Progress Series} 359:11-23.
#'
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#'
#' @examples
#' t1 <- cbind(1:100, sin(seq(0, 10 * 2 * pi, length.out = 100)))
#' t2 <- cbind(1:100, sin(seq(0, 10 * 2 * pi, length.out = 100) + 0.1 * pi))
#' t3 <- cbind(1:100, rnorm(100)) # white noise
#'
#' ## Compute wavelet spectra
#' wt.t1 <- wt(t1)
#' wt.t2 <- wt(t2)
#' wt.t3 <- wt(t3)
#'
#' ## Store all wavelet spectra into array
#' w.arr <- array(dim = c(3, NROW(wt.t1$wave), NCOL(wt.t1$wave)))
#' w.arr[1, , ] <- wt.t1$wave
#' w.arr[2, , ] <- wt.t2$wave
#' w.arr[3, , ] <- wt.t3$wave
#'
#' ## Compute dissimilarity and distance matrices
#' w.arr.dis <- wclust(w.arr)
#' plot(hclust(w.arr.dis$dist.mat, method = "ward.D"),
#'      sub = "", main = "", ylab = "Dissimilarity", hang = -1)
#'
#' @export
wclust <- function(w.arr, quiet = FALSE) {

  num_waves <- nrow(w.arr)
  dist.matrix <- matrix(NA, nrow = num_waves, ncol = num_waves)
  k <- 1

  if (!quiet) {
    prog.bar <- txtProgressBar(min = 0, num_waves ^ 2, style = 3)
  }

  for (n in seq_len(num_waves)) {
    for (j in seq_len(num_waves)) {
      dist.matrix[n,j] <- wdist(w.arr[n, ,], w.arr[j, ,])
      k <- k + 1

      if (!quiet) {
        setTxtProgressBar(prog.bar, k)
      }
    }
  }

  if (!quiet) {
    close(prog.bar)
  }

  list(diss.mat = dist.matrix, # square dissimilarity matrix
       dist.mat = as.dist(dist.matrix)) # lower triangular distance matrix
}
