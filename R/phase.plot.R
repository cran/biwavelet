phase.plot <-
function (x, y, phases, arrow.size=0.05, arrow.lwd=2) {
  g=meshgrid(x, y)
  x=g$x
  y=g$y
  dx=cos(phases)*arrow.size
  dy=sin(phases)*arrow.size
  suppressWarnings(
    arrows(x - dx, y - dy, x + dx, y + dy, length=arrow.size, lwd=arrow.lwd, 
           code=2))
}

