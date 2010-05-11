#  File library/LDdiag/R/Copas.R/linkpick.R
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


linkpick <- function(formula, data){
  fm <- as.name(deparse(substitute(formula)))
  dt2 <- as.name(deparse(substitute(data)))
  txt=paste("fit<-boxcox(lm(",fm,",",dt2,"),plotit=FALSE)", sep='')
  fit <- eval(parse(text=txt))
  ans <- fit$x[fit$y==max(fit$y)]
  cat("\n  Lambda =",ans)
  if (-1.5 < ans      & ans< -.5) {
    cat (". 'Inverse' link function is a reasonable option.")
    } else if (-2.5 <= ans & ans < -1.5) {
    cat(". 'inverse quadratic' link is a reasonable option.\n\n")
    } else if (-.5 <= ans & ans <.25) {
    cat(". 'log' link is a reasonable option.\n\n")
    } else if (.25 <= ans & ans <.75) {
    cat(". 'Square root' link is a reasonable option.\n\n")
    } else if (.75 <= ans & ans <1.5){
    cat(". 'Identity' link is a reasonable option.\n\n")
    } else if (1.5 <= ans & ans <2.5){
    cat(". 'Square' link is a reasonable option.\n\n")
    }
}


