#  File library/LDdiag/R/Copas.R/pcorr.R
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


pcorr <- function(object,digits=4) {
  scip <- options()$scipen
  options(scipen=6)
  mu <- fitted(object)
  res <- resid(object)
  crc <- cor(mu,res,method="pearson")
  n <- length(mu)

  trpm <- .5*log((1+crc)/(1-crc))
  tl <- trpm - 1.96*sqrt(1/(n-3))
  th <- trpm + 1.96*sqrt(1/(n-3))
  rl <- round((exp(2*tl)-1)/(exp(2*tl)+1),digits)
  rh <- round((exp(2*th)-1)/(exp(2*th)+1),digits)
  
  cat("\n  Pearson correlation coefficient:", paste(round(crc,digits),sep=""), "\n")
  cat("  The 95% confidence interval is", paste("(",rl,",",rh,")", sep=""),"\n")
  options(scipen=scip)
}



