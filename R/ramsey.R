#  File library/LDdiag/R/Copas.R/ramsey.R
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


ramsey <- function(object,power=4){
  if(power<3) stop("Minimum power is 3.")
  scip <- options()$scipen
  options(scipen=4)
  dname <- deparse(substitute(object))
  mf <- model.frame(object)
  X <- model.matrix(terms(object),mf)
  Y <- model.response(mf)
  Z <- t(sapply(fitted(object),"^",2:power))
  XZ <- cbind(X,Z)
  r1 <- lm.fit(X,Y)$residuals
  r2 <- lm.fit(XZ,Y)$residuals
  r1sq <- sum(r1^2) 
  r2sq <- sum(r2^2)
  k <- ncol(X)
  n <- nrow(X)
  q <- ncol(Z)
  df1 <- q; df2 <- n-(k+q)
  reset <- df2/df1 * (r1sq-r2sq)/r2sq
  names(reset) <- "RESET"
  dfs <- c(df1, df2)
  names(dfs) <- c("df1","df2")
  pvalue <- pf(reset,df1,df2,lower.tail=FALSE)
  rval <- list(statistic=reset, parameter=dfs,p.value=pvalue,data.name=dname, method="Ramsey RESET Test(Fiited Values)")
  class(rval) <- "htest"
  cat("\nThe Ramsey RESET Test p-value:",round(pvalue,4))
  if (pvalue>0.05) {
  cat(". This model does PASS this test.\n")
  }else{
  cat(". This model doesn't PASS this test.\n")
  }
  options(scipen=scip)
  return(rval)
}


