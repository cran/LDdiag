#  File library/LDdiag/R/Copas.R/pregibon.R
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


pregibon <- function(object,digits=4) UseMethod("pregibon")

pregibon.default <- function(object,digits=4) {
  scip <- options()$scipen
  options(scipen=4)
  if(as.character(object$call)[1]=="lm") pregibon.lm(object)
  else  if(as.character(object$call)[1]=="glm") pregibon.glm(object,digits)
  else cat("The function call doesn't support yet.\n")
  options(scipen=scip)
}

pregibon.lm <- function(object, digits=4) {
  dname <- deparse(substitute(object))
  mf <- model.frame(object)  #model frame
  Y <- model.response(mf) #from model frame
  hat <- fitted(object)
  hatsq <- hat^2
  coef_ <- cbind(hat,hatsq)
  k <- ncol(coef_)+1
  newfit <- lm(Y~coef_)
  ses <- sqrt(diag(vcov(newfit)))
  tvals <- coef(newfit)/ses
  tvals <- tvals[k]
  names(tvals) <- "T.value"
  defs <-newfit$df
  names(defs) <- "df"
  pvalue <- 2*pt(abs(tvals),defs,lower.tail=FALSE)
  rval <- list(statistic=tvals, parameter=defs,p.value=pvalue,data.name=dname, method="Pregibon's Link Test")
  class(rval) <- "htest"
  cat("\nThe Pregibon Test p-value:",round(pvalue,digits))
  if (pvalue>0.05) { cat(". This model DOES PASS this test.\n")
  }else{ cat(". This model DOES NOT PASS this test.\n") }
  return(rval)
}

pregibon.glm <- function(object, digits=4) {
  testf <- function(object) { #this function modified from summary.glm()
    df.r <- object$df.residual
   	dispersion <-
  	if(object$family$family %in% c("poisson", "binomial"))  1
  	else if(df.r > 0) {
    sum(object$residuals^2)/ df.r
  	}
    p <- object$rank
    p1 <- 1:p
    Qr <- object$qr
    coef.p <- object$coefficients[Qr$pivot[p1]]
    covmat.unscaled <- chol2inv(Qr$qr[p1,p1,drop=FALSE])
    covmat <- dispersion*covmat.unscaled
    var.cf <- diag(covmat)
      s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    pvalue <- 2*pt(-abs(tvalue), df.r)
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    return(coef.table)
    }
  dname <- deparse(substitute(object))
  fm <- object$'family'
  mf <- model.frame(object)
  X <- model.matrix(terms(object),mf)
  Y <- model.response(mf)
  xb <- fitted(object)
  xb2 <- xb^2
  X <- cbind(1,xb,xb2)
  fits <- glm.fit(X,Y, family=fm)
  ans <- testf(fits)
  defs <- fits$df.residual
  names(defs) <-"df"
  tvals <- ans[3,3]
  names(tvals) <- "T.value"
  pvalue <- ans[3,4]
  rval <- list(statistic=tvals, parameter=defs,p.value=pvalue,data.name=dname, method="Pregibon's Link Test")
  class(rval) <- "htest"

  cat("\nThe Pregibon Test p-value:",round(pvalue,digits))
  if (pvalue>0.05) { cat(". This model DOES PASS this test.\n")
  }else{ cat(". This model DOES NOT PASS this test.\n") }
  return(rval)
}



