#  File library/LDdiag/R/Copas.R/modifiedHL.R
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


modifiedHL <- function(object,digits=4){
  mf <- model.frame(object)
  mf$pred <- fitted(object)
  mf$Residuals <- resid(object)
  mf <- mf[order(mf$pred),]
  mf$seq <- 1:nrow(mf)
  xcuts <-  quantile(1:nrow(mf),prob=seq(0.1:1, by=0.1))
  mf$group <- cut(mf$seq,xcuts,labels=FALSE)
  mf$group <- ifelse(is.na(mf$group),0,mf$group)
  mf$group <- as.factor(mf$group)
  ans <- anova(lm(Residuals~group, mf))
  pvalue <-lapply(ans, "[[", 1)$"Pr(>F)"
  cat("\nThe modified Hosmer-Lemeshow Test p-value:",round(pvalue,digits))
  if (pvalue>0.05) {
  cat(". This model does PASS this test.\n")
  }else{
  cat(". This model doesn't PASS this test.\n")
  }

  cat("\n")
  print(ans)
  cat("\n")
}


