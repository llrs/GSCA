`metaDI` <-
function(corr1, corr2, GSdefList, nperm, permDI = FALSE){
  
  nGS <- length(GSdefList)
  GSdefList2 <- GSdefList
  for (i in 1:nGS){
    GSdefList2[[i]] <- as.numeric(match(GSdefList[[i]], rownames(corr1)))
  }
  
  dmatrix <- .Call("dEuc2perm", corr1, corr2, GSdefList2, nperm)
  rownames(dmatrix) <- names(GSdefList)
  colnames(dmatrix) <- paste("P", 1:nperm, sep = "")

  pvalue <- rep(NA, nGS)
  names(pvalue) <- names(GSdefList)
  for (i in 1:nGS){
    pvalue[i] <- sum(dmatrix[i, ] <= dmatrix[i, 1], na.rm = TRUE)/
      (nperm - sum(is.na(dmatrix[i, ])))
  }

  if (permDI){
    return(list(DI = dmatrix[, 1], pvalue = pvalue, permv = dmatrix))
  } else {
    return(list(DI = dmatrix[, 1], pvalue = pvalue, permv = NULL))
  }
}

