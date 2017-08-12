#' The function to run a meta-GSCA.
#'
#' This function can be used to run a meta-GSCA, described in Choi and
#' Kendziorski (2009). Unlike \code{singleDC}, the study-specific values (e.g.,
#' gene-gene pairwise correlations, difference in condition-specific
#' correlation signs) need to be pre-calculated and provided as input. For each
#' gene set defined in \code{GSdefList}, the dispersion index is calculated
#' between two given studies and returned.
#'
#' Gene pairs are randomly permuted across gene sets for \code{nperm} times.
#' Permutation-based p-values are calculated, based on the rank of observed
#' \code{DI} among permuted index values.
#'
#' Gene pairs are permuted across gene sets, as described in Choi and
#' Kendziorski (2009). For each permutation, annotated gene pairs (gene pairs
#' which belong to at least one gene set) are randomly re-assigned to gene
#' sets, and dispersion indices (DIs) are calculated based on those random gene
#' sets. As focus is on preservation, the p-value for each gene set is
#' calculated as:
#'
#' p = sum(permutation DIs <= observed DI) / nperm .
#'
#' @param corr1,corr2 Two symmetric matrices from two studies of interest. For
#' meta-GSCA, a difference matrix is used, of two condition-specific
#' correlation sign matrices within each study. The lower triangle is used.
#' \code{rownames(data)} is used to subset a sub-matrix from \code{data} for
#' each gene set. (Rows must be named by gene IDs used in \code{GSdefList}. For
#' example, if \code{GSdefList} defines gene sets in Entrez Gene IDs,
#' \code{rownames(data)} should be Entrez Gene IDs.
#' @param GSdefList A list of character vectors that define gene sets. Each
#' entry of this list is a gene set.
#' @param nperm The desired number of permutations.
#' @param permDI TRUE/FALSE. If set TRUE, dispersion index values from
#' permutation are saved and returned; if FALSE, permutation-based dispersion
#' index values are not returned. Default is FALSE.
#' @return \item{DI}{The dispersion index vector for each gene set.}
#' \item{pvalue}{The permutation-based p-value for each gene set.}
#' \item{permv}{The permutation-based DI matrix, of \code{nperm} columns. The
#' first column is identical to what is returned by \code{DI}.}
#' @note Currently, \code{metaDI} implements meta-analysis of gene-gene
#' pairwise correlations from two studies. In addition to meta-GSCA described
#' in Choi and Kendziorski (2009), which uses the sign difference, raw
#' correlations can be input to investigate preservation of them.
#' @author YounJeong Choi
#' @references Choi and Kendziorski, submitted
#' @keywords GSCA meta-GSCA
#' @examples
#'
#' data(LungCancer3)
#' GS <- LungCancer3$info$GSdef
#' GSdesc <- LungCancer3$info$Name
#'
#' data.grouped <- list(
#' Tumor = list(Harvard = LungCancer3$data$Harvard[, 1:139],
#' Michigan = LungCancer3$data$Michigan[, 1:86]),
#' Normal = list(Harvard = LungCancer3$data$Harvard[, 140:156],
#' Michigan = LungCancer3$data$Michigan[, 87:96]))
#'
#' corr.t <- lapply(lapply(data.grouped$Tumor, t), cor, use = "pairwise.complete.obs")
#' corr.n <- lapply(lapply(data.grouped$Normal, t), cor, use = "pairwise.complete.obs")
#' cor.diff <- list(Harvard = corr.t$Harvard - corr.n$Harvard,
#' Michigan = corr.t$Michigan - corr.n$Michigan)
#'
#' cor.diff.sign <- list(
#' Harvard = apply((cor.diff$Harvard > 0), 2, as.numeric) -
#' apply((cor.diff$Harvard < 0), 2, as.numeric),
#' Michigan = apply((cor.diff$Michigan > 0), 2, as.numeric) -
#' apply((cor.diff$Michigan < 0), 2, as.numeric))
#'
#' for (i in 1:length(cor.diff.sign)) {
#' rownames(cor.diff.sign[[i]]) <- colnames(cor.diff.sign[[i]])
#' }
#' dist.HM <- metaDI(cor.diff.sign$Harvard, cor.diff.sign$Michigan, GS, 3, permDI = TRUE)
#' @export
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

