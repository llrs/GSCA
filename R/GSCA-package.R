#' Gene Set Co-expression Analysis
#'
#'
#' \describe{
#' \item{plotMNW}{Draws a multi-edge network from correlation matrices from multiple studies.}
#' \item{plotNW}{Draws a network from a correlation matrix.}
#' \item{singleDC}{gene set differential co-expression analysis}
#' \item{metaDI}{Calculates the dispersion index of the gene list}}
#'
#' @name GSCA-package
#' @aliases GSCA
#' @docType package
#' @author YounJeong Choi
#'
#' Maintainer: YounJeong Choi <ychoi@@biostat.wisc.edu>
#' @references Choi and Kendziorski (2009)
#' @keywords GSCA
NULL

#' Three lung cancer microarray data sets and matched annotation
#'
#' The three human lung cancer microarray data sets, from Stanford (Garber et
#' al., 2001), Harvard (Bhattacharjee et al., 2001), and Michigan (Beer et al.,
#' 2002). The 3,924 common Entrez Gene IDs are represented, matched across the
#' three studies. The 3,649 common gene sets (GO categories and KEGG pathways)
#' are represented, defined in Entrez Gene ID.
#'
#' \code{LungCancer3$data} is a list of 3 R matrices named 'Harvard',
#' 'Stanford', and 'Michigan', respectively. Each matrix contains 3,924 rows as
#' the genes have been matched across three studies. Rows are named by Entrez
#' Gene IDs.
#'
#' The number of columns corresponds to the number of arrays used in each
#' study. Columns are named \code{Tumor1}, \code{Tumor2}, ..., and
#' \code{Normal1}, \code{Normal2}, ..., similarly for three matrices.
#' \code{LungCancer3$data$Harvard} (Harvard study) has 156 columns of 139 tumor
#' vs. 17 normal samples; \code{LungCancer3$data$Stanford} (Stanford study) has
#' 46 columns of 41 tumor vs. 5 normal samples; and
#' \code{LungCancer3$data$Michigan} (Michigan study) has 96 columns of 86 tumor
#' vs. 10 normal samples. For tumor samples, only adenocarcinoma samples have
#' been included for consistency.
#'
#' \code{LungCancer3$info} is again a list of two objects, named \code{GSdef}
#' and \code{Name}, respectively.
#'
#' \code{LungCancer3$info$GSdef} is a list of 3,649 gene set definitions for
#' the 3,924 genes by GO categories and KEGG pathways. This list is named by
#' the gene set IDs such as \code{"GO:0007169"} for GO categories and
#' \code{"00920"} for KEGG pathways. Each entry of this list is a character
#' vector of Entrez Gene IDs, which the gene set consists of. For example,
#' \code{LungCancer3$info$GSdef[[147]]} returns \code{c("348", "25", "27")}.
#'
#' \code{LungCancer3$info$Name} is a character vector of length 3,649
#' corresponding to the gene sets defined in \code{LungCancer3$info$GSdef},
#' named by the gene set IDs. \code{GSdef} and \code{Name} have the gene sets
#' in the same order. For example, \code{LungCancer3$info$GSdef[[123]]} and
#' \code{LungCancer3$info$Name[123]} are both \code{"GO:0032774"}.
#'
#' @name LungCancer3
#' @docType data
#' @format A list of two sub-lists named 'data' and 'info', respectively.
#' @references Beer, D. G. et al. (2002) Gene-expression profiles predict
#' survival of patients with lung adenocarcinoma. Nature Medicine, 8, 816-824.
#'
#' Bhattacharjee, A. et al. (2001) Classification of human lung carcinomas by
#' mRNA expression profiling reveals distinct adenocarcinoma subclasses. Proc.
#' Natl Acad. Sci., 98, 13790-13795.
#'
#' Garber, M. E. et al. (2001) Diversity of gene expression in adenocarcinoma
#' of the lung. Proc. Natl Acad. Sci., 98, 13784-13789.
#' @source Harvard data (Bhattacharjee et al.)
#' http://www.broad.mit.edu/mpr/lung/
#'
#' Stanford data (Garber et al.)
#' http://smd.stanford.edu/cgi-bin/data/viewDetails.pl?exptid=12827&viewSet=1
#'
#' Michigan data (Beer et al.)
#' http://dot.ped.med.umich.edu:2000/ourimage/pub/Lung/index.html
#' @keywords lung cancer, microarray meta-analysis
#' @examples
#' data(LungCancer3)
#' str(LungCancer3$data)
#' str(LungCancer3$info$GSdef[1:3])
#' LungCancer3$info$Name[1:3]
NULL



