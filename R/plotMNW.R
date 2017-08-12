#' Draw a multi-edge network from multiple studies.
#'
#' This function draws a network with multiple lines per edge, representing
#' multiple study-specific correlations for the corresponding gene pairs. Nodes
#' can be colored by DE information, and each line in edges can be colored
#' based on the correlation magnitude and direction.
#'
#' Edge lines are sorted first and put in order separately for every edge, for
#' effective display of study agreement.
#'
#' @param cormatrix.list A list of symmetric correlation matrices from multiple
#' studies.
#' @param genes Node names to be displayed in the network plot. If \code{NULL}
#' (default), \code{rownames(cormatrix.list[[1]])} is used.
#' @param mycolors Edge color palette. The default is \code{bluered(201)}
#' implemented by \code{gplots}, which sets c("blue", "white", "red") for the
#' correlation of c(-1, 0, 1).
#' @param ncolor Node color for selected genes specified by \code{node.de}.
#' @param node.de A vector of 0's and 1's. The "1" nodes are colored in
#' \code{ncolor}.
#' @param r The radius of the big circle where nodes are placed on
#' @param nr The radius of each node
#' @param jt Distance between adjacent lines in an edge
#' @param lwd The edge line thickness
#' @param xlab1 Sub-text to be place below the network
#' @param text.cex Node font size
#' @param mygrey Grey scale
#' @param color logical: Should the edges be colored ?
#' @param corv.decrease logical: Should the edges be ordered by decreasing size ?
#' @param range.cutoff Numeric value of the ranges or "sign"
#' @param range.agree Numeric value that affects the edges
#' @param sign.agree Numeric value that affects edges
#' @param \dots The rest of arguments are passed to \code{plot}.
#' @note Because edges tend to be thicker as the number of studies grow larger,
#' including too many genes (nodes) may result in a noninformative network
#' plot.
#' @author YounJeong Choi
#' @importFrom grDevices grey
#' @importFrom graphics lines plot symbols text
#' @importFrom methods new
#' @importFrom stats as.dist cor
#' @importFrom utils combn
#' @importFrom gplots bluered
#' @references Choi and Kendziorski (2009)
#' @keywords GSCA network multi-edge
#' @export
#' @examples
#' data(LungCancer3)
#' GS <- LungCancer3$info$GSdef
#' GSdesc <- LungCancer3$info$Name
#'
#' setid <- "GO:0008033"
#' gid <- GS[[setid]]
#' ss <- c("KARS", "SARS", "AARS", "SSB", "POP1", "RPP30")
#'
#' data.list <- list(Harvard = LungCancer3$data$Harvard[gid, 140:156],
#' Stanford = LungCancer3$data$Stanford[gid, 42:46],
#' Michigan = LungCancer3$data$Michigan[gid, 87:96])
#'
#' cormatrix.list <- lapply(lapply(data.list, t), cor,
#' use = "pairwise.complete.obs")
#'
#' plotMNW(cormatrix.list, genes = ss, mycolors = bluered(201), ncolor =
#' "yellow", node.de = c(rep(1, 3), rep(0, 3)), lwd = 5, jt = 0.3)
#'
`plotMNW` <-
function(cormatrix.list, genes = NULL,
  mycolors = bluered(201), mygrey = rev(grey(seq(0,1, by = .01))),
  ncolor = "white", node.de = NULL,
  r = 10, nr = 2, jt = .1, color = TRUE, lwd = 1, xlab1 = "",
  corv.decrease = FALSE, text.cex = 1, range.cutoff = 2,
  range.agree = NULL, sign.agree = 1, ...){

  make.edgename <- function(gpair){
    paste(gpair[1], " ~ ", gpair[2], sep = "")
  }

  nS <- length(cormatrix.list)
  hid <- rownames(cormatrix.list[[1]])
  ng <- length(hid)

  if (is.null(range.agree)) range.agree <- nS
  if (range.agree > nS) range.agree <- nS

  if (is.null(genes)){
    genes <- hid
  }

  edges <- combn(genes, 2, make.edgename)
  ne <- length(edges)

  corr.gs <- t(sapply(cormatrix.list, as.dist))
  colnames(corr.gs) = edges
  rownames(corr.gs) = names(cormatrix.list)

  theta = 0:(ng - 1) * 2 * pi / ng
  fs = r + nr

  cx = r * cos(theta)
  cy = r * sin(theta)
  names(cx) = genes
  names(cy) = genes

  gpairs = combn(genes, 2)
  plot(c(-1,1) * fs, c(-1,1) * fs, axes = F, xlab = xlab1, ylab = "",
       type = "n", ...)

  ned = 0 # number of edges drawn
  for ( i in 1:ne ){

    n1 = gpairs[1, i]
    n2 = gpairs[2, i]

    if (cx[n1] != cx[n2]){
      t = (cy[n2] - cy[n1]) / (cx[n2] - cx[n1])
    } else {
      t = 2.78588e+15 # tangent becomes Inf, so put in the biggest value
    }
    x1 = cx[n1] + (1:nS - mean(1:nS)) * jt * t / sqrt(1 + t^2) * (-1)
    x2 = cx[n2] + (1:nS - mean(1:nS)) * jt * t / sqrt(1 + t^2) * (-1)
    y1 = cy[n1] + (1:nS - mean(1:nS)) * jt / sqrt(1 + t^2)
    y2 = cy[n2] + (1:nS - mean(1:nS)) * jt / sqrt(1 + t^2)

    if (color){
      corv = sort(round(corr.gs[,edges[i]] * 100) + 101,
        decreasing = corv.decrease)
    } else {
      corv = sort(round(abs(corr.gs[,edges[i]] * 100)),
        decreasing = corv.decrease)
    }

    d = FALSE
    if (range.cutoff == "sign"){
      if (sum(corr.gs[,edges[i]] >= 0) >= sign.agree |
          sum(corr.gs[,edges[i]] <= 0) >= sign.agree){
        d = TRUE
      }

    } else {

      if (range.agree == nS){
        if (diff(range(corr.gs[,edges[i]])) <= range.cutoff){
          d = TRUE
        }
      } else if (range.agree < nS){
        range2 = apply(combn(corr.gs[,edges[i]], range.agree, range), 2, diff)
        if (min(range2) <= range.cutoff){
          d = TRUE
        }
      }
    }

    if (d){
      ned = ned + 1
      for ( j in 1:nS ){
        if (color){
          lines(c(x1[j], x2[j]), c(y1[j], y2[j]),
                lwd = lwd, col = mycolors[corv[j]])
        } else {
          lines(c(x1[j], x2[j]), c(y1[j], y2[j]),
                lwd = lwd, col = mygrey[corv[j]])
        }
      }
    }

  }

  if (is.null(node.de)){
    nc <- rep(ncolor, ng)
  } else if (is.numeric(node.de)){
    nc <- rep("white", ng)
    numc <- length(unique(node.de))
    if (numc > 1){
      for (j in 1:(numc - 1)) nc[which(node.de == j)] <- ncolor[j]
    }
  } else {
    index <- which(node.de)
    nc <- rep("white", ng)
    nc[index] <- ncolor
  }

  for (i in 1:ng) {
    symbols(cx[i], cy[i], circles = nr, add = T, inches = FALSE, bg = nc[i])
  }
  text(cx, cy, genes, cex = text.cex)
  if (range.cutoff < 2 | sign.agree > 1){
    message(ned, "out of", ne, "edges are drawn.\n")
    return(ned)
  }
}

