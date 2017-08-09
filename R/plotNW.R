#' The function to draw a network from a correlation matrix
#' 
#' This function is a custom wrapper of \code{plot.graph} implemented in
#' \code{Rgraphviz}. It draws a network from a given correlation matrix. Nodes
#' can be colored by DE information, and edges can be colored based on the
#' correlation magnitude and direction.
#' 
#' 
#' @param genes Node names to be displayed in the network plot. If \code{NULL}
#' (default), \code{rownames(cormatrix)} is used.
#' @param cormatrix A symmetric correlation matrix to draw a network from.
#' @param node.de A vector of 0's and 1's. The "1" nodes are colored in
#' \code{ncolor}.
#' @param ncolor Node color for selected genes specified by \code{node.de}.
#' @param ecolor Edge color palette. The default is \code{bluered(201)}
#' implemented by \code{gplots}, which sets c("blue", "white", "red") for the
#' correlation of c(-1, 0, 1).
#' @param ewidth Edge width. Default is set to 1.
#' @param gtype Graph type (\code{layout}) to be passed to \code{plot.graph}.
#' One of \code{dot}, \code{neato}, \code{twopi}, \code{circo}, and \code{fdp}.
#' The default is \code{neato}.
#' @param \dots The rest of arguments are passed to \code{plot.graph}.
#' @author YounJeong Choi
#' @seealso \code{\link{plot.graph}}, \code{\link{Rgraphviz}}
#' @references Choi and Kendziorski (2009)
#' @keywords GSCA network
#' @examples
#' 
#' data(LungCancer3)
#' GS <- LungCancer3$info$GSdef
#' GSdesc <- LungCancer3$info$Name
#' 
#' setid <- "GO:0019216"
#' gid <- GS[[setid]]
#' ss <- c("SERPINA3", "SOD1", "SCAP", "NPC2", "ADIPOQ", "PRKAA1", "AGT",
#' "PPARA", "BMP6", "BRCA1")
#' 
#' plotNW(genes = ss, cormatrix = cor(t(LungCancer3$data$Michigan[gid, 87:96]), use = "pairwise.complete.obs"), node.de = c(rep(1, 5), rep(0, 5)), ncolor = "yellow", ecolor = bluered(201), ewidth = 5, gtype = "circo")
#' 
`plotNW` <-
function(genes = NULL, cormatrix,
                    node.de = NULL, ncolor = "red", 
                    ecolor = NULL, ewidth = 1, 
                    gtype = "neato", ...){

  # build a graph
  ng <- nrow(cormatrix)
  if (is.null(genes)) genes <- rownames(cormatrix)
  g1 = new("graphNEL", nodes = genes)
  for ( k in 2:ng ){
    for ( m in 1:(k-1) ){
      g1 = addEdge(from = genes[k], to = genes[m], graph = g1)
    }
  }

  # DE info (color the node if DE)
  if (is.null(node.de)) node.de <- rep(0, ng)
  names(node.de) <- genes
  node.de[which(node.de == 0)] <- "white"
  node.de[which(node.de == 1)] <- ncolor

  # edges
  if (is.null(ecolor)) ecolor <- rep("black", 201)
  n <- floor(length(ecolor) / 2)

  # draw the graph
  attrs <- list(node = list(shape="ellipse", fixedsize = FALSE))
  nAttrs <- list(fillcolor = node.de)
  eAttrs <- list()
  lc <- NULL; lw <- NULL
  for ( k in 1:(ng-1) ){
    for ( m in (k+1):ng ){
      lw <- c(lw, ewidth)
      lc <- c(lc, ecolor[round(cormatrix[k,m] * n) + n + 1])
      names(lc)[length(lc)] <- paste(genes[k], "~", genes[m], sep = "")
      names(lw)[length(lw)] <- paste(genes[k], "~", genes[m], sep = "")
    }
  }

  eAttrs$color <- lc
  eAttrs$lwd <- lw
  plot(g1, gtype, attrs = attrs, nodeAttrs = nAttrs, edgeAttrs = eAttrs, ...)
}

