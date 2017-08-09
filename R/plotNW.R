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

