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
    cat(ned, "out of", ne, "edges are drawn.\n")
    return(ned)
  }
}

