\name{plotMNW}
\alias{plotMNW}
\title{The function to draw a multi-edge network from a list of
  correlation matrices from multiple studies.}
\description{
  This function draws a network with multiple lines per edge,
  representing multiple study-specific correlations for the
  corresponding gene pairs. Nodes can be colored by DE information, and
  each line in edges can be colored based on the correlation magnitude
  and direction.  
}
\usage{
plotMNW(cormatrix.list, genes = NULL, mycolors = bluered(201), ncolor = "white", node.de = NULL, r = 10, nr = 2, jt = 0.1, lwd = 1, xlab1 = "", text.cex = 1, ...)
}
\arguments{
  \item{cormatrix.list}{A list of symmetric correlation matrices from
    multiple studies.}
  \item{genes}{Node names to be displayed in the network plot. If
    \code{NULL} (default), \code{rownames(cormatrix.list[[1]])} is used.}
  \item{mycolors}{Edge color palette. The default is \code{bluered(201)}
    implemented by \code{gplots}, which sets c("blue", "white", "red")
    for the correlation of c(-1, 0, 1).}
  \item{ncolor}{Node color for selected genes specified by \code{node.de}.}
  \item{node.de}{A vector of 0's and 1's. The "1" nodes are colored in
    \code{ncolor}.}
  \item{r}{The radius of the big circle where nodes are placed on}
  \item{nr}{The radius of each node}
  \item{jt}{Distance between adjacent lines in an edge}
  \item{lwd}{The edge line thickness}
  \item{xlab1}{Sub-text to be place below the network}
  \item{text.cex}{Node font size}
  \item{\dots}{The rest of arguments are passed to \code{plot}.}
}
\details{
  Edge lines are sorted first and put in order separately for every
  edge, for effective display of study agreement.
}
\references{Choi and Kendziorski (2009)}
\author{YounJeong Choi}
\note{
  Because edges tend to be thicker as the number of studies grow
  larger, including too many genes (nodes) may result in a
  noninformative network plot.
}
\examples{
data(LungCancer3)
GS <- LungCancer3$info$GSdef
GSdesc <- LungCancer3$info$Name

setid <- "GO:0008033"
gid <- GS[[setid]]
ss <- c("KARS", "SARS", "AARS", "SSB", "POP1", "RPP30")

data.list <- list(Harvard = LungCancer3$data$Harvard[gid, 140:156],
Stanford = LungCancer3$data$Stanford[gid, 42:46],
Michigan = LungCancer3$data$Michigan[gid, 87:96])

cormatrix.list <- lapply(lapply(data.list, t), cor,
use = "pairwise.complete.obs")

plotMNW(cormatrix.list, genes = ss, mycolors = bluered(201), ncolor =
"yellow", node.de = c(rep(1, 3), rep(0, 3)), lwd = 5, jt = 0.3)
}
\keyword{GSCA}
\keyword{network}
\keyword{multi-edge}

