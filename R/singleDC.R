#' The function to run a single-study GSCA, differential co-expression (DC)
#' analysis
#' 
#' This function runs a single-study GSCA, differential co-expression (DC)
#' analysis, described in Choi and Kendziorski (2009).  The condition-specific
#' gene-gene pairwise correlations are first calculated; then for each gene set
#' defined in \code{GSdefList}, the dispersion index is calculated across
#' condition-specific correlations.
#' 
#' Samples are randomly permuted across conditions for \code{nperm} times.
#' Permutation-based p-values are calculated, based on the rank of observed
#' \code{DI} among permuted index values.
#' 
#' Samples (columns) are permuted across conditions. For each permutation,
#' condition-specific correlations are re-calculated based on permuted samples,
#' and dispersion indices (DIs) are calculated based on those permutation-based
#' correlations. As focus is on difference, the p-value for each gene set is
#' calculated as:
#' 
#' p = sum(permutation DIs >= observed DI) / nperm .
#' 
#' @param data A data matrix of rows representing genes and columns
#' representing arrays. \code{rownames(data)} is used to subset a sub-matrix
#' from \code{data} for each gene set. (Rows must be named by gene IDs used in
#' \code{GSdefList}. For example, if \code{GSdefList} defines gene sets in
#' Entrez Gene IDs, \code{rownames(data)} should be Entrez Gene IDs.
#' @param group A numeric vector that specifies the number of arrays (columns)
#' in each condition. For example, if \code{c(10, 5)} is provided, first 10
#' columns of the \code{data} matrix are used for one condition and the next 5
#' are used for the other condition.
#' @param GSdefList A list of character vectors that define gene sets. Each
#' entry of this list is a gene set.
#' @param nperm The desired number of permutations.
#' @return \item{DI}{The dispersion index vector for each gene set.}
#' \item{pvalue}{The permutation-based p-value for each gene set.}
#' \item{permv}{The permutation-based DI matrix, of \code{nperm} columns. The
#' first column is identical to what is returned by \code{DI}.}
#' @note Currently, \code{singleDC} implements DC analysis for two conditions
#' (e.g., tumor vs. normal) and three conditions (e.g., AA, AB, and BB
#' genotypes). For three conditions, pairwise DIs are first calculated and
#' averaged (internally).
#' @author YounJeong Choi
#' @references Choi and Kendziorski, submitted.
#' @keywords GSCA DC
#' @examples
#' 
#' data(LungCancer3)
#' GS <- LungCancer3$info$GSdef
#' GSdesc <- LungCancer3$info$Name
#' dc.M <- singleDC(data = LungCancer3$data$Michigan, group = c(86, 10),
#' GSdefList = GS, nperm = 3)
#' 
`singleDC` <- function(data, group, GSdefList, nperm){

   DI=rep(0,length(GSdefList))
   pvalue=rep(0,length(GSdefList))
   permv=matrix(0,length(GSdefList),nperm)

   for(j in 1:length(GSdefList)){
      permv[j,1]=compute.test.stat(data[GSdefList[[j]],],group)
   }
   
   print('Perm 1 completed')

   for(i in 2:nperm){
      data.perm=data[,sample(length(data[1,]))]

      for(j in 1:length(GSdefList)){
         permv[j,i]=compute.test.stat(data.perm[GSdefList[[j]],],group)
      }
   
      print(paste('Perm',i,'completed'))
   }

DI=permv[,1]
pvalue=rowSums(permv>=permv[,1])/length(permv[1,])

ret=NULL

names(DI)=names(GSdefList)
names(pvalue)=names(GSdefList)
dimnames(permv)[[1]]=names(GSdefList)
permnames=NULL
for(m in 1:nperm){
permnames[m]=paste('P',m,sep='')
}
dimnames(permv)[[2]]=permnames

ret$DI=DI
ret$pvalue=pvalue
ret$permv=permv

return(ret)
}


compute.test.stat = function(fixed.gs.data,group){

n.pairs<-choose(nrow(fixed.gs.data),2)
n.groups<-length(group)
csum.group<-cumsum(group)

corr.mats=NULL

if(n.groups==2){
corr.mats[[1]]=cor(t(fixed.gs.data[,1:(group[1])]))
corr.mats[[2]]=cor(t(fixed.gs.data[,(group[1]+1):(csum.group[2])]))}

if(n.groups==3){
corr.mats[[1]]=cor(t(fixed.gs.data[,1:(group[1])]))
corr.mats[[2]]=cor(t(fixed.gs.data[,(group[1]+1):(csum.group[2])]))
corr.mats[[3]]=cor(t(fixed.gs.data[,(csum.group[2]+1):(csum.group[3])])) }

if(n.groups>3){
  for(i in 1:n.groups){
 if(i==1){corr.mats[[i]]=cor(t(fixed.gs.data[,1:group[i]]))} 
else{corr.mats[[i]]=cor(t(fixed.gs.data[,(csum.group[i-1]+1):csum.group[i]]))}
  }}

vals=NULL
k=1
for(i in 1:(n.groups-1)){
   for(j in 2:n.groups){
      sq.diff.mat=(corr.mats[[i]]-corr.mats[[j]])^2
      vals[[k]]=sqrt((1/n.pairs)*sum(sq.diff.mat[upper.tri(sq.diff.mat)]))
      k=k+1
   }
}

test.stat=mean(vals)
return(test.stat)
}               


