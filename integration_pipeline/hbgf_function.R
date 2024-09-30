# https://www.math.chalmers.se/Stat/Grundutb/GU/MSA220/S18/ClusterEnsembles.pdf

# https://github.com/cran/IntClust/blob/master/R/Functions.R

#' @title Hybrid bipartite graph formulation
#' 
#' @description Hybrid Bipartite Graph Formulation (HBGF) is a graph-based consensus multi-source clustering technique. The method builds a 
#' bipartite graph in which the two types of vertices are represented by the objects on one hand and the clusters of the partitions on the 
#' other hand. An edge is only present between an object vertex and a cluster vertex indicating that the object belongs to that cluster. 
#' The graph can be partitioned with Spectral clustering \insertCite{Ng2000}{IntClust}.
#' 
#' Adapted from IntClust ! Added faster algorithms for svd and kmeans.
#' 
#' @param ClusterEnsembles a obs x clusterings data frame with rownames
#' @param nstart how often to restart the kmeans algorithm with random centroids. Only used if it falls back due to the precomputed centroids cusing an error
#' @param graphPartitioning A character string indicating the preferred graph partitioning algorithm. For now only spectral clustering ("Spec") is implemented. Defaults to "Spec".
#' @param optimalk An estimate of the final optimal number of clusters. Default is 7.
#' @param verbose boolean: whether to print messages about progress
#' @return The returned value is a list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @references \insertRef{Fern2004}{IntClust} 
#' 
#' 
HBGF<-function(ClusterEnsembles,nstart = 50,graphPartitioning="Spec",optimalk=7,verbose=TRUE,maxDropPct = 0.5,useMiniBatch=TRUE){

  if(verbose) message(" Running data preparation ")
  # construct a graph G from the cluster ensemble: based on matlab code from Brodley and Fern (http://web.engr.oregonstate.edu/~xfern/)

  for(i in 1:ncol(ClusterEnsembles)){
    colnames(ClusterEnsembles)[i]=paste("Solution_",i,sep="")
  }
  
  drop_pct = 1-length(which(apply(ClusterEnsembles,2,function(x){length(unique(x))}) > 1)) / ncol(ClusterEnsembles)
  if(verbose) message(" Dropping ",round(drop_pct*100,2),"% of clusterings.")
  if(drop_pct > maxDropPct){stop("Would drop more clustering columns (",drop_pct,") than maxDropPct (",maxDropPct,"). Stopping.")}
  ClusterEnsembles = ClusterEnsembles[,which(apply(ClusterEnsembles,2,function(x){length(unique(x))}) > 1)]
  
  A=lapply(seq(ncol(ClusterEnsembles)),function(i) model.matrix(~as.factor(ClusterEnsembles[,i])-1))
  # A=Reduce(cbind,A)
  A= do.call(cbind,A) # faster !
  
  rownames(A)=rownames(ClusterEnsembles)
  for(c in 1:ncol(A)){
    colnames(A)[c]=paste("Cluster ",c,sep="")
  }
  
  #Need W or just the connectivity matrix A? Ferns Algorthym uses A
  if(graphPartitioning=="Spec"){
    if(verbose) message(" Calculate D from bipartite adjacency matrix A: ",dim(A)[1]," x ",dim(A)[2])
    D<-diag(sqrt(colSums(A)))
    if(verbose) message(" Calculate L from A and D ")
    L<-A%*%solve(D)
    
    if(verbose) message(" Running svd on matrix L: ",dim(L)[1]," x ",dim(L)[2])
    #SingularValues= base::svd(L,nu=optimalk,nv=optimalk)
    message(" Tol: ",1," , maxit: ",1e+8) # need to 
    SingularValues = irlba::irlba(L,nu=optimalk,nv=optimalk,maxit=1e+8,tol=1,work = optimalk*2,verbose=F) # much faster ! but gives bad results if tol is not lowerred to 1 and maxit (possibly) increased)
    U=SingularValues$u
    V=SingularValues$v
    U=U/matrix(rep(sqrt(rowSums(U^2)),optimalk),nrow=nrow(A),ncol=optimalk,byrow=FALSE)
    V=V/matrix(rep(sqrt(rowSums(V^2)),optimalk),nrow=ncol(A),ncol=optimalk,byrow=FALSE)
    permutated=gtools::permute(1:nrow(A))
    centers=U[permutated[1],,drop=FALSE]
    c=rep(0,nrow(A))
    c[permutated[1]]=2*optimalk
    
    #finsih this last for loop
    if(verbose) message(" Running for loop ")
    for(j in 2:optimalk){
      
      c=c+abs(U%*%t(centers[j-1,,drop=FALSE]))
      m = which.min(c)
      centers = rbind(centers,U[m,])
      c[m] = 2*optimalk
      
    }
    
    #k-means clustering
    # using base R
    if(!"ClusterR" %in% rownames(installed.packages()) | !useMiniBatch){
      if(verbose) message(" Running kmeans using base R implemenation.")
      clusterid=try(kmeans(x=rbind(U,V),centers=centers,iter.max=200),silent=TRUE)
      if(class(clusterid)=="try-error"){
        if(verbose) message("Falling back to random initialization with ",nstart," starts and iter.max=200 for k= ",nrow(centers)," clusters.")
        clusterid=try(kmeans(x=rbind(U,V),centers=nrow(centers),iter.max=200,nstart=nstart),silent=TRUE)
      }
      Clusters=clusterid$cluster[1:nrow(A)]
      
    # or using MiniBatchKmeans
    }else{
      #library(ClusterR)
      if(verbose) message(" Running kmeans using ClusterR::MiniBatchKmeans.")
      clusterid=try(ClusterR::MiniBatchKmeans(data=rbind(U,V),clusters=nrow(centers),num_init =nstart ,max_iters=100,CENTROIDS=centers,
                                              batch_size = min(nrow(rbind(U,V)),2000), init_fraction = 0.2, initializer = 'kmeans++', early_stop_iter = 10,verbose = FALSE),
                    silent=TRUE)
      if(class(clusterid)[1]=="try-error"){
        if(verbose) message("Falling back to random initialization with ",nstart," starts and iter.max=200 for k= ",nrow(centers)," clusters.")
        clusterid=try(ClusterR::MiniBatchKmeans(data=rbind(U,V),clusters=nrow(centers),num_init =nstart ,max_iters=100,
                                                batch_size = min(nrow(rbind(U,V)),2000), init_fraction = 0.2, initializer = 'kmeans++', early_stop_iter = 10,verbose = FALSE),
                      silent=TRUE)
      }
      Clusters = predict(object = clusterid, newdata = rbind(U,V))[1:nrow(A)]
    }
    
    # prepare output
    names(Clusters)=rownames(A)
    clusters=unique(Clusters)
    order=c()
    for(j in clusters){
      order=c(order,which(Clusters==j))
    }
    order.lab=as.character(order)
  }
  if(verbose) message(" Complete ")
  Out=list(DistM=NULL,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
  attr(Out,"method")="Ensemble"
  return(Out)
  
}
