require(igraph)
require(pcalg)
require(gRbase)
generate_data <- function(n_ph, n_c, sp_cph, sp_ph, N) {
  # n_ph : Number of Variables 
  # n_c : Number of confounders 
  # sp_cph : Sparsity of grpah Confounders -> Variables 
  # sp_ph : Sparsity of graph of Variables
  
  n_variables = n_ph + n_c
  n_p <- 0
  n_v <- 0
  amat <- matrix(0,ncol=n_variables,nrow=n_variables)
  
  # amat(i,j) <=> i -> j
  names <- c()
  # Connect variables with each other. No cycles
  n_connections = floor(sp_ph * n_ph)
  for(i in 1:n_ph) {
    names[i] <- paste("V",i)
  }
  gNelDAG <- randomDAG(n=n_ph,prob=sp_ph,lB=-1,uB=1,V=names[(n_p+1):(n_p+n_ph)])
  dagamat <- get.adjacency(igraph.from.graphNEL(gNelDAG),attr="weight",sparse=FALSE)
  
  amat[(1):(n_ph),(1):(n_ph)] <- dagamat
  
  # Connect confounders to phenotypes
  n_connections=floor(sp_cph*n_ph)
  for(i in 1:(n_c)) {
    names[i+n_v+n_p+n_ph] <- paste("C",i)
    pathways <- sample(x=c(1:n_ph),size=n_connections,replace=FALSE)
    if(length(pathways)==0)
      next()
    for (j in pathways)
      amat[i+n_p+n_ph,j+n_p] <- runif(n=1,min=-1,max=1)
  }
  rownames(amat) <- colnames(amat) <- names
  #Now, build the corresponding DAG
  mydag <- graph.adjacency(amat,mode="directed",weighted=TRUE)
  reordering <- topological.sort(mydag,"out")
  cols <- (c(rep("orange",n_ph),rep("darkgray",n_c)))[reordering]
  amat <- get.adjacency(mydag,attr="weight",sparse=FALSE)
  amat <- amat[reordering,reordering] ** 2
  #Generate data according to the DAG
  data <- rmvDAG(n=N,as(amat, "graphNEL"),back.compatible = TRUE)
  
  #Identify Hidden/Obs nodes
  obs <- hid <- c()
  
  i = 1
  for(name in colnames(amat)) {
    sp <- strsplit(name,split=" ",fixed=TRUE)
    
    if(is.element(sp[[1]][1],c("V")))
      obs <- c(obs,i)
    else
      hid <- c(hid,i)
    i <- i + 1
  }
  
  covmatFull <- cov(data)
  covmatObs <- cov(data[obs,obs])
  
  # Moralise the graph.
  morGraph <- moralize(graph.adjacency(amat,mode="directed",weighted=TRUE))
  
  list(amat=amat, dag=graph.adjacency(amat,mode="directed",weighted=TRUE),cols=cols,ord=reordering,data=data,
       names=rownames(amat),obs=obs,ug=morGraph)
}