require(igraph)
require(pcalg)
require(gRbase)

generate_conditional_data <- function(n_v, n_ph, n_p, n_c, sp_pph, sp_cph, sp_ph, sp_v,N) {
  # n_v : Number of variants
  # n_ph : Number of phenotypes
  # n_p : Number of pathways
  # n_c : Number of confounders
  # sp_pph : Sparsity of graph Pathways -> Phenotypes
  # sp_cph : Sparsity of graph Confounders -> Phenotypes
  # sp_ph : Sparsity of graph of Phenotypes
  # sp_v : Sparsity of graph : Variants -> Pathways
  
  n_variables = n_v + n_ph + n_p + n_c
  amat <- matrix(0,ncol=n_variables,nrow=n_variables)
  
  # amat(i,j) <=> i -> j
  names <- c()
  # Start by connecting variants to pathways
  n_cpaths <- floor(sp_v * n_p)
  if(n_v > 0) {
    for (i in 1:n_v) {
      names[i] <- paste("V",i)
      #Pick pathways at random
      n_paths <- 1 #sample(x=c(1:n_cpaths), size=1)
      pathway <- sample(x=c(1:n_p),size=n_paths)
      amat[i,pathway+n_v] <- runif(n=n_paths,min=-1,max=1)
    }
  }
  # Connect pathways to phenotypes
  n_connections = floor(sp_pph * n_ph)
  if (n_p>0) {
    for (i in 1:n_p) {
      names[i+n_v] <- paste("P",i)
      pathways <- sample(x=c(1:n_ph),size=n_connections,replace=FALSE)
      if(length(pathways)==0)
        next()
      for (j in pathways)
        amat[i+n_v,j+n_p+n_v] <- runif(n=1,min=-1,max=1)
    }
  }
  # Connect pehnotypes with each other. No cycles
  n_connections = floor(sp_ph * n_ph)
  for(i in 1:n_ph) {
    names[i+n_v+n_p] <- paste("Ph",i)
  }
  gNelDAG <- randomDAG(n=n_ph,prob=sp_ph,lB=-1,uB=1,V=names[(n_v+n_p+1):(n_v+n_p+n_ph)])
  dagamat <- get.adjacency(igraph.from.graphNEL(gNelDAG),attr="weight",sparse=FALSE)
  
  amat[(1+n_v+n_p):(n_ph+n_v+n_p),(1+n_v+n_p):(n_ph+n_v+n_p)] <- dagamat
  
  # Connect confounders to phenotypes
  n_connections=floor(sp_cph*n_ph)
  for(i in 1:(n_c)) {
    names[i+n_v+n_p+n_ph] <- paste("C",i)
    pathways <- sample(x=c(1:n_ph),size=n_connections,replace=FALSE)
    if(length(pathways)==0)
      next()
    for (j in pathways)
      amat[i+n_v+n_p+n_ph,j+n_v+n_p] <- runif(n=1,min=-1,max=1)
  }
  rownames(amat) <- colnames(amat) <- names
  #Now, build the corresponding DAG
  mydag <- graph.adjacency(amat,mode="directed",weighted=TRUE)
  reordering <- topological.sort(mydag,"out")
  cols <- (c(rep("green",n_v),rep("pink",n_p),rep("orange",n_ph),rep("darkgray",n_c)))[reordering]
  amat <- get.adjacency(mydag,attr="weight",sparse=FALSE)
  amat <- amat[reordering,reordering] ** 2
  #mydag <- graph.adjacency(amat,mode="directed",weighted=TRUE)
  
  #Generate data according to the DAG
  data <- rmvDAG(n=N,as(amat, "graphNEL"), back.compatible = TRUE)
  
  #Identify Hidden/Obs nodes
  obs <- hid <- c()
  i = 1
  for(name in colnames(amat)) {
    sp <- strsplit(name,split=" ",fixed=TRUE)
    
    if(is.element(sp[[1]][1],c("Ph","V")))
      obs <- c(obs,i)
    else
      hid <- c(hid,i)
    i <- i + 1
  }
  
  covmatFull <- cov(data)
  covmatObs <- cov(data[obs,obs])
  # Moralise the graph.
  morGraph <- moralize(graph.adjacency(amat,mode="directed",weighted=TRUE))
  
  list(dag=graph.adjacency(amat,mode="directed",weighted=TRUE),cols=cols,data=data,
       names=rownames(amat),obs=obs,hid=hid, ug=morGraph)
}