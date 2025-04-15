TOPAS = function(network,
                 seeds,
                 expansion_steps = 2,
                 cores = 1){
  
  if (is.null(network) | ncol(network) < 2){
    stop('Network file is missing.', call. = FALSE)
  } 
  if (is.null(seeds)) {
    stop('Seeds are missing.', call. = FALSE)
  } 
  if (is.null(expansion_steps)) {
    expansion_steps = 2
  }
  if (is.null(cores)) {
    cores = parallel::detectCores()
  }
  if (class(cores) != 'numeric' | parallel::detectCores() < cores){
    stop('Please, introduce a numeric for cores.')
  }
  
  #########################################
  ########## STEP 1 - Full Seed Network
  message(paste0('STEP 1: Full Seed Network'))
  seeds_graph = igraph::simplify(igraph::graph_from_data_frame(network, directed = FALSE),
                                 remove.multiple = TRUE, remove.loops = TRUE)
  seeds = as.character(seeds)
  
  #########################################
  ########## STEP 2 - Largest Connected Module
  message(paste0('STEP 2: Largest Connected Module (using ',cores,ifelse(cores>1,' cores)', ' core)')))
  
  .extract_lcc = function(graph){
    ## This function takes a graph in input,
    ## computes all connected components, and 
    ## outputs the one with the largest number of seeds
    best_cc = 1
    max_seeds = 0
    graph_lcc = igraph::components(graph)
    for (i in 1:igraph::count_components(graph)){
      seeds.cc = sum(seeds %in% names(graph_lcc$membership[graph_lcc$membership==i]))
      if (seeds.cc > max_seeds){
        best_cc = i
        max_seeds = seeds.cc
      }
    }
    return(
      igraph::induced_subgraph(
        graph, 
        which(igraph::V(graph)$name %in%
                names(graph_lcc$membership[graph_lcc$membership==best_cc])))
    )
  }
  
  ## Extract the largest seed connected component
  graph_lcc = .extract_lcc(seeds_graph)
  ## Update seeds
  seeds = seeds[seeds %in% igraph::V(graph_lcc)$name]
  
  ## To speed up the process, extract first the distances between seeds.
  ## In the first steps, compute shortest paths between any two of seeds
  ## only if distance is lower than expansion step parameter
  d_m = 
    igraph::distances(
      graph_lcc,
      v = seeds,
      to = seeds,
      weights = NULL
    )
  
  .sp_compute = function(source_v,graph){
    ## This function computes all shortest paths between 
    ## a seed and all other seeds and returns only those
    ## which are shorter than input "expansion_steps"
    d = d_m[source_v,]
    dest_v = names(d[d>1 & d<=(expansion_steps+1)])
    
    sp = igraph::all_shortest_paths(
      graph,
      from = which(igraph::V(graph)$name == source_v),
      to = which(igraph::V(graph)$name %in% dest_v),
      weights = NULL)[['res']]
    
    return(
      unlist(lapply(sp,function(s){
        names(s[[2:(length(s)-1)]])}))
    ) 
  }
  
  ## Proceed with computation of shortest path
  ## Set up parallelization, either on 1 or more cores
  pb = utils::txtProgressBar(min=0, max=length(seeds), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  cl = parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  boot = foreach::foreach(i = seeds, .options.snow = opts)
  connectors = unique(
    unlist(
      foreach::`%dopar%`(boot, .sp_compute(i, graph_lcc)))
  )
  parallel::stopCluster(cl)
  
  ## Extract a subgraph composed of seeds and potential connectors
  seeds_graph = 
    igraph::induced_subgraph(
      graph_lcc, 
      unique(c(which(igraph::V(graph_lcc)$name %in% c(seeds,connectors)))))
  
  ## If no edge is retrieved, the output is NULL
  if(igraph::ecount(seeds_graph)<1){
    message('\nNo module found!')
    return(NULL)
  }
  
  ## If more than one component exists, take the largest seed component
  if(igraph::count_components(seeds_graph)>1){
    seeds_graph = .extract_lcc(seeds_graph)
  }
  
  ## Update seeds
  seeds = seeds[seeds %in% igraph::V(seeds_graph)$name]
  
  #########################################
  ########## STEP 3 - Pruned Module
  message(paste0('\nSTEP 3: Pruned Module'))
  
  ## Execute Random Walk With restart from seeds
  df = data.frame(v=igraph::V(seeds_graph)$name,val=0)
  df$val = ifelse(df$v %in% seeds, 1, 0)
  restart_list = as.matrix(df$val)
  rownames(restart_list) = df$v
  
  set.seed(1996)
  stationary_probability = 
    dnet::dRWR(seeds_graph,
               normalise = 'row',
               setSeeds = restart_list,
               restart = 0.75,
               verbose=F)
  
  df$stationary_probability = as.numeric(stationary_probability)
  
  set.seed(1996)
  ## Randomize the alphabetical order first, it adds nothing if all probabilities vary
  df = df[sample(nrow(df)),]
  ## order by probability afterwards
  df = df[order(df$stationary_probability),]; rownames(df) = NULL 
  
  ## Iteratively test connectors, from the one with lowest stationary probability to highest.
  for (v in df[df$val==0,'v']){
    if (!v %in% igraph::V(seeds_graph)$name) {
      next
    }
    
    ## Temporary delete a candidate connector from the module.
    tmp_graph = igraph::delete_vertices(
      seeds_graph,
      which(igraph::V(seeds_graph)$name==v))
    
    ## If removing a vertex generated more than 1 connected components,
    ## then check if seed coverage decreases. 
    ## Else remove vertex permanently.
    if (igraph::count_components(tmp_graph)>1){
      best_cc = 1
      max_seeds = 0
      tmp_graph_lcc = igraph::components(tmp_graph)
      
      for (i in 1:igraph::count_components(tmp_graph)){
        seeds.cc = sum(seeds %in% names(tmp_graph_lcc$membership[tmp_graph_lcc$membership==i]))
        if (seeds.cc > max_seeds){
          best_cc = i
          max_seeds = seeds.cc
        }
      }
      
      if ((length(seeds) - max_seeds)==0){
        tmp_graph_lcc = 
          igraph::induced_subgraph(
            tmp_graph, 
            which(igraph::V(tmp_graph)$name %in% names(tmp_graph_lcc$membership[tmp_graph_lcc$membership==best_cc])))
        
        seeds_graph = rlang::duplicate(tmp_graph_lcc)
      } 
    } else {
      seeds_graph = rlang::duplicate(tmp_graph)
    }
  }
  return(as.data.frame(igraph::get.edgelist(seeds_graph)))
}
