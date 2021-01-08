# - split whole sampleXotu table into plot level sampleXotu table ====
split_sampleXotu_plot <- function(sampleXotu = bac_sampleXotu, plot = sample_metadata$plot, occurNum = 0){
  res <- list()
  for(i in unique(plot)){
    
    temp_df <- sampleXotu[plot %in% i, ]
    temp_df <- temp_df[,colSums(temp_df) > 0]
    temp_df <- temp_df[, colSums(temp_df > 0) >= occurNum]
    res[[i]] <- temp_df
    
  }
  return(res)
}


# - split whole sampleXsoil table into plot level sampleXsoil table ====
split_sampleXsoil_plot <- function(sampleXsoil = soil9_stan, plot = sample_metadata$plot){
  res <- list()
  for(i in unique(plot)){
    
    temp_df <- sampleXsoil[plot %in% i, ]
    res[[i]] <- temp_df
  }
  return(res)
}


# - Extract coefficient matrix from SPIEC.easi result and convert it into network graph ====
coff_to_igraph <- function(spiec, sym_beta = F){
  if(!require(SpiecEasi)) stop("Please install package 'SpiecEasi'!")
  if(!require(igraph)) stop("Please install package 'igraph'!")
  if(sym_beta){
    adj_mat <- as.matrix(symBeta(getOptBeta(spiec)))
    rownames(adj_mat) <- colnames(adj_mat) <- colnames(spiec$select$est$data)
    adj_mat_temp <- adj_mat
    adj_mat_T <- t(adj_mat)
    # assign cooperation (+/+) 1, 
    adj_mat[(adj_mat_temp >0) & (adj_mat_T >0)] <- 1
    # competition (-/-) 2, 
    adj_mat[adj_mat_temp <0 & adj_mat_T <0] <- 2
    if(!isSymmetric.matrix(adj_mat)) stop("Error!!! Check your codes!")
    network <- graph.adjacency(adj_mat, mode = "upper", diag = F, weighted = T)
    # assign different color to different edge type
    E(network)[E(network)$weight == 1]$color <- "#346830" # COLOR: green
    E(network)[E(network)$weight == 2]$color <- "#ce3536" # COLOR: red
    
    return(network)
  } else {
    adj_mat <- as.matrix(getOptBeta(spiec))
    rownames(adj_mat) <- colnames(adj_mat) <- colnames(spiec$select$est$data)
    adj_mat_temp <- adj_mat
    adj_mat_T <- t(adj_mat)
    # assign cooperation (+/+) 1, 
    adj_mat[(adj_mat_temp >0) & (adj_mat_T >0)] <- 1
    # competition (-/-) 2, 
    adj_mat[adj_mat_temp <0 & adj_mat_T <0] <- 2
    # commensalism (+/0) 3, 
    adj_mat[adj_mat_temp>0 & adj_mat_T == 0] <- 3
    adj_mat[adj_mat_temp == 0 & adj_mat_T >0] <- 3
    # amensalism (-/0) 4,
    adj_mat[adj_mat_temp <0 & adj_mat_T == 0] <- 4
    adj_mat[adj_mat_temp == 0 & adj_mat_T < 0] <- 4
    # exploitation (+/-) 5
    adj_mat[adj_mat_temp < 0 & adj_mat_T > 0] <- 5
    adj_mat[adj_mat_temp > 0 & adj_mat_T < 0] <- 5
    
    if(!isSymmetric.matrix(adj_mat)) stop("Error!!! Check your codes!")
    network <- graph.adjacency(adj_mat, mode = "upper", diag = F, weighted = T)
    # assign different color to different edge type
    E(network)[E(network)$weight == 1]$color <- "#346830" # COLOR: green
    E(network)[E(network)$weight == 2]$color <- "#ce3536" # COLOR: red
    E(network)[E(network)$weight == 3]$color <- "#f3a001" # COLOR: yellow
    E(network)[E(network)$weight == 4]$color <- "#f07579" # COLOR: pink
    E(network)[E(network)$weight == 5]$color <- "#211a28" # COLOR: black
    
    return(network)
  }
}


# - Remove isolated vertex ====
rm_isolated_vertex <- function(igraph){
  bad_vs <- V(igraph)[degree(igraph) == 0]
  res <- delete.vertices(igraph, bad_vs)
  return(res)
}


# - Define a function to calculate network size and average degree ====
cal_net_property <- function(x) {
  network_size <- length(V(x))
  average_degree <- mean(igraph::degree(x))
  res <- data.frame(network_size, average_degree)
  return(res)
}


cal_net_property <- function(x) {
  network_size <- length(V(x))
  average_degree <- mean(igraph::degree(x))
  
  transi <- transitivity(x, type = "global")
  # average_cc <- mean(transitivity(x, type = "local", isolates = "zero"))
  average_gd <- mean_distance(x)
  
  modu <- NA
  try(modu <- modularity(cluster_fast_greedy(x, weights = NULL, merges = F, modularity = T)))
  
  res <- data.frame(network_size, average_degree, average_gd, modu, transi)
  return(res)
}


# - calculate plot level alpha diversity and beta diversity ====
cal_plot_alpha_beta_mean <- function(sampleXotu, col2rowname = F, method = "jaccard", tre = bac_tre){
  
  if(col2rowname){
    sampleXotu <- as.data.frame(sampleXotu)
    rownames(sampleXotu) <- sampleXotu[,1]
    sampleXotu <- sampleXotu[,-1]
  }
  
  sampleXotu <- sampleXotu[,colSums(sampleXotu) > 0]
  
  alpha <- rowSums(sampleXotu > 0)
  alpha_mean <- mean(alpha)
  alpha_sd <- sd(alpha)
  if(method == "jaccard"){
    beta <- vegan::vegdist(sampleXotu, method = "jaccard", binary = T)
    beta_mean <- mean(beta)
    beta_sd <- sd(beta)
  } else if(method == "bray"){
    beta <- vegan::vegdist(sampleXotu, method = "bray")
    beta_mean <- mean(beta)
    beta_sd <- sd(beta)
  } else if(method == "unifrac"){
    if(!require(PhyloMeasures)) stop("Please install package 'PhyloMeasures'!")
    beta <- unifrac.query(tre, sampleXotu)
    beta_mean <- mean(beta)
    beta_sd <- sd(beta)
  } else {
    stop("method must be in jaccard, bray, unifrac")
  }
  
  res <- data.frame(key = c("OTUnum", method), mean = c(alpha_mean, beta_mean), 
                    sd = c(alpha_sd, beta_sd))
  return(res)
}


# - calculate pairwise spearman correlation ====
cal_rp_from_sampleXotu <- function(sampleXotu){
  if(!require(WGCNA)) stop("Please install package 'WGCNA'!")
  
  res <- list()
  df <- corAndPvalue(sampleXotu, method = "spearman")
  res$r <- df$cor
  res$p <- df$p
  return(res)
}


# - convert correlation matrix into igraph object ====
corr_to_igraph <- function(rp, r = 0.6, p = 0.05){
  if(!require(igraph)) stop("Please install package 'igraph'!")
  rp$r[is.na(rp$r)] <- 0
  rp$p[is.na(rp$p)] <- 0
  adj_mat <- rp$r
  adj_mat[(abs(rp$r) < r) | (rp$p > p)] <- 0 
  
  network <- graph_from_adjacency_matrix(adj_mat, mode = "upper", diag = F, weighted = T)
  # assign different color to different edge type
  E(network)[E(network)$weight > 0]$color <- "#346830" # COLOR: green
  E(network)[E(network)$weight < 0]$color <- "#ce3536" # COLOR: red
  
  return(network)
}


# - Define a function to calculate network properties ====
cal_spearNet_property <- function(x) {
  if(!require(igraph)) stop("Please install package 'igraph'!")
  
  network_size <- length(V(x))
  total_edge <- length(E(x))
  
  # cooperation edge
  pp_edge <- sum(E(x)$weight > 0)
  # competition edge
  nn_edge <- sum(E(x)$weight < 0)
  # PN_ratio <- pp_edge / nn_edge
  density <- edge_density(x, loops = F)
  average_degree <- mean(igraph::degree(x))
  transi <- transitivity(x, type = "global")
  average_gd <- mean_distance(x)
  
  modu <- NA
  try(modu <- modularity(cluster_fast_greedy(x, weights = NULL, merges = F, modularity = T)))
  
  res <- data.frame(network_size, total_edge, pp_edge, nn_edge, density, average_degree, average_gd, modu, transi)
  
  return(res)
  
}


cal_spearNet_property_use_diffthreshold <- function(ls, r = c(0.6,0.7,0.8,0.9), p = 0.05){
  library(tidyverse)
  res <- c()
  
  for (i in r) {
    ig <- lapply(ls, corr_to_igraph, r = i, p = p)
    ig_rmIso <- lapply(ig, rm_isolated_vertex)
    ig_rmIso_property <- lapply(ig_rmIso, cal_spearNet_property)
    ig_rmIso_property <- do.call(rbind, ig_rmIso_property)
    ig_rmIso_2property <- ig_rmIso_property %>% rownames_to_column(var = "plot") %>%
      select(plot, network_size, average_degree) %>%
      gather(key = "key", value = "value", network_size, average_degree) %>%
      mutate(r = i, key = as_factor(key))
    res <- rbind(res, ig_rmIso_2property)
  }
  return(res)
}


# - calculate p value for SparCC ====
cal_p_for_sparcc <- function(sampleXotu, times = 99, cores = 10){
  if(!require(SpiecEasi)) stop("Please install package 'SpiecEasi'!")
  if(!require(abind)) stop("Please install package 'abind'!")
  if(!require(doParallel)) stop("Please install package 'doParallel'!")
  
  
  obs_sparcc <- sparcc(sampleXotu)$Cor
  
  cl <- makeCluster(cores, type = "SOCK")
  registerDoParallel(cl)
  
  null_sparcc <- foreach(i=1:times, .packages = "SpiecEasi") %dopar% {
    null_sampleXotu <- apply(sampleXotu, 2, sample, replace = T)
    sparcc(null_sampleXotu)$Cor
  }
  stopCluster(cl)
  
  null_sparcc <- abind(null_sparcc, along = 3)
  obs_null <- abind(obs_sparcc, null_sparcc, along = 3)
  
  obs_rank <- apply(obs_null, c(1,2), rank)[1,,]
  obs_rank <- obs_rank / times
  
  sparcc_p <- obs_rank
  sparcc_p[obs_rank <= 0.025 | obs_rank >= 0.975] <- 0
  sparcc_p[sparcc_p != 0] <- 1
  return(sparcc_p)
}


# - convert sparcc correlation matrix into igraph object ====
sparcc_to_igraph <- function(spcc, r = 0.6, p = 0.05){
  if(!require(igraph)) stop("Please install package 'igraph'!")
  spcc$Cor[is.na(spcc$Cor)] <- 0
  adj_mat <- spcc$Cor
  adj_mat[(abs(spcc$Cor) < r) | spcc$p > p] <- 0 
  
  network <- graph_from_adjacency_matrix(adj_mat, mode = "upper", diag = F, weighted = T)
  # assign different color to different edge type
  E(network)[E(network)$weight > 0]$color <- "#346830" # COLOR: green
  E(network)[E(network)$weight < 0]$color <- "#ce3536" # COLOR: red
  
  return(network)
}


cal_sparccNet_property_use_diffthreshold <- function(ls, r = c(0.6,0.7,0.8,0.9), p = 0.05){
  if(!require(tidyverse)) stop("Please install package 'tidyverse'!")
  res <- c()
  
  for (i in r) {
    ig <- lapply(ls, sparcc_to_igraph, r = i, p = p)
    ig_rmIso <- lapply(ig, rm_isolated_vertex)
    ig_rmIso_property <- lapply(ig_rmIso, cal_spearNet_property)
    ig_rmIso_property <- do.call(rbind, ig_rmIso_property)
    ig_rmIso_2property <- ig_rmIso_property %>% rownames_to_column(var = "plot") %>%
      select(plot, network_size, average_degree) %>%
      gather(key = "key", value = "value", network_size, average_degree) %>%
      mutate(r = i, key = as_factor(key))
    res <- rbind(res, ig_rmIso_2property)
  }
  return(res)
}


# - generate random network and calcuate network properties ====
collect_randNet_property <- function(ig, rand_num = 999){
  if(!require(igraph)) stop("Please install package 'igraph'!")
  if(!require(tidyverse)) stop("Please install package 'igraph'!")
  
  v_num <- length(V(ig))
  e_num <- length(E(ig))
  
  rand_ig_pro <- list()
  for(i in 1:rand_num){
    rand_ig_pro[[i]] <- cal_net_property(sample_gnm(v_num, e_num))
  }
  
  df_pro <-  rand_ig_pro %>%
    do.call(rbind, .)
  
  res <- list(mean = summarize_all(df_pro, mean), sd = summarise_all(df_pro, sd))
  return(res)
  
}


# - generate network for guild ====
generate_net_for_guild <- function(sampleXotu_l, guild2otu, r = 0.6, p = 0.05){
  
  res <- list()
  
  for (plot in names(sampleXotu_l)){
    res[[plot]] <- list()
    for(guild in names(guild2otu)){
      otus <- intersect(colnames(sampleXotu_l[[plot]]), guild2otu[[guild]])
      if(length(otus) < 5){
        net_pro <- NULL
      } else {
        df <- sampleXotu_l[[plot]][,otus]
        rp <- cal_rp_from_sampleXotu(df)
        ig <- corr_to_igraph(rp, r = r, p = p)
        obs_net <- cal_net_property(ig)
        rand_net <- collect_randNet_property(ig)
        rand_net <- do.call(rbind, rand_net)
        
        net_pro <- rbind(obs_net, rand_net)
      }
      res[[plot]][[guild]] <- net_pro
    }
  }
  return(res)
}


# - cal observed and random network properties ====
cal_obsNet_randNet_property <- function(ig){
  obs_net <- cal_net_property(ig)
  rand_net <- collect_randNet_property(ig, rand_num = 99)
  rand_net <- do.call(rbind, rand_net)
        
  net_pro <- rbind(obs_net, rand_net)
  return(net_pro)
}
