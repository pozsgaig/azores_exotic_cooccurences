# converts NAs to zeros
NA20<- function(x)
{
    x[is.na(x)]<-0
    x}

#### function to calculate number of positive and negative edges and the sum of those weights per node ####
poznetlinks<-function(net, sign = c(-1, 1), weighted=F, weight_col=NULL)
{
    net_df<-as_long_data_frame(net)
    all_nodes<-V(net)$name
    all_res<-rep(0, length(all_nodes))
    names(all_res)<-all_nodes
    selected_df<-net_df[net_df$sign == sign, ]
    # only continues if there are lines with the selected sign value
    if (nrow(selected_df)>0)
    {
        selected_nodes<-unique(c(selected_df$from_name, selected_df$to_name))
        if(isFALSE(weighted)){
            res<-sapply(selected_nodes, function(x) nrow(selected_df[selected_df$from_name == x | selected_df$to_name == x, ]))
            names(res)<-selected_nodes
        } else
        {
            if(is.null(weight_col)){
                cat("You must provide a column to sum!")
                stop()} else
                {res<-sapply(selected_nodes, function(x) sum(selected_df[selected_df$from_name == x | selected_df$to_name == x, weight_col]))
                names(res)<-selected_nodes
                }
        }
        
        all_res[names(res)]<-res}
    
    # returns or all zeroes or calculated values
    return(all_res)
}


##### Function to calculate a lot of network properties. Output is a list. ####
netprop<-function(net){
    
    require(igraph)
    
    diameter<-igraph::diameter
    degree<-igraph::degree
    closeness<-igraph::closeness
    betweenness<-igraph::betweenness
    evcent<-igraph::evcent
    
    ##### Network properties ####
    
    ##### Node-related ####
    
    # Number of nodes
    num_nodes<-gorder(net)
    
    # Number of unconnected nodes
    num_isolated<-length(degree(net)[degree(net)==0])
    
    # Proportion of unconnected nodes to all nodes
    prop_isolated<-ifelse(length(degree(net)[degree(net)==0])==0, 0, 
                          length(degree(net)[degree(net)==0])/gorder(net))
    
    # mean reach in two steps
    mean_reach_2 <- mean((ego_size(net, 2)-1)/(gorder(net)-1))
    
    # mean reach in three steps
    # mean_reach_3 <- mean((ego_size(net, 3)-1)/(vcount(net)-1))
    
    ##### Edge=related ####
    # Number of edges
    edgenum<-gsize(net)
    
    
    # Number of positive edges
    poz_edgenum<-length(E(net)[E(net)$sign==1])
    
    # Number of negative edges
    neg_edgenum<-length(E(net)[E(net)$sign==-1])
    
    # Proportion of positive edges to all edges
    poz_edge_ratio<-ifelse(length(E(net)[E(net)$sign==1])==0, 0,
                           length(E(net)[E(net)$sign==1])/gsize(net))
    
    # Proportion of negative edges to all edges
    neg_edge_ratio<-ifelse(length(E(net)[E(net)$sign==-1])==0, 0,
                           length(E(net)[E(net)$sign==-1])/gsize(net))
    
    # Average weights
    mean_edge_weight<-mean(E(net)$weight)
    sd_edge_weight<-sd(E(net)$weight)
    
    #Average negative weights
    mean_neg_edge_weight<-{xx<-E(net)[E(net)$sign==-1]$weight
    if (length(xx)>0)
    {mean(xx, na.rm=T)} else
    {0}}
    sd_neg_edge_weight<-ifelse(mean_neg_edge_weight>0, sd(E(net)[E(net)$sign==1]$weight), 0)
    
    # Average positive weights
    mean_poz_edge_weight<-{xx<-E(net)[E(net)$sign==1]$weight
    if (length(xx)>0)
    {mean(xx, na.rm=T)} else
    {0}}
    
    sd_poz_edge_weight<-ifelse(mean_poz_edge_weight>0, sd(E(net)[E(net)$sign==1]$weight), 0)
    
    # Ratio of the mean of positive weights to the mean of all weights
    poz_edge_weight_ratio<-ifelse(mean(E(net)[E(net)$sign==1]$weight)==0, 0,
                                  mean(E(net)[E(net)$sign==1]$weight)/mean(E(net)$weight))
    
    # Ratio of the mean of negative weights to the mean of all weights
    neg_edge_weight_ratio<-ifelse(neg_edgenum==0, 0,
                                  mean(E(net)[E(net)$sign==-1]$weight)/mean(E(net)$weight))
    
    
    ##### Centrality-related #####
    
    
    # Mean degree
    mean_degree<-mean(igraph::degree(net))
    
    # Mean of the realised degree (the proportion of realised links to all possible links)
    realised_degree<-mean(igraph::degree(net)/centr_degree_tmax(net))
    degree_cent<-centralization.degree (net, normalized = F, loop=F)$centralization 
    normalized_degree_cent<-centralization.degree(net, normalized = T, loop=F)$centralization 
    
    
    # centralities
    
    # Eigenvector centrality
    mean_eig <- mean(centralization.evcent(net, normalized = F, directed=F)$vector)                    
    eig_cent<-centralization.evcent(net, normalized = F, directed=F)$centralization
    normalized_eig_cent<-centralization.evcent(net, normalized = T, directed=F)$centralization 
    
    # "Hub" centrality
    mean_hubs <- mean(hub.score(net)$vector)  
    
    # "Authority" centrality
    mean_authorities <- mean(authority.score(net)$vector)   
    
    # Closeness centrality is meaningful only for connected graphs. 
    # In disconnected graphs, consider using the harmonic centrality with harmonic_centrality
    if(is_connected(net))
    {
        # Closeness centrality
        mean_closeness <- mean(centralization.closeness(net, normalized = F, mode = "all")$res, na.rm=T)                  
        closeness_cent<-centralization.closeness(net, normalized = F, mode = "all")$centralization 
        normalized_closeness_cent<-centralization.closeness(net, normalized = T, mode = "all")$centralization 
    } else
        
    {
        # Closeness centrality
        mean_closeness <- mean(harmonic_centrality (net, normalized = T, mode = "all"), na.rm=T)                  
        closeness_cent<-centralize(harmonic_centrality(net, normalized = T, mode = "all"), 
                                   centr_clo_tmax(net), normalized = F)
        normalized_closeness_cent<-centralize(harmonic_centrality(net, normalized = T, mode = "all"), 
                                              centr_clo_tmax(net), normalized = T)
    }
    
    
    # Vertex betweenness centrality
    mean_betweenness <- mean(centralization.betweenness(net, normalized = F, directed=F)$res)              
    betweenness_cent<-centralization.betweenness(net, normalized = F, directed=F)$centralization 
    normalized_betweenness_cent<-centralization.betweenness(net, normalized = T, directed=F)$centralization 
    
    
    
    # Diameter: Measures the length of the longest geodesic
    netw_diam<-diameter(net)
    
    # Girth: Measures the length of the shortest circle in the graph
    netw_girth<-girth(net)$girth
    
    # Radius: Measures the smallest eccentricity in the graph
    
    netw_rad<-radius(net)
    
    
    # Minimum cut: The minimum number of edges to remove in order to split the graph into two clusters
    netw_min_cut<-min_cut(net)
    
    # Minimum cut proportion: The percent of minimum number of edges to remove in order to split the graph into two clusters to all edges
    netw_min_cut_prop<-min_cut(net)/gsize(net)
    
    # Mean distance: Calculates the mean distance between all node pairs in the graph
    netw_mean_dist<-mean_distance(net)
    
    # graph_adhesion: Gives the minimum edge connectivity.
    netw_adhesion<-edge_connectivity(net)
    
    # graph_assortativity: Measures the propensity of similar nodes to be connected.
    netw_assort<-assortativity_degree(net)
    
    
    # connectedness, connectance*2 (?)
    netw_edge_dens<-edge_density(net, loops = FALSE)
    
    
    # modularity 
    
    # clustering
    
    # number of motifs
    motif_num<-count_motifs(net, 3)
    
    # clique count
    max_num_cliques<-count_max_cliques(net)
    
    # component count
    comp_num<-count_components(net)
    
    # size of the largest clique
    clique_num<-clique_num(net)
    
    ##### Node properties ####   
    
    # https://igraph.org/r/doc/vertex_connectivity.html - also in mean
    
    # https://igraph.org/r/doc/edge_connectivity.html- also in mean
    
    # centralities
    deg <- igraph::degree(net)                               # Degree centrality
    rel_deg <- igraph::degree(net)/(igraph::gorder(net)-1)   # Degree centrality relative to the overall number of nodes 
    eig <- igraph::evcent(net)$vector                        # Eigenvector centrality
    hubs <- hub.score(net)$vector                            # "Hub" centrality
    authorities <- authority.score(net)$vector               # "Authority" centrality
    closeness <- ifelse(deg==0, 0, igraph::closeness(net, normalized = F))                      # Closeness centrality
    norm_closeness <- ifelse(deg==0, 0, igraph::closeness(net, normalized = T)) # Normalized closeness centrality
    betweenness <- igraph::betweenness(net, normalized = F)                  # Vertex betweenness centrality
    norm_betweenness <- igraph::betweenness(net, normalized = T)  # Vertex betweenness centrality
    node_page_rank <- page_rank(net)$vector                  # Google page rank
    
    
    # Number of negative edges linked - vulnerability
    vulnerability<-poznetlinks(net, sign=-1, weighted = F)
    # Sum of negative weights
    weighted_vulnerability<-poznetlinks(net, sign=-1, weighted = T, weight_col = "weight")
    # Relative number of negative edges linked - vulnerability
    rel_vulnerability<-ifelse(rel_deg==0, 0,
                              (poznetlinks(net, sign=-1, weighted = F)/rel_deg))
    
    
    # Number of positive edges linked
    resilience<-poznetlinks(net, sign=1)
    # Sum of positive weights
    weighted_resilience<-poznetlinks(net, sign=1, weighted = T, weight_col = "weight")
    # Relative number of negative edges linked - vulnerability
    rel_resilience<-ifelse(rel_deg==0, 0, 
                           (poznetlinks(net, sign=1, weighted = F)/rel_deg))
    
    # Reach in two steps
    reach_2 <- (ego_size(net, 2)-1)/(vcount(net)-1)
    
    # mean reach in three steps
    # reach_3 <- (ego_size(net, 3)-1)/(vcount(net)-1)
    
    ##### Edge properties - there are no edge properties other than weight. Maybe between groupness could be calculated #####
    
    ##### Output ####       
    node_df<-data.frame(deg, rel_deg, eig, authorities, closeness, norm_closeness, 
                        betweenness, norm_betweenness, vulnerability, 
                        weighted_vulnerability, rel_vulnerability, resilience, 
                        weighted_resilience, rel_resilience, reach_2)
    rownames(node_df)<-V(net)$name
    
    node_vars<-c("Degree", "Relative degree", "Eigenvector centrality", "Authorities centrality", "Closeness centrality", 
                 "Normalized closeness centrality", "Betweenness centrality", "Normalized betweenness centrality", 
                 "Vulnerability", "Weighted vulnerability",  "Relative vulnerability", "Resilience", "Weighted resilience", 
                 "Relative resilience", "Reach (2)")
    
    netw_df<- c(num_nodes, 
                num_isolated, prop_isolated, 
                mean_reach_2, 
                edgenum, poz_edgenum, neg_edgenum, 
                poz_edge_ratio, neg_edge_ratio,
                mean_edge_weight, sd_edge_weight, 
                mean_neg_edge_weight, sd_neg_edge_weight, 
                mean_poz_edge_weight, sd_poz_edge_weight,
                poz_edge_weight_ratio, neg_edge_weight_ratio, 
                mean_degree, realised_degree, degree_cent, normalized_degree_cent, 
                mean_eig, eig_cent, normalized_eig_cent, 
                mean_hubs, 
                mean_authorities, 
                mean_closeness, closeness_cent, normalized_closeness_cent, 
                mean_betweenness, betweenness_cent, normalized_betweenness_cent,
                netw_diam, netw_girth, netw_rad, 
                netw_min_cut, netw_min_cut_prop, netw_mean_dist, 
                netw_adhesion, netw_assort, 
                netw_edge_dens, 
                motif_num, max_num_cliques, comp_num, clique_num)
    
    names(netw_df)<-c("num_nodes", 
                      "num_isolated", "prop_isolated", 
                      "mean_reach_2", 
                      "edgenum", "poz_edgenum", "neg_edgenum", 
                      "poz_edge_ratio", "neg_edge_ratio",
                      "mean_edge_weight", "sd_edge_weight", 
                      "mean_neg_edge_weight", "sd_neg_edge_weight", 
                      "mean_poz_edge_weight", "sd_poz_edge_weight",
                      "poz_edge_weight_ratio", "neg_edge_weight_ratio", 
                      "mean_degree", "realised_degree", "degree_cent", "normalized_degree_cent", 
                      "mean_eig", "eig_cent", "normalized_eig_cent", 
                      "mean_hubs", 
                      "mean_authorities", 
                      "mean_closeness", "closeness_cent", "normalized_closeness_cent",
                      "mean_betweenness", "betweenness_cent", "normalized_betweenness_cent",
                      "netw_diam", "netw_girth", "netw_rad", 
                      "netw_min_cut", "netw_min_cut_prop", "netw_mean_dist", 
                      "netw_adhesion", "netw_assort", 
                      "netw_edge_dens", 
                      "motif_num", "max_num_cliques", "comp_num", "clique_num")
    
    
    netw_vars<-c("Number of nodes", 
                 "Number of isolates", "Propoprtion of isolates", 
                 "Mean reach (2)", 
                 "Number of edges", "Number of positive edges", "Number of negative edges", 
                 "Positive edge ratio", "Negative edge ratio",
                 "Mean edge weight", "Edge weight SD", 
                 "Mean negative edge weight", "Negative edge weight SD", 
                 "Mean positive edge weight", "Positive edge weight SD",
                 "Positive edge weight ratio", "Negative edge weight ratio", 
                 "Mean degree", "Realised degree", "Degree centrality", "Normalized degree centrality", 
                 "Mean eigenvector centrality", "Eigenvector centrality", "Normalized eigenvector centrality",
                 "Mean hubs centrality", 
                 "Mean authorities centrality", 
                 "Mean closeness centrality", "Closeness centrality", "Normalized closeness centrality",
                 "Mean betweenness centrality", "Betweenness centrality", "Normalized betweenness centrality", 
                 "Diameter", "Girth", "Radius", 
                 "Minimum cut", "Minimum cut proportion", "Mean distance", 
                 "Adhesion", "Assortativity", 
                 "Edge density", 
                 "Number of motifs", "Maximum number of cliques", "Number of components", "Number of cliques")
    
    standardised_netw_vars<-c("realised_degree", "prop_isolated", "poz_edge_ratio", "neg_edge_ratio", 
                              "mean_edge_weight", "mean_neg_edge_weight", "mean_poz_edge_weight", 
                              "poz_edge_weight_ratio", "neg_edge_weight_ratio",
                              "normalized_degree_cent", "normalized_eig_cent",
                              "normalized_closeness_cent", "normalized_betweenness_cent",
                              "netw_min_cut_prop")
    
    standardised_node_vars<-c("rel_deg", "rel_vulnerability", "rel_resilience", "norm_closeness", "norm_betweenness")
    
    reslist<-list(node_prop = node_df, node_vars = node_vars, netw_prop = netw_df, netw_vars = netw_vars, 
                  standardised_netw_vars=standardised_netw_vars, standardised_node_vars = standardised_node_vars)
    
    return(reslist)
    
}



set_weights<-function(x){
    apply(x, 1, function(y) {
        if(is.na(y[1])){aa = y[2]}
        if(is.na(y[2])){aa = y[1]}
        if(!is.na(y[1]) & !is.na(y[2])) {aa = (sum(y[1], y[2]))/2}
        aa
    })
}


ordered_ggpirate<-function(data, x_ax, y_ax, colour2)
{
    
    ordered_varpart1<-with(data, aggregate(get(y_ax), list(get(x_ax)), mean))
    colnum<-nrow(ordered_varpart1)
    ordered_varpart<-ordered_varpart1[order(ordered_varpart1[,"x"]), "Group.1"]
    data[[x_ax]]<-factor(data[[x_ax]], levels = ordered_varpart)
    frm <- reformulate(glue("{y_ax}"),glue("{x_ax}"))
    frm <- eval(parse(text=glue("{y_ax}~{x_ax}")))
    
    # checks if all variables have enough variance/elements for analysis
    dattest<-sapply(unique(as.character(data[[x_ax]])), function(x) sd(data[data[[x_ax]]==x, y_ax], na.rm = T))
    dattest[is.na(dattest)]<-0
    
    if (any(dattest==0))
    {
        
        p<-ggplot(data, aes_string(x = x_ax, y = y_ax )) +
            geom_pirate(aes_string(colour = colour2), bars = FALSE,
                        points_params = list(shape = 19, alpha = 0.2),
                        lines_params = list(size = 0.8))+
            scale_colour_manual(values = brewer.pal(colnum, "RdYlBu")[order(ordered_varpart1[,"x"])])+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  plot.margin = margin(1, 1, 1, 1, "cm"),)
    } else
    {
        wkt<-pairwise_wilcox_test(frm, data = data, p.adjust.method = "fdr")
        kt<-kruskal_test(data = data,  frm)
        
        my_comparisons <- lapply(1:nrow(wkt), function(x) {if(wkt[x, "p.adj"]<=0.05)
            unlist(wkt[x, c("group1", "group2")])
        })
        my_comparisons <- my_comparisons %>% discard(is.null)
        textcol<-ifelse(kt$p<=0.05, "red", "black")
        p<-ggplot(data, aes_string(x = x_ax, y = y_ax )) +
            geom_pirate(aes_string(colour = colour2), bars = FALSE,
                        points_params = list(shape = 19, alpha = 0.2),
                        lines_params = list(size = 0.8))+
            scale_colour_manual(values = brewer.pal(colnum, "RdYlBu")[order(ordered_varpart1[,"x"])])+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  plot.margin = margin(1, 1, 1, 1, "cm"),) + 
            stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
            stat_compare_means(label.x.npc = "left", label.y = (min(data[[y_ax]])-(abs(max(data[[y_ax]])*0.1))), color = textcol)     # Add global p-value
    }
    p 
}


radarchart_func<-function(x){
  data<-as.data.frame(t(apply(island_measured, 2, max)))
  data[1,1]<-gorder(full_network_B1_B2) # all species number in the metanetwork
  data<-rbind(data, rep(0, ncol(island_measured)))
  colnames(data)<-colnames(island_measured)
  print(x)
  data<-rbind(data, t(island_measured[x,]))
  data2<-data
  # rescale numbers, leave proportions as they are
  data2[,c(1,4,5)]<-apply(data[,c(1,4,5)], 2, scales::rescale, c(0,1))
  data2[1, c(2,3,6)]<-1 # setting maximum for prop data
  svglite(filename = paste0(x, "_radar.svg"), width = 10, height = 10,
          bg = "white")
  par(oma=c(0,0,0,0), mar=c(0,0,2,0), cex=3)
  radarchart(as.data.frame(data2), axistype=1 , 
             
             #custom polygon
             pcol=island_cols[x] , pfcol=adjustcolor(island_cols[x], alpha.f = 0.4) , plwd=4 , 
             
             #custom the grid
             cglcol="dark grey", cglty=1, axislabcol="dark grey", caxislabels=seq(0,100,25), cglwd=3,
             
             #custom labels
             vlcex=0.8, title = x)
  dev.off()
}


ordered_ggpirate<-function(data, x_ax, y_ax, colour2, sidelab, title)
{
  ordered_varpart1<-with(data, aggregate(get(y_ax), list(get(x_ax)), mean))
  colnum<-nrow(ordered_varpart1)
  
  # reviewer asked for changing the order
  #ordered_varpart<-ordered_varpart1[order(ordered_varpart1[,"x"]), "Group.1"]
  ordered_varpart<-ordered_varpart1[, "Group.1"]
  # change ends here
  
  data[[x_ax]]<-factor(data[[x_ax]], levels = ordered_varpart)
  frm <- reformulate(glue("{y_ax}"),glue("{x_ax}"))
  frm <- eval(parse(text=glue("{y_ax}~{x_ax}")))
  
  # checks if all variables have enough variance/elements for analysis
  dattest<-sapply(unique(as.character(data[[x_ax]])), function(x) sd(data[data[[x_ax]]==x, y_ax], na.rm = T))
  dattest[is.na(dattest)]<-0
  
  if (any(dattest==0))
  {
    
    p<-ggplot(data, aes_string(x = x_ax, y = y_ax )) +
      geom_pirate(aes_string(colour = colour2), bars = FALSE,
                  points_params = list(shape = 19, alpha = 0.2),
                  lines_params = list(size = 0.8))+
      # scale_colour_manual(values = brewer.pal(colnum, "RdYlBu")[order(ordered_varpart1[,"x"])])+
      scale_colour_manual(values = brewer.pal(colnum, "RdYlBu"))+
      labs(tag = title, y = sidelab)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            plot.margin = margin(1, 1, 1, 1, "cm"),)
  } else
  {
    wkt<-pairwise_wilcox_test(frm, data = data, p.adjust.method = "fdr")
    kt<-kruskal_test(data = data,  frm)
    
    my_comparisons <- lapply(1:nrow(wkt), function(x) {if(wkt[x, "p.adj"]<=0.05)
      unlist(wkt[x, c("group1", "group2")])
    })
    my_comparisons <- my_comparisons %>% discard(is.null)
    textcol<-ifelse(kt$p<=0.05, "red", "black")
    p<-ggplot(data, aes_string(x = x_ax, y = y_ax )) +
      geom_pirate(aes_string(colour = colour2), bars = FALSE,
                  points_params = list(shape = 19, alpha = 0.2),
                  lines_params = list(size = 0.8))+
      # scale_colour_manual(values = brewer.pal(colnum, "RdYlBu")[order(ordered_varpart1[,"x"])])+
      scale_colour_manual(values = brewer.pal(colnum, "RdYlBu"))+
      labs(tag = title, y = sidelab)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            plot.margin = margin(1, 1, 1, 1, "cm"),) + 
      # stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
      stat_compare_means(label.x.npc = "left", label.y = (min(data[[y_ax]])-(abs(max(data[[y_ax]])*0.1))), color = textcol)     # Add global p-value
  }
  p 
}

### https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R


#Use the function pairwise.adonis() with following arguments
#
#x = community table
#
#factors = a column or vector with all factors to be tested pairwise
#
#sim.function = which function to calculate the similarity matrix. eg 'daisy' or 'vegdist' default is 'vegdist'
# NOTE that if you wnat to use daisy, you need to install library 'cluster'
#
#sim.method = similarity method from daisy or vegdist: default is 'bray' 
#
#p.adjust.m = the p.value correction method, one of the methods supported by p.adjust(); default is 'bonferroni'
#
#The function will return a table with the pairwise factors, F-values, R^2, p.value and adjusted p.value
#
# load in runnig R session with: 
# source('pairwise.adonis.txt')
#

# example:
# data(iris)
# pairwise.adonis(iris[,1:4],iris$Species)
#
#[1] "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
#                    pairs    F.Model        R2 p.value p.adjusted sig
#1    setosa vs versicolor  552.51266 0.8493496   0.001      0.003   *
#2     setosa vs virginica 1001.54509 0.9108722   0.001      0.003   *
#3 versicolor vs virginica   91.82959 0.4837475   0.001      0.003   *

# similarity euclidean from vegdist and holm correction
# pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')

# similarity manhattan from daisy and bonferroni correction
# pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='daisy',sim.method='manhattan',p.adjust.m='bonferroni')
##############################
pairwise.anosim <- function(x, factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  R_stat = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = anosim(x1, factors[factors %in% c(co[1,elem],co[2,elem])])
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    R_stat = c(R_stat,ad$statistic);
    p.value = c(p.value,ad$signif)
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.1] <-'.'
  sig[p.adjusted <= 0.05] <-'*'
  sig[p.adjusted <= 0.01] <-'**'
  sig[p.adjusted <= 0.001] <-'***'
  
  pairw.res = data.frame(pairs,R_stat,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 


