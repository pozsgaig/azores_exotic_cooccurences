#### 1. Setting the scene ####

##### 1.1. Loading necessary libraries #####

library(rstudioapi) # for selecting working directory
library(igraph) # for network analysis
library(SpiecEasi) # for making the networks
library(phyloseq) # for the phyloseq object
library(glue)
library(rstatix) # for stats in a tibble format
library(tidyverse)
library(Rsenal) # for merging data frame of node properties with igraph nodes (makeVertexAtt)
library(plotrix) # for calculating SE - std.error()
library(xlsx) #writing to xlsx
library(vegan)
library(reshape2)
library(psych) # for correlation test
library(MuMIn) # for regression

# for plotting
library(RColorBrewer)
library(pals) # for palettes with many colour categories
library(ggpirate)
library(fmsb) # for radar charts
library(ggord) # for ordination ggplots
library(corrplot) # for network property correlations
library(ggpubr)
library(wesanderson)
library(waffle) # for little boxes showing proportions
library(grid) #ggplots on one page
library(gridBase)
library(gridExtra)
library(svglite) # for svg export with text elements

##### 1.2. Selecting working directory #####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


##### 1.3. Loading necessary functions #####

source("Frontiers_exotic _species_cooccurrence_functions.R")

##### 1.4. Loading data #####
load("data_4_exotic_paper.RDA")

# listing natives and exotics 
NAT<-unlist(split(sp_taxonomy$Dev_MF, sp_taxonomy$IND_NON_IND)[2:3])
EXO<-unlist(split(sp_taxonomy$Dev_MF, sp_taxonomy$IND_NON_IND)[4:5])
END<-unlist(split(sp_taxonomy$Dev_MF, sp_taxonomy$N_E_I)[1])

# listing insects and arachnids
INSECT<-sp_taxonomy[sp_taxonomy$Class=="Insecta", "Dev_MF"]
ARACH<-sp_taxonomy[sp_taxonomy$Class=="Arachnida", "Dev_MF"]


#################################################################################
#################################################################################
#################################################################################



#### 2. Preparing networks ####
##### 2.1. Preparing BALA 1 datasets and network #####
B1_species_mat<-with(BALA1_NFF, tapply(Total, list(Dev_MF, Site_code), sum))
B1_species_mat<-NA20(B1_species_mat)

# applying filters for rare species
B1_species_mat<-B1_species_mat[rowSums(B1_species_mat)>10,]
B1_binary_sp_mat<-B1_species_mat
B1_binary_sp_mat[B1_binary_sp_mat>0]<-1
B1_species_mat<-B1_species_mat[rowSums(B1_binary_sp_mat)>3,colSums(B1_binary_sp_mat)>3]



B1_taxonomy_mat<-sp_taxonomy[sp_taxonomy$Dev_MF %in% rownames(B1_species_mat), c("Dev_MF", "Genus", "Family", "Order")]
rownames(B1_taxonomy_mat)<-B1_taxonomy_mat$Dev_MF
B1_taxonomy_mat$Dev_MF<-NULL
B1_taxonomy_mat<-as.matrix(B1_taxonomy_mat)
B1_taxonomy_mat<-B1_taxonomy_mat[rownames(B1_species_mat),]

B1_site_env_data<-site_env_data[site_env_data$Site.code %in% colnames(B1_species_mat), ]
rownames(B1_site_env_data)<-B1_site_env_data$Site.code
B1_site_env_data$Site.code<-NULL
B1_site_env_data<-B1_site_env_data[colnames(B1_species_mat),]

identical(rownames(B1_site_env_data), colnames(B1_species_mat))
identical(rownames(B1_taxonomy_mat), rownames(B1_species_mat))

B1_phyloseq_object = phyloseq(otu_table(B1_species_mat, taxa_are_rows = T), tax_table(B1_taxonomy_mat), 
                              sample_data(B1_site_env_data))

set.seed(123)
B1_net <- spiec.easi(B1_phyloseq_object, method='slr', r=5, lambda.min.ratio=1e-2,
                     nlambda=20, pulsar.params=list(rep.num=50))

# https://bytemeta.vip/repo/zdk123/SpiecEasi/issues/198
B1_se.icov <- B1_net$select$est$icov[[getOptInd(B1_net)]]
B1_se.cov <- solve(B1_se.icov)
B1_weighted.adj.mat <- B1_se.cov*getRefit(B1_net)
B1_ig_net <- adj2igraph(B1_weighted.adj.mat,  vertex.attr=list(name=taxa_names(B1_phyloseq_object)))



##### 2.2. Preparing BALA 2 datasets and network #####
B2_species_mat<-with(BALA2_NFF, tapply(Total, list(Dev_MF, Site_code), sum))
B2_species_mat<-NA20(B2_species_mat)

# applying filters for rare species
B2_species_mat<-B2_species_mat[rowSums(B2_species_mat)>10,]
B2_binary_sp_mat<-B2_species_mat
B2_binary_sp_mat[B2_binary_sp_mat>0]<-1
B2_species_mat<-B2_species_mat[rowSums(B2_binary_sp_mat)>3,colSums(B2_binary_sp_mat)>3]


B2_taxonomy_mat<-sp_taxonomy[sp_taxonomy$Dev_MF %in% rownames(B2_species_mat), c("Dev_MF", "Genus", "Family", "Order")]
rownames(B2_taxonomy_mat)<-B2_taxonomy_mat$Dev_MF
B2_taxonomy_mat$Dev_MF<-NULL
B2_taxonomy_mat<-as.matrix(B2_taxonomy_mat)
B2_taxonomy_mat<-B2_taxonomy_mat[rownames(B2_species_mat),]

B2_site_env_data<-site_env_data[site_env_data$Site.code %in% colnames(B2_species_mat), ]
rownames(B2_site_env_data)<-B2_site_env_data$Site.code
B2_site_env_data$Site.code<-NULL
B2_site_env_data<-B2_site_env_data[colnames(B2_species_mat),]

identical(rownames(B2_site_env_data), colnames(B2_species_mat))
identical(rownames(B2_taxonomy_mat), rownames(B2_species_mat))

B2_phyloseq_object = phyloseq(otu_table(B2_species_mat, taxa_are_rows = T), tax_table(B2_taxonomy_mat), 
                              sample_data(B2_site_env_data))

set.seed(123)
B2_net <- spiec.easi(B2_phyloseq_object, method='slr', r=5, lambda.min.ratio=1e-2,
                     nlambda=20, pulsar.params=list(rep.num=50))
# https://bytemeta.vip/repo/zdk123/SpiecEasi/issues/198
B2_se.icov <- B2_net$select$est$icov[[getOptInd(B2_net)]]
B2_se.cov <- solve(B2_se.icov)
B2_weighted.adj.mat <- B2_se.cov*getRefit(B2_net)
B2_ig_net <- adj2igraph(B2_weighted.adj.mat,  vertex.attr=list(name=taxa_names(B2_phyloseq_object)))


##### 2.3. Making full network and setting attributes #####
full_network_B1_B2<-igraph::union(B1_ig_net, B2_ig_net)


netweights_data<-as.data.frame(get.edge.attribute(full_network_B1_B2))
E(full_network_B1_B2)$weight<-set_weights(netweights_data)
full_network_B1_B2<-delete_edge_attr(full_network_B1_B2, "weight_1")
full_network_B1_B2<-delete_edge_attr(full_network_B1_B2, "weight_2")


E(full_network_B1_B2)$sign<-ifelse(E(full_network_B1_B2)$weight<0, -1, 1)
E(full_network_B1_B2)$weight<-abs(E(full_network_B1_B2)$weight)
E(full_network_B1_B2)$color<-ifelse(E(full_network_B1_B2)$sign==-1, "blue", "red")
V(full_network_B1_B2)$degree<-igraph::degree(full_network_B1_B2)


#################################################################################
#################################################################################
#################################################################################




#### 3. Calculating network properties and basic analysis with them ####
full_network_props<-netprop(full_network_B1_B2)

##### 3.1. Subset graph to island subnetworks #####

sp_list<- lapply(unique(BALA_NFF$Island), function(x) unique(BALA_NFF[BALA_NFF$Island==x, "Dev_MF"]))
names(sp_list)<-unique(BALA_NFF$Island)

island_subnetworks<-lapply(sp_list, function(x) 
    induced_subgraph(full_network_B1_B2, which(V(full_network_B1_B2)$name %in% x)))



##### 3.2. Subset graph to site subnetworks #####
sp_site_list<- lapply(unique(BALA_NFF$Site_code), function(x) unique(BALA_NFF[BALA_NFF$Site_code==x, "Dev_MF"]))
names(sp_site_list)<-unique(BALA_NFF$Site_code)

site_subnetworks<-lapply(sp_site_list, function(x) 
    induced_subgraph(full_network_B1_B2, which(V(full_network_B1_B2)$name %in% x)))
length(site_subnetworks)

# calculating network properties for each subnetwork
site_subnetworks_ntw_props<-lapply(site_subnetworks, netprop)
subnetw_prop_list<-lapply(names(site_subnetworks_ntw_props), 
                          function(x) data.frame(site = x, island = substr(x, 1, 3),
                                                 t(site_subnetworks_ntw_props[[x]][["netw_prop"]])))
site_subnetworks_props_df<-do.call("rbind", subnetw_prop_list)


site_subnetworks_modules<-lapply(site_subnetworks, 
                                 function(x)edge.betweenness.community(x) )                                

site_subnetworks_props_df$modularity<- 
    sapply(site_subnetworks_modules, function(x) modularity(x))

# calculate number of introduced per site, and the ratio of them to all species

site_subnetworks_props_df$Site_EXO_num<-sapply(site_subnetworks, 
                                               function(x) length(which(V(x)$name %in% EXO)))
site_subnetworks_props_df$Site_EXO_prop<-sapply(site_subnetworks, 
                                                function(x) length(which(V(x)$name %in% EXO))/gorder(x))


##### 3.3. Summary stats for networks #####

# Identify keystone, core and peripheral spp

all_node_props<-full_network_props$node_prop

degree_winners<-all_node_props[rev(order(all_node_props$deg)), ][1:5,]
degree_loosers<-all_node_props[order(all_node_props$deg), ][1:5,]

rel_degree_winners<-all_node_props[rev(order(all_node_props$rel_deg)),  ][1:5,]
rel_degree_losers<-all_node_props[order(all_node_props$rel_deg),  ][1:5,]

eig_winners<-all_node_props[rev(order(all_node_props$eig)), ][1:5,]
eig_losers<-all_node_props[order(all_node_props$eig), ][1:5,]

closeness_winners<-all_node_props[rev(order(all_node_props$closeness)), ][1:5,]
closeness_losers<-all_node_props[order(all_node_props$closeness), ][1:5,]


sp_taxonomy[sp_taxonomy$Dev_MF %in% rownames(degree_winners),]
sp_taxonomy[sp_taxonomy$Dev_MF %in% rownames(degree_loosers),]

sp_taxonomy[sp_taxonomy$Dev_MF %in% rownames(rel_degree_winners),]
sp_taxonomy[sp_taxonomy$Dev_MF %in% rownames(rel_degree_losers),]

sp_taxonomy[sp_taxonomy$Dev_MF %in% rownames(eig_winners),]
sp_taxonomy[sp_taxonomy$Dev_MF %in% rownames(eig_losers),]

sp_taxonomy[sp_taxonomy$Dev_MF %in% rownames(closeness_winners),]
sp_taxonomy[sp_taxonomy$Dev_MF %in% rownames(closeness_losers),]

# most negative links
vulnerability_losers<-all_node_props[rev(order(all_node_props$vulnerability)), ][1:5,]
sp_taxonomy[sp_taxonomy$Dev_MF %in% rownames(vulnerability_losers),]



#################################################################################
#################################################################################
#################################################################################



#### 4. Testing if there is a difference in the frequency of links between nativity (NEI) categories ####

##### 4.1. Data preparation #####
full_edgelist<-data.frame(get.edgelist(full_network_B1_B2), sign = edge_attr(full_network_B1_B2, "sign"))
colnames(full_edgelist)[1:2]<-c("from", "to")
aa<-merge(full_edgelist, sp_taxonomy[, c("Dev_MF", "N_E_I")], by.x = "from", 
          by.y = "Dev_MF", all.x = T, all.y = F, all=F)
colnames(aa)[4]<-"from_NEI"
full_edgelist<-merge(aa, sp_taxonomy[, c("Dev_MF", "N_E_I")], by.x = "to", 
                     by.y = "Dev_MF", all.x = T, all.y = F, all=F)
colnames(full_edgelist)[5]<-"to_NEI"

# All combinations of NEI categories
combs<-list(c("E", "E"), c("E", "I"), c("E", "N"), c("I", "I"), c("I", "N"), c("N", "N"))

##### 4.1 Calculating measured values and simulation process #####
# Calculating the proportions of negative links for each NEI combinations
NEI_percents<-sapply(combs, function (x) 
{aa<-full_edgelist[(full_edgelist$from_NEI == x[1] & full_edgelist$to_NEI == x[2]) |
                       (full_edgelist$from_NEI == x[2] & full_edgelist$to_NEI == x[1]), ]
(nrow(aa[aa$sign<0,])/nrow(aa))*100
})
names(NEI_percents)<-sapply(combs, function(x) paste(x, collapse = " - "))


# Calculating the proportions for each NEI combinations
NEI_prefs<-sapply(combs, function (x) 
{aa<-full_edgelist[(full_edgelist$from_NEI == x[1] & full_edgelist$to_NEI == x[2]) |
                       (full_edgelist$from_NEI == x[2] & full_edgelist$to_NEI == x[1]), ]
nrow(aa)/nrow(full_edgelist)*100
})
names(NEI_prefs)<-sapply(combs, function(x) paste(x, collapse = " - "))


# Simulation for preferential attachments of nodes
NEI_pref_simul<-data.frame(comb = names(NEI_prefs))
for (n in 1:1000)
{
    new_edgelist<-data.frame(full_edgelist[,1:3], 
                             from_NEI = sample(c("N", "E", "I"), nrow(full_edgelist), replace = T),
                             to_NEI = sample(c("N", "E", "I"), nrow(full_edgelist), replace = T))
    
    df<-data.frame(sapply(combs, function (x) 
    {aa<-new_edgelist[(new_edgelist$from_NEI == x[1] & new_edgelist$to_NEI == x[2]) |
                          (new_edgelist$from_NEI == x[2] & new_edgelist$to_NEI == x[1]), ]
    (nrow(aa)/nrow(full_edgelist))*100
    }))
    
    colnames(df)<-paste0("S", n)
    NEI_pref_simul<-cbind(NEI_pref_simul, df)
}
par(ask=T)
lapply(1:length(NEI_prefs), function(x) {hist(as.numeric(NEI_pref_simul[x, 2:ncol(NEI_pref_simul)]))
    abline(v = NEI_prefs[x], col="red")
    mtext(NEI_prefs[x], side=3)})


# test if proportions are different

sapply(1:length(NEI_prefs), function(x) t.test(as.numeric(NEI_pref_simul[x, 2:ncol(NEI_pref_simul)]),
                                               mu = NEI_prefs[x], alternative = "less")$p.value)

sapply(1:length(NEI_prefs), function(x) t.test(as.numeric(NEI_pref_simul[x, 2:ncol(NEI_pref_simul)]),
                                               mu = NEI_prefs[x], alternative = "greater")$p.value)

sapply(1:length(NEI_prefs), function(x) t.test(as.numeric(NEI_pref_simul[x, 2:ncol(NEI_pref_simul)]),
                                               mu = NEI_prefs[x])$p.value)


# Simulation for negative link numbers between NEI
# swap NEI values randomly 10000 times and compare real values

NEI_percent_simul<-data.frame(comb = names(NEI_percents))
for (n in 1:1000)
{
    new_edgelist<-data.frame(full_edgelist[,1:3], 
                             from_NEI = sample(c("N", "E", "I"), nrow(full_edgelist), replace = T),
                             to_NEI = sample(c("N", "E", "I"), nrow(full_edgelist), replace = T))
    
    df<-data.frame(sapply(combs, function (x) 
    {aa<-new_edgelist[(new_edgelist$from_NEI == x[1] & new_edgelist$to_NEI == x[2]) |
                          (new_edgelist$from_NEI == x[2] & new_edgelist$to_NEI == x[1]), ]
    (nrow(aa[aa$sign<0,])/nrow(aa))*100
    }))
    
    colnames(df)<-paste0("S", n)
    NEI_percent_simul<-cbind(NEI_percent_simul, df)
}

# test if proportions are different
sapply(1:length(NEI_percents), function(x) t.test(as.numeric(NEI_percent_simul[x, 2:ncol(NEI_percent_simul)]),
                                                  mu = NEI_percents[x], alternative = "less")$p.value)

sapply(1:length(NEI_percents), function(x) t.test(as.numeric(NEI_percent_simul[x, 2:ncol(NEI_percent_simul)]),
                                                  mu = NEI_percents[x], alternative = "greater")$p.value)

sapply(1:length(NEI_percents), function(x) t.test(as.numeric(NEI_percent_simul[x, 2:ncol(NEI_percent_simul)]),
                                                  mu = NEI_percents[x])$p.value)

##### 4.2. Statistical tests #####


## RECALCULATE NUMBERS!
# for I-N and N-N
# I - N + N - I = 17+19 = 36
prop.test(x = c(58, 45), n = c(nrow(full_edgelist), nrow(full_edgelist)),
          alternative = "less")

# for I-E and E-E
# I - E + E -I = 20+24 = 44
prop.test(x = c(48, 67), n = c(nrow(full_edgelist), nrow(full_edgelist)),
          alternative = "less")

# for I-E + I-N and N+N + N-E + E-E
# I-E + E-I + I-N + N-I= 20+24 + 17+19 = 80
# N-N + E-N + N-E + E-E= 66+98 + 62 + 88 = 314
prop.test(x = c(106, 209), n = c(nrow(full_edgelist), nrow(full_edgelist)),
          alternative = "less")



#################################################################################
#################################################################################
#################################################################################


#### 5. Preparing dataset for camparing NEI networks #####


NEI_subnetworks<-list(
  NAT = induced_subgraph(full_network_B1_B2, which(V(full_network_B1_B2)$name %in% NAT)),
  EXO = induced_subgraph(full_network_B1_B2, which(V(full_network_B1_B2)$name %in% EXO))
)  


#### Subset graph to native/exotic per island and calculating network properties for them ####

# making subnetworks
NEI_island_subnetworks<-lapply(island_subnetworks, function(x)
{subntwk<-list(
  NAT = induced_subgraph(x, which(V(x)$name %in% NAT)),
  EXO = induced_subgraph(x, which(V(x)$name %in% EXO))
)})  
NEI_island_subnetworks<-do.call(list, unlist(NEI_island_subnetworks, recursive=FALSE))

# calculating network properties
NEI_island_subnetworks_ntw_props<-lapply(NEI_island_subnetworks, netprop)
subnetw_prop_list<-lapply(names(NEI_island_subnetworks_ntw_props), 
                          function(x) data.frame(NEI = substr(x, 5, 7), island = substr(x, 1, 3),
                                                 t(NEI_island_subnetworks_ntw_props[[x]][["netw_prop"]])))
NEI_island_subnetworks_ntw_props_df<-do.call("rbind", subnetw_prop_list)

# calcualting modularity
NEI_island_subnetworks_modules<-lapply(NEI_island_subnetworks, 
                                       function(x)edge.betweenness.community(x) )                                


NEI_island_subnetworks_ntw_props_df$modularity<- 
  sapply(NEI_island_subnetworks_modules, function(x) modularity(x))


#### Subset graph to native/exotic per site and calculating network properties for them ####

# making subnetworks
NEI_site_subnetworks<-lapply(site_subnetworks, function(x)
{subntwk<-list(
  NAT = induced_subgraph(x, which(V(x)$name %in% NAT)),
  EXO = induced_subgraph(x, which(V(x)$name %in% EXO))
)})  
NEI_site_subnetworks<-do.call(list, unlist(NEI_site_subnetworks, recursive=FALSE))

# calculating network properties
NEI_site_subnetworks_ntw_props<-lapply(NEI_site_subnetworks, netprop)
subnetw_prop_list<-lapply(names(NEI_site_subnetworks_ntw_props), 
                          function(x) data.frame(NEI = substr(x, nchar(x)-2, nchar(x)),
                                                 site = substr(x, 4, nchar(x)-4),
                                                 t(NEI_site_subnetworks_ntw_props[[x]][["netw_prop"]])))
NEI_site_subnetworks_ntw_props_df<-do.call("rbind", subnetw_prop_list)

# calcualting modularity
NEI_site_subnetworks_modules<-lapply(NEI_site_subnetworks, 
                                     function(x)edge.betweenness.community(x) )                                


NEI_site_subnetworks_ntw_props_df$modularity<- 
  sapply(NEI_site_subnetworks_modules, function(x) modularity(x))

