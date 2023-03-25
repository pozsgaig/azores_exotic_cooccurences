unique(paste(BALA_NFF$Site_code, BALA_NFF$Project))
length(which(!is.na(match(unique(BALA_NFF[BALA_NFF$Project=="BALA", "Site_code"]), unique(BALA_NFF[BALA_NFF$Project=="BALA2", "Site_code"])))))


# transect spp matrices data for paper 
# all BALA
species_mat<-with(BALA_NFF, tapply(Total, list(Dev_MF, Site_code), sum))
species_mat<-NA20(species_mat)
ncol(species_mat)
nrow(species_mat)

# BALA 1
species_mat<-with(BALA1_NFF, tapply(Total, list(Dev_MF, Site_code), sum))
species_mat<-NA20(species_mat)
species_mat<-species_mat[rowSums(species_mat)>10,]
ncol(species_mat)
nrow(species_mat)

binary_sp_mat<-species_mat
binary_sp_mat[binary_sp_mat>0]<-1
binary_sp_mat<-binary_sp_mat[rowSums(binary_sp_mat)>3,colSums(binary_sp_mat)>3]
ncol(binary_sp_mat)
nrow(binary_sp_mat)


# BALA2
species_mat<-with(BALA2_NFF, tapply(Total, list(Dev_MF, Site_code), sum))
species_mat<-NA20(species_mat)
species_mat<-species_mat[rowSums(species_mat)>10,]
ncol(species_mat)
nrow(species_mat)

binary_sp_mat<-species_mat
binary_sp_mat[binary_sp_mat>0]<-1
binary_sp_mat<-binary_sp_mat[rowSums(binary_sp_mat)>3,colSums(binary_sp_mat)>3]
ncol(binary_sp_mat)
nrow(binary_sp_mat)

nrow(all_node_props[which(V(full_network_B1_B2)$name %in% EXO), ])
nrow(all_node_props[which(V(full_network_B1_B2)$name %in% END), ])
nrow(all_node_props[which(V(full_network_B1_B2)$name %in% NAT), ])

full_network_B1_B2

length(E(full_network_B1_B2)[E(full_network_B1_B2)$sign==-1])
length(E(full_network_B1_B2)[E(full_network_B1_B2)$sign==1])

# exotic spp abundances
exo_spp<-BALA_NFF[BALA_NFF$Dev_MF %in% V(full_network_B1_B2)$name & BALA_NFF$N_E_I=="I",]
exo_spp[order(exo_spp$Total),]
# SM 2
sm2<-full_network_props$node_prop

colnames(sm2)
nrow(sm2)
sm2<-sm2[, c("deg", "rel_deg", "vulnerability", "rel_vulnerability", 
             "norm_closeness", "norm_betweenness")]
colnames(sp_taxonomy)
sm2<-merge(sm2, sp_taxonomy[, c("Dev_MF","Species_name", "Order", "Family","N_E_I")], by.x= "row.names", by.y = "Dev_MF")

write.table(sm2, file = "sm2a.csv", sep=",", row.names = F, col.names = T)


gorder(full_network_B1_B2)
gsize(full_network_B1_B2)
edge_density(full_network_B1_B2)
sort(igraph::degree(full_network_B1_B2))
length(E(full_network_B1_B2)[E(full_network_B1_B2)$sign==-1])
length(E(full_network_B1_B2)[E(full_network_B1_B2)$sign==1])
full_network_props

names(site_subnetworks)
V(full_network_B1_B2)$NEI<-makeVertexAtt(full_network_B1_B2, sp_taxonomy, vname= "N_E_I", by.g = "name", by.df = "Dev_MF")
V(full_network_B1_B2)$NEI
table(V(full_network_B1_B2)$NEI)



# Figure 2
# setting colour
V(full_network_B1_B2)$NEI_color<-c("seagreen4", "brown4", "seagreen2", "light grey")[as.numeric(as.factor(V(full_network_B1_B2)$NEI))]

data.frame(nei = V(full_network_B1_B2)$NEI, col = V(full_network_B1_B2)$NEI_color)
# highest number of negative and positive links should be checked

length(unique(sp_taxonomy$N_E_I))

# setting colors for each order 
V(full_network_B1_B2)$Order<-makeVertexAtt(full_network_B1_B2, df=sp_taxonomy, vname='Order', 
                                           by.df='Dev_MF', by.g='name')

length(unique(V(full_network_B1_B2)$Order))

V(full_network_B1_B2)$Order_color<-cols25(18)[as.numeric(as.factor(V(full_network_B1_B2)$Order))]


# setting colors for each class 
V(full_network_B1_B2)$Class<-makeVertexAtt(full_network_B1_B2, df=sp_taxonomy, vname='Class', 
                                           by.df='Dev_MF', by.g='name')

V(full_network_B1_B2)$Class_color<-
  wes_palette("GrandBudapest2",4)[as.numeric(as.factor(V(full_network_B1_B2)$Class))]


sp_list<- lapply(unique(BALA_NFF$Island), function(x) unique(BALA_NFF[BALA_NFF$Island==x, "Dev_MF"]))
names(sp_list)<-unique(BALA_NFF$Island)

island_subnetworks<-lapply(sp_list, function(x) 
  induced_subgraph(full_network_B1_B2, which(V(full_network_B1_B2)$name %in% x)))

# one large network
net_Isolated<-V(full_network_B1_B2)[igraph::degree(full_network_B1_B2)==0]$name
connected_net<-delete_vertices(full_network_B1_B2, net_Isolated)


mat<-as.matrix(as_adjacency_matrix(connected_net, names = T))

# setting the color matrix for the links
col_mat<-as.matrix(as_adjacency_matrix(connected_net, names = T, attr = "sign"))
col_mat[col_mat==1]<-adjustcolor("red", alpha.f = 0.2)
col_mat[col_mat==-1]<-adjustcolor("blue", alpha.f = 0.2)
col_mat[col_mat==0]<-NA 

# setting up groups and node colours by order and NEI, respectively
v_groups<-structure(V(connected_net)$Order, names=V(connected_net)$name)
v_col<-structure(V(connected_net)$NEI_color, names=V(connected_net)$name)
v_order_df<-data.frame(Cl = V(connected_net)$Class, 
                       Ord = V(connected_net)$Order, 
                       NEI = V(connected_net)$NEI,
                       SP = V(connected_net)$name)

v_order<-v_order_df[order(v_order_df$Cl, v_order_df$Ord, v_order_df$NEI), "SP"]

pdf(file = "fig2.pdf", width = 12, height = 10)
# png has to be saved manually from the pdf

graphics::layout(mat = matrix(c(2, 1), 
                              nrow = 2, 
                              ncol = 1),
                 heights = c(1, 2))


par(mar = c(0,0,0,15), oma=c(0,0,0,0), xpd=T)

circos.clear()
chordDiagram(mat, grid.col= v_col, col=col_mat, group = v_groups, symmetric = T,
             annotationTrack = "grid",
             annotationTrackHeight = 0.04,
             preAllocateTracks = list(
               track.height = mm_h(6),
               track.margin = c(mm_h(1), 0))
)

lapply(unique(V(connected_net)$Order), function(x){
  sector_2_highlight<-V(connected_net)[V(connected_net)$Order==x]$name
  sector_col<-V(connected_net)[V(connected_net)$Order==x]$Order_color[1]
  highlight.sector(sector_2_highlight, track.index = 1, col = sector_col,niceFacing = TRUE)
})


legend("right", inset=c(-0.25,0), pch = 15, col = cols25(18), 
       legend = sort(unique(V(full_network_B1_B2)$Order)), cex = 1.5)

isolated_props<-table(sp_taxonomy[sp_taxonomy$Dev_MF %in% net_Isolated, "N_E_I"])
isolated_props<-isolated_props[c(1,3,2)]
names(isolated_props)<-c("Endemic", "Native", "Exotic")

plot.new()
vps <- baseViewports()
pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
vp1 <-plotViewport(c(0,0,0,1))
p<-waffle(isolated_props, colors = c("seagreen4", "seagreen2", "brown4"),
          rows = 3, title = "Unconnected species")


print(p, vp = vp1)

dev.off()


# Network figure size of the circle can be the relative abundance of the species in the community 
# and the gradient colour can be the closeness/betweenness
colnames(site_subnetworks_props_df)


island_summaries<-site_subnetworks_props_df %>% 
  select(c(island, num_nodes, edgenum, netw_edge_dens, modularity)) %>% 
  group_by(island) %>% 
  summarise_all(c(mean = mean, std.error = std.error)) %>% 
  dplyr::arrange(island)

island_cols<-structure(brewer.pal(7, "RdYlBu"), names=sort(names(island_subnetworks)))

island_cols[c("SJG", "SMG")]<-c("#fafa57", "#b0ecfc") # yellow and blue were too bright


island_measured<-sapply(island_subnetworks, function(x) {
  negrat<-ifelse(length(E(x)[E(x)$sign==1])==0, 0,
                 length(E(x)[E(x)$sign==1])/gsize(x))
  norm_closen<-centralize(harmonic_centrality(x, normalized = T, mode = "all"), 
                          centr_clo_tmax(x), normalized = T)
  EXO_prop<-length(which(V(x)$name %in% EXO))/gorder(x)
  
  aa<-c(gorder(x), edge_density(x), negrat, norm_closen)
  bb<-edge.betweenness.community(x)
  aa<-c(aa, modularity(bb), EXO_prop)})
# rownames(island_measured)<-c("num_nodes", "netw_edge_dens", "neg_edge_ratio", 
#                              "normalized_closeness", "modularity", "EXO_prop") 
rownames(island_measured)<-c("Number of nodes", "Connectance", "Negative edge\nratio", 
                             "Normalised closeness centrality", "Modularity", 
                             "\nProportion of\nexotic species") 

island_measured<-t(island_measured)

lapply(rownames(island_measured), radarchart_func)

dev.off()

island_measured<-island_measured[order(rownames(island_measured)),]
island_topologies<-cbind(island_summaries, island_measured)

write.xlsx2(island_topologies, file = "summaries_4_paper1.xlsx", sheetName = "Island")




#######################################################################
#### making the plots and pairwise comparison plots for both dbRDAs ####
# Fig. 3.
colnames(site_subnetworks_props_df)
ntw_variables_of_interest<-colnames(site_subnetworks_props_df)[c(3, 7, 5, 11, 20, 43, 23, 26, 31, 34, 35, 44, 48)]

env_variables_of_interest<-colnames(island_NFF_area)[c(2, 3, 35, 36, 37, 4, 10, 28, 20, 31, 15)]


alldat<-merge(site_subnetworks_props_df, island_NFF_area, by.x = "island", by.y = "Island")
nrow(site_subnetworks_props_df)
nrow(alldat)
colnames(alldat)

comm_mat<-alldat[, ntw_variables_of_interest]
comm_mat_scaled<-apply(comm_mat, 2,  scale) # z-score, can be used with Euclidean dist
comm_mat_scaled<-as.data.frame(comm_mat_scaled)
rownames(comm_mat_scaled)<-alldat$site

env_mat<-alldat[, env_variables_of_interest]
env_mat_scaled<-apply(env_mat, 2,  scale)
env_mat_scaled<-as.data.frame(env_mat_scaled)
env_mat_scaled$Island<-alldat$island
rownames(env_mat_scaled)<-alldat$site


comm_mat_scaled<-na.omit(comm_mat_scaled)
env_mat_scaled<-env_mat_scaled[rownames(env_mat_scaled) %in% rownames(comm_mat_scaled),]

# dbrda with Eucledian distances is the same as RDA
mod0 <- rda(comm_mat_scaled ~ 1, env_mat_scaled[, -ncol(env_mat_scaled)])  # Model with intercept only
mod1 <- rda(comm_mat_scaled ~ ., env_mat_scaled[, -ncol(env_mat_scaled)])  # Model with all explanatory variables
summary(mod1)

set.seed(666)

step.res <- ordistep(mod1, perm.max = 100000, direction = "both")


##### manual edit
# repeating the same as in paper
step.res<-rda(formula = comm_mat_scaled ~ NFF_prop + EXO_num + EXO_prop + Altitude_Average + Temp.summer_Median
              + Rh_Average, data = env_mat_scaled[, -ncol(env_mat_scaled)])

# manual edit ends  


step.res$anova  # Summary table
summary(step.res)

anova(step.res)


set.seed(123)
ano_res<-anosim(comm_mat_scaled, env_mat_scaled$Island, distance = "euc",
                permutations=9999)

p_anosim<-pairwise.anosim(x = comm_mat_scaled,  factors = env_mat_scaled$Island,  sim.method = 'euc', p.adjust.m ='fdr')


# see what properties cause the differences
sigpairs<-gsub(" vs ", "_", p_anosim[p_anosim$p.adjusted<0.05, "pairs"])

(sim <- with(env_mat_scaled, simper(comm_mat_scaled, Island)))
sumsim<-summary(sim)

sapply(sumsim[sigpairs], function(x) rownames(x[1:4,]))


ano_text<-paste("ANOSIM R =", round(ano_res$statistic, 3), 
                ifelse(round(ano_res$signif, 3)<=0.001, "P < 0.001", paste("P =", round(ano_res$signif, 3))))

P_ord_1<-ggord(step.res, env_mat_scaled$Island, cols = island_cols, polylntyp = "solid", poly = T, 
               alpha_el=0.4, veclsz = 1, veccol = "red", arrow=0.3, size=1.5, xlims = c(-1.2,1.2), ylims = c(-1.5, 1),
               grp_title = "Islands", coord_fix = F)+
  annotate(geom="text", Inf,-Inf, label=ano_text, hjust = 1.2, vjust = -0.3)


p_anosim$Pair1<-sapply(strsplit(p_anosim$pairs, " vs "), "[", 1)
p_anosim$Pair2<-sapply(strsplit(p_anosim$pairs, " vs "), "[", 2)
p_anosim$R_stat<-ifelse(p_anosim$p.adjusted<0.05, p_anosim$R_stat, NA)

P_pairs_1<-ggplot(p_anosim, ggplot2::aes(y = Pair1, 
                                         x = Pair2, fill = R_stat)) + 
  geom_tile(color = "white", width = 0.97, height = 0.97) + 
  # geom_text(ggplot2::aes(label = "Islands"), na.rm = TRUE) + # use this if you want to show text in the tiles
  scale_fill_gradient(low = "mistyrose2", high = "brown1", na.value = "grey90",
                      limits = c(min(p_anosim$R_stat), max(p_anosim$R_stat)), name = NULL) +
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "right") + 
  labs(x = "", y = "", title = "Pairwise island comparisons") + 
  coord_cartesian(expand = FALSE)

####* dbRDA to show island differences based on species pairs and relationship with environmental variables *####


comm_mat<-t(sapply(site_subnetworks, function(x){
  full_EL<-paste(get.edgelist(full_network_B1_B2)[,1], get.edgelist(full_network_B1_B2)[,2], sep="--")
  sub_EL<-paste(get.edgelist(x)[,1], get.edgelist(x)[,2], sep="--")
  as.integer(full_EL %in% sub_EL)
}))

colnames(comm_mat)<-paste(get.edgelist(full_network_B1_B2)[,1], get.edgelist(full_network_B1_B2)[,2], sep="--")
rownames(comm_mat)

comm_mat<-merge(comm_mat, site_subnetworks_props_df[,1:2], by.x="row.names", by.y="site")
nrow(comm_mat)
env_variables_of_interest<-colnames(island_NFF_area)[c(2, 3, 35, 36, 37, 4, 10, 28, 20, 31, 15)]


alldat<-merge(comm_mat, island_NFF_area, by.x = "island", by.y = "Island")
nrow(alldat)
colnames(alldat)
rownames(alldat)<-alldat$Row.names
comm_mat<-alldat[, 3:332]
rownames(comm_mat)

comm_mat<-comm_mat[rowSums(comm_mat)>0,]

env_mat<-alldat[, env_variables_of_interest]
env_mat_scaled<-apply(env_mat, 2,  scale)
env_mat_scaled<-as.data.frame(env_mat_scaled)
env_mat_scaled$Island<-alldat$island
rownames(env_mat_scaled)<-rownames(alldat)

comm_mat<-na.omit(comm_mat)
env_mat_scaled<-env_mat_scaled[rownames(env_mat_scaled) %in% rownames(comm_mat),]

# dbrda with Eucledian distances is the same as RDA

mod0 <- dbrda(comm_mat ~ 1, env_mat_scaled[, -ncol(env_mat_scaled)], distance = "jaccard")  # Model with intercept only
mod1 <- dbrda(comm_mat ~ ., env_mat_scaled[, -ncol(env_mat_scaled)], distance = "jaccard")  # Model with all explanatory variables
summary(mod1)


set.seed(123)

step.res <- ordistep(mod1, perm.max = 100000, direction = "both")
step.res$anova  # Summary table
summary(step.res)
anova.cca(step.res)

set.seed(123)
ano_res<-anosim(comm_mat, env_mat_scaled$Island, distance = "jaccard",
                permutations=9999)

p_anosim<-pairwise.anosim(x = comm_mat,  factors = env_mat_scaled$Island,  sim.method = 'jaccard', p.adjust.m ='fdr')

(sim <- with(env_mat_scaled, simper(comm_mat, Island)))
summary(sim)

ano_text<-paste("ANOSIM R =", round(ano_res$statistic, 3), 
                ifelse(round(ano_res$signif, 3)<=0.001, "P < 0.001", paste("P =", round(ano_res$signif, 3))))

P_ord_2<-ggord(step.res, env_mat_scaled$Island, cols = island_cols, polylntyp = "solid", poly = T, 
               alpha_el=0.4, veclsz = 1, veccol = "red", arrow=0.3, size=1.5, xlims = c(-1.2,1.5), ylims = c(-1.0, 0.7),
               grp_title = "Islands", coord_fix = F)+
  annotate(geom="text", Inf,-Inf, label=ano_text, hjust = 1.2, vjust = -0.3)


p_anosim$Pair1<-sapply(strsplit(p_anosim$pairs, " vs "), "[", 1)
p_anosim$Pair2<-sapply(strsplit(p_anosim$pairs, " vs "), "[", 2)
p_anosim$R_stat<-ifelse(p_anosim$p.adjusted<0.05, p_anosim$R_stat, NA)

P_pairs_2<-ggplot(p_anosim, ggplot2::aes(y = Pair1, 
                                         x = Pair2, fill = R_stat)) + 
  geom_tile(color = "white", width = 0.97, height = 0.97) + 
  # geom_text(ggplot2::aes(label = "Islands"), na.rm = TRUE) + # use this if you want to show text in the tiles
  scale_fill_gradient(low = "mistyrose2", high = "brown1", na.value = "grey90",
                      limits = c(min(p_anosim$R_stat), max(p_anosim$R_stat)), name = NULL) +
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "right") + 
  labs(x = "", y = "", title = "Pairwise island comparisons") + 
  coord_cartesian(expand = FALSE)


ggsave(plot = ggarrange(plotlist =  list(P_ord_1, P_pairs_1, P_ord_2, P_pairs_2), widths = c(2, 1), 
                        heights = c(1, 1), labels = LETTERS[1:4]),
       filename="Fig_3.pdf",  
       device = "pdf", height = 10, width = 15,
       units = "in")


########

# all links

site_NEI_all_link_props<-lapply(site_subnetworks, function(m) {
  full_edgelist<-data.frame(get.edgelist(m), sign = edge_attr(m, "sign"))
  colnames(full_edgelist)[1:2]<-c("from", "to")
  aa<-merge(full_edgelist, sp_taxonomy[, c("Dev_MF", "N_E_I")], by.x = "from", 
            by.y = "Dev_MF", all.x = T, all.y = F, all=F)
  colnames(aa)[4]<-"from_NEI"
  full_edgelist<-merge(aa, sp_taxonomy[, c("Dev_MF", "N_E_I")], by.x = "to", 
                       by.y = "Dev_MF", all.x = T, all.y = F, all=F)
  colnames(full_edgelist)[5]<-"to_NEI"
  NEI_prefs<-sapply(combs, function (x) 
  {aa<-full_edgelist[(full_edgelist$from_NEI == x[1] & full_edgelist$to_NEI == x[2]) |
                       (full_edgelist$from_NEI == x[2] & full_edgelist$to_NEI == x[1]), ]
  (nrow(aa)/nrow(full_edgelist))*100
  })
  names(NEI_prefs)<-sapply(combs, function(x) paste(x, collapse = " - "))
  NEI_prefs
})
dat<-melt(do.call("rbind", site_NEI_all_link_props))
dat<-na.omit(dat)


pA <- ggboxplot(dat, x = "Var2", y = "value",
                color = "Var2", palette = "jco",
                add = "jitter", add.params = list(alpha=0.1))+
  stat_compare_means(label.y = 50)+
  ylab(label="Proportion of links\nto all link")+
  theme(legend.position = "none")


wilcox_dat<-tidy(pairwise.wilcox.test(dat$value, dat$Var2, p.adjust.method = "fdr"))


wilcox_dat$Significance<-"n.s"
wilcox_dat[wilcox_dat$p.value<=0.05, "Significance"]<-"p < 0.05"
wilcox_dat[wilcox_dat$p.value<=0.01, "Significance"]<- "p < 0.01"
wilcox_dat[wilcox_dat$p.value<=0.001, "Significance"]<- "p < 0.001"

wilcox_dat$Significance<-as.factor(wilcox_dat$Significance)
levels(wilcox_dat$Significance)<-c(levels(wilcox_dat$Significance), "p < 0.05")

pB<-ggplot(wilcox_dat, ggplot2::aes(y = group1 , 
                                    x = group2, fill = Significance)) + 
  geom_tile(color = "white", width = 0.97, height = 0.97) + 
  scale_fill_manual(values = c("white", "red", "orange", "yellow"))+
  # geom_text(ggplot2::aes(label = "Islands"), na.rm = TRUE) + # use this if you want to show text in the tiles
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "right") + 
  labs(x = "", y = "", title = "Pairwise comparisons of linking preferences")+
  coord_cartesian(expand = FALSE)


# Calculating the proportions of negative links for each NEI combinations 
site_NEI_neg_link_props<-lapply(site_subnetworks, function(m) {
  full_edgelist<-data.frame(get.edgelist(m), sign = edge_attr(m, "sign"))
  colnames(full_edgelist)[1:2]<-c("from", "to")
  aa<-merge(full_edgelist, sp_taxonomy[, c("Dev_MF", "N_E_I")], by.x = "from", 
            by.y = "Dev_MF", all.x = T, all.y = F, all=F)
  colnames(aa)[4]<-"from_NEI"
  full_edgelist<-merge(aa, sp_taxonomy[, c("Dev_MF", "N_E_I")], by.x = "to", 
                       by.y = "Dev_MF", all.x = T, all.y = F, all=F)
  colnames(full_edgelist)[5]<-"to_NEI"
  NEI_percents<-sapply(combs, function (x) 
  {bb<-full_edgelist[(full_edgelist$from_NEI == x[1] & full_edgelist$to_NEI == x[2]) |
                       (full_edgelist$from_NEI == x[2] & full_edgelist$to_NEI == x[1]), ]
  (nrow(bb[bb$sign<0,])/nrow(bb))*100
  })
  names(NEI_percents)<-sapply(combs, function(x) paste(x, collapse = " - "))
  NEI_percents
})
dat<-melt(do.call("rbind", site_NEI_neg_link_props))


pC <- ggboxplot(dat, x = "Var2", y = "value",
                color = "Var2", palette = "jco",
                add = "jitter", add.params = list(alpha=0.1))+
  #stat_compare_means(comparisons = my_comparisons, hide.ns = T)+
  stat_compare_means(label.y = 50)+
  ylab(label="Proportion of negative links\nto all links")+
  theme(legend.position = "none")


wilcox_dat<-tidy(pairwise.wilcox.test(dat$value, dat$Var2, p.adjust.method = "fdr"))

wilcox_dat$Significance<-"n.s"
wilcox_dat[wilcox_dat$p.value<=0.05, "Significance"]<-"p < 0.05"
wilcox_dat[wilcox_dat$p.value<=0.01, "Significance"]<- "p < 0.01"
wilcox_dat[wilcox_dat$p.value<=0.001, "Significance"]<- "p < 0.001"

wilcox_dat$Significance<-as.factor(wilcox_dat$Significance)
levels(wilcox_dat$Significance)<-c(levels(wilcox_dat$Significance), "p < 0.05")

pD<-ggplot(wilcox_dat, ggplot2::aes(y = group1 , 
                                    x = group2, fill = Significance)) + 
  geom_tile(color = "white", width = 0.97, height = 0.97) + 
  scale_fill_manual(values = c("white", "red", "orange", "yellow"))+
  # geom_text(ggplot2::aes(label = "Islands"), na.rm = TRUE) + # use this if you want to show text in the tiles
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "right") + 
  labs(x = "", y = "", title = "Pairwise comparisons of linking preferences")+
  coord_cartesian(expand = FALSE)


ggsave(plot = ggarrange(plotlist =  list(pA, pB, pC, pD), widths = c(2, 1), 
                        heights = c(1, 1), labels = LETTERS[1:4]),
       filename="fig_4.pdf",  
       device = "pdf", height = 10, width = 15,
       units = "in")


#### Comparing native/exotic networks ####
colnames(NEI_site_subnetworks_ntw_props_df)
NEI_site_subnetworks_ntw_props_df$NEI[NEI_site_subnetworks_ntw_props_df$NEI=="NAT"]<-"IND"


NEI_site_subnetworks_ntw_props_df$normalized_betweenness_cent<-
  log(NEI_site_subnetworks_ntw_props_df$normalized_betweenness_cent+0.001)

NEI_site_subnetworks_ntw_props_df$normalized_closeness_cent<-
  log(NEI_site_subnetworks_ntw_props_df$normalized_closeness_cent+0.001)

NEI_site_subnetworks_ntw_props_df$realised_degree <-
  log(NEI_site_subnetworks_ntw_props_df$realised_degree+0.0001)

NEI_site_subnetworks_ntw_props_df$edgenum<-
  log(NEI_site_subnetworks_ntw_props_df$edgenum+1)


vars_to_plot<-c("num_nodes", "edgenum", "netw_edge_dens", "prop_isolated", 
                "neg_edge_ratio", "modularity")
varnames<-c("Number of nodes", "Number of edges\n(log-transformed)",
            "Connectance\n(log-transformed)",
            "Proportion of isolated nodes", 
            "Proportion of negative edges","Modularity")


kruskal_test(NEI_site_subnetworks_ntw_props_df, modularity~NEI)
pairwise.wilcox.test(NEI_site_subnetworks_ntw_props_df$modularity,
                     NEI_site_subnetworks_ntw_props_df$NEI, p.adjust.method = "fdr")

NEI_network_props<-lapply(1:length(vars_to_plot), 
                          function(x) {print(x)
                            ordered_ggpirate(NEI_site_subnetworks_ntw_props_df, 
                                             x_ax = "NEI", y_ax = vars_to_plot[x], 
                                             colour2 = "NEI",
                                             sidelab = varnames[x],
                                             title = LETTERS[x])+
                              xlab(label = "Nativity")})

ggsave(plot = marrangeGrob(NEI_network_props, nrow=2, ncol=3, top = ""),
       filename="Fig_5.pdf",  
       device = "pdf", height = 12, width = 18,
       units = "in")



#### Comparing island networks ####
colnames(site_subnetworks_props_df)
site_subnetworks_props_df$realised_degreelog<-log(0.001+site_subnetworks_props_df$netw_edge_dens)



ind_network_props<-lapply(1:length(vars_to_plot), 
                          function(x) {print(x)
                            ordered_ggpirate(site_subnetworks_props_df, x_ax = "island", y_ax = vars_to_plot[x], 
                                             colour2 = "island", sidelab = varnames[x], title = LETTERS[x])})


ggsave(plot = marrangeGrob(ind_network_props, nrow=2, ncol=2, 
                           layout_matrix = matrix(1:4, 2, 2, TRUE)),
       filename="SM3_island_boxplots.pdf",  
       device = "pdf", height = 12, width = 12,
       units = "in")



#### Compare node properties between native and exotics ####

all_node_props$NEI<-NA
all_node_props[which(V(full_network_B1_B2)$name %in% NAT), "NEI"]<-"IND"
all_node_props[which(V(full_network_B1_B2)$name %in% EXO), "NEI"]<-"EXO"
# NA NEIs exist
V(full_network_B1_B2)$name[which(!V(full_network_B1_B2)$name %in% EXO & !V(full_network_B1_B2)$name %in% NAT)]
# unknown NEI removed
NEI_dat<-all_node_props[!is.na(all_node_props$NEI),]

# non-isolated nodes are kept!
# NEI_dat<-NEI_dat[NEI_dat$deg>0,]
vars_to_plot<-colnames(NEI_dat)[c(1,2,6,8,9, 11)]
varnames<-c("Degree", "Relative degree", "Normalised closeness centrality", 
            "Normalised betweenness centrality", "Vulnerability", "Relative vulnerability")



NEI_node_props<-lapply(1:length(vars_to_plot), 
                       function(x) {print(x)
                         ordered_ggpirate(NEI_dat, x_ax = "NEI", 
                                          y_ax = vars_to_plot[x], colour2 = "NEI",
                                          sidelab = varnames[x], title = LETTERS[x])+
                           xlab(label="Nativity")})

ggsave(plot = marrangeGrob(NEI_node_props, nrow=2, ncol=3),
       filename="SM_nodes.pdf",  
       device = "pdf", height = 12, width = 16,
       units = "in")



#### Island network propeties' relationship with island environmental variables ####


colnames(site_subnetworks_props_df)

site_subnetworks_props_df$realised_cliques<-
  site_subnetworks_props_df$clique_num/site_subnetworks_props_df$max_num_cliques

ntw_variables_of_interest<-site_subnetworks_props_df[, c("island","num_nodes", 
                                                         "edgenum", 
                                                         "prop_isolated", 
                                                         "neg_edge_ratio", 
                                                         "netw_edge_dens", 
                                                         "normalized_closeness_cent",
                                                         "normalized_betweenness_cent", 
                                                         "modularity")]


rownames(ntw_variables_of_interest)<-site_subnetworks_props_df$site
nrow(ntw_variables_of_interest)

colnames(island_NFF_area)
env_vars<-island_NFF_area[,c(1:4, 6,7,9,10, 14, 15, 19, 20, 24, 25, 27, 28, 30, 32, 33,35,36, 37)]

lm_dat<-merge(ntw_variables_of_interest, env_vars, by.x = "island", by.y = "Island",
              all.x = T, all.y = F)

nrow(lm_dat)
sapply(lm_dat, class)

# removing meaningless columns
keep<-!is.na(apply(lm_dat[, 2:ncol(lm_dat)], 2, sd, na.rm=T)) & 
  apply(lm_dat[, 2:ncol(lm_dat)], 2, sd, na.rm=T)>0
lm_dat<-lm_dat[, c(T, keep)]

# scaling variables
lm_dat[,2:ncol(lm_dat)]<-apply(lm_dat[,2:ncol(lm_dat)], 2, scale)



all_corr<-corr.test(lm_dat[, 2:ncol(lm_dat)], method = "pearson", adjust = "none")
all_corr$r[all_corr$p >0.05]<-0
all_corr_mat<-all_corr$r

varnames<-
  c("Number of nodes", "Number of edges", 
    "Proportion of\nisolated nodes",
    "Negative edge\nratio", "Connectance", 
    "Normalised closeness\ncentrality",
    "Normalised betweenness\ncentrality",
    "Modularity", "Area of native forest",
    "Island area", "Average altitude", 
    "Maximum altitude", "Average slope",
    "Maximum slope", "Mean precipitation",
    "Precipitation range", 
    "Mean relative\nhumidity", "Humidity range",
    "Mean temperature", "Temperature range",
    "Mean summer\nprecipitation", 
    "Maximum summer\nprecipitation", 
    "Minimum summer\nprecipitation",
    "Mean summer\ntemperature", 
    "Maximum summer\ntemperature", 
    "Minimum summer\ntemperature",
    "Area proportion of\nnative forest",
    "Number of\nexotic species",
    "Proportion of\nexotic species"
  )
names(varnames)<-dimnames(all_corr_mat)[[1]]
# since there may be the same negative and positive values in the same column,
# the absolute value of each number has to be taken
absColSums<-function(dat){apply(dat, 2, function(x)sum(abs(x)))}
absRowSums<-function(dat){apply(dat, 1, function(x)sum(abs(x)))}

# removing variables showing no correlation at all
removed_corrs<-all_corr_mat[absRowSums(all_corr_mat)==1, absColSums(all_corr_mat)==1]
all_confirmed_corr<-all_corr_mat[!absRowSums(all_corr_mat)==1, !absColSums(all_corr_mat)==1]

ntw_conf_corrs<-all_confirmed_corr[rownames(all_confirmed_corr) %in% colnames(site_subnetworks_props_df), 
                                   colnames(all_confirmed_corr) %in% colnames(site_subnetworks_props_df)]

env_conf_corrs<-all_confirmed_corr[rownames(all_confirmed_corr) %in% colnames(island_NFF_area), 
                                   colnames(all_confirmed_corr) %in% colnames(island_NFF_area)]

mixed_conf_corrs<-all_confirmed_corr[rownames(all_confirmed_corr) %in% colnames(site_subnetworks_props_df), 
                                     colnames(all_confirmed_corr) %in% colnames(island_NFF_area)]

ntw_ord<-rev(order(absColSums(ntw_conf_corrs)))
env_ord<-rev(order(absColSums(env_conf_corrs)))
mixed_col_ord<-order(absColSums(mixed_conf_corrs))
mixed_row_ord<-order(absRowSums(mixed_conf_corrs))

corr_col<-brewer.pal(n = 10, name = "RdBu")[10:1]
dimnames(ntw_conf_corrs)[[1]]<-varnames[dimnames(ntw_conf_corrs)[[1]]]
dimnames(ntw_conf_corrs)[[2]]<-varnames[dimnames(ntw_conf_corrs)[[2]]]

dimnames(env_conf_corrs)[[1]]<-varnames[dimnames(env_conf_corrs)[[1]]]
dimnames(env_conf_corrs)[[2]]<-varnames[dimnames(env_conf_corrs)[[2]]]

dimnames(mixed_conf_corrs)[[1]]<-varnames[dimnames(mixed_conf_corrs)[[1]]]
dimnames(mixed_conf_corrs)[[2]]<-varnames[dimnames(mixed_conf_corrs)[[2]]]

pdf(file = "SM_Correlations.pdf", width = 12, height = 10)
corrplot(ntw_conf_corrs[,ntw_ord], tl.cex=.8, col = corr_col)
corrplot(env_conf_corrs[,env_ord], tl.cex=.8, col = corr_col)
corrplot(mixed_conf_corrs[,mixed_col_ord], tl.cex=.8, is.corr=FALSE, col = corr_col)
dev.off()






##### Citations ####
citation("lme4")
citation("lmerTest")
citation("MuMIn")
citation("ggplot2")
citation("reshape2")
citation("phyloseq")
