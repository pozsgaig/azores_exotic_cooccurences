#### Test if network properties are in relationship with the number and proportion of exotics. ####
  ## This is done separately for the full network, for native networks only, and for the residuals of ##
  ## the regressions with the native networks. Plots and data frames are stored/saved. ##

site_subnetworks_props_df$Site_NAT_num<-
  site_subnetworks_props_df$num_nodes-site_subnetworks_props_df$Site_EXO_num

# ntw_variables_of_interest<-colnames(site_subnetworks_props_df)[c(3, 7, 5, 11, 20, 21, 23, 26, 31, 34, 35, 44, 48, 51)]

ntw_variables_of_interest<-colnames(site_subnetworks_props_df)[c(3, 7, 5, 11, 20, 43, 31, 34, 48)]
names(ntw_variables_of_interest)<-c("Number of nodes", "Number of edges", "Proportion of isolated nodes", 
            "Proportion of negative edges", "Mean degree", "Connectance",  
            "Normalised closeness centrality", 
            "Normalised betweenness centrality",
            "Modularity")

####* 1. Regression on the full network against EXO variables *####


  #### ** 1.1. Preparing data ** ####

tempdf<-site_subnetworks_props_df

# square root transform proportions
tempdf$Site_EXO_prop<-sqrt(tempdf$Site_EXO_prop)
tempdf$prop_isolated<-sqrt(tempdf$prop_isolated)
tempdf$prop_isolated<-sqrt(tempdf$prop_isolated)
tempdf$realised_degree<-sqrt(tempdf$realised_degree)
tempdf$neg_edge_ratio<-sqrt(tempdf$neg_edge_ratio)

# the only zero value is removed
tempdf[tempdf$Site_EXO_prop==0, ]
tempdf<-tempdf[!tempdf$Site_EXO_prop==0, ]


  #### ** 1.2. Regression against exotic NUMBERS ** ####

      #### *** 1.2.1. Generating and saving plots *** ####

fn_EXO_num_plots<-lapply(names(ntw_variables_of_interest), function(x){
  reg_formula <- as.formula(paste(ntw_variables_of_interest[x], "~", "Site_EXO_num", "+", "(1|island)"))
  mod<-lme4::lmer(reg_formula, data = tempdf)
  print(anova(mod))
  print(car::Anova(mod))
  mod_pval<-ifelse(car::Anova(mod)[3]<0.001, "p<0.001", 
                   round(car::Anova(mod)[3], 3))
  text_col<-unlist(ifelse(car::Anova(mod)[3]<0.05, "red", "black"))
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  
  effects_EXO_num <- as.data.frame(effects::effect(term= "Site_EXO_num", mod= mod))
  dat<-tempdf[, c("site", "island", "Site_EXO_num", ntw_variables_of_interest[x])]
  
  print(paste("P =", mod_pval, ", R2 =", mod_R2))
  ggplot() + 
    geom_point(data=dat, aes_string("Site_EXO_num", ntw_variables_of_interest[x]), color="light grey") + 
    geom_point(data=effects_EXO_num, aes(x=Site_EXO_num, y=fit)) +
    geom_line(data=effects_EXO_num, aes(x=Site_EXO_num, y=fit), color=text_col[1]) +
    geom_ribbon(data= effects_EXO_num, aes(x=Site_EXO_num, ymin=lower, ymax=upper), alpha= 0.3, fill=text_col[1])+
    ggplot2::annotate("text", x = max(dat[,"Site_EXO_num"])*0.9, 
             y = max(dat[,ntw_variables_of_interest[x]])*0.9, 
             label = paste("P =", mod_pval, ", R2 =", mod_R2), color=text_col[1])+
    ylab(label=x)+
    xlab(label="Number of exotic species")+
    theme_classic()
  
  
}
)
ggsave(plot = marrangeGrob(fn_EXO_num_plots, nrow=3, ncol=3),
       filename="SM5.pdf", 
       device = "pdf", height = 21, width = 21,
       units = "in")

      #### *** 1.2.2. Storing df of R2 and significance values *** ####
fn_EXO_num_lm_results<-t(sapply(ntw_variables_of_interest, function(x){
  formula <- as.formula(paste(x, "~", "Site_EXO_num", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  print(car::Anova(mod))
  mod_pval<-round(car::Anova(mod)[3], 3)
  F_val<-round(anova(mod)$`F value`, 3)
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  structure(c(mod_R2, F_val, mod_pval), names = c("R_squered", "F_Value", "P_value"))
}
))

write.xlsx2(fn_EXO_num_lm_results, file = "lmer_out.xlsx",
            sheetName = "fn_EXO_num", append = F)


  #### ** 1.3. Regression against exotic PROPORTIONS ** ####

      #### *** 1.3.1. Generating and saving plots *** ####
fn_EXO_prop_plots<-lapply(names(ntw_variables_of_interest), function(x){
  formula <- as.formula(paste(ntw_variables_of_interest[x], "~", "Site_EXO_prop", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  mod_pval<-ifelse(car::Anova(mod)[3]<0.001, "p<0.001", 
                   round(car::Anova(mod)[3], 3))
  text_col<-ifelse(car::Anova(mod)[3]<0.05, "red", "black")
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  
  effects_EXO_num <- as.data.frame(effects::effect(term= "Site_EXO_prop", mod= mod))
  dat<-tempdf[, c("site", "island", "Site_EXO_prop", ntw_variables_of_interest[x])]
  
  ggplot() + 
    geom_point(data=dat, aes_string("Site_EXO_prop", ntw_variables_of_interest[x]), color="light grey") + 
    geom_point(data=effects_EXO_num, aes(x=Site_EXO_prop, y=fit)) +
    geom_line(data=effects_EXO_num, aes(x=Site_EXO_prop, y=fit), color=text_col[1]) +
    geom_ribbon(data= effects_EXO_num, aes(x=Site_EXO_prop, ymin=lower, ymax=upper), alpha= 0.3, fill=text_col[1])+
    annotate("text", x = max(dat[,"Site_EXO_prop"])*0.9, 
             y = max(dat[,ntw_variables_of_interest[x]])*0.9, 
             label = paste("P =", mod_pval, ", R2 =", mod_R2), color=text_col[1])+
    ylab(label=x)+
    xlab(label="Proportion of exotic species")+
    theme_classic()
  
  
}
)
ggsave(plot = marrangeGrob(fn_EXO_prop_plots, nrow=3, ncol=3),
       filename="SM6.pdf", 
       device = "pdf", height = 21, width = 21,
       units = "in")



      #### *** 1.3.2. Storing df of R2 and significance values *** ####
fn_EXO_prop_lm_results<-t(sapply(names(ntw_variables_of_interest), function(x){
  formula <- as.formula(paste(ntw_variables_of_interest[x], "~", "Site_EXO_prop", "+", "(1|island)"))
  print(formula)
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  mod_pval<-round(car::Anova(mod)[3], 3)
  F_val<-round(anova(mod)$`F value`, 3)
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  structure(c(mod_R2, F_val, mod_pval), names = c("R_squered", "F_Value", "P_value"))
}
))

write.xlsx2(fn_EXO_prop_lm_results, file = "lmer_out.xlsx",
            sheetName = "fn_EXO_prop", append = T)


  #### ** 1.4. Regression against native NUMBERS ** ####

# Since there is a correlation between the number of natives and the number of nodes, 
# we double check if the correlations between number of nodes and number/proportion of exotics still
# exist if we remove the effect of natives from the system (i.e. run model on the residuals)
# residuals are stored here

    #### *** 1.4.1. Generating and saving plots *** ####
NAT_num_plots<-lapply(ntw_variables_of_interest, function(x){
  formula <- as.formula(paste(x, "~", "Site_NAT_num", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  mod_pval<-ifelse(car::Anova(mod)[3]<0.001, "p<0.001", 
                   round(car::Anova(mod)[3], 3))
  text_col<-ifelse(car::Anova(mod)[3]<0.05, "red", "black")
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  
  effects_NAT_num <- as.data.frame(effects::effect(term= "Site_NAT_num", mod= mod))
  dat<-tempdf[, c("site", "island", "Site_NAT_num", x)]
  
  print(paste("P =", mod_pval, ", R2 =", mod_R2))
  ggplot() + 
    geom_point(data=dat, aes_string("Site_NAT_num", x), color="light grey") + 
    geom_point(data=effects_NAT_num, aes(x=Site_NAT_num, y=fit)) +
    geom_line(data=effects_NAT_num, aes(x=Site_NAT_num, y=fit), color=text_col[1]) +
    geom_ribbon(data= effects_NAT_num, aes(x=Site_NAT_num, ymin=lower, ymax=upper), alpha= 0.3, fill=text_col[1])+
    annotate("text", x = max(dat[,"Site_NAT_num"])*0.9, y = max(dat[,x])*0.9, label = paste("P =", mod_pval, ", R2 =", mod_R2), color=text_col[1])+
    theme_classic()
  
  
}
)

      


# adding the one extra plot


mod<-lme4::lmer(Site_EXO_num ~Site_NAT_num + (1|island), dat=tempdf)

anova(mod)
effects_NAT_num <- as.data.frame(effects::effect(term= "Site_NAT_num", mod= mod))
mod_pval<-ifelse(car::Anova(mod)[3]<0.001, "p<0.001", 
                 round(car::Anova(mod)[3], 3))
text_col<-ifelse(car::Anova(mod)[3]<0.05, "red", "black")
mod_R2<-round(r.squaredGLMM(mod)[2], 3)
NAT_num_plots[["EXO_NAT"]]<-ggplot() + 
  geom_point(data=tempdf, aes(Site_NAT_num, Site_EXO_num), color="light grey") + 
  geom_point(data=effects_NAT_num, aes(x=Site_NAT_num, y=fit)) +
  geom_line(data=effects_NAT_num, aes(x=Site_NAT_num, y=fit), color=text_col[1]) +
  geom_ribbon(data= effects_NAT_num, aes(x=Site_NAT_num, ymin=lower, ymax=upper), alpha= 0.3, fill=text_col[1])+
  annotate("text", x = max(tempdf[,"Site_NAT_num"])*0.9, y = max(tempdf[,"Site_EXO_num"])*0.9, 
           label = paste("P =", mod_pval, ", R2 =", mod_R2), color=text_col[1])+
  theme_classic()




ggsave(plot = marrangeGrob(NAT_num_plots, nrow=2, ncol=2),
       filename="fn_B1_B2_NAT_num_plots_one_netw.pdf",  
       device = "pdf", height = 12, width = 12,
       units = "in")


    #### *** 1.4.2. Storing df of R2 and significance values, as well as a list of the residuals for later analysis *** ####
fn_NAT_num_lm_residuals<-as.data.frame(sapply(ntw_variables_of_interest, function(x){
  formula <- as.formula(paste(x, "~", "Site_NAT_num", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  resid(mod)
}
))

fn_NAT_num_lm_residuals$site<-tempdf$site
fn_NAT_num_lm_residuals$island<-tempdf$island
fn_NAT_num_lm_residuals$Site_EXO_num<-tempdf$Site_EXO_num
fn_NAT_num_lm_residuals$Site_EXO_prop<-sqrt(tempdf$Site_EXO_prop)






#### * 2. Regression for native subnetworks * ####



  #### ** 2.1. Preparing data ** ####
# making subnetworks
# splitting each site subnetwork to native and exotic networks
NEI_site_subnetworks<-lapply(site_subnetworks, function(x)
{subntwk<-list(
  NAT = induced_subgraph(x, which(V(x)$name %in% NAT)),
  EXO = induced_subgraph(x, which(V(x)$name %in% EXO))
)})

# liftin up one list level
NEI_site_subnetworks<-do.call(list, unlist(NEI_site_subnetworks, recursive=FALSE))



NEI_site_subnetworks_ntw_props<-lapply(NEI_site_subnetworks, netprop)
NEI_site_subnetw_prop_list<-lapply(names(NEI_site_subnetworks_ntw_props), 
                          function(x) data.frame(NEI = substr(x, (nchar(x)-2), nchar(x)), site = substr(x, 1, (nchar(x)-4)),
                                                 t(NEI_site_subnetworks_ntw_props[[x]][["netw_prop"]])))
                          
NEI_site_subnetworks_ntw_props_df<-do.call("rbind", NEI_site_subnetw_prop_list)

# calcualting modularity
NEI_site_subnetworks_modules<-lapply(NEI_site_subnetworks, 
                                       function(x)edge.betweenness.community(x) )                                


NEI_site_subnetworks_ntw_props_df$modularity<- 
  sapply(NEI_site_subnetworks_modules, function(x) modularity(x))

NEI_site_subnetworks_ntw_props_df$realised_cliques<-
  NEI_site_subnetworks_ntw_props_df$clique_num/NEI_site_subnetworks_ntw_props_df$max_num_cliques

# selecting NAT networks only and adding exoric numbers and proportion from the FULL SITE SUBNETWORK

NAT_site_props<-NEI_site_subnetworks_ntw_props_df[NEI_site_subnetworks_ntw_props_df$NEI=="NAT", ]
identical(NAT_site_props$site, site_subnetworks_props_df$site)

NAT_site_props$Site_EXO_num<-site_subnetworks_props_df$Site_EXO_num
NAT_site_props$Site_EXO_prop<-site_subnetworks_props_df$Site_EXO_prop
NAT_site_props$island<-substr(NAT_site_props$site, 1, 3)
NAT_site_props$NEI<-NULL


# colnames(ntw_variables_of_interest)[-1]
ntw_variables_of_interest<-colnames(site_subnetworks_props_df)[c(3, 7, 5, 11, 20, 43, 23, 26, 31, 34, 35, 44, 48)]
names(ntw_variables_of_interest)<-c("Number of nodes", "Number of edges", "Proportion of isolated nodes",
                                    "Proportion of negative edges", "Mean degree", "Connectance",
                                    "Normalised degree centrality",
                                    "Normalised eigenvector centrality",
                                    "Normalised closeness centrality",
                                    "Normalised betweenness centrality",
                                    "Network diameter",
                                    "Motif number",
                                    "Modularity")


tempdf<-NAT_site_props

colnames(tempdf)
# tempdf<-as.data.frame(sapply(colnames(tempdf), function(x) {print(x)
#   if (grepl("prop", x)) 
#   {sqrt(tempdf[,x])} else
#   {tempdf[,x]}}))
# tempdf[,ntw_variables_of_interest]<-sapply(ntw_variables_of_interest,  
#                                            function(x) as.numeric(tempdf[,x]))

# tempdf$Site_EXO_num<-as.numeric(tempdf$Site_EXO_num)
# tempdf$Site_EXO_prop<-as.numeric(tempdf$Site_EXO_prop)


# square root transform proportions
tempdf$Site_EXO_prop<-sqrt(tempdf$Site_EXO_prop)
tempdf$prop_isolated<-sqrt(tempdf$prop_isolated)
tempdf$prop_isolated<-sqrt(tempdf$prop_isolated)
tempdf$realised_degree<-sqrt(tempdf$realised_degree)
tempdf$neg_edge_ratio<-sqrt(tempdf$neg_edge_ratio)

# the only zero value is removed
tempdf<-tempdf[!tempdf$Site_EXO_prop==0, ]

# adding the number of natives as well
tempdf$Site_NAT_num<-tempdf$num_nodes-tempdf$Site_EXO_num
# tempdf<-na.omit(tempdf)


  #### ** 2.2. Regression against exotic NUMBERS ** ####

colnames(tempdf)
      #### *** 2.2.1. Generating and saving plots *** ####
nn_EXO_num_plots<-lapply(names(ntw_variables_of_interest), function(x){
  xx<-ntw_variables_of_interest[x]
  formula <- as.formula(paste(xx, "~", "Site_EXO_num", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  mod_pval<-ifelse(car::Anova(mod)[3]<0.001, "p<0.001", 
                   round(car::Anova(mod)[3], 3))
  text_col<-ifelse(car::Anova(mod)[3]<0.05, "red", "black")
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  
  effects_EXO_num <- as.data.frame(effects::effect(term= "Site_EXO_num", mod= mod))
  dat<-tempdf[, c("site", "island", "Site_EXO_num", xx)]
  
  print(paste("P =", mod_pval, ", R2 =", mod_R2))
  ggplot() + 
    geom_point(data=dat, aes_string("Site_EXO_num", xx), color="light grey") + 
    geom_point(data=effects_EXO_num, aes(x=Site_EXO_num, y=fit)) +
    geom_line(data=effects_EXO_num, aes(x=Site_EXO_num, y=fit), color=text_col[1]) +
    geom_ribbon(data= effects_EXO_num, aes(x=Site_EXO_num, ymin=lower, ymax=upper), alpha= 0.3, fill=text_col[1])+
    annotate("text", x = max(dat[,"Site_EXO_num"])*0.9, y = max(dat[,xx])*0.9, 
             label = paste("P =", mod_pval, ", R2 =", mod_R2), color=text_col[1])+
    theme_classic()
  
  
}
)

ggsave(plot = marrangeGrob(nn_EXO_num_plots, nrow=2, ncol=2),
       filename="nn_B1_B2_NAT_netw_EXO_num_plots_one_netw.pdf",  
       device = "pdf", height = 12, width = 12,
       units = "in")

ggarrange(plotlist = nn_EXO_num_plots, nrow=3, ncol=5)

      #### *** 2.2.2. Storing df of R2 and significance values *** ####

nn_EXO_num_lm_results<-t(sapply(names(ntw_variables_of_interest), function(x){
  xx<-ntw_variables_of_interest[x]
  formula <- as.formula(paste(xx, "~", "Site_EXO_num", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  mod_pval<-round(car::Anova(mod)[3], 3)
  F_val<-round(anova(mod)$`F value`, 3)
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  structure(c(mod_R2, F_val, mod_pval), names = c("R_squered", "F_Value", "P_value"))
}
))

write.xlsx2(nn_EXO_num_lm_results, file = "lmer_out.xlsx",
            sheetName = "nn_EXO_num", append = T)

  #### ** 2.3. Regression against exotic PROPORTIONS ** ####
      #### *** 2.3.1. Generating and saving plots *** ####



nn_EXO_prop_plots<-lapply(names(ntw_variables_of_interest), function(x){
  xx<-ntw_variables_of_interest[x]
  formula <- as.formula(paste(xx, "~", "Site_EXO_prop", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  mod_pval<-ifelse(car::Anova(mod)[3]<0.001, "p<0.001", 
                   round(car::Anova(mod)[3], 3))
  text_col<-ifelse(car::Anova(mod)[3]<0.05, "red", "black")
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  
  effects_EXO_num <- as.data.frame(effects::effect(term= "Site_EXO_prop", mod= mod))
  dat<-tempdf[, c("site", "island", "Site_EXO_prop", xx)]
  ggplot() + 
    geom_point(data=dat, aes_string("Site_EXO_prop", xx), color="light grey") + 
    geom_point(data=effects_EXO_num, aes(x=Site_EXO_prop, y=fit)) +
    geom_line(data=effects_EXO_num, aes(x=Site_EXO_prop, y=fit), color=text_col[1]) +
    geom_ribbon(data= effects_EXO_num, aes(x=Site_EXO_prop, ymin=lower, ymax=upper), alpha= 0.3, fill=text_col[1])+
    annotate("text", x = max(dat[,"Site_EXO_prop"])*0.9, 
             y = max(dat[,xx])*0.9, 
             label = paste("P =", mod_pval, ", R2 =", mod_R2), color=text_col[1])+
    ylab(label=x)+
    xlab(label="Proportion of exotic species")+
    theme_classic()
  
  
}
)
ggsave(plot = marrangeGrob(nn_EXO_prop_plots, nrow=3, ncol=3),
       filename="SM8.pdf", 
       device = "pdf", height = 21, width = 21,
       units = "in")

ggarrange(plotlist = nn_EXO_prop_plots, nrow=3, ncol=5)



      #### *** 2.3.2. Storing df of R2 and significance values *** ####

nn_EXO_prop_lm_results<-t(sapply(ntw_variables_of_interest, function(x){
  formula <- as.formula(paste(x, "~", "Site_EXO_prop", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  mod_pval<-round(car::Anova(mod)[3], 3)
  F_val<-round(anova(mod)$`F value`, 3)
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  structure(c(mod_R2, F_val, mod_pval), names = c("R_squered", "F_Value", "P_value"))
}
))

write.xlsx2(nn_EXO_prop_lm_results, file = "lmer_out.xlsx",
            sheetName = "nn_EXO_prop", append = T)








#### ** 2.4. Regression against native NUMBERS ** ####

# Since there is a correlation between the number of natives and the number of nodes, 
# we double check if the correlations between number of nodes and number/proportion of exotics still
# exist if we remove the effect of natives from the system (i.e. run model on the residuals)
# residuals are stored here

#### *** 2.4.1. Generating and saving plots *** ####
nn_NAT_num_plots<-lapply(names(ntw_variables_of_interest), function(x){
  print(x)
  xx<-ntw_variables_of_interest[x]
  formula <- as.formula(paste(xx, "~", "Site_NAT_num", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf, REML = T)
  print(anova(mod))
  mod_pval<-ifelse(car::Anova(mod)[3]<0.001, "p<0.001", 
                   round(car::Anova(mod)[3], 3))
  text_col<-ifelse(car::Anova(mod)[3]<0.05, "red", "black")
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  
  effects_NAT_num <- as.data.frame(effects::effect(term= "Site_NAT_num", mod= mod))
  dat<-tempdf[, c("site", "island", "Site_NAT_num", xx)]
  
  print(paste("P =", mod_pval, ", R2 =", mod_R2))
  ggplot() + 
    geom_point(data=dat, aes_string("Site_NAT_num", xx), color="light grey") + 
    geom_point(data=effects_NAT_num, aes(x=Site_NAT_num, y=fit)) +
    geom_line(data=effects_NAT_num, aes(x=Site_NAT_num, y=fit), color=text_col[1]) +
    geom_ribbon(data= effects_NAT_num, aes(x=Site_NAT_num, ymin=lower, ymax=upper), alpha= 0.3, fill=text_col[1])+
    annotate("text", x = max(dat[,"Site_NAT_num"])*0.9, y = max(dat[,xx])*0.9, 
             label = paste("P =", mod_pval, ", R2 =", mod_R2), color=text_col[1])+
    theme_classic()
  
  
}
)




# adding the one extra plot


mod<-lme4::lmer(Site_EXO_num ~Site_NAT_num + (1|island), dat=tempdf)

anova(mod)
effects_NAT_num <- as.data.frame(effects::effect(term= "Site_NAT_num", mod= mod))
mod_pval<-ifelse(car::Anova(mod)[3]<0.001, "p<0.001", 
                 round(car::Anova(mod)[3], 3))
text_col<-ifelse(car::Anova(mod)[3]<0.05, "red", "black")
mod_R2<-round(r.squaredGLMM(mod)[2], 3)
NAT_num_plots[["EXO_NAT"]]<-ggplot() + 
  geom_point(data=tempdf, aes(Site_NAT_num, Site_EXO_num), color="light grey") + 
  geom_point(data=effects_NAT_num, aes(x=Site_NAT_num, y=fit)) +
  geom_line(data=effects_NAT_num, aes(x=Site_NAT_num, y=fit), color=text_col[1]) +
  geom_ribbon(data= effects_NAT_num, aes(x=Site_NAT_num, ymin=lower, ymax=upper), alpha= 0.3, fill=text_col[1])+
  annotate("text", x = max(tempdf[,"Site_NAT_num"])*0.9, y = max(tempdf[,"Site_EXO_num"])*0.9, 
           label = paste("P =", mod_pval, ", R2 =", mod_R2), color=text_col[1])+
  theme_classic()




ggsave(plot = marrangeGrob(nn_NAT_num_plots, nrow=2, ncol=2),
       filename="nn_B1_B2_NAT_num_plots_one_netw.pdf",  
       device = "pdf", height = 12, width = 12,
       units = "in")


#### *** 2.4.2. Storing df of R2 and significance values, as well as a list of the residuals for later analysis *** ####
nn_NAT_num_lm_residuals<-as.data.frame(sapply(names(ntw_variables_of_interest), function(x){
  xx<-ntw_variables_of_interest[x]
  formula <- as.formula(paste(xx, "~", "Site_NAT_num", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  resid(mod)
}
))

nn_NAT_num_lm_residuals$site<-tempdf$site
nn_NAT_num_lm_residuals$island<-tempdf$island
nn_NAT_num_lm_residuals$Site_EXO_num<-tempdf$Site_EXO_num
nn_NAT_num_lm_residuals$Site_EXO_prop<-sqrt(tempdf$Site_EXO_prop)





#### * 3. Regression on the residuals from the regression on full network * ####

    #### ** 3.1. Preparing data ** ####
tempdf<-fn_NAT_num_lm_residuals
colnames(tempdf)[colnames(tempdf) %in% names(ntw_variables_of_interest)]<-ntw_variables_of_interest[names(ntw_variables_of_interest) %in% colnames(tempdf)]

    #### ** 3.2. Regression against exotic NUMBERS ** ####

        #### *** 3.2.1. Generating and saving plots *** ####
res_EXO_num_plots<-lapply(ntw_variables_of_interest[-c(7,8,11,12,14)], function(x){
  # x = "num_nodes"
  # rm(x)
  formula <- as.formula(paste(x, "~", "Site_EXO_num", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  mod_pval<-ifelse(car::Anova(mod)[3]<0.001, "p<0.001", 
                   round(car::Anova(mod)[3], 3))
  text_col<-ifelse(car::Anova(mod)[3]<0.05, "red", "black")
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  
  effects_EXO_num <- as.data.frame(effects::effect(term= "Site_EXO_num", mod= mod))
  dat<-tempdf[, c("site", "island", "Site_EXO_num", x)]
  print(paste("P =", mod_pval, ", R2 =", mod_R2))
  ggplot() + 
    geom_point(data=dat, aes_string("Site_EXO_num", x), color="light grey") + 
    geom_point(data=effects_EXO_num, aes(x=Site_EXO_num, y=fit)) +
    geom_line(data=effects_EXO_num, aes(x=Site_EXO_num, y=fit), color=text_col[1]) +
    geom_ribbon(data= effects_EXO_num, aes(x=Site_EXO_num, ymin=lower, ymax=upper), alpha= 0.3, fill=text_col[1])+
    annotate("text", x = max(dat[,"Site_EXO_num"])*0.9, y = max(dat[,x])*0.9, 
             label = paste("P =", mod_pval, ", R2 =", mod_R2), color=text_col[1])+
    theme_classic()
  
  
}
)
ggsave(plot = marrangeGrob(res_EXO_num_plots, nrow=2, ncol=2),
       filename="res_B1_B2_EXO_num_plots_one_netw.pdf",  
       device = "pdf", height = 12, width = 12,
       units = "in")

        #### *** 3.2.2. Storing df of R2 and significance values *** ####
fn_EXO_num_lm_results<-t(sapply(ntw_variables_of_interest[-c(7,8,11,12,14)], function(x){
  formula <- as.formula(paste(x, "~", "Site_EXO_num", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  mod_pval<-round(car::Anova(mod)[3], 3)
  F_val<-round(anova(mod)$`F value`, 3)
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  structure(c(mod_R2, F_val, mod_pval), names = c("R_squered", "F_Value", "P_value"))
}
))

write.xlsx2(fn_EXO_num_lm_results, file = "lmer_out.xlsx",
            sheetName = "fn_res_EXO_num", append = T)






#### * 4. Regression on the residuals from the regression on native network * ####

#### ** 4.1. Preparing data ** ####
tempdf<-nn_NAT_num_lm_residuals
colnames(tempdf)[colnames(tempdf) %in% names(ntw_variables_of_interest)]<-ntw_variables_of_interest[names(ntw_variables_of_interest) %in% colnames(tempdf)]


#### ** 4.2. Regression against exotic NUMBERS ** ####

#### *** 4.2.1. Generating and saving plots *** ####
res_EXO_num_plots<-lapply(ntw_variables_of_interest, function(x){
  print(x)
  # x = "num_nodes"
  # rm(x)
  formula <- as.formula(paste(x, "~", "Site_EXO_num", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  mod_pval<-ifelse(car::Anova(mod)[3]<0.001, "p<0.001", 
                   round(car::Anova(mod)[3], 3))
  text_col<-ifelse(car::Anova(mod)[3]<0.05, "red", "black")
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  
  effects_EXO_num <- as.data.frame(effects::effect(term= "Site_EXO_num", mod= mod))
  dat<-tempdf[, c("site", "island", "Site_EXO_num", x)]
  print(paste("P =", mod_pval, ", R2 =", mod_R2))
  ggplot() + 
    geom_point(data=dat, aes_string("Site_EXO_num", x), color="light grey") + 
    geom_point(data=effects_EXO_num, aes(x=Site_EXO_num, y=fit)) +
    geom_line(data=effects_EXO_num, aes(x=Site_EXO_num, y=fit), color=text_col[1]) +
    geom_ribbon(data= effects_EXO_num, aes(x=Site_EXO_num, ymin=lower, ymax=upper), alpha= 0.3, fill=text_col[1])+
    annotate("text", x = max(dat[,"Site_EXO_num"])*0.9, y = max(dat[,x])*0.9, 
             label = paste("P =", mod_pval, ", R2 =", mod_R2), color=text_col[1])+
    theme_classic()
  
  
}
)
ggsave(plot = marrangeGrob(res_EXO_num_plots, nrow=2, ncol=2),
       filename="res_B1_B2_EXO_num_plots_one_netw.pdf",  
       device = "pdf", height = 12, width = 12,
       units = "in")

#### *** 4.2.2. Storing df of R2 and significance values *** ####
nn_EXO_num_lm_results<-t(sapply(ntw_variables_of_interest, function(x){
  formula <- as.formula(paste(x, "~", "Site_EXO_num", "+", "(1|island)"))
  mod<-lme4::lmer(formula, data = tempdf)
  print(anova(mod))
  mod_pval<-round(car::Anova(mod)[3], 3)
  F_val<-round(anova(mod)$`F value`, 3)
  mod_R2<-round(r.squaredGLMM(mod)[2], 3)
  structure(c(mod_R2, F_val, mod_pval), names = c("R_squered", "F_Value", "P_value"))
}
))




    


