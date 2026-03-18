library(tidyverse)
library(vegan)
library(patchwork)
library(nlme)


meta<-read_csv("C:/Users/jnh23xdx/OneDrive - Bangor University/Work/Madog Trip/Madog Results/pre_cruise_pressure.csv")
meta<-meta %>%
  mutate(Station=as.factor(ID),
         Folk=as.factor(Folk),
         Folk_simplified=as.factor(Folk_simplified),
         Pressure_raw = Pressure^2)
# meta<-meta[,c(7,2:6,8)]
trawl_meta <-read_csv("C:/Users/jnh23xdx/OneDrive - Bangor University/Work/Madog Trip/Madog Results/Beam Trawl/Trawl_id_pressures.csv")
trawl_meta <- trawl_meta %>%
  rename(trawl_id = Trawl_id) %>%
  mutate(trawl_id = as.factor(trawl_id))



benthos<-read_csv("C:/Users/jnh23xdx/OneDrive - Bangor University/Work/Madog Trip/Madog Results/Beam Trawl/Species data/Beam Trawl Benthos_edited for NA and zero.csv")
ben_p <- benthos %>%
  mutate(Station = as.factor(Station),
         Number = ifelse(Number =="P", 0, Number), ## Set to 0 or 1?
         Number = as.numeric(Number),
         Number = Number/Subsample_perc,
         Weight = Weight/Subsample_perc,
         trawl_id = as.factor(trawl_id),
         t_length = trawl_meta$length[match(trawl_id, trawl_meta$trawl_id)],
         length_scale = 500/t_length,
         Number = Number*length_scale,
         Weight = Weight*length_scale)
ben_p_2 <- ben_p %>%
  group_by(trawl_id, Species) %>% #Station or trawl_id
  summarise(count = sum(Number))
ben_p_2.1 <- ben_p_2 %>%
  mutate(Station = ben_p$Station[match(trawl_id, ben_p$trawl_id)],
                Station = as.factor(Station))
ben_p_3<-pivot_wider(ben_p_2, names_from = Species, values_from = count, values_fill = 0)
ben_station<-as.character(ben_p_3$trawl_id)
ben_p_3<-ben_p_3[,2:71]
rownames(ben_p_3) <- ben_station

ben_m_2 <- ben_p %>%
  group_by(trawl_id, Species) %>% #Station or trawl_id
  summarise(biomass = sum(Weight))
ben_mass<-pivot_wider(ben_m_2, names_from = Species, values_from = biomass, values_fill = 0)
ben_mass<-ben_mass[,2:71]
rownames(ben_mass) <- ben_station


#### NMDS ####
epi_trans<-ben_p_3[rownames(ben_p_3) %in% div_meta$trawl_id, ] %>%
  mutate(across(c(1:70), ~.x^0.25)) ## ben_mass for biomass, ben_p_3 for abundance
benMDS<-metaMDS(epi_trans, try = 9999, trymax=9999, k = 2)
meta_sites = as.character(div_meta$trawl_id) ## Requires running some code further down first

nmds_SiteScores <-
  # get nmds scores 
  as.data.frame(scores(benMDS)$sites) %>%
  # change rownames (site) to a column 
  cbind(., meta_sites) %>%
  mutate(trawl_id = meta_sites,
         Station = ben_p$Station[match(trawl_id, ben_p$trawl_id)]
  ) %>%
  # join our habitat type (grouping variable) to each site 
  inner_join(., meta, by = "Station") %>%
  inner_join(., trawl_meta, by = "trawl_id") %>%
  mutate(Anchored = ifelse(Avg_pressure == 0, "no", "yes")) %>% 
  subset(., Folk_simplified != "Mud to Muddy Sand")
# Extract species scores
nmds_SpeciesScores <- 
  as.data.frame(scores(benMDS, "species"))
nmds_SpeciesScores$species <- rownames(nmds_SpeciesScores)  # create a column of species, from the rownames of species.scores
nmds_SpeciesScores <- cbind(nmds_SpeciesScores, pval.simp = simper.anchor$Yes_No$p)
sig.spp.scrs <- subset(nmds_SpeciesScores, pval.simp<=0.05)


# get centroid 
Habitat_Centroid <- 
  nmds_SiteScores %>% 
  group_by(Anchored) %>% 
  summarise(axis1 = mean(NMDS1),
            axis2 = mean(NMDS2)) %>% 
  ungroup()

# extract convex hull
habitat.hull <- 
  nmds_SiteScores %>% 
  group_by(Anchored) %>%
  slice(chull(NMDS1, NMDS2))

nmds_stress <- benMDS$stress

# use ggplot to plot 
ggplot() + 
  
  # add site scores
  geom_point(data = nmds_SiteScores, 
             aes(x=NMDS1, y=NMDS2, colour = Anchored, shape = Folk_simplified)) + ##size = Avg_pressure
  
  # add species scores
  # geom_text(data = nmds_SpeciesScores,
  # aes(x=NMDS1, y=NMDS2, label = species)) +
  # # 
  # add centroid
  # geom_point(data = Habitat_Centroid,
  #            aes(x = axis1, y = axis2, color = Anchored),
  #            size = 5, shape = 17) +
  #
  # Add arrows for species
  # geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2),
  #              arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  # ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = species),
  #                          cex = 3, direction = "both", segment.size = 0.25) + #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

  # 
  # # add convex hull
# geom_polygon(data = habitat.hull,
#              aes(x = NMDS1, y = NMDS2, fill = Anchored, group = Anchored),
#              alpha = 0.30) +
# # 
# add stress value
# annotate("text", x = 0.35, y = 0.75,
#          label = paste("2d stress =", round(nmds_stress, 3))) +
  # edit theme
  labs(x = "NMDS1", y = "NMDS2") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = .5),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5))

ordiplot(benMDS, display = "si", type = 'n')
points(benMDS)
orditorp(benMDS, display = 'sp')
fit<-envfit(benMDS ~ Folk_simplified + Avg_pressure, div_meta, perm = 999)
with(div_meta, ordispider(benMDS, Anchored, col="skyblue"))
plot(fit)
epi.spp.fit<-envfit(benMDS, epi_trans)
plot(epi.spp.fit, p.max = 0.01, col = "black", cex = 0.7)



perm_epi<-adonis2(epi_trans ~ Anchored+Folk_simplified, data = div_meta, permutations = 9999, by = "margin") # Significant but very small R2
perm_epi
summary(perm_epi)

simper.anchor <- with(div_meta, simper(epi_trans, Anchored))
summary(simper.anchor)

indspec.epi <- multipatt(epi_trans, cluster = div_meta$Anchored, 
                           control = how(nperm=9999), duleg = TRUE)
summary(indspec.epi)

ggplot(subset(epi_species, Species == "Buccinum undatum"), aes(x = Avg_pressure.x, y = log(biomass+1))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() 


#### Diversity metrics ####

H <- diversity(ben_p_3) ## Calculate indices
simp <- diversity(ben_p_3, index = "simpson")
S <- specnumber(ben_p_3)
S
J <- H/log(S)
J
diversities<-cbind(H,simp,S,J,ben_station) ## COmbine into one df
diversities<-as.data.frame(diversities[,c(5,1:4)]) %>% ## Make sure columns all correct
  mutate(trawl_id=as.factor(ben_station),
         H = as.numeric(H),
         simp = as.numeric(simp),
         J = as.numeric(J),
         S = as.numeric(S),
         Station = ben_p$Station[match(trawl_id, ben_p$trawl_id)])
diversities<-as.data.frame(diversities[,c(6,2:5,7)])
div_meta_epi<-inner_join(diversities, meta, by = "Station")
div_meta_epi<-inner_join(div_meta, trawl_meta, by = "trawl_id") %>% ## Create df with meta
  mutate(Anchored = as.factor(Anchored))

ben_mass <- ben_p[,c(1,2,3,6,7)] ## Same process with biomass
ben_mass <- ben_mass %>%
  mutate(trawl_id=as.factor(trawl_id)) %>%
  group_by(trawl_id, Species) %>% #Station or trawl_id
  summarise(biomass = sum(Weight), abundance = sum(Number))
Ben_mass_group <- ben_mass %>%
  group_by(trawl_id) %>%
  summarise(totbiomass = sum(biomass), tot.abund = sum(abundance))
complete_meta <- inner_join(Ben_mass_group, div_meta, by = "trawl_id") %>% ## combined df with meta, biomass, abundance and indices
  mutate(Anchored = ifelse(Avg_pressure == 0, "no", "yes"))

ggplot(subset(complete_meta, Folk != "(g)mS"), aes(x=log(Avg_pressure+1), y=S)) +
  geom_point(aes(colour = Folk_simplified)) +
  geom_smooth(method = "lm", aes(colour = Folk_simplified)) +
  # facet_wrap(~Folk_simplified,
  #          scales = "free") +
  theme_classic() +
  labs(x = "ln(Hours))",
       y = "ln(Abundance)",
       colour = "Sediment Type")

ggplot(subset(complete_meta, Folk != "(g)mS"), aes(y = S, fill = Anchored, x = Folk_simplified)) +
  geom_boxplot() +
  theme_classic()

beam_glm <- glm.nb(S~ Avg_pressure*Folk_simplified, data = subset(complete_meta, Folk != "(g)mS"))
summary.glm(beam_glm)

beam_glm_mass <-glm(log(totbiomass+1) ~ Avg_pressure, data = subset(complete_meta, Folk_simplified == "Coarse Sediment"))
summary(beam_glm_mass)

beam_glm_S <- glm(S~ log(Avg_pressure+1)*Folk_simplified, data = subset(complete_meta, Folk != "(g)mS"), family = poisson())
summary.glm(beam_glm_S)


tallied_epi <- ben_mass %>%
  group_by(Species) %>%
  summarize(frequency = n()) %>% 
  ungroup()

epi_species <- benthos %>%
  mutate(trawl_id = as.factor(trawl_id),
         Species = as.factor(Species)) %>%
  group_by(trawl_id, Species, .drop = FALSE) %>% #Station or trawl_id ()
  summarise(count = sum(as.numeric(Number)),
            biomass = sum(as.numeric(Weight)), 
            .groups = "drop") %>% ###
  inner_join(div_meta_epi, by = "trawl_id") %>%
  mutate(Anchored = ifelse(Avg_pressure.x == 0, "no", "yes"),
         count = as.numeric(count))




ggplot(subset(epi_species, Species == "Buccinum undatum" & Folk_simplified == "Coarse Sediment"), aes(x = Folk_simplified, y = log(biomass+1), fill = Anchored)) +
  geom_boxplot() +
  theme_classic()

ggplot(subset(epi_species, Species == "Pagurus bernhardus"), aes(x = Avg_pressure.x, y = log(biomass+1))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() 

hist(log(epi_species$biomass[epi_species$Species == "Buccinum undatum" & epi_species$Folk_simplified == "Coarse Sediment"]+1))

tallied_epi_Sand <- subset(epi_species, Folk_simplified == "Sand" & count >0) %>%
  group_by(Species) %>%
  summarize(frequency = n()) %>% 
  ungroup()
tally_epi_sand = slice_max(tallied_epi_Sand, n=6, order_by = frequency)
common_epi_df_sand = epi_species[epi_species$Species %in% tally_epi_sand$Species, ]

ggplot(subset(common_epi_df_sand, Folk_simplified == "Sand"), aes(y = log(biomass +1), x = Anchored)) +
  geom_boxplot() +
  facet_wrap(~Species,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic()

ggplot(subset(common_epi_df_sand, Folk_simplified == "Sand"), aes(y = log(count +1), x = Avg_pressure.x)) +
  geom_point() +
  facet_wrap(~Species, scale="free_y") +
  geom_smooth(method = "lm")




gam_epi<-mgcv::gam(S ~ s(Avg_pressure) +
            Folk_simplified, family = poisson(), data = complete_meta)
summary(gam_epi)
plot(gam_epi)

#### Classes ####

# benthos_class<- ben_p %>%
#   mutate(trawl_id=as.factor(trawl_id),
#          Class = as.factor(Class)) %>%
#   group_by(trawl_id, Class, .drop = FALSE) %>% #Station or trawl_id
#   summarise(biomass = sum(Weight), abundance = sum(Number), .groups = "drop") %>%
#   mutate(Station = ben_p$Station[match(trawl_id, ben_p$trawl_id)])
# 
# benthos_class_meta <- inner_join(benthos_class, meta, by = "Station")
# 
# ggplot(subset(benthos_class_meta, Folk_simplified == "Coarse Sediment"), aes(x=Anchored, y = abundance, fill = Class)) +
#   geom_bar(stat = "identity", position = "fill")
# 
# ggplot(subset(benthos_class_meta, Folk_simplified == "Sand"), aes(x=Anchored, y = abundance, fill = Class)) +
#   geom_bar(stat = "identity", position = "fill")
# 
# 
# ggplot(data = subset(benthos_class_meta, Folk != "(g)mS"), aes(y=log(biomass+1), fill = Anchored, x = Folk_simplified)) +
#   geom_boxplot() +
#   scale_y_continuous(limits = c(0, NA)) +
#   facet_wrap(~Class, scale="free_y")

#### Sessile only ####

benthos_sessile <- read_csv("C:/Users/jnh23xdx/OneDrive - Bangor University/Work/Madog Trip/Madog Results/Beam Trawl/Species data/Beam Trawl Benthos_edited for NA and zero.csv")

ben_s <- benthos %>%
  mutate(Station = as.factor(Station),
         Number = ifelse(Number =="P", 0, Number), ## Set to 0 or 1?
         Number = as.numeric(Number),
         Number = Number/Subsample_perc,
         Weight = Weight/Subsample_perc,
         trawl_id = as.factor(trawl_id))
ben_s_2 <- ben_s %>%
  group_by(trawl_id, Species) %>% #Station or trawl_id
  summarise(count = sum(Number))
ben_s_2.1 <- ben_s_2 %>%
  mutate(Station = ben_s$Station[match(trawl_id, ben_s$trawl_id)],
         Station = as.factor(Station))
ben_s_3<-pivot_wider(ben_s_2, names_from = Species, values_from = count, values_fill = 0)
ben_station<-as.character(ben_s_3$trawl_id)
ben_s_3<-ben_s_3[,2:71]
rownames(ben_s_3) <- ben_station

benthos_sessile_2 <- filter(benthos_sessile, mob_Sessile>0) %>%
  mutate(Station = as.factor(Station),
         trawl_id=as.factor(trawl_id),
         Class = as.factor(Class),
         Number = ifelse(Number=="P", 0, Number),
         Number = as.numeric(Number),
         Number = Number/Subsample_perc,
         Weight = Weight/Subsample_perc) %>%
  group_by(trawl_id, .drop = FALSE) %>% #Station or trawl_id
  summarise(biomass = sum(Weight), abundance = sum(Number), .groups = "drop") %>%
mutate(Station = benthos_sessile$Station[match(trawl_id, benthos_sessile$trawl_id)],
               Station = as.factor(Station))
benthos_sessile_meta <- inner_join(benthos_sessile_2, meta, by = "Station")

ggplot(subset(benthos_sessile_meta, Folk != "(g)mS"), aes(y = abundance, fill = Anchored, x = Folk_simplified)) +
  geom_boxplot()

ggplot(data = benthos_sessile_meta, aes(x=Pressure, y=log(biomass+1), color = Folk_simplified)) +
  geom_point() +
  geom_smooth(method = lm) +
  # geom_smooth(method = lm, inherit.aes = F, aes(x=Pressure, y=abundance), colour = "BLACK") +
  scale_y_continuous(limits = c(0, NA)) #+
  facet_wrap(~Class, scale="free_y")

sessile_grouped <- filter(benthos_sessile, mob_Sessile>0) %>%
  mutate(Station = as.factor(Station),
         trawl_id=as.factor(trawl_id),
         Class = as.factor(Class),
         Number = ifelse(Number=="P", 1, Number),
         Number = as.numeric(Number),
         Number = Number/Subsample_perc,
         Weight = Weight/Subsample_perc) %>%
  group_by(Station, .drop = FALSE) %>%
  summarise(totbiomass = sum(Weight), abundance = sum(Number), .groups = "drop") %>%
  mutate(Station = as.factor(Station)) %>%
  inner_join(., meta, by = "Station")

####Traits stuff ####

benthos_traits <- ben_p[1:60] %>%
  rowwise() %>%
  mutate(across(13:18, ~ . / sum(c_across(13:18))),
         across(19:24, ~ . / sum(c_across(19:24))),
         across(25:28, ~ . / sum(c_across(25:28))),
         across(29:32, ~ . / sum(c_across(29:32))),
         across(33:35, ~ . / sum(c_across(33:35))),
         across(36:41, ~ . / sum(c_across(36:41))),
         across(42:45, ~ . / sum(c_across(42:45))),
         across(46:51, ~ . / sum(c_across(46:51))),
         across(52:55, ~ . / sum(c_across(52:55))),
         across(56:60, ~ . / sum(c_across(56:60))))

epi_traits_2<-benthos_traits %>%
  rowwise() %>%
  mutate(across(13:18, ~ . * (Weight)),
         across(19:24, ~ . * (Weight)),
         across(25:28, ~ . * (Weight)),
         across(29:32, ~ . * (Weight)),
         across(33:35, ~ . * (Weight)),
         across(36:41, ~ . * (Weight)),
         across(42:45, ~ . * (Weight)),
         across(46:51, ~ . * (Weight)),
         across(52:55, ~ . * (Weight)),
         across(56:60, ~ . * (Weight)))

epi_traits_3 <- epi_traits_2 %>% 
  pivot_longer(cols = !c(1:12), names_to = "trait", values_to = "score") %>%
  mutate(Station = as.factor(Station),
         trawl_id=as.factor(trawl_id)) %>%
  inner_join(., meta, by = "Station")

epi_traits_4 <- epi_traits_3 %>%
  group_by(trawl_id, trait) %>%
  summarise(g_score = sum(score, na.rm = T)) %>% # mean or sum?
  inner_join(., div_meta, by = "trawl_id") %>% 
  mutate(t_group = regmatches(trait, regexpr("^.+?[_ ]", trait)))


dfs_epi = split(subset(epi_traits_4, Folk_simplified != "Mud to Muddy Sand"), f = subset(epi_traits_4, Folk_simplified != "Mud to Muddy Sand")$t_group)
# apply ggplot function and write to list
gg_l_epi = lapply(dfs_epi, function(x) {
  ggplot(x, aes(x = log(Avg_pressure+1), y = log(g_score+1), colour = trait), ylab = '') + 
    # geom_bar(position = "fill", stat = "identity") + 
    geom_point() + 
    geom_smooth(method = "lm") +    
    facet_wrap(vars(t_group), strip.position = "left", 
               labeller = as_labeller(
                 c(b_ = "Bioturbation", 
                   ed_ = "Egg Development", 
                   f_ = "Feeding Mechanism", 
                   sr_ = "Maximum Size", 
                   m_ = "Morphology", 
                   l_ = "Lifespan", 
                   ld_ = "Larva Development", 
                   lh_ = " Living Habit",
                   mob_ = "Mobility",
                   sp_ = "Sediment Position"),
                 label_wrap_gen(width = 10)),
               axes = "all") +
    theme(strip.background = element_blank(), 
          strip.placement = "outside", 
          axis.title.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.ticks.length.x = unit(0, "pt"),
          plot.margin = margin(0,0,4,0,"pt"),
          panel.background = element_rect(colour = "black", fill = NA),
          element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.height= unit(8, 'pt'),
          legend.key.width= unit(8, 'pt')) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(y=NULL)
})
# patchwork
wrap_plots(gg_l_epi) + plot_layout(axis_titles = "collect", guides = "auto", axes = "collect") 

gg_l_epi$f_


ggplot(data = subset(epi_traits_4, trait == "f_Scavenger"), aes(x= Avg_pressure, y=g_score)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  theme_classic() +
  geom_function(fun = function(x) exp(4.7435) * x^0.1132)  
  


dfs_epi_coarse = split(subset(epi_traits_4, Folk_simplified == "Coarse Sediment"), f = subset(epi_traits_4, Folk_simplified == "Coarse Sediment")$t_group)

gg_l_epi_point = lapply(dfs_epi_coarse, function(x) {
  ggplot(x, aes(x = log(Avg_pressure+1), y = log(g_score+1), colour = trait), ylab = '') + 
    geom_point() + 
    geom_smooth(method = "lm") +
    facet_wrap(vars(t_group), strip.position = "left", 
               labeller = as_labeller(
                 c(b_ = "Bioturbation", 
                   ed_ = "Egg Development", 
                   f_ = "Feeding Mechanism", 
                   sr_ = "Maximum Size", 
                   m_ = "Morphology", 
                   l_ = "Lifespan", 
                   ld_ = "Larva Development", 
                   lh_ = " Living Habit",
                   mob_ = "Mobility",
                   sp_ = "Sediment Position"),
                 label_wrap_gen(width = 10)),
               axes = "all") +
    theme(strip.background = element_blank(), 
          strip.placement = "outside", 
          axis.title.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.ticks.length.x = unit(0, "pt"),
          plot.margin = margin(0,0,4,0,"pt"),
          panel.background = element_rect(colour = "black", fill = NA),
          element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.height= unit(8, 'pt'),
          legend.key.width= unit(8, 'pt')) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(y=NULL)
})

wrap_plots(gg_l_epi_point) + plot_layout(axis_titles = "collect", guides = "auto", axes = "keep") 

wrap_plots(gg_l_epi_point$f_ + gg_l_epi_point$lh_ + gg_l_epi_point$m_ + gg_l_epi_point$mob_) + plot_layout(axis_titles = "collect", guides = "auto", axes = "keep") 

gg_l_epi_point_bar = lapply(dfs_epi_coarse, function(x) {
  ggplot(x, aes(x = Anchored, y = log(g_score+1), fill = trait), ylab = '') + 
    geom_bar(position = "fill", stat = "identity") + 
    facet_wrap(vars(t_group), strip.position = "left", 
               labeller = as_labeller(
                 c(b_ = "Bioturbation", 
                   ed_ = "Egg Development", 
                   f_ = "Feeding Mechanism", 
                   sr_ = "Maximum Size", 
                   m_ = "Morphology", 
                   l_ = "Lifespan", 
                   ld_ = "Larva Development", 
                   lh_ = " Living Habit",
                   mob_ = "Mobility",
                   sp_ = "Sediment Position"),
                 label_wrap_gen(width = 10)),
               axes = "all") +
    theme(strip.background = element_blank(), 
          strip.placement = "outside", 
          axis.title.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.ticks.length.x = unit(0, "pt"),
          plot.margin = margin(0,0,4,0,"pt"),
          panel.background = element_rect(colour = "black", fill = NA),
          element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.height= unit(8, 'pt'),
          legend.key.width= unit(8, 'pt')) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(y=NULL)
})

wrap_plots(gg_l_epi_point_bar$f_ + gg_l_epi_point_bar$lh_ + gg_l_epi_point_bar$m_ + gg_l_epi_point_bar$mob_) + plot_layout(axis_titles = "collect", guides = "auto", axes = "keep") 


dfs_epi_sand = split(subset(epi_traits_4, Folk_simplified == "Sand"), f = subset(epi_traits_4, Folk_simplified == "Sand")$t_group)

gg_l_epi_point_sand = lapply(dfs_epi_sand, function(x) {
  ggplot(x, aes(x = log(Avg_pressure+1), y = log(g_score+1), colour = trait), ylab = '') + 
    geom_point() + 
    geom_smooth(method = "lm") +
    facet_wrap(vars(t_group), strip.position = "left", 
               labeller = as_labeller(
                 c(b_ = "Bioturbation", 
                   ed_ = "Egg Development", 
                   f_ = "Feeding Mechanism", 
                   sr_ = "Maximum Size", 
                   m_ = "Morphology", 
                   l_ = "Lifespan", 
                   ld_ = "Larva Development", 
                   lh_ = " Living Habit",
                   mob_ = "Mobility",
                   sp_ = "Sediment Position"),
                 label_wrap_gen(width = 10)),
               axes = "all") +
    theme(strip.background = element_blank(), 
          strip.placement = "outside", 
          axis.title.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.ticks.length.x = unit(0, "pt"),
          plot.margin = margin(0,0,4,0,"pt"),
          panel.background = element_rect(colour = "black", fill = NA),
          element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.height= unit(8, 'pt'),
          legend.key.width= unit(8, 'pt')) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(y=NULL)
})

wrap_plots(gg_l_epi_point_sand) + plot_layout(axis_titles = "collect", guides = "auto", axes = "keep") 

traitglm <- glm(log(g_score+1) ~ Anchored*trait, data = subset(epi_traits_4, Folk_simplified == "Coarse Sediment"))
summary(traitglm)
traitglmsand <- glm(log(g_score+1) ~ log(Avg_pressure+1)*trait, data = subset(epi_traits_4, Folk_simplified == "Sand"))
summary(traitglmsand)
hist(log(epi_traits_4$g_score[epi_traits_4$Folk_simplified == "Coarse Sediment"]+1))
traitglmtot <- glm(log(g_score+1) ~ Avg_pressure*trait, data = epi_traits_4)
summary(traitglmtot)

epi.ses.glm <- glm(log(g_score+1) ~ Avg_pressure, data = subset(epi_traits_4, trait == "f_Predator" & Folk_simplified == "Coarse Sediment"))
summary(epi.ses.glm)

epi.trait.glm <- subset(epi_traits_4, Folk_simplified != "Mud to Muddy Sand") %>%
  nest(data = -trait) %>%
  mutate(fit = map(data, ~glm(log(g_score+1) ~ log(Avg_pressure+1), data = .x)),
         tidied = map(fit, broom::tidy)) %>%
  unnest(tidied) %>%
  mutate(p.value = round(p.value, digits = 3))  %>% 
  mutate(glance = map(fit, broom::glance)) %>% 
  unnest(glance)

epi.trait.glm.cs <- subset(epi_traits_4, Folk_simplified == "Coarse Sediment") %>%
  nest(data = -trait) %>%
  mutate(fit = map(data, ~glm(log(g_score+1) ~ log(Avg_pressure+1), data = .x)),
         tidied = map(fit, broom::tidy)) %>%
  unnest(tidied) %>%
  mutate(p.value = round(p.value, digits = 3))  %>% 
  mutate(glance = map(fit, broom::glance)) %>% 
  unnest(glance)

hist(log(epi_traits_4$Avg_pressure+1))

epi_trait_grouped <- epi_traits_4 %>%
  ungroup() %>%
  group_by(trawl_id, t_group) %>%
  mutate(prop_score = g_score/sum(g_score) * 100)

dfs_epi_g = split(subset(epi_trait_grouped, Folk_simplified != "Mud to Muddy Sand"), 
                  f = subset(epi_trait_grouped, Folk_simplified != "Mud to Muddy Sand")$t_group)

gg_l_epi = lapply(dfs_epi_g, function(x) {
  ggplot(x, aes(x = log(Avg_pressure+1), y = prop_score, colour = trait), ylab = '') + 
    # geom_bar(position = "fill", stat = "identity") +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(vars(t_group), strip.position = "left", 
               labeller = as_labeller(
                 c(b_ = "Bioturbation", 
                   ed_ = "Egg Development", 
                   f_ = "Feeding Mechanism", 
                   sr_ = "Maximum Size", 
                   m_ = "Morphology", 
                   l_ = "Lifespan", 
                   ld_ = "Larva Development", 
                   lh_ = " Living Habit",
                   mob_ = "Mobility",
                   sp_ = "Sediment Position"),
                 label_wrap_gen(width = 10)),
               axes = "all") +
    theme(strip.background = element_blank(), 
          strip.placement = "outside", 
          axis.title.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.ticks.length.x = unit(0, "pt"),
          plot.margin = margin(0,0,4,0,"pt"),
          panel.background = element_rect(colour = "black", fill = NA),
          element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.height= unit(8, 'pt'),
          legend.key.width= unit(8, 'pt')) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(y=NULL)
})
wrap_plots(gg_l_epi) + plot_layout(axis_titles = "collect", guides = "auto", axes = "collect") 



epi_nmds_df<-epi_traits_3 %>% 
  group_by(trawl_id, trait) %>%
  summarise(g_score = sum(score, na.rm = T)) %>%
  group_by(trawl_id) %>%
  pivot_wider(., names_from = trait, values_from = g_score, values_fill = 0) %>%
  ungroup()
meta_sites = as.character(div_meta$trawl_id)

epi_nmds_df<-epi_nmds_df[,2:49]
rownames(epi_nmds_df) <- meta_sites

epi_nmds_df <- epi_nmds_df %>%
  mutate(across(c(1:48), ~.x^0.25))
epi.trait.MDS<-metaMDS(epi_nmds_df, try = 9999, trymax=9999, k = 2)

nmds_SiteScores <-
  # get nmds scores 
  as.data.frame(scores(epi.trait.MDS)$sites) %>%
  # change rownames (site) to a column 
  cbind(., meta_sites) %>%
  mutate(trawl_id = meta_sites,
         Station = ben_p$Station[match(trawl_id, ben_p$trawl_id)]
  ) %>%
  # join our habitat type (grouping variable) to each site 
  inner_join(., meta, by = "Station") %>%
  inner_join(., trawl_meta, by = "trawl_id") %>%
  mutate(Anchored = ifelse(Avg_pressure == 0, "no", "yes")) %>% 
  subset(., Folk_simplified != "Mud to Muddy Sand")
# Extract species scores
nmds_SpeciesScores <- 
  as.data.frame(scores(benMDS, "species"))
nmds_SpeciesScores$species <- rownames(nmds_SpeciesScores)  # create a column of species, from the rownames of species.scores

# get centroid 
Habitat_Centroid <- 
  nmds_SiteScores %>% 
  group_by(Anchored) %>% 
  summarise(axis1 = mean(NMDS1),
            axis2 = mean(NMDS2)) %>% 
  ungroup()

# extract convex hull
habitat.hull <- 
  nmds_SiteScores %>% 
  group_by(Anchored) %>%
  slice(chull(NMDS1, NMDS2))

nmds_stress <- benMDS$stress

# use ggplot to plot 
ggplot() + 
  
  # add site scores
  geom_point(data = nmds_SiteScores, 
             aes(x=NMDS1, y=NMDS2, colour = Anchored)) + ##size = Avg_pressure
  
  # add species scores
  # geom_text(data = nmds_SpeciesScores,
  # aes(x=NMDS1, y=NMDS2, label = species)) +
  # # 
  # add centroid
  geom_point(data = Habitat_Centroid,
             aes(x = axis1, y = axis2, color = Anchored),
             size = 5, shape = 17) +
  # 
  # # add convex hull
  # geom_polygon(data = habitat.hull,
  #              aes(x = NMDS1, y = NMDS2, fill = Anchored, group = Anchored),
  #              alpha = 0.30) +
  # # 
  # add stress value
  annotate("text", x = 0.35, y = 0.75,
           label = paste("2d stress =", round(nmds_stress, 3))) +
  # edit theme
  labs(x = "NMDS1", y = "NMDS2") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = .5),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5))




