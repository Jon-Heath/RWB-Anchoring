inf_div_trans<-inf_div %>%
  mutate(across(c(1:370), ~.x^0.25))
infMDS<-metaMDS(inf_div_trans, try = 9999, trymax=9999, k = 2)
nmds_SiteScores <-
  # get nmds scores 
  as.data.frame(scores(infMDS)$sites) %>%
  # change rownames (site) to a column 
  cbind(., inf_station) %>%
  mutate(Station_id = inf_station,
         Station = as.factor(infauna$Station[match(Station_id, infauna$Station_id)]))%>%
  # join our habitat type (grouping variable) to each site 
  inner_join(., meta, by = "Station") %>% ### Meta from Neater Beam Trawl file currently
  inner_join(., grab_press, by ="Station_id") %>%
  mutate(Anchored = ifelse(Pressure.y == 0, "no", "yes"))


nmds_SpeciesScores <- 
  as.data.frame(scores(infMDS, "species"))
nmds_SpeciesScores$species <- rownames(nmds_SpeciesScores)
nmds_SpeciesScores <- cbind(nmds_SpeciesScores, pval.simp = simper.anchor.inf$yes_no$p)
sig.spp.scrs.inf <- subset(nmds_SpeciesScores, pval.simp<=0.05)

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
ggplot() + 
  
  # add site scores
  geom_point(data = nmds_SiteScores, 
             aes(x=NMDS1, y=NMDS2, colour = Anchored, size = Pressure.y))  +
  
  # geom_text(data = nmds_SpeciesScores,
  #                       aes(x=NMDS1, y=NMDS2, label = species),
  #                       size = 2) +
  # ## Add centroid
  # geom_point(data = Habitat_Centroid,
  #            aes(x = axis1, y = axis2, color = Anchored),
  #            size = 5, shape = 17) +
  # # # add convex hull
  # geom_polygon(data = habitat.hull,
  #              aes(x = NMDS1, y = NMDS2, fill = Anchored, group = Anchored),
  #              alpha = 0.30) +
  # annotate("text", x = 0.8, y = 0.75,
  #          label = paste("2d stress =", round(nmds_stress, 3))) +
  # # Add species arrows
  geom_segment(data = sig.spp.scrs.inf, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs.inf, aes(x=NMDS1, y=NMDS2, label = species),  
                           cex = 3, direction = "both", segment.size = 0.25) +
                           
  # edit theme
  labs(x = "NMDS1", y = "NMDS2") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = .5),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.2))         
perm_inf<-adonis2(inf_div_trans ~ Pressure.x*Folk_simplified, data = infdiv_meta, permutations = 9999) # Significant but very small R2
perm_inf
summary(perm_epi)

simper.anchor.inf <- with(infdiv_meta, simper(inf_div_trans, Anchored))
summary(simper.anchor.inf)

library(indicspecies)
indspec.inf <- multipatt(inf_div_trans, cluster = infdiv_meta$Anchored, duleg = TRUE)
summary(indspec.inf)

#### coarse #####

inf_div_coarse <- inf_div_trans[c(1:3,7:9,13:15,22:24,31:36,52:60),]
inf_coarse_press <- infdiv_meta[c(1:3,7:9,13:15,22:24,31:36,52:60),]
c.inf.stat <- inf_station[c(1:3,7:9,13:15,22:24,31:36,52:60)]

infMDS.coarse<-metaMDS(inf_div_coarse, try = 9999, trymax=9999, k = 2)
nmds_SiteScores <-
  # get nmds scores 
  as.data.frame(scores(infMDS.coarse)$sites) %>%
  # change rownames (site) to a column 
  cbind(., inf_station[c(1:3,7:9,13:15,22:24,31:36,52:60)]) %>%
  mutate(Station_id = inf_station[c(1:3,7:9,13:15,22:24,31:36,52:60)],
         Station = as.factor(infauna$Station[match(Station_id, infauna$Station_id)]))%>%
  # join our habitat type (grouping variable) to each site 
  inner_join(., grab_press, by ="Station_id") %>%
  mutate(Anchored = ifelse(Pressure == 0, "no", "yes"))


nmds_SpeciesScores <- 
  as.data.frame(scores(infMDS.coarse, "species"))
nmds_SpeciesScores$species <- rownames(nmds_SpeciesScores)

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
ggplot() + 
  
  # add site scores
  geom_point(data = nmds_SiteScores, 
             aes(x=NMDS1, y=NMDS2, colour = Anchored, size = Pressure))  +
  
  # geom_text(data = nmds_SpeciesScores,
  #                       aes(x=NMDS1, y=NMDS2, label = species),
  #                       size = 2) +
  # ## Add centroid
  geom_point(data = Habitat_Centroid,
             aes(x = axis1, y = axis2, color = Anchored),
             size = 5, shape = 17) +
  # # add convex hull
  geom_polygon(data = habitat.hull,
               aes(x = NMDS1, y = NMDS2, fill = Anchored, group = Anchored),
               alpha = 0.30) +
  annotate("text", x = 0.8, y = 0.75,
           label = paste("2d stress =", round(nmds_stress, 3))) +
  # edit theme
  labs(x = "NMDS1", y = "NMDS2") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = .5),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.2))    

perm_inf_c<-adonis2(inf_div_coarse ~ Anchored, data = inf_coarse_press, permutations = 9999) # Significant but very small R2
perm_inf_c
summary(perm_inf_c)

simper.anchor.inf <- with(inf_coarse_press, simper(inf_div_coarse, Anchored))
summary(simper.anchor.inf)

library(indicspecies)
indspec.inf.c <- multipatt(inf_div_coarse, cluster = inf_coarse_press$Anchored, 
                           control = how(nperm=9999), duleg = TRUE)
summary(indspec.inf.c)

#### Sand ####

inf_div_sand <- inf_div_trans[c(4:6,10:12,16:21,25:30,37:51),]
inf_sand_press <- infdiv_meta[c(4:6,10:12,16:21,25:30,37:51),]
s.inf.stat <- inf_station[c(4:6,10:12,16:21,25:30,37:51)]

infMDS.sand<-metaMDS(inf_div_sand, try = 9999, trymax=9999, k = 2)
nmds_SiteScores <-
  # get nmds scores 
  as.data.frame(scores(infMDS.sand)$sites) %>%
  # change rownames (site) to a column 
  cbind(., inf_station[c(4:6,10:12,16:21,25:30,37:51)]) %>%
  mutate(Station_id = inf_station[c(4:6,10:12,16:21,25:30,37:51)],
         Station = as.factor(infauna$Station[match(Station_id, infauna$Station_id)]))%>%
  # join our habitat type (grouping variable) to each site 
  inner_join(., grab_press, by ="Station_id") %>%
  mutate(Anchored = ifelse(Pressure == 0, "no", "yes"))


nmds_SpeciesScores <- 
  as.data.frame(scores(infMDS.sand, "species"))
nmds_SpeciesScores$species <- rownames(nmds_SpeciesScores)

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
ggplot() + 
  
  # add site scores
  geom_point(data = nmds_SiteScores, 
             aes(x=NMDS1, y=NMDS2, colour = Anchored, size = Pressure))  +
  
  # geom_text(data = nmds_SpeciesScores,
  #                       aes(x=NMDS1, y=NMDS2, label = species),
  #                       size = 2) +
  # ## Add centroid
  geom_point(data = Habitat_Centroid,
             aes(x = axis1, y = axis2, color = Anchored),
             size = 5, shape = 17) +
  # # add convex hull
  geom_polygon(data = habitat.hull,
               aes(x = NMDS1, y = NMDS2, fill = Anchored, group = Anchored),
               alpha = 0.30) +
  annotate("text", x = 0.8, y = 0.75,
           label = paste("2d stress =", round(nmds_stress, 3))) +
  # edit theme
  labs(x = "NMDS1", y = "NMDS2") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = .5),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.2))  

perm_inf_s<-adonis2(inf_div_sand ~ Anchored, data = inf_sand_press, permutations = 9999) # Significant but very small R2
perm_inf_s
summary(perm_inf_s)

simper.anchor.inf.s <- with(inf_sand_press, simper(inf_div_sand, Anchored))
summary(simper.anchor.inf.s)

library(indicspecies)
indspec.inf.s <- multipatt(inf_div_sand, cluster = inf_sand_press$Anchored, 
                           control = how(nperm=9999), duleg = TRUE)
summary(indspec.inf.s)

