library(tidyverse)
library(vegan)
library(MASS)


grab_press<-read_csv("C:/Users/jnh23xdx/OneDrive - Bangor University/Work/Madog Trip/Madog Results/Grab Samples/Grab Sample pressures.csv")
colnames(grab_press)[11]<-"Pressure"
colnames(grab_press)[4]<-"Station"
grab_press <- grab_press %>% 
  mutate(raw_pressure = Pressure^2,
       Station = as.factor(Station))

infauna<-read_csv("C:/Users/jnh23xdx/OneDrive - Bangor University/Work/Madog Trip/Madog Results/Grab Samples/Species data/Benthic infauna RWB anchoring traits.csv")
infauna<- infauna %>%
  mutate(Station_id = as.factor(Station_id),
         ScientificName_accepted = as.factor(ScientificName_accepted),
         Number = as.numeric(Number))

inf_2<-infauna %>%
  group_by(Station_id, ScientificName_accepted, .drop = FALSE) %>% #Station or trawl_id
  summarise(count = sum(Number))

inf_div<-pivot_wider(na.omit(inf_2), names_from = ScientificName_accepted, values_from = count, values_fill = 0)
inf_station<-as.character(inf_div$Station_id)
inf_div <- inf_div[,2:371]
rownames(inf_div) <- inf_station

#### NMDS ####
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
             aes(x=NMDS1, y=NMDS2, colour = Anchored, shape = Folk_simplified))  + ###, size = Pressure.y
  
  # geom_text(data = nmds_SpeciesScores,
  #                       aes(x=NMDS1, y=NMDS2, label = species),
  #                       size = 2) +
  # ## Add centroid
  # geom_point(data = Habitat_Centroid,
  #            aes(x = axis1, y = axis2, color = Anchored),
  #            size = 5, shape = 17) +
  # # add convex hull
  # geom_polygon(data = habitat.hull,
  #              aes(x = NMDS1, y = NMDS2, fill = Anchored, group = Anchored),
  #              alpha = 0.30) +
  # annotate("text", x = 0.8, y = 0.75,
  #        label = paste("2d stress =", round(nmds_stress, 3))) +
  # edit theme
  labs(x = "NMDS1", y = "NMDS2") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = .5),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.2))         
perm_inf<-adonis2(inf_div_trans ~ Anchored + Folk_simplified, data = infdiv_meta, permutations = 9999, by = "margin") # Significant but very small R2
perm_inf
summary(perm_epi)

simper.anchor.inf <- with(infdiv_meta, simper(inf_div_trans, Anchored))
summary(simper.anchor.inf)

library(indicspecies)
indspec.inf <- multipatt(inf_div_trans, cluster = infdiv_meta$Anchored, duleg = TRUE)
summary(indspec.inf)
#####

H <- diversity(inf_div)
simp <- diversity(inf_div, index = "simpson")
S <- specnumber(inf_div)
S
J <- H/log(S)
J
inf_diversities<-cbind(H,simp,S,J,inf_station)

inf_diversities<-as.data.frame(inf_diversities[,c(5,1:4)]) %>%
  mutate(Station_id=as.factor(inf_station),
         H = as.numeric(H),
         simp = as.numeric(simp),
         J = as.numeric(J),
         S = as.numeric(S),
         Station = as.factor(infauna$Station[match(Station_id, infauna$Station_id)]))
infdiv_meta<-inner_join(inf_diversities, meta, by = "Station")
infdiv_meta<-inner_join(infdiv_meta, grab_press, by = "Station_id")%>%
  mutate(Anchored = ifelse(Pressure.y == 0, "no", "yes"))

tallied <- infauna %>% group_by(ScientificName_accepted) %>%
  summarize(frequency = n()) %>% 
  ungroup()

tallied_class <- infauna %>% group_by(Class)%>%
  summarize(frequency = n()) %>% 
  ungroup()
  
inf_3 <- infauna %>%
  mutate(Station_id = as.factor(Station_id),
         ScientificName_accepted = as.factor(ScientificName_accepted)) %>%
  group_by(Station_id, ScientificName_accepted, .drop = FALSE) %>% #Station or trawl_id ()
  summarise(count = sum(Number), .groups = "drop") %>% ###()
  inner_join(grab_press, by = "Station_id") %>%
  inner_join(subset(meta, select = -c(Pressure, Pressure_raw)), by = "Station")%>%
  mutate(Anchored = ifelse(Pressure == 0, "no", "yes"))

inf_abundance <- subset(infauna, !is.na(Number)) %>%
  mutate(Station_id = as.factor(Station_id)) %>%
  group_by(Station_id) %>%
  summarise(tot_abundance = sum(Number))

infdiv_meta<-inner_join(infdiv_meta, inf_abundance, by = "Station_id")

ggplot(infdiv_meta, aes(y = S, x = Folk_simplified, fill = Anchored)) + ## Different to result of glm on anchoring pressure
  geom_boxplot() +
  # facet_wrap(~Folk_simplified)
  theme_classic()

ggplot(data = infdiv_meta, aes(x = Pressure.y, y = S)) +
  geom_point() +
  geom_function(data = subset(infdiv_meta, Folk_simplified == "Sand"), fun = function(x) exp((3.9858371 - 0.6988028 - 0.0004274*x + 0.001623*x))) +  
  geom_function(data = subset(infdiv_meta, Folk_simplified == "Coarse Sediment"), fun = function(x) exp(3.9858371 + (-0.0004274*x))) +
  theme_classic() +
  facet_wrap(~Folk_simplified, 
             scale="free", 
             axis = "all") +
  labs(x = "Anchor Hours",
       y = "Species Richness") 

ggplot(data = infdiv_meta, aes(x = Pressure.y, y = S)) +
  geom_point() +
  theme_classic() +
  labs(x = "Anchor Hours",
       y = "Species Richness") +
  geom_smooth(method = "lm") #+
  facet_wrap(~Folk_simplified, scale="free_x") 

ggplot(data = infdiv_meta, aes(x = Pressure.y, y = tot_abundance)) +
  geom_point() +
  facet_wrap(~Folk_simplified, 
             scale="free", 
             axis = "all") +
  theme_classic() +
  geom_function(data = subset(infdiv_meta, Folk_simplified == "Coarse Sediment"), fun = function(x) exp(5.5292361 + (-0.0005614*x))) +
  geom_function(data = subset(infdiv_meta, Folk_simplified == "Sand"), fun = function(x) exp((5.5292361 - 0.8003071  - 0.0005614*x + 0.0029111*x)))
  


ggplot(subset(inf_3, ScientificName_accepted == "Urothoe elegans"), aes(x = Anchored, y = count)) +
  geom_boxplot()

#### Traits ####

inf_5 <- inner_join(inf_3, subset(infauna, select = c(10:68)), by = "ScientificName_accepted", multiple = "first") %>%
  group_by(Station_id, f_Surface_deposit) %>%
  summarise(n = sum(count)) %>%
  mutate(Station_id = Station_id.x) %>%###()
  inner_join(grab_press, by = "Station_id") %>%
  inner_join(subset(meta, select = -c(Pressure, Pressure_raw)), by = "Station")%>%
  mutate(Anchored = ifelse(Pressure == 0, "no", "yes")) %>%
  na.omit(.)

ggplot(inf_5, aes(x=Anchored, y = n, fill = as.factor(f_Surface_deposit))) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Folk_simplified)

ggplot(subset(inf_5, f_Subsurface_deposit != "#N/A"), aes(x=Pressure, y = n)) +
  geom_pie()
  facet_wrap(~as.factor(f_Subsurface_deposit), scale = "free")
  
  
traits_v2 <- infauna %>% 
    pivot_longer(cols = !c(1:19, 68), names_to = "trait", values_to = "score") %>%
    mutate(Station = as.factor(Station),
           Station_id=as.factor(Station_id)) %>%
    inner_join(., meta, by = "Station")
  
traits_v3 <- traits_v2 %>%
  group_by(Station_id, trait) %>%
  mutate(score = as.numeric(score)) %>%
  summarise(g_score = mean(score, na.omit = "T")) %>%
  inner_join(., infdiv_meta, by = "Station_id")

ggplot(traits_v3, aes(x = Folk_simplified, fill = Anchored, y = g_score)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(vars(as.factor(trait)))

traits_v3 = traits_v3 %>% 
  mutate(t_group = regmatches(trait, regexpr("^.+?[_ ]", trait)))

ggplot(traits_v3, aes(x = Anchored, fill = dummy, y = g_score)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(col = vars(Folk_simplified), rows = vars(as.factor(t_group)))
  

install.packages(patchwork)
library(patchwork)
         
dfs = split(traits_v3, f = traits_v3$t_group)
# apply ggplot function and write to list
gg_l = lapply(dfs, function(x) {
  ggplot(x, aes(x = interaction(Anchored, Folk_simplified), y = g_score, fill = trait), ylab = '') + 
    geom_bar(position = "dodge", stat = "identity") + 
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
          plot.margin = margin(0,0,0,0,"cm"),
          panel.background = element_rect(colour = "black", fill = NA),
          element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(y=NULL)
})
# patchwork
wrap_plots(gg_l, ncol = 1) + plot_layout(axis_titles = "collect", guides = "collect", axes = "collect") 




gg_l

ggplot(traits_v3, aes(x = Anchored, y = g_score, fill = trait), ylab = '') + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_grid(t_group ~ Folk_simplified) +
  theme(strip.background = element_blank(), strip.placement = "outside", axis.title.y = element_blank()) +
  theme_bw() 

ggplot(data = infdiv_meta, aes(x = Pressure.y, y = S, colour = Folk_simplified)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

#### Models ####

perm_inf<-adonis2(inf_div_trans ~ Anchored*Folk_simplified, data = infdiv_meta, permutations = 9999, by = "margin") # Significant but very small R2
perm_inf

hist(infdiv_meta$S)
# inf_glm_richness <- glm(formula = S ~ Pressure.y,
#                                 data = infdiv_meta,
#                        family = quasipoisson())
inf_glm_richness_nb <-MASS::glm.nb(formula = S ~ Pressure.y*Folk_simplified,  ## Decrease in coarse, despite boxplot showing potential increase in both?
                data = infdiv_meta)

summary.glm(inf_glm_richness_nb)

infdiv_meta$Folk_simplified <- droplevels(infdiv_meta$Folk_simplified)

newdata2 <- data.frame(
  Pressure.y = rep(seq(from = min(infdiv_meta$Pressure.y), to = max(infdiv_meta$Pressure.y), length.out = 100), 2),
  Folk_simplified = factor(rep(1:2, each = 100), levels = 1:2, labels =
                  levels(infdiv_meta$Folk_simplified)))

newdata2 <- cbind(newdata2, predict(inf_glm_richness_nb, newdata2, type = "link", se.fit=TRUE))
newdata2 <- within(newdata2, {
  SpeciesRichness <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata2, aes(Pressure.y, SpeciesRichness)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = Folk_simplified), alpha = .25) +
  geom_line(aes(colour = Folk_simplified), size = 2) +
  labs(x = "Math Score", y = "Predicted Days Absent")





par(mfrow = c(2,2))
plot(inf_glm_richness_nb)
plot(infdiv_meta$Folk_simplified, resid(inf_glm_richness, xlab = "Folk"))
pscl::odTest(inf_glm_richness_nb) ## Significant so use nb model

inf_glm_abundance <-glm.nb(formula = tot_abundance ~ Pressure.y*Folk_simplified, 
                       data = infdiv_meta)
summary.glm(inf_glm_abundance)
# inf_glm_abundance_poi <- glm(formula = tot_abundance ~ Pressure.y*Folk_simplified, 
#                         data = infdiv_meta,
#                         family = quasipoisson())
# summary(inf_glm_abundance_poi)
pscl::odTest(inf_glm_abundance) ## Significant so use nb model

inf_glm_diversity <- glm(formula = H ~ Pressure.y*Folk_simplified,
                         data = infdiv_meta)
summary(inf_glm_diversity)

inf_glm_simpson <- glm(formula = simp ~ Pressure.y*Folk_simplified,
                         data = infdiv_meta)
summary(inf_glm_simpson)

ggplot(data = subset(infdiv_meta, Folk_simplified == "Coarse Sediment"), aes(x = Pressure.y, y = S)) +
  geom_point() +
  geom_function(fun = function(x) exp(3.9858371 + (-0.0004274*x))) +
  geom_smooth() +
  theme_classic()

ggplot(data = subset(infdiv_meta, Folk_simplified == "Sand"), aes(x = Pressure.y, y = S)) +
  geom_point() +
  geom_function(fun = function(x) exp((3.9858371 - 0.6988028 - 0.0004274*x + 0.001623*x))) +
  geom_smooth() +
  theme_classic()
par(mfrow = c(1,1))  


# GAM who knows tbh
# press_cs <- infdiv_meta$Pressure.y[infdiv_meta$Folk_simplified == "Coarse Sediment"] 
# s_cs <- infdiv_meta$S[infdiv_meta$Folk_simplified == "Coarse Sediment"]
# plot(press_cs,s_cs)
# m1<- gam(s_cs ~ lo(press_cs, span = 0.3))
# plot(m1, se = T)
# library(mgcv)
# m2<-gam(s_cs ~ s(press_cs, fx = F, k = -1, bs = "cr"))
# plot(m2, se = F)
# m3 <- gam(infdiv_meta$S ~ s(infdiv_meta$Pressure.y) + factor(infdiv_meta$Folk_simplified))
# summary(m3)
m4 <- mgcv::gam(S ~ s(Pressure.y, by = Folk_simplified, bs = "cr") +
            Folk_simplified, family = poisson(), data = infdiv_meta)
summary(m4)
plot(m4, xlim = c(0,300), pages = 1)
infdiv_meta$Folk_simplified <- as.character(infdiv_meta$Folk_simplified)
infdiv_meta$Folk_simplified <- as.factor(infdiv_meta$Folk_simplified)

m5 <- mgcv::gam(S ~ s(Pressure.y, k = 2), family = poisson(), data = subset(infdiv_meta, Folk_simplified == "Coarse Sediment"), method = 'REML')
summary(m5)
plot.gam(m5)
gam.check(m5)
m6 <- mgcv::gam(S ~ s(Pressure.y), family = poisson(), data = subset(infdiv_meta, Folk_simplified == "Sand"))
summary(m6)
plot.gam(m6)

