## Needs code from Neater Beam Trawl First

trawl_meta <-read_csv("C:/Users/jnh23xdx/OneDrive - Bangor University/Work/Madog Trip/Madog Results/Beam Trawl/Trawl_id_pressures.csv")
trawl_meta <- trawl_meta %>%
  rename(trawl_id = Trawl_id) %>%
  mutate(trawl_id = as.factor(trawl_id))

trawl_all_meta<-inner_join(diversities, trawl_meta, by = "trawl_id") %>% 
  inner_join(., Ben_mass_group, by = "trawl_id") %>%
  mutate(Anchored = ifelse(Avg_pressure < 1, "no", "yes")) %>% ## Create df with meta
  mutate(Anchored = as.factor(Anchored))

##Biomass plots
ggplot(data = trawl_all_meta, aes(x = Anchored, y = log(totbiomass))) + 
  geom_boxplot()

ggplot(data = trawl_all_meta, aes(x = Avg_pressure, y = log(totbiomass))) +
  geom_point() +
  geom_smooth(method = "lm")

##Abundance plots
ggplot(data = trawl_all_meta, aes(x = Anchored, y = log(tot.abund))) +
  geom_boxplot()

ggplot(data = trawl_all_meta, aes(x = Avg_pressure, y = log(tot.abund))) +
  geom_point() +
  geom_smooth(method = "lm")
 

##Indices plots
ggplot(data = trawl_all_meta, aes(x = Anchored, y = S)) +
  geom_boxplot()

ggplot(data = trawl_all_meta, aes(x = Avg_pressure, y = S)) +
  geom_point() +
  geom_smooth(method = "lm")

epi_species_all <- benthos %>%
  mutate(trawl_id = as.factor(trawl_id),
         Species = as.factor(Species)) %>%
  group_by(trawl_id, Species, .drop = FALSE) %>% #Station or trawl_id ()
  summarise(count = sum(as.numeric(Number)),
            biomass = sum(as.numeric(Weight)), 
            .groups = "drop") %>% ###
  inner_join(trawl_meta, by = "trawl_id") %>%
  mutate(Anchored = ifelse(Avg_pressure == 0, "no", "yes"),
         count = as.numeric(count))

ggplot(subset(epi_species_all, Species == "Pagurus bernhardus"), aes(x = Anchored, y = count)) +
  geom_boxplot() +
  theme_classic()

ggplot(subset(epi_species_all, Species == "Pagurus bernhardus"), aes(x = Avg_pressure, y = count)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() 


tally_epi_common = slice_max(tallied_epi, n=9, order_by = frequency)
common_epi_df_all = epi_species_all[epi_species_all$Species %in% tally_epi_common$Species, ]

ggplot(common_epi_df_all, aes(y = log(biomass+1), x = Anchored, fill = Folk_simplified)) +
  geom_boxplot() +
  facet_wrap(~Species,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic() 

ggplot(common_epi_df_all, aes(y = biomass, x = Avg_pressure)) +
  geom_point() +
  facet_wrap(~Species, scale="free_y") +
  geom_smooth(method = "lm")






trawl_all_meta$tot.abund = round(trawl_all_meta$tot.abund, digits = 0)
beam_glm_all <- glm.nb(S ~ Avg_pressure, data = trawl_all_meta)
summary(beam_glm_all)

beam_glm_mass_all <-glm(totbiomass ~ Avg_pressure, data = trawl_all_meta, family = Gamma())
summary(beam_glm_mass_all)

