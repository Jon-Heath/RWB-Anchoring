library(gapminder)

# Species stuff excluding extra trawls #

tallied_epi <- ben_mass %>% ## Tally of no. of trawls containing each species
  group_by(Species) %>%
  summarize(frequency = n()) %>% 
  ungroup()

tallied_epi_mass <- ben_mass %>% ## Tally of no. of trawls containing each species
  group_by(Species) %>%
  summarize(frequency = sum(biomass)) %>% 
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

tally_epi_common = slice_max(tallied_epi, n=16, order_by = frequency)
common_epi_df = epi_species[epi_species$Species %in% tally_epi_common$Species, ]

ggplot(subset(common_epi_df, Folk_simplified != "Mud to Muddy Sand"), aes(y = log(biomass+1), x = Folk_simplified, fill = Anchored)) +
  geom_boxplot() +
  facet_wrap(~Species,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic() +
  labs(fill = "Sediment Type")

ggplot(subset(common_epi_df, Folk_simplified != "Mud to Muddy Sand"), aes(y = log(biomass+1), x = log(Avg_pressure.x+1))) +
  geom_point() +
  facet_wrap(~Species,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic() +
  geom_smooth(method = "lm")

## Just Sand ##

tallied_epi_Sand <- subset(epi_species, Folk_simplified == "Sand" & count >0) %>%
  group_by(Species) %>%
  summarize(frequency = n()) %>% 
  ungroup()
tally_epi_sand = slice_max(tallied_epi_Sand, n=16, order_by = frequency)
common_epi_df_sand = epi_species[epi_species$Species %in% tally_epi_sand$Species, ]

ggplot(subset(common_epi_df_sand, Folk_simplified == "Sand"), aes(y = log(biomass +1), x = Anchored)) +
  geom_boxplot() +
  facet_wrap(~Species,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic()

ggplot(subset(common_epi_df_sand, Folk_simplified == "Sand"), aes(y = log(biomass +1), x = Avg_pressure.x)) +
  geom_point() +
  facet_wrap(~Species, scale="free_y") +
  geom_smooth(method = "lm") +
  theme_classic()

### Coarse ###

tallied_epi_coarse <- subset(epi_species, Folk_simplified == "Coarse Sediment" & count >0) %>%
  group_by(Species) %>%
  summarize(frequency = n()) %>% 
  ungroup()
tally_epi_coarse = slice_max(tallied_epi_coarse, n=16, order_by = frequency)
common_epi_df_coarse = epi_species[epi_species$Species %in% tally_epi_coarse$Species, ]

ggplot(subset(common_epi_df_coarse, Folk_simplified == "Coarse Sediment"), aes(y = log(biomass +1), x = Anchored)) +
  geom_boxplot() +
  facet_wrap(~Species,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic()

ggplot(subset(common_epi_df_coarse, Folk_simplified == "Coarse Sediment"), aes(y = log(biomass +1), x = sqrt(Avg_pressure.x))) +
  geom_point() +
  facet_wrap(~Species,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic() +
  geom_smooth(method = "lm")

coarse.ep.glm<- glm(log(biomass+1) ~ Avg_pressure.x*Species, data = common_epi_df)
summary(coarse.ep.glm)
anova(coarse.ep.glm, test = "F")


p.ber.coarseglm <- glm(log(biomass+1) ~ Avg_pressure.x, data = subset(common_epi_df, Species == "Alcyonidium diaphanum" & Folk_simplified != "Coarse Sediment")) ### GLM model for each species?
summary(p.ber.coarseglm)

ggplot(subset(epi_species, Species == "Psammechinus miliaris" & Folk_simplified != "Mud to Muddy Sand"), aes(x = log(Avg_pressure.x+1), y = log(biomass+1), colour = Folk_simplified)) +
  geom_point() +
  geom_smooth(method = "lm") +
  #geom_function(fun = function(x) exp(2.9378) * x^0.1565) +  
  theme_classic() #+
  # geom_line(data = data.frame(x = p.ber.coarseglm$model$`log(Avg_pressure.x + 1)`,
  #                             y = predict(p.ber.coarseglm)),
  #           aes(x = x, y = y))

plot(p.ber.coarseglm)
hist(log(epi_species$biomass[epi_species$Species == "Psammechinus miliaris"]))
ggplot(subset(epi_species, Species == "Pagurus bernhardus" & Folk_simplified != "Mud to Muddy Sand"), aes(x = Folk_simplified, y = log(biomass+1), fill = Anchored)) +
  geom_boxplot() +
  theme_classic()


## Model for each species hopefully

pressure_model <- function(df) {
  glm(log(biomass+1) ~ Avg_pressure.x, data = df)
}
options(scipen=999)

common_epi_df_2 <- subset(common_epi_df_coarse) %>%
  nest(data = -Species) %>%
  mutate(model = map(data, pressure_model),
         fit = map(data, ~glm(log(biomass+1) ~ Avg_pressure.x, data = .x)),
         tidied = map(fit, broom::tidy)) %>%
  unnest(tidied) %>%
  mutate(p.value = round(p.value, digits = 3))  %>% 
  mutate(glance = map(model, broom::glance)) %>% 
  unnest(glance)
common_epi_df_2



common_epi_df_2 <- common_epi_df_2 
resids <- unnest(common_epi_df_2, resids)
resids
resids %>% 
  ggplot(aes(Avg_pressure.x, resid)) +
  geom_line(aes(colour = Species), alpha = 1 / 3) + 
  geom_smooth(se = FALSE)

glance <- common_epi_df_2 %>% 
  mutate(glance = map(model, broom::glance)) %>% 
  unnest(glance, .drop = TRUE)
glance
glance %>%
  arrange(r.squared)

## Including extra trawls ##

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

tally_epi_common = slice_max(tallied_epi, n=9, order_by = frequency)
common_epi_df_all = epi_species_all[epi_species_all$Species %in% tally_epi_common$Species, ]

ggplot(common_epi_df_all, aes(y = log(biomass+1), x = Anchored)) +
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

#### NMDS Stuff ####

comb_nmds <- cbind(epi_trans, inf_div_trans)



