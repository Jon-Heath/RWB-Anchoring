tallied_inf_Sand <- subset(inf_3, Folk_simplified == "Sand" & count >0) %>%
  group_by(ScientificName_accepted) %>%
  summarize(frequency = n()) %>% 
  ungroup()
tally_inf_sand = slice_max(tallied_inf_Sand, n=9, order_by = frequency)
common_inf_df_sand = inf_3[inf_3$ScientificName_accepted %in% tally_inf_sand$ScientificName_accepted, ]

ggplot(subset(common_inf_df_sand, Folk_simplified == "Sand"), aes(y = log(count+1), x = Anchored)) +
  geom_boxplot() +
  facet_wrap(~ScientificName_accepted,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic()

ggplot(subset(common_inf_df_sand, Folk_simplified == "Sand"), aes(y = log(count+1), x = Pressure)) +
  geom_point() +
  facet_wrap(~ScientificName_accepted,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  geom_smooth(method = "lm") +
  theme_classic()

### Coarse ###

tallied_inf_coarse <- subset(inf_3, Folk_simplified == "Coarse Sediment" & count >0) %>%
  group_by(ScientificName_accepted) %>%
  summarize(frequency = n()) %>% 
  ungroup()
tally_inf_coarse = slice_max(tallied_inf_coarse, n=12, order_by = frequency)
common_inf_df_coarse = inf_3[inf_3$ScientificName_accepted %in% tally_inf_coarse$ScientificName_accepted, ]

ggplot(subset(common_inf_df_coarse, Folk_simplified == "Coarse Sediment"), aes(y = log(count +1), x = Anchored)) +
  geom_boxplot() +
  facet_wrap(~ScientificName_accepted,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic()

ggplot(subset(common_inf_df_coarse, Folk_simplified == "Coarse Sediment"), aes(y = count, x = Pressure)) +
  geom_point() +
  facet_wrap(~ScientificName_accepted,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic() +
  geom_smooth(method = "lm")

## Combined ##

tally_inf_common = slice_max(na.omit(tallied), n=16, order_by = frequency)
common_inf_df = inf_3[inf_3$ScientificName_accepted %in% tally_inf_common$ScientificName_accepted, ]

ggplot(common_inf_df, aes(y = log(count+1), x = Folk_simplified, fill = Anchored)) +
  geom_boxplot() +
  facet_wrap(~ScientificName_accepted,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic()

ggplot(common_inf_df, aes(y = log(count+1), x = Pressure)) +
  geom_point() +
  facet_wrap(~ScientificName_accepted,
             axes = "all",
             axis.labels = "all_y",
             scales = "free_y") +
  theme_classic() +
  geom_smooth(method = "lm")
