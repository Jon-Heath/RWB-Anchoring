library(dplyr)
library(tidyverse)
library(ggplot2)
library(sp)
library(terra)
library(sf)
library(lwgeom)
library(units)
library(circular)
library(broom)
library(data.table)
library(raster)
library(zoo)
library(ncdf4)
library(raster)
library(lattice)
library(CFtime)
library(RColorBrewer)
library(fasterize)
library(spatialEco)
library(scales)
library(gapminder)
library(patchwork)
library(trip)
library(ggbreak)

##### Metadata #####
EquipNumb <- read.csv("EN_Table.csv")

raw_metadata <-read_csv("self_reported_vessel_length_width_final.csv",
                        col_names = TRUE,
                        col_types = cols(
                          MMSI = readr::col_factor(),
                          shipname = readr::col_factor(),
                          imo = readr::col_factor(),
                          length_m = col_double(),
                          width_m = col_double(),
                          class = readr::col_factor()
                        ))

tankers<-subset(raw_metadata, class =="tanker") %>% 
  mutate(Al = ((length_m)^2)*((-0.0001245*length_m) + 0.07626),
         C = 1.131 * width_m - 0.2479*length_m +1.33,
         Hc = 0.0318 * length_m + 3.278,
         Af = length_m^2 * (-0.00003303*length_m + 0.02094),
         LB = length_m*width_m,
         AlLB = Al/LB,
         CL = C/length_m,
         Cx = -0.922 + 0.507*AlLB + 1.162*CL,
         DWT = -4167.36 + 6.04 * (length_m*width_m) + 0.00044*(length_m*width_m)^2 - 0.00000000213 * (length_m * width_m)^3,
         disp = (DWT/83)*100,
         EN = disp^(2/3) + 2*Af + 0.1*Al,
         WSA = 9.56*DWT^0.63) ### Calculate current force acting on a ship (wetted surface area) ### Calculate current force acting on a ship (wetted surface area))

tankers<-subset(tankers, DWT > 0) ####### FIX WITH MORE ACCURATE DWT REGRESSIONS????

for (i in 1:nrow(tankers)) {
  if (tankers$EN[i] <= 240) {
    tankers$EN_type[i] = 1}
  else if (tankers$EN[i] <= 280 & tankers$EN[i] > 240) {
    tankers$EN_type[i] = 2}
  else if (tankers$EN[i] <= 320 & tankers$EN[i] > 280) {
    tankers$EN_type[i] = 3}
  else if (tankers$EN[i] <= 360 & tankers$EN[i] > 320) {
    tankers$EN_type[i] = 4}
  else if (tankers$EN[i] <= 400 & tankers$EN[i] > 360) {
    tankers$EN_type[i] = 5}
  else if (tankers$EN[i] <= 450 & tankers$EN[i] > 400) {
    tankers$EN_type[i] = 6}
  else if (tankers$EN[i] <= 500 & tankers$EN[i] > 450) {
    tankers$EN_type[i] = 7}
  else if (tankers$EN[i] <= 550 & tankers$EN[i] > 500) {
    tankers$EN_type[i] = 8}
  else if (tankers$EN[i] <= 600 & tankers$EN[i] > 550) {
    tankers$EN_type[i] = 9}
  else if (tankers$EN[i] <= 660 & tankers$EN[i] > 600) {
    tankers$EN_type[i] = 10}
  else if (tankers$EN[i] <= 720 & tankers$EN[i] > 660) {
    tankers$EN_type[i] = 11}
  else if (tankers$EN[i] <= 780 & tankers$EN[i] > 720) {
    tankers$EN_type[i] = 12}
  else if (tankers$EN[i] <= 840 & tankers$EN[i] > 780) {
    tankers$EN_type[i] = 13}
  else if (tankers$EN[i] <= 910 & tankers$EN[i] > 840) {
    tankers$EN_type[i] = 14}
  else if (tankers$EN[i] <= 980 & tankers$EN[i] > 910) {
    tankers$EN_type[i] = 15}
  else if (tankers$EN[i] <= 1060 & tankers$EN[i] > 980) {
    tankers$EN_type[i] = 16}
  else if (tankers$EN[i] <= 1140 & tankers$EN[i] > 1060) {
    tankers$EN_type[i] = 17}
  else if (tankers$EN[i] <= 1220 & tankers$EN[i] > 1140) {
    tankers$EN_type[i] = 18}
  else if (tankers$EN[i] < 1300 & tankers$EN[i] > 1220) {
    tankers$EN_type[i] = 19}
  else if (tankers$EN[i] < 1390 & tankers$EN[i] > 1300) {
    tankers$EN_type[i] = 20}
  else if (tankers$EN[i] < 1480 & tankers$EN[i] > 1390) {
    tankers$EN_type[i] = 21}
  else if (tankers$EN[i] < 1570 & tankers$EN[i] > 1480) {
    tankers$EN_type[i] = 22}
  else if (tankers$EN[i] < 1670 & tankers$EN[i] > 1570) {
    tankers$EN_type[i] = 23}
  else if (tankers$EN[i] < 1790 & tankers$EN[i] > 1670) {
    tankers$EN_type[i] = 24}
  else if (tankers$EN[i] < 1930 & tankers$EN[i] > 1790) {
    tankers$EN_type[i] = 25}
  else if (tankers$EN[i] < 2080 & tankers$EN[i] > 1930) {
    tankers$EN_type[i] = 26}
  else if (tankers$EN[i] < 2230 & tankers$EN[i] > 2080) {
    tankers$EN_type[i] = 27}
  else if (tankers$EN[i] < 2380 & tankers$EN[i] > 2230) {
    tankers$EN_type[i] = 28}
  else if (tankers$EN[i] < 2530 & tankers$EN[i] > 2380) {
    tankers$EN_type[i] = 29}
  else if (tankers$EN[i] < 2700 & tankers$EN[i] > 2530) {
    tankers$EN_type[i] = 30}
  else if (tankers$EN[i] < 2870 & tankers$EN[i] > 2700) {
    tankers$EN_type[i] = 31}
  else if (tankers$EN[i] < 3040 & tankers$EN[i] > 2870) {
    tankers$EN_type[i] = 32}
  else if (tankers$EN[i] < 3210 & tankers$EN[i] > 3040) {
    tankers$EN_type[i] = 33}
  else if (tankers$EN[i] < 3400 & tankers$EN[i] > 3210) {
    tankers$EN_type[i] = 34}
  else if (tankers$EN[i] < 3600 & tankers$EN[i] > 3400) {
    tankers$EN_type[i] = 35}
  else if (tankers$EN[i] < 3800 & tankers$EN[i] > 3600) {
    tankers$EN_type[i] = 36}
  else if (tankers$EN[i] < 4000 & tankers$EN[i] > 3800) {
    tankers$EN_type[i] = 37}
  else if (tankers$EN[i] < 4200 & tankers$EN[i] > 4000) {
    tankers$EN_type[i] = 38}
  else if (tankers$EN[i] < 4400 & tankers$EN[i] > 4200) {
    tankers$EN_type[i] = 39}
  else if (tankers$EN[i] < 4600 & tankers$EN[i] > 4400) {
    tankers$EN_type[i] = 40}
  else if (tankers$EN[i] < 4800 & tankers$EN[i] > 4600) {
    tankers$EN_type[i] = 41}
  else if (tankers$EN[i] < 5000 & tankers$EN[i] > 4800) {
    tankers$EN_type[i] = 42}
  else if (tankers$EN[i] < 5200 & tankers$EN[i] > 5000) {
    tankers$EN_type[i] = 43}
  else if (tankers$EN[i] < 5500 & tankers$EN[i] > 5200) {
    tankers$EN_type[i] = 44}
  else if (tankers$EN[i] < 5800 & tankers$EN[i] > 5500) {
    tankers$EN_type[i] = 45}
  else if (tankers$EN[i] < 6100 & tankers$EN[i] > 5800) {
    tankers$EN_type[i] = 46}
  else if (tankers$EN[i] < 6500 & tankers$EN[i] > 6100) {
    tankers$EN_type[i] = 47}
  else if (tankers$EN[i] < 6900 & tankers$EN[i] > 6500) {
    tankers$EN_type[i] = 48}
  else if (tankers$EN[i] < 7400 & tankers$EN[i] > 6900) {
    tankers$EN_type[i] = 49}
  else if (tankers$EN[i] < 7900 & tankers$EN[i] > 7400) {
    tankers$EN_type[i] = 50}
  else if (tankers$EN[i] < 8400 & tankers$EN[i] > 7900) {
    tankers$EN_type[i] = 51}
  else if (tankers$EN[i] < 8900 & tankers$EN[i] > 8400) {
    tankers$EN_type[i] = 52}
  else if (tankers$EN[i] < 9400 & tankers$EN[i] > 8900) {
    tankers$EN_type[i] = 53}
  else if (tankers$EN[i] < 10000 & tankers$EN[i] > 9400) {
    tankers$EN_type[i] = 54}
  else if (tankers$EN[i] < 10700 & tankers$EN[i] > 10000) {
    tankers$EN_type[i] = 55}
  else if (tankers$EN[i] < 11500 & tankers$EN[i] > 10700) {
    tankers$EN_type[i] = 56}
  else if (tankers$EN[i] < 12400 & tankers$EN[i] > 11500) {
    tankers$EN_type[i] = 57}
  else if (tankers$EN[i] < 13400 & tankers$EN[i] > 12400) {
    tankers$EN_type[i] = 58}
  else if (tankers$EN[i] < 14600 & tankers$EN[i] > 13400) {
    tankers$EN_type[i] = 59}
  else if (tankers$EN[i] < 16000 & tankers$EN[i] > 14600) {
    tankers$EN_type[i] = 60}
} ### Loop to add EN index to tankers dataset (can be copied and edited for cargo ships)

tankers<-tankers %>%
  mutate(Anchor_Weight = EquipNumb$Anchor.Mass..kg.[match(EN_type, EquipNumb$EN_type)],
         Chain_Weight = EquipNumb$Weight..KG.m.[match(EN_type, EquipNumb$EN_type)],
         Anchor_Weight_water = Anchor_Weight * 0.87,
         Chain_Weight_water = Chain_Weight * 0.87)

cargos <- subset(raw_metadata, class == "cargo") %>%
  mutate(Al = length_m*((0.06116*length_m) + 4.634),
         C = 0.06344*length_m - 18.26,
         Hc = 0.02892 * length_m +3.996,
         Af = length_m * (0.007522*length_m + 2.259),
         LB = length_m*width_m,
         AlLB = Al/LB,
         CL = C/length_m,
         Cx = -0.922 + 0.507*AlLB + 1.162*CL,
         DWT = 0.034 * (length_m * width_m)^1.62,
         disp = (DWT/70)*100,
         EN = disp^(2/3) + 2*Af + 0.1*Al,
         WSA = 14.24*DWT^0.596)

for (i in 1:nrow(cargos)) {
  if (cargos$EN[i] < 240) {
    cargos$EN_type[i] = 1}
  else if (cargos$EN[i] > 240 & cargos$EN[i] < 280) {
    cargos$EN_type[i] = 2}
  else if (cargos$EN[i] > 280 & cargos$EN[i] < 320) {
    cargos$EN_type[i] = 3}
  else if (cargos$EN[i] > 320 & cargos$EN[i] < 360) {
    cargos$EN_type[i] = 4}
  else if (cargos$EN[i] < 400 & cargos$EN[i] > 360) {
    cargos$EN_type[i] = 5}
  else if (cargos$EN[i] < 450 & cargos$EN[i] > 400) {
    cargos$EN_type[i] = 6}
  else if (cargos$EN[i] < 500 & cargos$EN[i] > 450) {
    cargos$EN_type[i] = 7}
  else if (cargos$EN[i] < 550 & cargos$EN[i] > 500) {
    cargos$EN_type[i] = 8}
  else if (cargos$EN[i] < 600 & cargos$EN[i] > 550) {
    cargos$EN_type[i] = 9}
  else if (cargos$EN[i] < 660 & cargos$EN[i] > 600) {
    cargos$EN_type[i] = 10}
  else if (cargos$EN[i] < 720 & cargos$EN[i] > 660) {
    cargos$EN_type[i] = 11}
  else if (cargos$EN[i] < 780 & cargos$EN[i] > 720) {
    cargos$EN_type[i] = 12}
  else if (cargos$EN[i] < 840 & cargos$EN[i] > 780) {
    cargos$EN_type[i] = 13}
  else if (cargos$EN[i] < 910 & cargos$EN[i] > 840) {
    cargos$EN_type[i] = 14}
  else if (cargos$EN[i] < 980 & cargos$EN[i] > 910) {
    cargos$EN_type[i] = 15}
  else if (cargos$EN[i] < 1060 & cargos$EN[i] > 980) {
    cargos$EN_type[i] = 16}
  else if (cargos$EN[i] < 1140 & cargos$EN[i] > 1060) {
    cargos$EN_type[i] = 17}
  else if (cargos$EN[i] < 1220 & cargos$EN[i] > 1140) {
    cargos$EN_type[i] = 18}
  else if (cargos$EN[i] < 1300 & cargos$EN[i] > 1220) {
    cargos$EN_type[i] = 19}
  else if (cargos$EN[i] < 1390 & cargos$EN[i] > 1300) {
    cargos$EN_type[i] = 20}
  else if (cargos$EN[i] < 1480 & cargos$EN[i] > 1390) {
    cargos$EN_type[i] = 21}
  else if (cargos$EN[i] < 1570 & cargos$EN[i] > 1480) {
    cargos$EN_type[i] = 22}
  else if (cargos$EN[i] < 1670 & cargos$EN[i] > 1570) {
    cargos$EN_type[i] = 23}
  else if (cargos$EN[i] < 1790 & cargos$EN[i] > 1670) {
    cargos$EN_type[i] = 24}
  else if (cargos$EN[i] < 1930 & cargos$EN[i] > 1790) {
    cargos$EN_type[i] = 25}
  else if (cargos$EN[i] < 2080 & cargos$EN[i] > 1930) {
    cargos$EN_type[i] = 26}
  else if (cargos$EN[i] < 2230 & cargos$EN[i] > 2080) {
    cargos$EN_type[i] = 27}
  else if (cargos$EN[i] < 2380 & cargos$EN[i] > 2230) {
    cargos$EN_type[i] = 28}
  else if (cargos$EN[i] < 2530 & cargos$EN[i] > 2380) {
    cargos$EN_type[i] = 29}
  else if (cargos$EN[i] < 2700 & cargos$EN[i] > 2530) {
    cargos$EN_type[i] = 30}
  else if (cargos$EN[i] < 2870 & cargos$EN[i] > 2700) {
    cargos$EN_type[i] = 31}
  else if (cargos$EN[i] < 3040 & cargos$EN[i] > 2870) {
    cargos$EN_type[i] = 32}
  else if (cargos$EN[i] < 3210 & cargos$EN[i] > 3040) {
    cargos$EN_type[i] = 33}
  else if (cargos$EN[i] < 3400 & cargos$EN[i] > 3210) {
    cargos$EN_type[i] = 34}
  else if (cargos$EN[i] < 3600 & cargos$EN[i] > 3400) {
    cargos$EN_type[i] = 35}
  else if (cargos$EN[i] < 3800 & cargos$EN[i] > 3600) {
    cargos$EN_type[i] = 36}
  else if (cargos$EN[i] < 4000 & cargos$EN[i] > 3800) {
    cargos$EN_type[i] = 37}
  else if (cargos$EN[i] < 4200 & cargos$EN[i] > 4000) {
    cargos$EN_type[i] = 38}
  else if (cargos$EN[i] < 4400 & cargos$EN[i] > 4200) {
    cargos$EN_type[i] = 39}
  else if (cargos$EN[i] < 4600 & cargos$EN[i] > 4400) {
    cargos$EN_type[i] = 40}
  else if (cargos$EN[i] < 4800 & cargos$EN[i] > 4600) {
    cargos$EN_type[i] = 41}
  else if (cargos$EN[i] < 5000 & cargos$EN[i] > 4800) {
    cargos$EN_type[i] = 42}
  else if (cargos$EN[i] < 5200 & cargos$EN[i] > 5000) {
    cargos$EN_type[i] = 43}
  else if (cargos$EN[i] < 5500 & cargos$EN[i] > 5200) {
    cargos$EN_type[i] = 44}
  else if (cargos$EN[i] < 5800 & cargos$EN[i] > 5500) {
    cargos$EN_type[i] = 45}
  else if (cargos$EN[i] < 6100 & cargos$EN[i] > 5800) {
    cargos$EN_type[i] = 46}
  else if (cargos$EN[i] < 6500 & cargos$EN[i] > 6100) {
    cargos$EN_type[i] = 47}
  else if (cargos$EN[i] < 6900 & cargos$EN[i] > 6500) {
    cargos$EN_type[i] = 48}
  else if (cargos$EN[i] < 7400 & cargos$EN[i] > 6900) {
    cargos$EN_type[i] = 49}
  else if (cargos$EN[i] < 7900 & cargos$EN[i] > 7400) {
    cargos$EN_type[i] = 50}
  else if (cargos$EN[i] < 8400 & cargos$EN[i] > 7900) {
    cargos$EN_type[i] = 51}
  else if (cargos$EN[i] < 8900 & cargos$EN[i] > 8400) {
    cargos$EN_type[i] = 52}
  else if (cargos$EN[i] < 9400 & cargos$EN[i] > 8900) {
    cargos$EN_type[i] = 53}
  else if (cargos$EN[i] < 10000 & cargos$EN[i] > 9400) {
    cargos$EN_type[i] = 54}
  else if (cargos$EN[i] < 10700 & cargos$EN[i] > 10000) {
    cargos$EN_type[i] = 55}
  else if (cargos$EN[i] < 11500 & cargos$EN[i] > 10700) {
    cargos$EN_type[i] = 56}
  else if (cargos$EN[i] < 12400 & cargos$EN[i] > 11500) {
    cargos$EN_type[i] = 57}
  else if (cargos$EN[i] < 13400 & cargos$EN[i] > 12400) {
    cargos$EN_type[i] = 58}
  else if (cargos$EN[i] < 14600 & cargos$EN[i] > 13400) {
    cargos$EN_type[i] = 59}
  else if (cargos$EN[i] < 16000 & cargos$EN[i] > 14600) {
    cargos$EN_type[i] = 60}
}

cargos<-subset(cargos, DWT > 0) ####### FIX WITH MORE ACCURATE DWT REGRESSIONS????

cargos<-cargos %>%
  mutate(Anchor_Weight = EquipNumb$Anchor.Mass..kg.[match(EN_type, EquipNumb$EN_type)],
         Chain_Weight = EquipNumb$Weight..KG.m.[match(EN_type, EquipNumb$EN_type)],
         Anchor_Weight_water = Anchor_Weight * 0.87,
         Chain_Weight_water = Chain_Weight * 0.87)

largeshipmeta<-rbind(tankers, cargos)

#### Initial filtering ####

redwharfalltotal<-read_csv("all_positions.csv")
redwharfalltotal<-subset(redwharfalltotal, select = -c(shipname, distance_from_port_m, distance_from_shore_m, vessel_id, seg_id, nnet_score, prod_shiptype, prod_geartype)) %>%
  mutate(year2 = year(redwharfalltotal$timestamp)) ## approx 15 mill
redwharfalltotal<-subset(redwharfalltotal, registry_vessel_class == "tanker" | registry_vessel_class == "cargo")
redwharfalltotal$registry_vessel_class<-as.factor(redwharfalltotal$registry_vessel_class)
redwharfalltotal<- st_as_sf(redwharfalltotal, coords = c("lon", "lat"), crs = 4326)

completeyrs<-redwharfalltotal %>%
  group_by(MMSI) %>%
  arrange(MMSI, timestamp, .by_group = TRUE) %>%##Group by ship ID
  distinct(geometry, timestamp, .keep_all = TRUE) %>%
  mutate(dist = st_distance(geometry, lag(geometry), by_element = T)) %>% ###Calculate distance between coordinates within groups
  subset(., as.integer(dist) !=0) %>%
  mutate(timediff = timestamp-lag(timestamp), ##Calculate time difference
         diff_hour = as.numeric(timediff, units = 'hours'), ## Convert into hours
         diff_hour = ifelse(diff_hour == 0, 0.0003, diff_hour), ### Change incorrect time differences of 0 into 1 second
         calc_heading = c(lwgeom::st_geod_azimuth(geometry),
                          set_units(NA, "degrees")), ### Calculate true heading
         calc_heading = set_units(calc_heading, "degrees"),
         calc_heading = drop_units(calc_heading), ## Sort out units
         calc_heading = ifelse(calc_heading < 0, 360 - abs(calc_heading), calc_heading), ## Calculate heading difference to next point
         heading_chng = abs(calc_heading-heading), #### Calculate difference between true direction and ship orientation
         #heading_chng = ifelse(heading_chng>180, 360-heading_chng, heading_chng), ### Sort out degree measurements
         type = "A", ## Set generic type to be overwritten later
         calc_heading = replace_na(calc_heading, 0)) %>% ### replace na with 0
  mutate(calc_sp = as.numeric(dist)/as.numeric(diff_hour), ## Calculate speed using calc time and distance
         sp_kn = calc_sp/1852,) %>% ### Calculate speed in knots ?Not sure why it's not already in knots?
  ungroup() %>%
  mutate(lastMMSI = dplyr::lag(MMSI))

completeyrs<-completeyrs %>% #[-c(1,887713:887714),] %>% ### redundant?
  group_by(MMSI) %>%
  mutate(calc_heading = replace_na(calc_heading, 0),
         heading_chng = replace_na(heading_chng, 9999),
         dist = replace_na(dist, set_units(9999, m)),
         timediff = replace_na(as.numeric(timediff), 0),
         diff_hour = replace_na(diff_hour, 9999),
         calc_sp = replace_na(calc_sp, 0),
         sp_kn = replace_na(sp_kn, 9999)) %>% ##### So that initial arrival points are correctly classed as motoring %>%
  ungroup() 

completeyrs<-completeyrs%>%
  mutate(type = ifelse(completeyrs$sp_kn < 1 & completeyrs$diff_hour < 1.5 & as.numeric(completeyrs$dist) < 500 & completeyrs$heading_chng > 10, ### Original 1.5
                       "stop",
                       "go")) %>% ### Determine points of suspected anchoring and motoring
  mutate(idv2=as.factor(cumsum(as.integer(dist)>2000| is.na(dist) == TRUE | diff_hour > 2 | dplyr::lag(MMSI) != MMSI | is.na(lastMMSI)) + 1),
         idv3 = as.factor(cumsum(dplyr::lag(type) != type | is.na(dplyr::lag(type))))) # %>% ### Create column of everytime behaviour changes

completeyrsstop<-subset(completeyrs, completeyrs$type=="stop") %>%
  group_by(MMSI) %>%
  mutate(dist = st_distance(geometry, lag(geometry), by_element = T), ###Calculate distance between coordinates within groups
         timediff = timestamp-lag(timestamp), ##Calculate time difference
         diff_hour = as.numeric(timediff, units = 'hours'),
         diff_hour = ifelse(diff_hour == 0, 0.0003, diff_hour)) %>%
  ungroup() %>%
  group_by(idv3) %>%
  mutate(count2=n()) %>%
  ungroup()### Dataframe of only stopped points

completeyrsstopdf <- completeyrsstop %>%
  st_coordinates() %>%
  as.data.frame() 
completeyrsstop <-cbind(completeyrsstop, completeyrsstopdf) ### Add actual coordinates as a column

#stopcompleteyrssml <- completeyrsstop[,c(1,2,5,16,21,23:26)]  ####### Dataset of just important columns atm
stopcompleteyrsavg <- completeyrsstop %>% ### Average coordinates
  group_by(idv3) %>%
  summarise(X_avg = mean(X), Y_avg = mean(Y))
stopcompleteyrsavg<-stopcompleteyrsavg %>% ## Create geometry column of average points
  st_as_sf(coords = c("X_avg", "Y_avg")) %>%
  mutate(id_lead = lead(idv3)) ### Column of next idv3 for use when calculating distances

completeyrsstop<-completeyrsstop %>%
  mutate(x_avg = stopcompleteyrsavg$X_avg[match(idv3, stopcompleteyrsavg$idv3)],
         y_avg = stopcompleteyrsavg$Y_avg[match(idv3, stopcompleteyrsavg$idv3)]) 

completeyrsstop<-st_transform(completeyrsstop,4326)
stopcompleteyrshead <- completeyrsstop %>% ## code for calculating the distances and time differences between anchoring groups to determine if they should be combined
  group_by(idv3) %>%
  filter(row_number()==1 | row_number() == n()) %>%
  mutate(merge_dist = as.numeric(st_distance(geometry, stopcompleteyrsavg$geometry[match(idv3, stopcompleteyrsavg$id_lead)], by_element = T))) %>%
  ungroup() %>%
  group_by(MMSI) %>%
  mutate(merge_time = as.numeric((timestamp-lag(timestamp)), units = 'hours')) %>%
  ungroup() %>%
  group_by(idv3) %>%
  filter(row_number()==1) %>%
  mutate(merge = ifelse(merge_dist<=500 & merge_time<1, "yes", "no")) ### Original = 3 hours

completeyrsstop <- completeyrsstop %>%
  ungroup() %>%
  mutate(merge = stopcompleteyrshead$merge[match(idv3, stopcompleteyrshead$idv3)],
         merge = as.factor(merge),
         idv3 = as.numeric(idv3),
         merge = replace_na(merge, "no")) %>%
  mutate(group2 = cumsum(merge == "no" & idv3 != lag(idv3, default = 1))) %>%
  group_by(group2) %>%##### Seems to work to combine anchoring events in the correct way  group_by(group2) %>%
  mutate(count = n())

completeyrsstop<-st_transform(completeyrsstop,32630)

my_points<-completeyrsstop %>% 
  group_by(group2) %>%
  summarise( geometry = st_combine( geometry ) ) %>%
  st_as_sf() %>%
  # circlesmin<-my_points %>%
  #   group_by(group2) %>%
  st_minimum_bounding_circle() %>%
  ungroup() %>%
  mutate(perimeter = st_perimeter(.),
         diameter = perimeter/pi) #### Calculate the maximum distance between points by using minimum bounding circle

stopcompleteyrssmllr<- completeyrsstop %>% # Calculate total time and distance between first and last point for each anchoring event
  group_by(group2) %>%
  filter(row_number()==1 | row_number() == n()) %>%
  mutate(tot_time = as.numeric((timestamp-lag(timestamp)), units = 'hours'),
         tot_dist = st_distance(geometry, lag(geometry), by_element = T),
         diameter = my_points$diameter[match(group2, my_points$group2)]) %>% #### True total distance
  subset(., !is.na(tot_time))

completeyrsstop<-completeyrsstop %>%
  group_by(group2) %>%
  mutate(tot_time = stopcompleteyrssmllr$tot_time[match(group2, stopcompleteyrssmllr$group2)],
         tot_dist = stopcompleteyrssmllr$tot_dist[match(group2, stopcompleteyrssmllr$group2)],
         diameter = stopcompleteyrssmllr$diameter[match(group2, stopcompleteyrssmllr$group2)],
         start_time = dplyr::first(timestamp),
         end_time = dplyr::last(timestamp))

completeyrsstop<-subset(completeyrsstop, count > 15 & tot_time > 1.5 & as.numeric(diameter) < 1100) ### Filter out likely anomalies or slow steaming

write.csv(completeyrs, file = "completeyrs.csv")

##### Wind NetCDF Reading in #####

gc()

#### Code from https://pjbartlein.github.io/REarthSysSci/netCDF.html ###Read
ncpath <- "~/"
ncname <- "Wind speed only"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")

ncin<- nc_open(ncfname)
print(ncin)
dname <- "wind_speed"

lon <- ncvar_get(ncin,"longitude")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin,"latitude")
nlat <- dim(lat)
head(lat)

print(c(nlon, nlat))

time <- ncvar_get(ncin,"time")
time
tunits <- ncatt_get(ncin,"time","units")
tunits
nt <- dim(time)
nt

# get depth
w_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(w_array)

# decode time
cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
cf
timestamps <- CFtimestamp(cf) # get character-string times
timestamps
time_cf <- CFparse(cf, timestamps) # parse the string into date components
time_cf
# replace netCDF fill values with NA's
w_east_array[w_array==fillvalue$value] <- NA
#length(na.omit(as.vector(w_array[,1])))

# get a single slice or layer (January)
# m <- 1
# w_slice <- w_array[,m]

# levelplot of the slice
# grid <- expand.grid(lon=lon, lat=lat)
# cutpts <- c(-5.6,-5.5,-5.4,-5.3,-5.2,-5.1)
# levelplot(w_east_slice ~ lon * lat, data=grid, cuts=11, pretty=T, 
#           col.regions=(rev(brewer.pal(10,"RdBu"))))

# reshape the array into vector for whole dataset
w_vec_long <- as.vector(w_array)
length(w_vec_long)
# reshape the vector into a matrix
w_mat <- matrix(w_vec_long, nrow=nlon*nlat, ncol=nt)
dim(w_mat)
head(na.omit(w_mat))
# create a dataframe
lonlat <- as.matrix(expand.grid(lon,lat))
windcomplete <- data.frame(cbind(lonlat,w_mat))
names(windcomplete) <- c("lon","lat",
                         "2018Jan","2018Feb","2018Mar","2018Apr","2018May","2018Jun","2018Jul","2018Aug","2018Sep","2018Oct","2018Nov","2018Dec",
                         "2019Jan","2019Feb","2019Mar","2019Apr","2019May","2019Jun","2019Jul","2019Aug","2019Sep","2019Oct","2019Nov","2019Dec",
                         "2020Jan","2020Feb","2020Mar","2020Apr","2020May","2020Jun","2020Jul","2020Aug","2020Sep","2020Oct","2020Nov","2020Dec",
                         "2021Jan","2021Feb","2021Mar","2021Apr","2021May","2021Jun","2021Jul","2021Aug","2021Sep","2021Oct","2021Nov","2021Dec",
                         "2022Jan","2022Feb","2022Mar","2022Apr","2022May","2022Jun","2022Jul","2022Aug","2022Sep","2022Oct","2022Nov","2022Dec",
                         "2023Jan","2023Feb","2023Mar","2023Apr","2023May","2023Jun","2023Jul","2023Aug","2023Sep","2023Oct","2023Nov","2023Dec",
                         "2024Jan")
# options(width=96)
#head(na.omit(tmp_df02, 20))

#wind2020<-windcomplete[,c(1:2,87:98)]
windcompleteyrsr<-rasterFromXYZ(windcomplete)

#### Wind ####

bathy <- raster("~/Merged admiralty bathy.tif") ### Bathymetry data
completeyrsstop$depth<-raster::extract(bathy, completeyrsstop)
completeyrsstop<-completeyrsstop %>%
  group_by(., group2) %>%
  mutate(depthavg = mean(depth, na.rm = T))

completeyrsstop<-st_transform(completeyrsstop, 4326)
completeyrsstop$month <-month(completeyrsstop$timestamp) # get months as number
completeyrsstop$month <-month.abb[completeyrsstop$month]
completeyrsstop$year <-year(completeyrsstop$timestamp)
completeyrsstop$match <- paste("X",completeyrsstop$year,completeyrsstop$month)
completeyrsstop$match <- as.factor(gsub(" ", "", completeyrsstop$match, fixed = TRUE))
windmatcompleteyrs<-raster::extract(windcompleteyrsr, completeyrsstop) ## Create matrix of wind at each coordinate
windmatcompleteyrs<-na.aggregate(windmatcompleteyrs) ### Convert any NA into column average value
colnames(windmatcompleteyrs)<- c(1:73)
thing<-cbind(row_number(completeyrsstop$MMSI), completeyrsstop$match) ### Create a matrix of row number and month values
completeyrsstop$wind<-windmatcompleteyrs[thing]  ### New column with wind speed by month
completeyrsstop<-completeyrsstop %>%
  group_by(., group2) %>%
  mutate(windavg = mean(wind, na.rm = T))

rm(windmatcompleteyrs)

#### Read in current NetCDF ####

ncpath <- "~/"
ncname <- "Current high res"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")

ncin<- nc_open(ncfname)
print(ncin)
dname <- "vo"

lon <- ncvar_get(ncin,"longitude")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin,"latitude")
nlat <- dim(lat)
head(lat)

print(c(nlon, nlat))

time <- ncvar_get(ncin,"time")
time
tunits <- ncatt_get(ncin,"time","units")
tunits
nt <- dim(time)
nt

# get depth
vo_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(vo_array)

# decode time
cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
cf
timestamps <- CFtimestamp(cf) # get character-string times
timestamps
time_cf <- CFparse(cf, timestamps) # parse the string into date components
time_cf
# replace netCDF fill values with NA's
vo_array[vo_array==fillvalue$value] <- NA
length(na.omit(as.vector(vo_array[,,1])))

# get a single slice or layer (January)
m <- 1
vo_slice <- vo_array[,,m]
#levelplot of the slice
grid <- expand.grid(lon=lon, lat=lat)
cutpts <- c(-5.6,-5.5,-5.4,-5.3,-5.2,-5.1)
# levelplot(si10_slice ~ lon * lat, data=grid, cuts=11, pretty=T,
#           col.regions=(rev(brewer.pal(10,"RdBu"))))
# reshape the array into vector
vo_vec_long <- as.vector(vo_array)
length(vo_vec_long)
# reshape the vector into a matrix
vo_mat <- matrix(vo_vec_long, nrow=nlon*nlat, ncol=nt)
dim(vo_mat)
# create a dataframe
lonlat <- as.matrix(expand.grid(lon,lat))
vocomplete <- data.frame(cbind(lonlat,vo_mat))
names(vocomplete) <- c("lon","lat",
                       "2018Jan","2018Feb","2018Mar","2018Apr","2018May","2018Jun","2018Jul","2018Aug","2018Sep","2018Oct","2018Nov","2018Dec",
                       "2019Jan","2019Feb","2019Mar","2019Apr","2019May","2019Jun","2019Jul","2019Aug","2019Sep","2019Oct","2019Nov","2019Dec",
                       "2020Jan","2020Feb","2020Mar","2020Apr","2020May","2020Jun","2020Jul","2020Aug","2020Sep","2020Oct","2020Nov","2020Dec",
                       "2021Jan","2021Feb","2021Mar","2021Apr","2021May","2021Jun","2021Jul","2021Aug","2021Sep","2021Oct","2021Nov","2021Dec",
                       "2022Jan","2022Feb","2022Mar","2022Apr","2022May","2022Jun","2022Jul","2022Aug","2022Sep","2022Oct","2022Nov","2022Dec",
                       "2023Jan","2023Feb","2023Mar","2023Apr","2023May","2023Jun","2023Jul","2023Aug","2023Sep","2023Oct","2023Nov","2023Dec",
                       "2024Jan")


vocompleteyrsr<-rasterFromXYZ(vocomplete)

#### uo ###
ncin<- nc_open(ncfname)
print(ncin)
dname <- "uo"

lon <- ncvar_get(ncin,"longitude")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin,"latitude")
nlat <- dim(lat)
head(lat)

print(c(nlon, nlat))

time <- ncvar_get(ncin,"time")
time
tunits <- ncatt_get(ncin,"time","units")
tunits
nt <- dim(time)
nt

# get depth
uo_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(uo_array)

# decode time
cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
cf
timestamps <- CFtimestamp(cf) # get character-string times
timestamps
time_cf <- CFparse(cf, timestamps) # parse the string into date components
time_cf
# replace netCDF fill values with NA's
uo_array[uo_array==fillvalue$value] <- NA
length(na.omit(as.vector(uo_array[,,1])))

# reshape the array into vector
uo_vec_long <- as.vector(uo_array)
length(uo_vec_long)
# reshape the vector into a matrix
uo_mat <- matrix(uo_vec_long, nrow=nlon*nlat, ncol=nt)
dim(uo_mat)
# create a dataframe
lonlat <- as.matrix(expand.grid(lon,lat))
uocomplete <- data.frame(cbind(lonlat,uo_mat))
names(uocomplete) <- c("lon","lat",
                       "2018Jan","2018Feb","2018Mar","2018Apr","2018May","2018Jun","2018Jul","2018Aug","2018Sep","2018Oct","2018Nov","2018Dec",
                       "2019Jan","2019Feb","2019Mar","2019Apr","2019May","2019Jun","2019Jul","2019Aug","2019Sep","2019Oct","2019Nov","2019Dec",
                       "2020Jan","2020Feb","2020Mar","2020Apr","2020May","2020Jun","2020Jul","2020Aug","2020Sep","2020Oct","2020Nov","2020Dec",
                       "2021Jan","2021Feb","2021Mar","2021Apr","2021May","2021Jun","2021Jul","2021Aug","2021Sep","2021Oct","2021Nov","2021Dec",
                       "2022Jan","2022Feb","2022Mar","2022Apr","2022May","2022Jun","2022Jul","2022Aug","2022Sep","2022Oct","2022Nov","2022Dec",
                       "2023Jan","2023Feb","2023Mar","2023Apr","2023May","2023Jun","2023Jul","2023Aug","2023Sep","2023Oct","2023Nov","2023Dec",
                       "2024Jan")


uocompleteyrsr<-rasterFromXYZ(uocomplete)
plot(uocompleteyrsr)
st_crs(vocompleteyrsr)
st_crs(completeyrsstop)

#### Current #####
completeyrsstop<-st_transform(completeyrsstop, 4326) ### Make sure projection is same as in rasters
vomatcompleteyrs<-raster::extract(vocompleteyrsr, completeyrsstop) ## Create matrix of current at each coordinate
colnames(vomatcompleteyrs)<- c(1:73) ### Change columns to numbers not months 
vomatcompleteyrs<-na.aggregate(vomatcompleteyrs) ### Convert any NA into column average value
vothing<-cbind(row_number(completeyrsstop$MMSI), completeyrsstop$match) ### Create a matrix of row number and month values
completeyrsstop$vo<-vomatcompleteyrs[vothing]  ### New column with wind speed by month

uomatcompleteyrs<-raster::extract(uocompleteyrsr, completeyrsstop) ## Create matrix of current at each coordinate
colnames(uomatcompleteyrs)<- c(1:73) ### Change columns to numbers not months 
uomatcompleteyrs<-na.aggregate(uomatcompleteyrs) ### Convert any NA into column average value
uothing<-cbind(row_number(completeyrsstop$MMSI), completeyrsstop$match) ### Create a matrix of row number and month values
completeyrsstop$uo<-uomatcompleteyrs[uothing]

completeyrsstop<-completeyrsstop %>%
  mutate(current_v = sqrt(uo^2 * vo^2)) %>%
  group_by(., group2) %>%
  mutate(velavg = mean(current_v, na.rm = T))

rm(uomatcompleteyrs, vomatcompleteyrs, thing)

#### Circles ####
completeyrsstop<-st_transform(completeyrsstop,32630) ## Convert to projected crs for radius calculations in m

stopdfcompleteyrs<-completeyrsstop %>%
  st_coordinates() %>%
  as.data.frame() %>%
  cbind(., as.data.frame(completeyrsstop)) ## Create dataset without geometry
groupnm<-unique(stopdfcompleteyrs$group2)

stopcentroids <- completeyrsstop %>%
  group_by(group2) %>%
  summarise( geometry = st_combine( geometry ) ) %>%
  st_as_sf() 
centroid <- st_centroid(stopcentroids) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(group2 = groupnm)

centres = NULL
names <- c("radius", "X", "Y","group2")
groupnm<-unique(stopdfcompleteyrs$group2)
names(stopdfcompleteyrs)[1] <- "X_coord"
names(stopdfcompleteyrs)[2] <- "Y_coord"

for (i in groupnm) {
  model = lsfit.circle(x = stopdfcompleteyrs$X_coord[stopdfcompleteyrs$group2 == paste(i)], 
                       y = stopdfcompleteyrs$Y_coord[stopdfcompleteyrs$group2 == paste(i)])
  centres = bind_rows(centres,model$coefficients)
  #print(model$coefficients)
  #print(i)
}
centres<-cbind(centres, groupnm)
colnames(centres) <-names
# centregeom<- st_as_sf(centres, coords = c("X", "Y"))
shrinkage <-NULL

shrinkage<-centres %>%
  mutate(Xnew = ifelse(centres$radius <600, centres$X[match(group2, centres$group2)], centroid$X[match(group2, centroid$group2)]),
         Ynew = ifelse(centres$radius <600, centres$Y[match(group2, centres$group2)], centroid$Y[match(group2, centroid$group2)])) %>%
  st_as_sf(., coords = c("Xnew", "Ynew"), crs = 32630) ### Create dataframe of points to shrink from (centroid or circle centre)

centres<-subset(centres, centres$radius<600) %>% #### Create dataset with large radii excluded (600m)
  mutate(MMSI = completeyrsstop$MMSI[match(group2, completeyrsstop$group2)],
         Name = completeyrsstop$Name[match(group2, completeyrsstop$group2)],
         group2 = group2)

centres<- st_as_sf(centres, coords = c("X", "Y"))
st_crs(centres) <-32630
#centres<-centres[,-c(5)]
centres<- centres[, c("MMSI", "group2", "geometry")] ## Rearrange columns

centres<-rbind(centres, completeyrsstop[,c("MMSI","group2","geometry")]) ### Combine filtered centrepoints with stopped dataset
centredf.2 <-centres %>%
  mutate(tot_time = stopcompleteyrssmllr$tot_time[match(group2, stopcompleteyrssmllr$group2)],
         windavg = completeyrsstop$windavg[match(group2, completeyrsstop$group2)],
         velavg = completeyrsstop$velavg[match(group2, completeyrsstop$group2)],
         depthavg = completeyrsstop$depthavg[match(group2, completeyrsstop$group2)],
         depthavg = abs(depthavg))

#### Depth and habitat####

largeshipcompleteyrs <- centredf.2 %>% #### Big df for 20completeyrs all cargo and tankers
  mutate(class = largeshipmeta$class[match(MMSI, largeshipmeta$MMSI)]) %>%
  group_by(group2) %>%
  mutate(## Metadata in Metadata worksheet
    IMO = largeshipmeta$imo[match(MMSI, largeshipmeta$MMSI)],
    length = largeshipmeta$length_m[match(MMSI, largeshipmeta$MMSI)],
    DWT = largeshipmeta$DWT[match(MMSI, largeshipmeta$MMSI)],
    breadth = largeshipmeta$width_m[match(MMSI, largeshipmeta$MMSI)], 
    Al = largeshipmeta$Al[match(MMSI, largeshipmeta$MMSI)],
    C = largeshipmeta$C[match(MMSI, largeshipmeta$MMSI)],
    Af = largeshipmeta$Af[match(MMSI, largeshipmeta$MMSI)],
    LB = largeshipmeta$LB[match(MMSI, largeshipmeta$MMSI)],
    AlLB = largeshipmeta$AlLB[match(MMSI, largeshipmeta$MMSI)],
    CL = largeshipmeta$CL[match(MMSI, largeshipmeta$MMSI)],
    Cx = largeshipmeta$Cx[match(MMSI, largeshipmeta$MMSI)],
    disp = largeshipmeta$disp[match(MMSI, largeshipmeta$MMSI)],
    EN = largeshipmeta$EN[match(MMSI, largeshipmeta$MMSI)],
    EN_type = largeshipmeta$EN_type[match(MMSI, largeshipmeta$MMSI)],
    anchor_weight = largeshipmeta$Anchor_Weight[match(MMSI, largeshipmeta$MMSI)],
    chain_weight = largeshipmeta$Chain_Weight[match(MMSI, largeshipmeta$MMSI)],
    anchor_weight_water = largeshipmeta$Anchor_Weight_water[match(MMSI, largeshipmeta$MMSI)],
    chain_weight_water = largeshipmeta$Chain_Weight_water[match(MMSI, largeshipmeta$MMSI)],
    chain_weight_water_KN = chain_weight_water * 0.009807,
    WSA = largeshipmeta$WSA[match(MMSI, largeshipmeta$MMSI)], #### Match dimension data to df
    Xw0 = abs(1/2 * 1.293 * Cx * Af * windavg^2), ## Wind load
    Xw0_KN = Xw0/1000, ## Wind load in kn
    Fw = 9.81 * 0.17 * velavg^1.83 * WSA,
    Fw_KN = Fw/1000,
    a = Xw0_KN +Fw_KN/chain_weight_water_KN, ### catenary variable calulation
    Hc = largeshipmeta$Hc[match(MMSI, largeshipmeta$MMSI)],  ### Hawsepipe height
    height = Hc+depthavg, ### Midpoint (hawsepipe) plus mean depth
    s = sqrt(height*(height + (2 * a))), ### Catenary length
    d = a*acosh((height/a) + 1), ### Horizontal catenary length
    k_tot_len = 39 * sqrt(depthavg), ### total chain length 39*sqrt(depth)
    l_on_seabed = k_tot_len - s, ### Length on seabed
    percent_seabed = l_on_seabed/k_tot_len * 100, ## mean = 83.63, median = 83.84,
    EnHr = tot_time*EN) %>%
  subset(., !is.na(depthavg) & !is.na(length))

#### Shrinking work ####

cxhullscompleteyrs <- largeshipcompleteyrs %>%
  group_by(group2) %>%
  summarise( geometry = st_combine( geometry ) ) %>%
  st_as_sf() %>%
  st_convex_hull()
st_crs(cxhullscompleteyrs)<-32630

cxhulls_scale<-extract.vertices(cxhullscompleteyrs) %>%
  cbind(., st_coordinates(.))

#cxhulls_old <- cxhulls_scale

cxhulls_scale <- cxhulls_scale %>% 
  arrange(group2) %>%
  group_by(group2) %>%
  mutate(X0 = shrinkage$X[match(group2, shrinkage$group2)],
         Y0 = shrinkage$Y[match(group2, shrinkage$group2)],
         b1 = X - X0,
         a1 = Y - Y0,
         l1 = sqrt(a1^2 + b1^2),
         l = largeshipcompleteyrs$l_on_seabed[match(group2, largeshipcompleteyrs$group2)],
         k = largeshipcompleteyrs$percent_seabed[match(group2, largeshipcompleteyrs$group2)],
         d = largeshipcompleteyrs$d[match(group2, largeshipcompleteyrs$group2)],
         tot_lngth = l+d,
         # red_length = d + largeshipcompleteyrs$length[match(group2, largeshipcompleteyrs$group2)],
         # l_cat_red = l1 - red_length,#### Maybe to include length in shrinking factor
         lngth_ratio = l/tot_lngth,
         d2 = ifelse(l1 - d < 0, l1*(l/l+d), l1-d),
         ratio = l/l1,
         l = ifelse(l>l1, l1*lngth_ratio, l),
         theta = ifelse(a1 != 0, atan(a1/b1), 0),
         a2 = (sin(abs(theta))) * (l),
         #a2 = ifelse(a1 < 0, -a2, a2),
         a2 = ifelse(a1*a2 > 0, a2, -a2),
         b2 = ifelse(theta != 0, (cos(theta)) * (l), 0),
         b2 = ifelse(b1*b2 > 0, b2, -b2),
         #b2 = ifelse(b1 < 0, -b2, b2),
         X2 = X0 + b2,
         Y2 = Y0 + a2) %>%
  st_drop_geometry() %>%
  as.data.frame(.) %>%
  st_as_sf(., coords = c("X2", "Y2"), crs = 32630)

cxhulls_scale <- cxhulls_scale %>%
  group_by(group2) %>%
  summarise( geometry = st_combine( geometry ) ) %>%
  st_as_sf() %>%
  st_convex_hull() %>% ### Looks like it works?
  mutate(EnHr = largeshipcompleteyrs$EnHr[match(group2, largeshipcompleteyrs$group2)],
         EN = largeshipcompleteyrs$EN[match(group2, largeshipcompleteyrs$group2)],
         s = largeshipcompleteyrs$s[match(group2, largeshipcompleteyrs$group2)],
         d = largeshipcompleteyrs$d[match(group2, largeshipcompleteyrs$group2)],
         percent_seabed = largeshipcompleteyrs$percent_seabed[match(group2, largeshipcompleteyrs$group2)],
         tot_time = largeshipcompleteyrs$tot_time[match(group2, largeshipcompleteyrs$group2)],
         start_time = completeyrsstop$start_time[match(group2, completeyrsstop$group2)],
         end_time = completeyrsstop$end_time[match(group2, completeyrsstop$group2)]) 
st_crs(cxhulls_scale)<-32630

cxhulls_scale <-cxhulls_scale %>%
  mutate(Depth = largeshipcompleteyrs$depthavg[match(group2, largeshipcompleteyrs$group2)],
         Year = completeyrsstop$year[match(group2, largeshipcompleteyrs$group2)],
         length = largeshipcompleteyrs$length[match(group2, largeshipcompleteyrs$group2)],
         windavg = largeshipcompleteyrs$windavg[match(group2, largeshipcompleteyrs$group2)],
         velavg = largeshipcompleteyrs$velavg[match(group2, largeshipcompleteyrs$group2)])

st_write(cxhulls_scale, "C:/Users/jnh23xdx/OneDrive - Bangor University/Work/Madog Trip/Madog Planning/Updated R outputs/Polygons/RWB Shrunk Polygons v2.shp", append = FALSE) ### Move to QGIS to shrink by %
st_write(cxhullscompleteyrs, "~/RWB Unshrunk Polygons.shp", append = FALSE) ### Move to QGIS to shrink by %

totalgrid<-read_sf(dsn = "~/Grid", layer = "Grid") ### Read in Grid of area

secty_pct<-st_intersection(totalgrid,cxhulls_scale) %>%
  mutate(intersect_area = st_area(.),
         intersect_pct = as.numeric(intersect_area/400),
         EnHr_Prop = as.numeric(intersect_pct*EnHr),
         start_time = completeyrsstop$start_time[match(group2, completeyrsstop$group2)],
         end_time = completeyrsstop$end_time[match(group2, completeyrsstop$group2)]) %>%
  dplyr::select(group2, id, intersect_area, intersect_pct, EnHr, EnHr_Prop, start_time, end_time) %>%
  st_drop_geometry()

time_intersectcompleteyrs<-secty_pct %>%
  group_by(id) %>%
  arrange(start_time, .by_group = TRUE) %>%
  mutate(time_diff = start_time - dplyr::lag(end_time),
         diff_hour = as.numeric(time_diff, units = 'hours')
  ) %>%
  subset(., !is.na(time_diff)) %>%
  ungroup(.) %>%
  mutate(time_se = sqrt(sum((diff_hour-mean(diff_hour))^2/(length(diff_hour)-1)))/sqrt(length(diff_hour)),
         avg_time = mean(diff_hour))

sumintersect20completeyrs<-secty_pct %>%
  group_by(id) %>%
  arrange(start_time, .by_group = TRUE) %>%
  mutate(time_diff = start_time - dplyr::lag(end_time)) %>%
  summarise(., A_EnHr = sum(EnHr_Prop), number = n(), avg_time = mean(time_diff, na.rm =TRUE)) %>%
  ungroup() %>%
  mutate(normal_A_EnHr = ((A_EnHr - min(A_EnHr))/(max(A_EnHr) - min(A_EnHr))),
         standard_A_AOIWDA = (A_EnHr - mean(A_EnHr))/sd(A_EnHr),
         avg_time = as.numeric(avg_time, units = 'hours'))

#totalgrid<-merge(totalgrid[, -c(2:5)], secty_pct, by = "id", all.x = TRUE)
totalgrid<-merge(totalgrid, sumintersect20completeyrs, by = "id", all.x = TRUE)

## Convert to raster ##

gc()

ext<-extent(totalgrid)
rr<-raster(ext, res = 20)
anchcompleteyrs<-fasterize(totalgrid, rr, field = sqrt(totalgrid$A_EnHr))
anchcompleteyrssqrt<-fasterize(totalgrid, rr, field = sqrt(totalgrid$A_EnHr))
anchcompleteyrs_n<-fasterize(totalgrid, rr, field = totalgrid$number)
anchcompleteyrs_time<-fasterize(totalgrid, rr, field = totalgrid$avg_time)


