
library(mizer)
library(tidyverse)
library(patchwork)


rm(list =ls()) # clear memory

#Biomass size data cleaning
bio_size <- read.csv("SOLACE_biomass_filtered.csv") 
bio_mean_wc <- read.csv("Biomass_av_wc.csv")
w_subsample <- read.csv("W_subsample.csv")

bio_size <- bio_size %>%
  mutate(Taxa_updated = if_else(Taxa_updated == "Thysannoessa macura", "Thysanoessa macrura", Taxa_updated))

w_subsample <- w_subsample %>%
  mutate(Taxa_updated = if_else(Taxa_updated == "Thysannoessa macura", "Thysanoessa macrura", Taxa_updated))

bio_mean_wc <- bio_mean_wc %>%
  mutate(Taxa_updated = if_else(Taxa_updated == "Thysannoessa macura", "Thysanoessa macrura", Taxa_updated))


w_subsample <- w_subsample %>% 
  # Rename Species values
  mutate(
    Taxa_updated = case_when(
    Taxa_updated == "Bathylagus spp" ~ "Bathylagus spp.",
    Taxa_updated == "Sergestes spp" ~ "Sergestes spp.",
    Taxa_updated == "Cyclothone spp" ~ "Cyclothone spp.",
    TRUE ~ Taxa_updated))



max_weight <- w_subsample %>%
  group_by(Taxa_updated) %>%
  summarise(max_g = max(weight_g, na.rm = TRUE))

write.csv(max_weight, "max_weight.csv")


w_subsample <- left_join(w_subsample, bio_mean_wc, by = "Taxa_updated")



bio_size <- w_subsample %>%
  group_by(Taxa_updated, weight_total_gm2_sp) %>%
  # Count frequency of each size
  count(weight_g, name = "freq") %>%
  # Calculate proportion of each size within species
  mutate(
    prop = freq / sum(freq),
    biomass_bin = prop * weight_total_gm2_sp
  ) %>%
  ungroup()



# 
# 
# 
# bio_size <- left_join(bio_size, bio_mean_wc, by = "Taxa_updated")

bio_size$species<-as.factor(bio_size$Taxa_updated)
levels(bio_size$species)
# 
# bio_size <- bio_size %>%
#   mutate(
#     Abd_Std_nos_m2 = case_when(
#       Depth_Stratum == "E"  ~ Abd_Std_nosm3_sp * 200,
#       Depth_Stratum == "UM" ~ Abd_Std_nosm3_sp * 200,
#       Depth_Stratum == "LM" ~ Abd_Std_nosm3_sp * 600,
#       TRUE ~ NA_real_  # optional fallback for unexpected values
#     )
#   )
# 
# 
# bio_size <- bio_size %>%
#   mutate(
#     Weight_Std_gm2 = case_when(
#       Depth_Stratum == "E"  ~ Weight_Std_gm3_sp * 200,
#       Depth_Stratum == "UM" ~ Weight_Std_gm3_sp * 200,
#       Depth_Stratum == "LM" ~ Weight_Std_gm3_sp * 600,
#       TRUE ~ NA_real_  # optional fallback for unexpected values
#     )
#   )


# # get average biomass across time. & space, then sum into weightclass bins
# biomass_mean<-bio_size %>%
#   group_by(Site, daylight, Depth_Stratum, Deployment, species) %>% 
#   summarise(StdWt = sum(Weight_Std_gm2)) %>% #one weight per species per deployment
#   ungroup() %>%
#   group_by(daylight, species, Depth_Stratum) %>%  
#   summarise(weight_mean_depth = mean(StdWt)) %>% #one weight per species per day/night & depth (average deployments)
#   ungroup() %>% 
#   group_by(daylight, species) %>%
#   summarise(
#     Weight_wc_gm2_sp = sum(weight_mean_depth, na.rm = TRUE) #one weight per species per whole water column (add each layer)
#   ) %>% 
#   ungroup() %>% 
#   group_by(species) %>% 
#   summarise(
#     weight_mean_gm2 = mean(Weight_wc_gm2_sp) #average day and night values
#   )
# 
# 
# bio_size <- left_join(bio_size, biomass_mean, by = "species")
# 







# read params

params_v1<-readRDS("params_v1.rds") 


species_params(params_v1)$species


# create weight bins from estimated  individual weights
bio_size$weightclass<- as.factor(cut(bio_size$weight_g,breaks=10^seq(-1.9,2.7,1)))



#check range of bins - looks OK
# 
# upp<-10^seq(-1.8,2.7,1)
# low<-10^seq(-1.9,2.6,1)
upp<-10^seq(-1.8,2.7,1.5)
low<-10^seq(-1.9,2.6,1.5)
widths<-upp-low
mid<-(upp + low)/2 # midpoint of weightclass

bio_size$upper<-bio_size$weightclass
levels(bio_size$upper)<-upp
bio_size$upper<-as.numeric(as.character(bio_size$upper))
bio_size$lower<-bio_size$weightclass
levels(bio_size$lower)<-low
bio_size$lower<-as.numeric(as.character(bio_size$lower))
bio_size$widths<-bio_size$upper-bio_size$lower

bio_size <- bio_size |> mutate(mid = exp((log(upper) + log(lower)) / 2))

# get average biomass across time. & space, then sum into weightclass bins
biomass_size_spectra<-bio_size %>%
  # group_by(Site, daylight, Depth_Stratum, Deployment,species,weightclass, mid,widths)%>%
  # summarise(total_biomass=sum(weight_mean_gm2))%>% # mean water column biomass in g m-2 - used to calibrate model
  # ungroup() %>%
  group_by(species,weightclass,mid,widths)%>%
  summarise(mean_total_biomass=mean(biomass_bin))%>%
  rename(Species=species)%>%
  ungroup()
  

biomass_size_spectra <- biomass_size_spectra %>%
  mutate(Species = if_else(Species == "Thysannoessa macura", "Thysanoessa macrura", Species))




# include bin widths in data
#biomass_size_spectra$widths<-diff(10^seq(-1.9,2.7,0.1))


#normalise biomass y-axis
biomass_size_spectra$norm_biomass<-biomass_size_spectra$mean_total_biomass/as.numeric(as.character(biomass_size_spectra$widths))

biomass_size_spectra$mid<-as.numeric(as.character(biomass_size_spectra$mid))

levels(biomass_size_spectra$Species)<-c("Atolla wyvillei","Bathylagus spp","Chaetognath","Cranchidae",          
"Cyclothone spp","Euphausia similis","Gymnosomata","Lampanyctus australis","Nannobrachium achirus", "Pyrosoma atlanticum","Sergestes spp","Themisto spp",         
 "Thysannoessa macura")

plotSpectra(params_v1,power=0,resource = F)  + geom_point(data=biomass_size_spectra,aes(x=mid,y=norm_biomass,colour = Species)) + facet_wrap(~Species) + theme_minimal()

mean_biomass_size_spectra <- ggplot()+ geom_point(data=biomass_size_spectra,aes(x=mid,y=norm_biomass,colour = Species)) + facet_wrap(~Species) + scale_y_log10() + scale_x_log10()



mod_spectra <- plotSpectra(params_v1,power=0,resource = F) + facet_wrap(~Species) + theme_minimal()



#percentage depth plot

percentage_depth <- read.csv("Percentage_depth.csv") #from Biomass_SOLACE.rmd

percentage_depth <- percentage_depth %>%
  mutate(Taxa_updated = if_else(Taxa_updated == "Thysannoessa macura", "Thysanoessa macrura", Taxa_updated))


# Make sure Depth_Stratum is ordered correctly (Epi > Umeso > Lmeso)
Percentage_depth <- percentage_depth %>%
  mutate(Depth_Stratum = factor(Depth_Stratum, levels = c("E", "UM", "LM")))

Percentage_depth$Species<-as.factor(Percentage_depth$Taxa_updated)

# Plot
ggplot(Percentage_depth, aes(x = percent_depth, y = Depth_Stratum)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  facet_wrap(~ Taxa_updated) +
  scale_y_discrete(limits = rev(levels(Percentage_depth$Depth_Stratum))) +  # Flip y-axis order
  labs(x = "Percent of Depth-Integrated Biomass", y = "Depth Stratum") +
  theme_bw() +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))



#  Clean names in Percentage_depth ---
Percentage_depth <- Percentage_depth %>%
  mutate(
    Species = case_when(
      Species == "Bathylagus spp." ~ "Bathylagus spp",
      Species == "Cyclothone spp." ~ "Cyclothone spp",
      Species == "Sergestes spp." ~ "Sergestes spp",
      TRUE ~ Species
    ),
    Depth_Stratum = factor(Depth_Stratum, levels = c("E", "UM", "LM")),
    Depth_Stratum = as.character(Depth_Stratum)  # Avoid list-columns
  )

mod_spectra <- plotSpectra(params_v1, power = 0, resource = FALSE, return_data = TRUE) 


#  Loop to create individual plots per species ---
# Identify species that exist in all 3 datasets
species_list <- intersect(
  intersect(unique(biomass_size_spectra$Species),
            unique(mod_spectra$Species)),
  unique(Percentage_depth$Species)
)

combined_plots <- list()

for (sp in species_list) {
  # Filter species data
  obs_sp <- filter(biomass_size_spectra, Species == sp)
  mod_sp <- filter(mod_spectra, Species == sp)
  depth_sp_data <- filter(Percentage_depth, Species == sp)
  
  if (nrow(depth_sp_data) == 0) next  # Skip if no depth data
  
  # Ensure all 3 depth layers are present (for consistent inset layout)
  depth_levels <- c("E", "UM", "LM")
  depth_sp_data <- depth_sp_data %>%
    complete(Depth_Stratum = depth_levels, fill = list(percent_depth = 0))
  
  x_limits <- range(c(mod_spectra$w, biomass_size_spectra$mid), na.rm = TRUE)
  y_limits <- range(c(mod_spectra$value, biomass_size_spectra$norm_biomass), na.rm = TRUE)
  
  
  # Create main spectra plot (modeled + observed)
  main_sp <- ggplot() +
    geom_line(data = mod_sp, aes(x = w, y = value), color = "black") +
    geom_point(data = obs_sp, aes(x=mid,y=norm_biomass), color = "darkblue") +
    scale_x_log10(limits = x_limits) +
    scale_y_log10(limits = y_limits) +
    theme_bw() +
    labs(x = NULL, y = NULL) +
    ggtitle(sp) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),  # Title font
      axis.text = element_text(size = 10),
      axis.title = element_blank()
    )
  
  # Create inset bar plot of depth contribution
  inset_sp <- ggplot(depth_sp_data,
                     aes(x = percent_depth,
                         y = factor(Depth_Stratum, levels = rev(depth_levels)),
                         fill = Depth_Stratum)) +
    geom_bar(stat = "identity", color = "black") +  # <- black outline around bars
    theme_void() +
    scale_fill_manual(values = c(
      "E" = "#a6cee3",    # light blue
      "UM" = "#1f78b4",  # medium blue
      "LM" = "#08306b"   # dark blue
    )) +
    theme(
      legend.position = "none"
    )
  
  combined <- cowplot::ggdraw(main_sp) +
    cowplot::draw_plot(inset_sp, x = 0.20, y = 0.15, width = 0.2, height = 0.25)
  
  
  combined_plots[[sp]] <- combined
}

#  Assemble all plots into final figure ---
final_plot <- wrap_plots(combined_plots, ncol = 4)

# Display the final combined plot
print(final_plot)


ggsave("Combined_Spectra_with_DepthInset.png", final_plot, width = 12, height = 8, dpi = 300)





