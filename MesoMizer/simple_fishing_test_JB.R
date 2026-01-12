##### Effects of different types of fishing on DIC pools sequestered for 100 years


## read in params object
params_v1<-readRDS("~/Library/CloudStorage/OneDrive-UniversityofTasmania/MesoMizer/MesoMizer/params_v1.rds") 

#check yield curves
plotYieldVsF(params_v1,species="Nannobrachium achirus",F_max=0.1)

## test 1. add fishing to mesopelagic fish (mostly lower mesopelagic)
params_fished1<-params_v1
nspp<-dim(params_v1@given_species_params)[1]
initial_effort(params_fished1)<-1 # using catchabillty to set fishing rather tahn effort
gear_params(params_fished1)$catchability<-rep(0,length=nspp)
gear_params(params_fished1)[c("Nannobrachium achirus, knife_edge_gear","Cyclothone spp, knife_edge_gear", "Bathylagus spp, knife_edge_gear"),"catchability"]<-0.2
params_fished1<-projectToSteady(params_fished1,t_max = 500)
plot(params_fished1)

## test 2. add fishing to euphasiids (mid mesopelagic)
params_fished2<-params_v1
initial_effort(params_fished2)<-0.2
gear_params(params_fished2)$catchability<-rep(0,length=nspp)
gear_params(params_fished2)[c("Euphausia similis, knife_edge_gear","Thysannoessa macura, knife_edge_gear"),"catchability"]<-0.2
params_fished2<-projectToSteady(params_fished2,t_max = 500)
plot(params_fished2)

## Plots of size spectra compared to unfished

plotlySpectraRelative(params_v1,params_fished1)

plotlySpectraRelative(params_v1,params_fished2)

## How do I remake plots of DIC pools? 

getBiomass(params_fished1)
getBiomass(params_fished2)


# Barplot

library(ggbreak)

data1<-data.frame(species=names(getBiomass(params_fished1)),perc_biomass_diff=(getBiomass(params_fished1)-getBiomass(params_v1))/getBiomass(params_v1)*100)
data1$fishing_scenario<-rep("mesopelagics")
data2<-data.frame(species=names(getBiomass(params_fished2)),perc_biomass_diff=(getBiomass(params_fished2)-getBiomass(params_v1))/getBiomass(params_v1)*100)
data2$fishing_scenario<-rep("krill")

data<-rbind(data1,data2)

# Barplot

# Create the side-by-side barplot
ggplot(data, aes(x = species, y = perc_biomass_diff, fill = fishing_scenario)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() + coord_flip() +  scale_y_break(c(100, 300)) +  scale_y_break(c(360, 1000)) 

# Break between 100 and 300

ggplot(data2, aes(x = species, y = biomass_diff, fill = fishing_scenario)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() + coord_flip()

ggplot(data1, aes(x=species, y=biomass_diff)) + 
  geom_bar(stat = "identity") + coord_flip()

ggplot(data2, aes(x=species, y=biomass_diff)) + 
  geom_bar(stat = "identity") + coord_flip()

barplot(getBiomass(params_fished2)-getBiomass(params_v1))

## save params & do carbion calculations
saveRDS(params_fished1,"~/Library/CloudStorage/OneDrive-UniversityofTasmania/MesoMizer/MesoMizer/params_fished1.rds")
saveRDS(params_fished2,"~/Library/CloudStorage/OneDrive-UniversityofTasmania/MesoMizer/MesoMizer/params_fished2.rds")
