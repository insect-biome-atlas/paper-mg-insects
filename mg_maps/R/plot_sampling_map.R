library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(MetBrewer)
library(terra)
library(tidyterra)
library(geodata)
library(scales)
library(ggrepel)
library(gridExtra)

# Sweden samples ----------------------------------------------------------

meta_mg <- read_tsv("iba_data/sites_metadata_MG.tsv")

# Encoding on metadata files is broken so names need to be cleaned manually.
# Strings to remove
rm_str <- paste(
          c("RŽserve" ,"de" , "Parc National" , 
            "SpŽciale" , "Ressources Naturelles",
            "la For\u0090t Naturelle", "Concession.*?(?=Kirindy)" ,
            "PaysageHarmonieuxduComplexeZonesHumis", "d'", "PaysageHarmonieuxProtŽgŽLoky",
            ".*(?=Betampona)", ".*(?=Ambatovy)" , ".*(?=Anjozorobe)" , "Mite", "Cap", 
            "\\(CNFEREF\\)" , ".*(?=Mahavavy)", ".*(?=Loky)", ".*(?=Oronjia)" , " " 
            ),
            collapse = "|")

meta_mg <- meta_mg %>% 
        mutate(parkID = str_remove_all(parkID , rm_str) , fixed = TRUE,
               parkNumber = as.integer(factor(parkID)))

# same but for madagascar -------------------------------------------------

# Get elevation, slope, and aspect data
mad_elev        <- elevation_30s("MDG" , path = ".") 
mad_slope       <- terrain(mad_elev, "slope", unit = "radians")
mad_aspect      <- terrain(mad_elev, "aspect", unit = "radians")
pal_greys <- hcl.colors(1000, "Grays")


# Decide on which slopes to shade
mad_hshade        <- shade(mad_slope, mad_aspect, 30, 270)
names(mad_hshade) <- "shades" # add names to index variable

# Get an index of which values to shade and to what degree
index <- mad_hshade %>%
          mutate(index_col = scales::rescale(shades, to = c(1, length(pal_greys)))) %>% 
          mutate(index_col = round(index_col)) %>%
          pull(index_col)

# Get a palette of grey colours to act as the "shade"
shade_cols <- pal_greys[index]

# Make the plot
mad_hshade_plot <- ggplot() +
  geom_spatraster(data = mad_hshade, fill = shade_cols, maxcell = Inf,alpha = 1)   # Avoid resampling with maxcell

# Scale min max values of elevation raster
elev_limits <- minmax(mad_elev) %>% as.vector() # Get min max
elev_limits <- c(floor(elev_limits[1] / 500), ceiling(elev_limits[2] / 500)) * 500 # Rounded to lower and upper 500
elev_limits <- pmax(elev_limits, 0) # Set min to 0


# save the hillshaded map to use in other scripts
saveRDS(mad_hshade_plot  , "mg_maps/R/mdg_hs_map.rds")

# plots ---------------------------------------------------------------------------------------

# Format meta-data
meta_mg <- meta_mg  |>
            mutate(trap_habitat=recode(trap_habitat , 
                                       "Dry_Forest"          = "Dry Forest",
                                        "Montane_Rainforest" = "Montane Forest",
                                        "Rainforest"         = "Wet Forest"),
                              malaise_trap_type = recode(malaise_trap_type , "Single_trap" = "Single trap"))


p1 <- mad_hshade_plot + # PLot hillshaded map
  geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
  ggspatial::annotation_scale(location = 'tl',width_hint = .6,text_cex = 2)+
  scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.4,direction = 1)+
  new_scale_fill()+
  geom_point(data=meta_mg , aes(longitude_WGS84 , latitude_WGS84 , 
                                fill = trap_habitat, shape = trap_habitat),size=5)+ 
  scale_fill_viridis_d(option="mako", end=.8)+
  scale_shape_manual(values = c(24, 21, 25))+
  theme_linedraw(base_size = 25)+
  labs(x = "Longitude" ,
       y = "Latitude" , 
       fill = "Elevation",
       shape = "Trap-type",
       colour = "Trap habitat")+
  theme(legend.key = element_rect(fill = "white"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(1.5, "cm"))+
  guides(fill = guide_legend(override.aes = list(size = 0)),
          color=guide_legend(override.aes=list(fill=NA,size=6)),
         shape=guide_legend(override.aes=list(fill=NA,size=6)))


meta_mg_us <- meta_mg %>% 
              select(parkID , siteID , longitude_WGS84 , latitude_WGS84,trap_habitat) %>% 
              slice(1, .by = parkID) %>% 
              mutate(key = paste(siteID , parkID)) %>% 
  mutate(key = factor(key, levels = key[order(as.numeric(siteID))]))

p2 <- mad_hshade_plot + # PLot hillshaded map
  ggspatial::annotation_scale(location = 'tl',width_hint = .6,text_cex = 2)+
  geom_point(data = meta_mg_us, aes(longitude_WGS84, latitude_WGS84, colour = key),
             size = 0, show.legend = TRUE)+
  geom_text_repel(data=meta_mg_us , aes(longitude_WGS84 , latitude_WGS84 , 
                                label = siteID),
             size=7, , fontface = "bold",
             arrow = arrow(length = unit(0.25, 'cm'), type = 'closed'))+
 theme_linedraw(base_size = 25)+
  scale_colour_viridis_d(option="mako", end=.8)+
  labs(x = "Longitude" ,
       y = "Latitude")+
  theme(
        legend.position = "bottom",
        legend.justification = c(0.1, 0),
        legend.title = element_blank(),
        legend.text = element_text (size = 16))+
  guides(
         color=guide_legend(override.aes=list(fill=NA,size=0)),
         shape=guide_legend(override.aes=list(fill=NA,size=0)))


# plot ----------------------------------------------------------------------------------------
tiff("mg_maps/figs/sample_sites.tiff", width = 1400, height = 1000, compression = "lzw")
(p2 + p1) 
dev.off()


browseURL("mg_maps/figs/sample_sites.tiff")

