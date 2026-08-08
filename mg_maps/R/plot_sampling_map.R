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
library(ggspatial)
library(ggplot2)
library(ggnewscale)

# Madagascar sites metadata -----------------------------------------------

meta_mg <- read_tsv("../../iba_data/sites_metadata_MG.tsv", show_col_types=FALSE)

# We want one coordinate for each site [avoid multiple overshadowing symbols for multitrap sites]
meta_mg <- meta_mg[!duplicated(meta_mg$siteID),]

# Fix park names; encoding is broken and names too long
parkNames <- c("Oronjia",
               "Montagne d'Ambre",
               "Analamerana",
               "Ankarana",
               "Loky Manambato",
               "Lokobe",
               "Marojejy",
               "Masoala",
               "Mahavavy Kinkony",
               "Baie de Baly",
               "Marotandrano",
               "Mananara-Nord",
               "Ankarafantsika",
               "Namoroka",
               "Analalava",
               "Betampona",
               "Ambohitantely",
               "Anjozorobe Angavo",
               "Ambatovy",
               "Maromizaha",
               "Tsingy de Bemaraha",
               "Kirindy (Conc. Forest.)",
               "Kirindy Mite",
               "Ranomafana",
               "Isalo",
               "Mikea",
               "Zombitse-Vohibasia",
               "Agnalazaha",
               "Beza-Mahafaly",
               "Tsimanampesotse",
               "Tsitongambarika",
               "Andohahela",
               "Cap de Sainte Marie")
meta_mg$parkID <- parkNames
meta_mg$parkNumber <- as.integer(factor(parkNames))


# Underlying map ----------------------------------------------------------

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


# plots ---------------------------------------------------------------------------------------

# Format meta-data
meta_mg <- meta_mg  |>
            mutate(trap_habitat=recode(trap_habitat , 
                                       "Dry_Forest"          = "Dry Forest",
                                       "Montane_Rainforest"  = "Montane forest",
                                       "Rainforest"          = "Wet forest"),
                   malaise_trap_type = recode(malaise_trap_type , "Single_trap" = "Single trap"))

# Make plot with trap type and habitat
p1 <- mad_hshade_plot + # PLot hillshaded map
  geom_spatraster(data = mad_elev, maxcell = Inf , show.legend = FALSE) +
#  ggspatial::annotation_scale(location = 'tl',width_hint = .6,text_cex = 2) +
  scale_fill_hypso_tint_c(limits = elev_limits , palette = "dem_poster",alpha = 0.4,direction = 1) +
  new_scale_fill() +
  geom_point(data=meta_mg , aes(longitude_WGS84 , latitude_WGS84 , 
                                fill = factor(malaise_trap_type), shape = factor(trap_habitat)),size=7) + 
  scale_fill_manual(values = c("Single trap" = "white", "Multitrap" = "grey60"),
                    breaks = c("Single trap", "Multitrap"),
                    labels = c("Single trap", "Multitrap")) +
#  scale_fill_viridis_d(option="mako", end=.8, guide=guide_legend()) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme_linedraw(base_size = 25) +
  labs(x = NULL, # "Longitude" ,
       y = NULL, # "Latitude" , 
       fill = "Traps",
       shape = "Habitat") +
  theme(legend.key = element_rect(fill = "white"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        axis.text.x = element_blank(), # element_text(size=17), 
        axis.text.y = element_blank()) + # element_text(size=17)) +
  guides(
         shape = guide_legend(override.aes=list(fill=NA, size=7), order=1),
         fill = guide_legend(override.aes = list(shape = 21, size = 7), order=2))


# Make plot with site names
# -------------------------

# Select unique parkIDs (already done above...)
meta_mg_us <- meta_mg %>% 
              select(parkID , siteID , longitude_WGS84 , latitude_WGS84,trap_habitat) %>% 
              slice(1, .by = parkID) %>% 
              mutate(key = paste(siteID , parkID)) %>% 
  mutate(key = factor(key, levels = key[order(as.numeric(siteID))]))

p2 <- mad_hshade_plot + # PLot hillshaded map
  ggspatial::annotation_scale(location = 'tl',width_hint = .6,text_cex = 2, pad_y=unit(0.5,"cm")) +
  geom_point(data = meta_mg_us, aes(longitude_WGS84, latitude_WGS84, colour = key),
             size = 1, show.legend = TRUE) +
  geom_text_repel(data=meta_mg_us , aes(longitude_WGS84 , latitude_WGS84 , 
                                label = siteID),
             size=7, , fontface = "bold",
             arrow = arrow(length = unit(0.25, 'cm'), type = 'closed')) +
  theme_linedraw(base_size = 25) +
  scale_colour_viridis_d(option="mako", end=.1) +      # Essentially the same colour...
  labs(x = "Longitude" ,
       y = "Latitude") +
  theme(
        legend.position = "bottom",
        legend.justification = c(0.1, 0),
        legend.title = element_blank(),
        legend.text = element_text (size = 18),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17)) +
  guides(
         color=guide_legend(override.aes=list(fill=NA,size=0, color="white")),
         shape=guide_legend(override.aes=list(fill=NA,size=0)))


# plot ----------------------------------------------------------------------------------------
ggsave(file="../figs/Fig_sample_sites.jpg",
       width = 21.0,
       height = 15.0,
       plot = (p2 + p1) + plot_layout(axis_titles="collect"))

