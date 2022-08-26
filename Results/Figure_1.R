#Figure 1. Location and elevation (in unit of meter) of weather stations
#used in the study over the East Asian monsoon (EAM) region.

source("H:\\Toolkit\\Mirror\\Comparison\\codes\\Functions.R", encoding = 'UTF-8')

df_meta <- read.csv("H:\\Toolkit\\Mirror\\Comparison\\EAM_alt.csv")

easm_meta <- subset(df_meta, lon <= 146 & lon >= 110 & lat <= 55 & lat >= 20)
easm_meta <- st_as_sf(easm_meta, coords = c('lon', 'lat'))
easm_meta <- st_set_crs(easm_meta, "+proj=longlat +datum=WGS84 +no_defs")

ter_seq9 <- c('#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014',
              '#cc4c02','#993404','#662506')


map <- ggplot() +
  borders(fill = "grey90") +
  geom_sf(aes(color = alt), data = easm_meta, size = 1.5) +
  scale_color_gradientn(name = "Elevation (m)", colors = terrain.colors(4), 
                        aesthetics = 'color', lim = c(0, 2500), oob = squish) +
  borders(fill = NA) +
  coord_sf(xlim = c(109.5, 146), ylim = c(19.5, 55), expand = FALSE) + 
  scalebar(easm_meta, dist = 500, dist_unit = "km", model = "WGS84",
           st.size = 2.5, st.bottom = TRUE, transform = TRUE, 
           location = "bottomright", 
           anchor = c(x = 143.5, y = 21), border.size = 0.5) +
  labs(x="Longitude (°E)", y="Latitude (°N)") + 
  theme_bw() +
  guides(color = guide_colourbar(barheight = 12, nbin = 30, 
                                 raster = FALSE,
                                 frame.colour = 'black', 
                                 ticks.colour = 'black')) +
  theme(panel.background = element_rect(fill = "azure"),
        legend.title = element_text(size = 9)) +
  geom_text(data = c.label, aes(x = long, y = lat, label = country), 
            inherit.aes = FALSE, color = "blue", fontface = "bold.italic")

map

ggsave("G:/내 드라이브/ICP Lab/Comparison/Figures/figure_1.png",
       plot = map, width = 5.5, height = 5.5, dpi = 500, units = "in")