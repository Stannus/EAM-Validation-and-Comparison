
source("H:\\Toolkit\\Mirror\\Comparison\\codes\\Functions.R", encoding = 'UTF-8')
eam <- read.csv("H:\\Toolkit\\Mirror\\Comparison\\EAM\\v3\\eam_rmse.csv")

eam[,'method'] <- factor(eam[,'method'], levels = c("rmse"))
eam[,'data'] <- factor(eam[,'data'], levels = c("ERA5", "JRA55", "NCEP2", "CRU"))
eam[,'significance'] <- factor(eam[,'significance'], levels = c(1, 0))

eam <- st_as_sf(eam, coords = c('lon', 'lat'))
eam <- st_make_valid(eam)
eam <- st_set_crs(eam, "+proj=longlat +datum=WGS84 +no_defs")

st_data <- list()
stn <- list(eam$var, eam$season)
st_data <- split(eam, stn)

for(i in 1:length(names(st_data))){
  assign(names(st_data)[i], anl_map(st_data[[i]]))
}

t <- avg_temp.ANN/avg_temp.DJF/avg_temp.MAM/avg_temp.JJA/avg_temp.SON + 
  plot_layout(guides = 'collect') + 
  plot_annotation(title = '(a) Temperature', theme = theme(plot.title = element_text(size = 20, face = 'bold')))

p <- pre.ANN/pre.DJF/pre.MAM/pre.SON/pre.JJA + 
  plot_layout(guides = 'collect') + 
  plot_annotation(title = '(b) Precipitation', 
                  theme = theme(plot.title = element_text(size = 20, face = 'bold'),
                                legend.spacing = unit(6.5, "inch"), 
                                legend.justification = "bottom"))

ggsave("G:/내 드라이브/ICP Lab/Comparison/Figures/figure_5a.png", device = 'png',
       plot = t, width = 14, height = 21, dpi = "retina", units = 'in')
ggsave("G:/내 드라이브/ICP Lab/Comparison/Figures/figure_5b.png", device = 'png',
       plot = p, width = 14, height = 21, dpi = "retina", units = 'in')