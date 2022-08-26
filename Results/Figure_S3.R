source("C:\\Users\\Minseok\\projects\\Comparison\\Functions.R", encoding = 'UTF-8')
eam <- read.csv("C:\\Users\\Minseok\\projects\\Comparison\\Revision\\eam_best_inc_CRU.csv")

eam[,'season'] <- factor(eam[,'season'], 
                             levels = c("ANN", "DJF", "MAM", "JJA", "SON"))
eam <- subset(eam, eam$method %in% c("Pearson", "RMSE", "Bias"))
eam[,'method'] <- factor(eam[,'method'], 
                             levels = c("Pearson", "RMSE", "Bias"))
eam[,'data'] <- factor(eam[,'data'], levels = c("ERA5", "JRA55", "NCEP2", "CRU"))

eam <- st_as_sf(eam, coords = c('lon', 'lat'))
eam <- st_make_valid(eam)
eam <- st_set_crs(eam,"+proj=longlat +datum=WGS84 +no_defs")

colorset <- c("#1b9e77","#d95f02", "#7570b3", "#ffd400") 

st_data <- list()
stn <- list(eam$var, eam$method)
st_data <- split(eam, stn)

for(i in 1:length(names(st_data))){
  assign(names(st_data)[i], best_map(st_data[[i]], subplot = TRUE))
}

avg_temp.Pearson <- avg_temp.Pearson + 
  theme(strip.text.x = element_text(hjust = 0.5, face = "bold", size = 14))
pre.Pearson <- pre.Pearson + 
  theme(strip.text.x = element_text(hjust = 0.5, face = "bold", size = 14))

avg_temp.Bias <- avg_temp.Bias + 
  theme(axis.text.x = element_text(), axis.ticks.x = element_line())
pre.Bias <- pre.Bias + 
  theme(axis.text.x = element_text(), axis.ticks.x = element_line())
  
t <- avg_temp.Pearson/avg_temp.RMSE/avg_temp.Bias + 
  plot_layout(guides = 'collect') + 
  plot_annotation(title = '(a) Temperature', 
                  theme = theme(plot.title = element_text(size = 20, face = 'bold'),
                                legend.position = 'bottom'))

p <- pre.Pearson/pre.RMSE/pre.Bias + 
  plot_layout(guides = 'collect') + 
  plot_annotation(title = '(b) Precipitation', 
                  theme = theme(plot.title = element_text(size = 20, face = 'bold'),
                                legend.position = 'bottom'))

ggsave("G:/내 드라이브/ICP Lab/Comparison/Figures/figure_s3a.png", plot = t, width = 14, height = 11.5, dpi = "retina", units = "in")
ggsave("G:/내 드라이브/ICP Lab/Comparison/Figures/figure_s3b.png", plot = p, width = 14, height = 11.5, dpi = "retina", units = "in")