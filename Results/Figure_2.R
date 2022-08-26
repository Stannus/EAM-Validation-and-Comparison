# Figure 2. Climatology of (a) temperature and (b) precipitation
# from observation data at each station in the EAM region
# using the annual and seasonal (DJF, MAM, JJA, and SON) averages
# over the forty years of 1981-2020 for North Korea, South Korea, and Japan
# and the thirty years of 1981-2010 for China.
# The scale bar of temperature is in the range of -25 to 30 °C
# and the midpoint of scale is at 2.5 °C.

source("H:\\Toolkit\\Mirror\\Comparison\\codes\\Functions.R", encoding = 'UTF-8')
easm_avg <- read.csv("H:\\Toolkit\\Mirror\\Comparison\\EAM\\v2\\eam_mean.csv")

season <- c("ANN", "DJF", "MAM", "JJA", "SON")
data <- c("STN", "ERA-5", "JRA-55", "NCEP2", "CRU_TS")

easm_avg[, "season"] <- factor(easm_avg[, "season"], levels = season)
easm_avg[, "data"] <- factor(easm_avg[, "data"], levels = data)

easm_sf <- st_as_sf(easm_avg, coords = c("lon", "lat"))
easm_sf <- st_set_crs(easm_sf, "+proj=longlat +datum=WGS84 +no_defs")

st_data <- list()
stn <- list(easm_sf$var, easm_sf$data)
st_data <- split(easm_sf, stn)

t <- mean_map(st_data$temp.STN)
p <- mean_map(st_data$pre.STN)

t
p

ggsave("G:/내 드라이브/ICP Lab/Comparison/Figures/figure_2a.png",
       plot = t, width = 14, height = 4, dpi = "retina", units = "in"
       )
ggsave("G:/내 드라이브/ICP Lab/Comparison/Figures/figure_2b.png",
       plot = p, width = 14, height = 4, dpi = "retina", units = "in"
       )
