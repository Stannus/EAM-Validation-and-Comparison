eam <- read.csv("H:\\Toolkit\\Mirror\\Comparison\\EAM\\eam_clm.csv")
eam <- eam[,-1]

season <- c("ANN", "DJF", "MAM", "JJA", "SON")
data <- c("ERA5", "JRA55", "NCEP2", "CRU")
#country <- c("asos", "amedas", "cma", "gts")

eam[,'season'] <- factor(eam[,'season'], levels = season)

t.ref <- eam[, c("ID", "season", "avg_TEMP")]
p.ref <- eam[, c("ID", "season", "avg_PRE")]

t.ref <- unique(t.ref)
p.ref <- unique(p.ref)

t.mod <- eam[, c("ID", "season", "data", "avg_t2m")]
p.mod <- eam[, c("ID", "season", "data", "avg_tp")]

colnames(t.mod)[4] <- "temp"
colnames(t.ref)[3] <- "temp"
colnames(p.mod)[4] <- "pre"
colnames(p.ref)[3] <- "pre"

t.mod$data <- factor(t.mod$data, levels = data)
p.mod$data <- factor(p.mod$data, levels = data)

#seperate vectors by data & seasons

t.ref.split <- split(t.ref[,c(1,3)], t.ref$season)
p.ref.split <- split(p.ref[,c(1,3)], p.ref$season)
t.mod.split <- split(t.mod[,c(1,4)], list(t.mod$data, t.mod$season))
p.mod.split <- split(p.mod[,c(1,4)], list(p.mod$data, p.mod$season))

colors = c("#0000FF", "#FF0000", "#FF00FF", "#00FF00")

corSeason <- 1:5
lab <- c("1", "2", "3", "4", "5")

png("G:/내 드라이브/ICP Lab/Comparison/Figures/figure_9.png", width = 12, height = 6, res = 320, units = "in")

par(mfrow = c(1,2), pty = "s")

for (i in 1:5){ #i: season
  for (j in 1:4){ #j: data
    if(i == 1 & j == 1){
      taylor.diagram2(t.ref.split[[i]][,2], t.mod.split[[((i-1)*4)+j]][,2], 
                     main = "", xlab = "", ylab = "Standard Deviations (Normalized)",
                     ref.sd = TRUE, normalize = TRUE, #sd.arcs = TRUE, 
                     pch = 16, col = colors[j], lwd = 2, text = corSeason[i])
      title("(a) Temperature", adj = 0)
    }else{
      taylor.diagram2(t.ref.split[[i]][,2], t.mod.split[[((i-1)*4)+j]][,2], 
                     add = TRUE, normalize = TRUE,
                     pch = 16, col = colors[j], lwd = 2, text = corSeason[i])}}}

for (i in 1:5){ #i: season
  for (j in 1:4){ #j: data
    if(i == 1 & j == 1){
      taylor.diagram2(p.ref.split[[i]][,2], p.mod.split[[((i-1)*4)+j]][,2], 
                     main = "", xlab = "", ylab = "Standard Deviation (Normalized)",
                     ref.sd = TRUE, normalize = TRUE, #sd.arcs = TRUE, 
                     pch = 16, col = colors[j], lwd = 2, text = corSeason[i])
      title("(b) Precipitation", adj = 0)
    }else{
      taylor.diagram2(p.ref.split[[i]][,2], p.mod.split[[((i-1)*4)+j]][,2], 
                     add = TRUE, normalize = TRUE,
                     pch = 16, col = colors[j], lwd = 2, text = corSeason[i])}}}

legend(1.4, 1.7, legend=season, pch=lab, border="white", box.lty=0, cex=0.85)
legend(0.95, 1.7, legend=data, fill=colors, border="white", box.lty=0, cex=0.85)

dev.off()