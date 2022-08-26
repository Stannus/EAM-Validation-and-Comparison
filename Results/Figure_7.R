eam <- read.csv("H:\\Toolkit\\Mirror\\Comparison\\EAM\\v2\\eam_clm.csv")
eam <- eam[,-1]

eam[,'season'] <- factor(eam[,'season'], 
                              levels = c("ANN", "DJF", "MAM", "JJA", "SON"))
eam[,'data'] <- factor(eam[,'data'], 
                            levels = c("ERA5", "JRA55", "NCEP2", "CRU"))


stn.t.avg <- eam[,c(-1, -3, -7, -8, -11, -12)]
stn.p.avg <- eam[,c(-1, -3, -5, -6, -9, -10)]

ssn.stn.t.avg <- split(stn.t.avg, stn.t.avg$season)
ssn.stn.p.avg <- split(stn.p.avg, stn.p.avg$season)

for(i in 1:length(names(ssn.stn.t.avg))){
  assign(paste0("t.", names(ssn.stn.t.avg)[i]), get_multi(ssn.stn.t.avg[[i]]))
}

for(i in 1:length(names(ssn.stn.p.avg))){
  assign(paste0("p.", names(ssn.stn.p.avg)[i]), get_multi2(ssn.stn.p.avg[[i]]))
}

t.ANN <- t.ANN + theme(strip.text.x = element_text(hjust = 0.5,
                                                   face = "bold", size = 14))
t.MAM <- t.MAM + theme(axis.title.y = element_text(face = "bold", size = 16))
t.SON <- t.SON + theme(axis.title.x = element_text(face = "bold", size = 16))

t <- t.ANN/t.DJF/t.MAM/t.JJA/t.SON +
  plot_layout(guides = 'collect') + 
  plot_annotation(title = '(a) Temperature', 
                  theme = theme(plot.title = element_text(size = 20, face = 'bold'),
                                legend.position = 'bottom'))

p.ANN <- p.ANN + theme(strip.text.x = element_text(hjust = 0.5, 
                                                   face = "bold", size = 14))
p.MAM <- p.MAM + theme(axis.title.y = element_text(face = "bold", size = 16))
p.SON <- p.SON + theme(axis.title.x = element_text(face = "bold", size = 16))

p <- p.ANN/p.DJF/p.MAM/p.JJA/p.SON +
  plot_layout(guides = 'collect') + 
  plot_annotation(title = '(b) Precipitation', 
                  theme = theme(plot.title = element_text(size = 20, face = 'bold'),
                                legend.position = 'bottom'))

ggsave("G:/내 드라이브/ICP Lab/Comparison/Figures/figure_7a.png", 
       plot = t, width = 14, height = 20, dpi = "retina", units = "in")
ggsave("G:/내 드라이브/ICP Lab/Comparison/Figures/figure_7b.png", 
       plot = p, width = 14, height = 20, dpi = "retina", units = "in")