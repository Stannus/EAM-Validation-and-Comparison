library("hydroGOF")

source("C:\\Users\\Minseok\\projects\\Comparison\\Functions.R", encoding = 'UTF-8')

#
meta_source <- c("C:\\Users\\Minseok\\projects\\Comparison\\Revision\\meta_revised.csv",
               "C:\\Users\\Minseok\\projects\\Comparison\\Revision\\metadata_NKO_station.csv",
               "C:\\Users\\Minseok\\projects\\Comparison\\Revision\\meta_lonlatel.csv", 
               "C:\\Users\\Minseok\\projects\\Comparison\\Revision\\meta_cma.csv")

meta.list <- list()
num <- c('station', 'ID', 'id', 'id')

for (i in 1:4){
  meta.list[[i]] <- read.csv(meta_source[i])
  meta.list[[i]] <- meta.list[[i]][,c(num[i], 'lon', 'lat')]
  colnames(meta.list[[i]]) <- c('station', 'lon', 'lat')
}
#

# asos data ---------------------------------------------------------------
dir.list <- c("eras_asos", "jra_asos", "ncep_asos", "cru_asos")
d.list <- c("ERA-5", 'JRA-55', 'NCEP2', 'CRU_TS')

t.result <- list()
p.result <- list()

station <- 
  read.csv("H:\\Toolkit\\Mirror\\Comparison\\Station\\all_month_asos.csv")
#meta <- read.csv("\\\\climatelab\\data\\ASOS\\meta_revised.csv", header = TRUE)

meta <- meta.list[[1]]

for (i in 1:4){
  t.result[[i]] <- read.csv(paste0("H:\\Toolkit\\Mirror\\Comparison\\", 
                                   dir.list[i], "\\", d.list[i],
                                   "_monthly mean(made)_station_c.csv"))
  p.result[[i]] <- read.csv(paste0("H:\\Toolkit\\Mirror\\Comparison\\",
                                   dir.list[i], "\\", d.list[i],
                                   "_monthly tp(made)_station_mm.csv"))
  t.result[[i]] <- t.result[[i]][,-1]
  p.result[[i]] <- p.result[[i]][,-1]
  colnames(t.result[[i]]) <- c(meta$station, "date", "year", "month")
  colnames(p.result[[i]]) <- c(meta$station, "date", "year", "month")
}

station <- station[,c('station','avg_temp','avg_prep','year','month')]
colnames(station) <- c("ID", "AVG_TEMP", "MON_PRE", "year", "month")

for(k in 1:4){
    
  data <- levels(factor(meta$station)) #오름차순 정렬
  l <- length(data)
    
  # 기온 
  re <- t.result[[k]] 
  re.list <- list()
    
    for (i in 1:l) {
      re.list[[i]] <- re[, c(l + 2, l + 3, i)]
      re.list[[i]] <- transform(re.list[[i]], "ID" = names(re)[i])
      re.list[[i]] <- re.list[[i]][, c(4, 1, 2, 3)]
      names(re.list[[i]]) <- c("ID", "year", "month", "avgtemp")
    }
    # 강수
    re <- p.result[[k]] 
    re.list2 <- list()
    
    for(i in 1:l){
      re.list2[[i]] <- re[,c(l+2,l+3,i)]
      re.list2[[i]] <- transform(re.list2[[i]], "ID" = names(re)[i])
      re.list2[[i]] <- re.list2[[i]][,c(4,1,2,3)]
      names(re.list2[[i]]) <- c("ID", "year", "month", "tp")
    }
    
    datalist <- list()
    datalist <- merge.rs(re.list, re.list2, station)
    
    name <- unique(station$ID) # levels(factor()) -> 자동으로 오름차순으로 정렬
    # name의 id 순서와 t.result의 행순서가 일치하는지 항상 확인
    
    #tmp <- avgvar(datalist, 1981, 2020) 
    
    NSE
    
    
    t.anl <- temp_anl(name, datalist, 1981, 2020)
    
    colnames(t.anl) <- c("station", "variable", "method", "coefficient", "p","df") 
    
    #
    
    coord.cor <- merge(t.anl, meta, by = "station", all = T) 
    
    coord.cor$coefficient <- as.numeric(coord.cor$coefficient)
    coord.cor$p <- as.numeric(coord.cor$p)
    
    method <- as.factor(coord.cor$method)
    variable <- strsplit(coord.cor$variable, "[.]")
    
    ssn <- vector()
    var <- vector()
    
    for (i in 1:length(variable)){
      ssn[i] <- variable[[i]][2]
      var[i] <- variable[[i]][1]
    }
    
    coord.cor$season <- ssn
    coord.cor$var <- var
    
    coord.cor$season <- factor(coord.cor$season, levels = c("DJF", "MAM", "JJA", 
                                                            "SON", "ANN", "MON"))
    coord.cor$method <- factor(coord.cor$method, 
                               levels = c("pearson", "spearman", "kendall", 
                                          'bias', 'rmse'))
    
    coord.cor <- transform(coord.cor, 
                           "significance" = ifelse(p < 0.01, "sig", "not sig"))
    
    coord.cor$significance <- factor(coord.cor$significance, 
                                     levels = c("sig", "not sig"))
    
    write.csv(coord.cor, paste0("H:\\Toolkit\\Mirror\\Comparison\\", 
                                dir.list[k], "\\all_analysis_v4_station.csv"))

