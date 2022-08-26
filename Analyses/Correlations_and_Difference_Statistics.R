# function
all_analysis <- function(meta) {
  for (k in 1:4) {
    #
    re <- t.result[[k]]
    data <- levels(factor(meta$station)) # 오름차순 정렬
    l <- length(data)

    re.list <- list()

    for (i in 1:l) {
      re.list[[i]] <- re[, c(l + 2, l + 3, i)]
      re.list[[i]] <- transform(re.list[[i]], "ID" = names(re)[i])
      re.list[[i]] <- re.list[[i]][, c(4, 1, 2, 3)]
      names(re.list[[i]]) <- c("ID", "year", "month", "avgtemp")
    }
    #
    re <- p.result[[k]]
    re.list2 <- list()

    for (i in 1:l) {
      re.list2[[i]] <- re[, c(l + 2, l + 3, i)]
      re.list2[[i]] <- transform(re.list2[[i]], "ID" = names(re)[i])
      re.list2[[i]] <- re.list2[[i]][, c(4, 1, 2, 3)]
      names(re.list2[[i]]) <- c("ID", "year", "month", "tp")
    }

    datalist <- list()
    datalist <- merge.rs(re.list, re.list2, station)

    name <- unique(station$ID) # levels(factor()) -> 자동으로 오름차순으로 정렬
    # name의 id 순서와 t.result의 행순서가 일치하는지 항상 확인

    t.anl <- temp_anl(name, datalist, 1981, 2020)

    colnames(t.anl) <- c("station", "variable", "method", "coefficient", "p", "df")

    coord.cor <- merge(t.anl, meta, by = "station", all = T)

    coord.cor$coefficient <- as.numeric(coord.cor$coefficient)
    coord.cor$p <- as.numeric(coord.cor$p)

    method <- as.factor(coord.cor$method)
    variable <- strsplit(coord.cor$variable, "[.]")

    ssn <- vector()
    var <- vector()

    for (i in 1:length(variable)) {
      ssn[i] <- variable[[i]][2]
      var[i] <- variable[[i]][1]
    }

    coord.cor$season <- ssn
    coord.cor$var <- var

    coord.cor$season <- factor(coord.cor$season, levels = c(
      "DJF", "MAM", "JJA",
      "SON", "ANN", "MON"
    ))
    coord.cor$method <- factor(coord.cor$method,
      levels = c(
        "pearson", "spearman", "kendall",
        "bias", "rmse"
      )
    )

    coord.cor <- transform(coord.cor,
      "significance" = ifelse(p < 0.01, "sig", "not sig")
    )

    coord.cor$significance <- factor(coord.cor$significance,
      levels = c("sig", "not sig")
    )

    write.csv(coord.cor, paste0(
      "\\\\climatelab\\Minseok\\Comparison\\",
      dir.list[k], "\\all_analysis_v4_station.csv"
    ))

    print(paste0(k, "th task has done"))
  }
}

#
meta_source <- c(
  "\\\\climatelab\\data\\ASOS\\meta_revised.csv",
  "\\\\climatelab\\data\\Northkorea_GTS\\metadata_NKO_station.csv",
  "\\\\climatelab\\data\\JMA\\meta_lonlatel.csv",
  "\\\\climatelab\\data\\CMA\\meta_cma.csv"
)

meta.list <- list()
num <- c("station", "ID", "id", "id")

for (i in 1:4) {
  meta.list[[i]] <- read.csv(meta_source[i])
  meta.list[[i]] <- meta.list[[i]][, c(num[i], "lon", "lat")]
  colnames(meta.list[[i]]) <- c("station", "lon", "lat")
}
#

# asos data ---------------------------------------------------------------
dir.list <- c("eras_asos", "jra_asos", "ncep_asos", "cru_asos")
d.list <- c("ERA-5", "JRA-55", "NCEP2", "CRU_TS")

t.result <- list()
p.result <- list()

station <-
  read.csv("\\\\climatelab\\Minseok\\Comparison\\Station\\all_month_asos.csv")
# meta <- read.csv("\\\\climatelab\\data\\ASOS\\meta_revised.csv", header = TRUE)

meta <- meta.list[[1]]

for (i in 1:4) {
  t.result[[i]] <- read.csv(paste0(
    "\\\\climatelab\\Minseok\\Comparison\\",
    dir.list[i], "\\", d.list[i],
    "_monthly mean(made)_station_c.csv"
  ))
  p.result[[i]] <- read.csv(paste0(
    "\\\\climatelab\\Minseok\\Comparison\\",
    dir.list[i], "\\", d.list[i],
    "_monthly tp(made)_station_mm.csv"
  ))
  t.result[[i]] <- t.result[[i]][, -1]
  p.result[[i]] <- p.result[[i]][, -1]
  colnames(t.result[[i]]) <- c(meta$station, "date", "year", "month")
  colnames(p.result[[i]]) <- c(meta$station, "date", "year", "month")
}

station <- station[, c("station", "avg_temp", "avg_prep", "year", "month")]
colnames(station) <- c("ID", "AVG_TEMP", "MON_PRE", "year", "month")

all_analysis(meta)

# amedas data ---------------------------------------------------------------

#
dir.list <- c("eras_amedas", "jra_amedas", "ncep_amedas", "cru_amedas")

t.result <- list()
p.result <- list()
#

for (i in 1:4) {
  t.result[[i]] <- read.csv(paste0(
    "\\\\climatelab\\Minseok\\Comparison\\",
    "AMeDAS_Reanalysis_Data\\", d.list[i],
    "_temperature.csv"
  ))
  p.result[[i]] <- read.csv(paste0(
    "\\\\climatelab\\Minseok\\Comparison\\",
    "AMeDAS_Reanalysis_Data\\", d.list[i],
    "_precipitation.csv"
  ))
  t.result[[i]] <- t.result[[i]][, -1]
  p.result[[i]] <- p.result[[i]][, -1]
  l <- 143
  colnames(t.result[[i]])[1:l] <-
    as.numeric(substr(colnames(t.result[[i]])[1:l], 2, 6))
  colnames(p.result[[i]])[1:l] <-
    as.numeric(substr(colnames(p.result[[i]])[1:l], 2, 6))
}

missing <- c(44226, 88971, 91107, 91146, 49251)
# CRU에서 나타나지 않는 섬 4개와 월 강수기록이 없는 후지산

meta <- meta.list[[3]]

for (i in 1:4) {
  t.result[[i]] <- t.result[[i]][, -(which(meta$station %in% missing))]
  p.result[[i]] <- p.result[[i]][, -(which(meta$station %in% missing))]
}

#
station <- read.csv(paste0(
  "\\\\climatelab\\Minseok\\Comparison\\Station\\",
  "all_month_amedas.csv"
))

station <- station[, c("id", "avg_temp", "avg_pre", "year", "month")]
colnames(station) <- c("ID", "AVG_TEMP", "MON_PRE", "year", "month")
station <- subset(station, !station$ID %in% missing)

meta <- subset(meta, !meta$station %in% missing) # 138개

all_analysis(meta)

# GTS(NK) data ---------------------------------------------------------------

#
dir.list <- c("eras_gts", "jra_gts", "ncep_gts", "cru_gts")

#
station <- read.csv("\\\\climatelab\\lab meeting documents\\북한\\NK\\r_result\\all_month_gts.csv")
meta <- meta.list[[2]]

rcn <- c(meta$station, "date", "year", "month")
for (i in 1:4) {
  t.result[[i]] <- read.csv(paste0(
    "\\\\climatelab\\Minseok\\Comparison\\",
    "GTS_Reanalysis_Data\\", d.list[i],
    "_temperature.csv"
  ))
  p.result[[i]] <- read.csv(paste0(
    "\\\\climatelab\\Minseok\\Comparison\\",
    "GTS_Reanalysis_Data\\", d.list[i],
    "_precipitation.csv"
  ))
  t.result[[i]] <- t.result[[i]][, -1]
  p.result[[i]] <- p.result[[i]][, -1]
  colnames(t.result[[i]]) <- rcn
  colnames(p.result[[i]]) <- rcn
}

#
station <- station[, c(2, 5, 8, 3, 4)]
colnames(station) <- c("ID", "AVG_TEMP", "MON_PRE", "year", "month")

all_analysis(meta)

# cma data ---------------------------------------------------------------

dir.list <- c("eras_cma", "jra_cma", "ncep_cma", "cru_cma")

meta <- meta.list[[4]]

rcn <- c(meta$station, "date", "year", "month")

for (i in 1:4) {
  t.result[[i]] <- read.csv(paste0(
    "\\\\climatelab\\Minseok\\Comparison\\CMA_Reanalysis_Data\\", d.list[i],
    "_temperature.csv"
  ))
  p.result[[i]] <- read.csv(paste0(
    "\\\\climatelab\\Minseok\\Comparison\\CMA_Reanalysis_Data\\", d.list[i],
    "_precipitation.csv"
  ))
  names(t.result[[i]])[2:518] <- as.numeric(substr(names(t.result[[i]])[2:518], 2, 6))
  names(p.result[[i]])[2:518] <- as.numeric(substr(names(p.result[[i]])[2:518], 2, 6))
  if (i == 4) {
    oob <- vector()
    for (k in 1:length(t.result[[4]])) {
      if (is.na(t.result[[4]][1, k])) {
        oob <- c(oob, names(t.result[[4]][k]))
        print(k) # 131, 132, 511번째 행
      }
    }
  }
}

for (i in 1:4) {
  t.result[[i]] <- t.result[[i]][, c(-1, -132, -133, -512)]
  p.result[[i]] <- p.result[[i]][, c(-1, -132, -133, -512)]
}

station <- read.csv("\\\\climatelab\\data\\CMA\\all_month_cma.csv")
station <- station[, c("ID", "AVG_TEMP", "AVG_PRE", "year", "month")]
colnames(station)[3] <- "MON_PRE"
station <- subset(station, !station$ID %in% oob)

meta <- subset(meta, !meta$station %in% oob) # 514개

all_analysis(meta)
