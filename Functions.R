# Function for correlation analysis and preprocessing

library("raster")
library("ggmap")
library("sp")
library("rgeos")
library("maptools")
library("ggrepel")
library("scales")
library("Metrics")
library("ggpubr")
library("sf")
library("dplyr")
library("spData")
library("ggplot2")
library("lwgeom")
library("classInt")
library("maptools")
library("ggsn")
library("ggpmisc")
library("gridExtra")
library("maps")
library("lemon")
library("patchwork")
library("plotrix")
library("hydroGOF")

avg.ssn <- function(data_ann, ssn) {
  avg_ann <- data.frame()
  a <- data_ann
  avg_ann[1, 1] <- ssn
  avg_ann[1, 2] <- a[1, "ID"] # ID
  avg_ann[1, 3] <- mean(a[, "t2m"], na.rm = TRUE)
  avg_ann[1, 4] <- mean(a[, "AVG_TEMP"], na.rm = TRUE)
  avg_ann[1, 5] <- mean(a[, "tp"], na.rm = TRUE)
  avg_ann[1, 6] <- mean(a[, "MON_PRE"], na.rm = TRUE)
  avg_ann[1, 7] <- sum(!is.na(a[, "t2m"]))
  avg_ann[1, 8] <- sum(!is.na(a[, "AVG_TEMP"]))
  avg_ann[1, 9] <- sum(!is.na(a[, "tp"]))
  avg_ann[1, 10] <- sum(!is.na(a[, "MON_PRE"]))

  return(avg_ann)
}

# required package :
#  library("raster")
#  library("ncdf4")
#
# input: f1 : raster, e : spatial coverage, station : station coordinate data frame
# output: dataframe of timeseries grid data at each station

extcell <- function(f1, e, station) {
  r1 <- stack(f1)

  result <- data.frame()
  crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"

  r1 <- crop(r1, e)

  for (j in 1:nlayers(r1)) { # length of month
    for (i in 1:nrow(station)) { # number of station
      point <- cbind(station[i, "lon"], station[i, "lat"]) # x가 경도, y가 위도
      result[j, i] <- extract(r1[[j]], point)
    }
  }

  print(paste(f1, "extraction is done."))

  return(result)
}

###

# required package :
#  library("raster")
#  library("ncdf4")
#  library("lubridate")
#  library("ncdf4.helpers")
# input: f1 : raster, result : extcelled dataframe
# output: dataframe + date, year, month in col

markdate <- function(result, date) {
  l <- ncol(result)

  for (i in 1:length(date)) {
    result[i, (l + 1)] <- date[i]
    result[i, (l + 2)] <- year(date[i])
    result[i, (l + 3)] <- month(date[i])
  }
  return(result)
}

###  test

###

# required package:
# input: asos : station data, era : gridded data, var : variation name
# output: gridded data list

tolist <- function(asos, era, var) {
  l <- length(levels(factor(asos$ID)))

  eralis <- list()

  for (i in 1:l) {
    eralis[[i]] <- era[, c(l + 2, l + 3, i)]
    eralis[[i]] <- transform(eralis[[i]], "ID" = names(era)[i])
    eralis[[i]] <- eralis[[i]][, c(4, 1, 2, 3)]
    names(eralis[[i]]) <- c("ID", "year", "month", var) # avgtemp or tp
  }

  return(eralis)
}
###


# required package:
# input:
# output:

merge.rs <- function(eralis, eralis2, asos) {
  datalist <- list()

  l <- length(levels(factor(asos$ID)))

  for (i in 1:l) {
    d <- as.data.frame(eralis[[i]])
    d2 <- as.data.frame(eralis2[[i]])
    g <- subset(asos, ID == d[1, 1])
    g2 <- subset(asos, ID == d2[1, 1])

    if (g[1, 1] == g2[1, 1]) {
      data <- merge(d, g, by = c("ID", "year", "month"))

      data2 <- merge(d2, g2, by = c("ID", "year", "month"))

      data3 <- merge(data, data2, by = colnames(data)[c(1:3, 5:6)])

      colnames(data3)[6] <- "t2m"

      datalist[[i]] <- data3
    }
  }
  return(datalist)
}

### test

# required package:
# input:
# output:

avgvar <- function(datalist, S, E) {
  mon <- list()
  ann <- list()
  djf <- list()
  mam <- list()
  jja <- list()
  son <- list()

  for (k in 1:l) {
    # 겨울의 경우, 전년 12월이 들어가므로, 계산의 편리를 위해 전년을 다음해로 바꾸기
    # 전년도 12월, 1월, 2월이 한쌍이 되도록 해야함... => 12월은 연도를 +1 로 하자

    temp <- datalist[[k]]

    temp <- subset(temp, year >= S & year <= E)

    ann_t1 <- aggregate(t2m ~ ID + year, data = temp, mean) # 평균기온
    ann_t2 <- aggregate(AVG_TEMP ~ ID + year, data = temp, mean)
    ann_t <- merge(ann_t1, ann_t2, by = c("ID", "year"), all = TRUE)

    ann_p1 <- aggregate(tp ~ ID + year, data = temp, mean) # 강수량
    ann_p2 <- aggregate(MON_PRE ~ ID + year, data = temp, mean)
    ann_p <- merge(ann_p1, ann_p2, by = c("ID", "year"), all = TRUE)

    ann[[k]] <- merge(ann_t, ann_p, by = c("ID", "year"))

    mon[[k]] <- temp

    temp[which(temp$month == 12), 2] <- temp[which(temp$month == 12), 2] + 1

    wi <- subset(temp, month >= 12 | month <= 2)
    sp <- subset(temp, month >= 3 & month <= 5)
    su <- subset(temp, month >= 6 & month <= 8)
    au <- subset(temp, month >= 9 & month <= 11)

    # 연도별 계절 평균값 구하기, 만약 한 달이라도 na가 있으면 그 연도는 계산 안됨.
    for (j in 0:3) {
      for (i in 1:nrow(wi)) {
        if (is.na(wi[i, 4 + j])) {
          wi[which(wi[, 2] == wi[i, 2]), 4 + j] <- NA
        }
      }
    }

    for (j in 0:3) {
      for (i in 1:nrow(sp)) {
        if (is.na(sp[i, 4 + j])) {
          sp[which(sp[, 2] == sp[i, 2]), 4 + j] <- NA
        }
      }
    }

    for (j in 0:3) {
      for (i in 1:nrow(su)) {
        if (is.na(su[i, 4 + j])) {
          su[which(su[, 2] == su[i, 2]), 4 + j] <- NA
        }
      }
    }


    for (j in 0:3) {
      for (i in 1:nrow(au)) {
        if (is.na(au[i, 4 + j])) {
          au[which(au[, 2] == au[i, 2]), 4 + j] <- NA
        }
      }
    }

    # 평균기온(겨울)
    djf_t1 <- aggregate(t2m ~ ID + year, data = wi, mean)
    djf_t2 <- aggregate(AVG_TEMP ~ ID + year, data = wi, mean)
    djf_t <- merge(djf_t1, djf_t2, by = c("ID", "year"), all = TRUE)

    # 강수량(겨울)
    djf_p1 <- aggregate(tp ~ ID + year, data = wi, mean)
    djf_p2 <- aggregate(MON_PRE ~ ID + year, data = wi, mean)
    djf_p <- merge(djf_p1, djf_p2, by = c("ID", "year"), all = TRUE)

    # 평균기온(봄)
    mam_t1 <- aggregate(t2m ~ ID + year, data = sp, mean)
    mam_t2 <- aggregate(AVG_TEMP ~ ID + year, data = sp, mean)
    mam_t <- merge(mam_t1, mam_t2, by = c("ID", "year"), all = TRUE)

    # 강수량(봄)
    mam_p1 <- aggregate(tp ~ ID + year, data = sp, mean)
    mam_p2 <- aggregate(MON_PRE ~ ID + year, data = sp, mean)
    mam_p <- merge(mam_p1, mam_p2, by = c("ID", "year"), all = TRUE)

    # 평균기온(여름)
    jja_t1 <- aggregate(t2m ~ ID + year, data = su, mean)
    jja_t2 <- aggregate(AVG_TEMP ~ ID + year, data = su, mean)
    jja_t <- merge(jja_t1, jja_t2, by = c("ID", "year"), all = TRUE)

    # 강수량(여름)
    jja_p1 <- aggregate(tp ~ ID + year, data = su, mean)
    jja_p2 <- aggregate(MON_PRE ~ ID + year, data = su, mean)
    jja_p <- merge(jja_p1, jja_p2, by = c("ID", "year"), all = TRUE)

    # 평균기온(가을)
    son_t1 <- aggregate(t2m ~ ID + year, data = au, mean)
    son_t2 <- aggregate(AVG_TEMP ~ ID + year, data = au, mean)
    son_t <- merge(son_t1, son_t2, by = c("ID", "year"), all = TRUE)

    # 강수량(가을)
    son_p1 <- aggregate(tp ~ ID + year, data = au, mean)
    son_p2 <- aggregate(MON_PRE ~ ID + year, data = au, mean)
    son_p <- merge(son_p1, son_p2, by = c("ID", "year"), all = TRUE)

    djf[[k]] <- merge(djf_t, djf_p, by = c("ID", "year"))
    mam[[k]] <- merge(mam_t, mam_p, by = c("ID", "year"))
    jja[[k]] <- merge(jja_t, jja_p, by = c("ID", "year"))
    son[[k]] <- merge(son_t, son_p, by = c("ID", "year"))
  }

  avg_ann <- data.frame()

  for (i in 1:l) {
    tmp1 <- avg.ssn(ann[[i]], "ANN")
    tmp2 <- avg.ssn(djf[[i]], "DJF")
    tmp3 <- avg.ssn(mam[[i]], "MAM")
    tmp4 <- avg.ssn(jja[[i]], "JJA")
    tmp5 <- avg.ssn(son[[i]], "SON")
    tmp6 <- avg.ssn(mon[[i]], "MON")
    avg_ann <- rbind(avg_ann, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)
  }

  colnames(avg_ann) <- c("season", "ID", "avg_t2m", "avg_TEMP", "avg_tp", "avg_PRE", "freq_t2m", "freq_TEMP", "freq_tp", "freq_PRE")

  return(avg_ann)
}

# required package:
# input:
# output:

tempcor <- function(name, datalist, S, E) { # <- 두개로 나누기(평균값 구하기, 상관분석 하기)

  # 상관분석(평균기온): 지점별, 월별. 켄달, 피어슨, 스피어만 상관분석

  resultli <- list()

  l <- length(name)

  for (i in 1:l) {
    temp <- datalist[[i]]
    temp <- subset(temp, year >= S & year <= E) # 분석기간 설정

    temp <- temp[complete.cases(temp[, c("AVG_TEMP", "MON_PRE")]), ]

    cor_t1 <- cor.test(temp$t2m, temp$AVG_TEMP, method = "pearson")
    cor_t2 <- cor.test(temp$t2m, temp$AVG_TEMP, method = "spearman")
    cor_t3 <- cor.test(temp$t2m, temp$AVG_TEMP, method = "kendall")

    cor_p1 <- cor.test(temp$tp, temp$MON_PRE, method = "pearson")
    cor_p2 <- cor.test(temp$tp, temp$MON_PRE, method = "spearman")
    cor_p3 <- cor.test(temp$tp, temp$MON_PRE, method = "kendall")

    result <- data.frame("variable", "cor.test", "coefficient", "p-value", "df")

    result[c(1:3), 1] <- "avg_temp.MON"
    result[1, 3] <- cor_t1[4] # correlation coefficient(pearson)
    result[1, 4] <- cor_t1[3] # p-value
    result[1, 5] <- cor_t1[2] # df
    result[2, 3] <- cor_t2[4] # correlation coefficient(spearman)
    result[2, 4] <- cor_t2[3] # p-value
    result[3, 3] <- cor_t3[4] # correlation coefficient(kendall)
    result[3, 4] <- cor_t3[3] # p-value

    result[c(4:6), 1] <- "pre.MON"
    result[4, 3] <- cor_p1[4] # correlation coefficient(pearson)
    result[4, 4] <- cor_p1[3] # p-value
    result[4, 5] <- cor_p1[2] # df
    result[5, 3] <- cor_p2[4] # correlation coefficient(spearman)
    result[5, 4] <- cor_p2[3] # p-value
    result[6, 3] <- cor_p3[4] # correlation coefficient(kendall)
    result[6, 4] <- cor_p3[3] # p-value

    result[, 2] <- c("pearson", "spearman", "kendall")
    resultli[[i]] <- result
  }

  #-------------------------------------------------------------------------------
  # 지점별 계절별

  for (k in 1:l) {
    # 겨울의 경우, 전년 12월이 들어가므로, 계산의 편리를 위해 전년을 다음해로 바꾸기
    # 전년도 12월, 1월, 2월이 한쌍이 되도록 해야함... => 12월은 연도를 +1 로 하자

    temp <- datalist[[k]]

    temp <- subset(temp, year >= S & year <= E)

    temp[which(temp$month == 12), 2] <- temp[which(temp$month == 12), 2] + 1

    wi <- subset(temp, month >= 12 | month <= 2)
    sp <- subset(temp, month >= 3 & month <= 5)
    su <- subset(temp, month >= 6 & month <= 8)
    au <- subset(temp, month >= 9 & month <= 11)

    # 연도별 계절 평균값 구하기, 만약 한 달이라도 na가 있으면 그 연도는 계산 안됨.
    for (j in 0:3) {
      for (i in 1:nrow(wi)) {
        if (is.na(wi[i, 4 + j])) {
          wi[which(wi[, 2] == wi[i, 2]), 4 + j] <- NA
        }
      }
    }

    for (j in 0:3) {
      for (i in 1:nrow(sp)) {
        if (is.na(sp[i, 4 + j])) {
          sp[which(sp[, 2] == sp[i, 2]), 4 + j] <- NA
        }
      }
    }

    for (j in 0:3) {
      for (i in 1:nrow(su)) {
        if (is.na(su[i, 4 + j])) {
          su[which(su[, 2] == su[i, 2]), 4 + j] <- NA
        }
      }
    }


    for (j in 0:3) {
      for (i in 1:nrow(au)) {
        if (is.na(au[i, 4 + j])) {
          au[which(au[, 2] == au[i, 2]), 4 + j] <- NA
        }
      }
    }

    # 평균기온(겨울)
    djf_t1 <- aggregate(t2m ~ ID + year, data = wi, mean)
    djf_t2 <- aggregate(AVG_TEMP ~ ID + year, data = wi, mean)
    djf_t <- merge(djf_t1, djf_t2, by = c("ID", "year"), all = TRUE)

    # 강수량(겨울)
    djf_p1 <- aggregate(tp ~ ID + year, data = wi, mean)
    djf_p2 <- aggregate(MON_PRE ~ ID + year, data = wi, mean)
    djf_p <- merge(djf_p1, djf_p2, by = c("ID", "year"), all = TRUE)

    # 평균기온(봄)
    mam_t1 <- aggregate(t2m ~ ID + year, data = sp, mean)
    mam_t2 <- aggregate(AVG_TEMP ~ ID + year, data = sp, mean)
    mam_t <- merge(mam_t1, mam_t2, by = c("ID", "year"), all = TRUE)

    # 강수량(봄)
    mam_p1 <- aggregate(tp ~ ID + year, data = sp, mean)
    mam_p2 <- aggregate(MON_PRE ~ ID + year, data = sp, mean)
    mam_p <- merge(mam_p1, mam_p2, by = c("ID", "year"), all = TRUE)

    # 평균기온(여름)
    jja_t1 <- aggregate(t2m ~ ID + year, data = su, mean)
    jja_t2 <- aggregate(AVG_TEMP ~ ID + year, data = su, mean)
    jja_t <- merge(jja_t1, jja_t2, by = c("ID", "year"), all = TRUE)

    # 강수량(여름)
    jja_p1 <- aggregate(tp ~ ID + year, data = su, mean)
    jja_p2 <- aggregate(MON_PRE ~ ID + year, data = su, mean)
    jja_p <- merge(jja_p1, jja_p2, by = c("ID", "year"), all = TRUE)

    # 평균기온(가을)
    son_t1 <- aggregate(t2m ~ ID + year, data = au, mean)
    son_t2 <- aggregate(AVG_TEMP ~ ID + year, data = au, mean)
    son_t <- merge(son_t1, son_t2, by = c("ID", "year"), all = TRUE)

    # 강수량(가을)
    son_p1 <- aggregate(tp ~ ID + year, data = au, mean)
    son_p2 <- aggregate(MON_PRE ~ ID + year, data = au, mean)
    son_p <- merge(son_p1, son_p2, by = c("ID", "year"), all = TRUE)

    # 결측값 제거
    djf_t <- djf_t[complete.cases(djf_t[, c("AVG_TEMP")]), ]
    mam_t <- mam_t[complete.cases(mam_t[, c("AVG_TEMP")]), ]
    jja_t <- jja_t[complete.cases(jja_t[, c("AVG_TEMP")]), ]
    son_t <- son_t[complete.cases(son_t[, c("AVG_TEMP")]), ]

    djf_p <- djf_p[complete.cases(djf_p[, c("MON_PRE")]), ]
    mam_p <- mam_p[complete.cases(mam_p[, c("MON_PRE")]), ]
    jja_p <- jja_p[complete.cases(jja_p[, c("MON_PRE")]), ]
    son_p <- son_p[complete.cases(son_p[, c("MON_PRE")]), ]

    # 겨울
    season <- c("DJF", "MAM", "JJA", "SON")

    cor_wi_t1 <- cor.test(djf_t$t2m, djf_t$AVG_TEMP, method = "pearson")
    cor_wi_t2 <- cor.test(djf_t$t2m, djf_t$AVG_TEMP, method = "spearman")
    cor_wi_t3 <- cor.test(djf_t$t2m, djf_t$AVG_TEMP, method = "kendall", exact = FALSE)

    cor_wi_p1 <- cor.test(djf_p$tp, djf_p$MON_PRE, method = "pearson")
    cor_wi_p2 <- cor.test(djf_p$tp, djf_p$MON_PRE, method = "spearman")
    cor_wi_p3 <- cor.test(djf_p$tp, djf_p$MON_PRE, method = "kendall", exact = FALSE)

    resultli[[k]][c(7:9), 1] <- paste0("avg_temp.", season[1])
    resultli[[k]][7, 3] <- cor_wi_t1[4] # correlation coefficient
    resultli[[k]][7, 4] <- cor_wi_t1[3] # p-value
    resultli[[k]][7, 5] <- cor_wi_t1[2] # df
    resultli[[k]][8, 3] <- cor_wi_t2[4] # correlation coefficient
    resultli[[k]][8, 4] <- cor_wi_t2[3] # p-value
    resultli[[k]][9, 3] <- cor_wi_t3[4] # correlation coefficient
    resultli[[k]][9, 4] <- cor_wi_t3[3] # p-value

    resultli[[k]][c(10:12), 1] <- paste0("pre.", season[1])
    resultli[[k]][10, 3] <- cor_wi_p1[4] # correlation coefficient
    resultli[[k]][10, 4] <- cor_wi_p1[3] # p-value
    resultli[[k]][10, 5] <- cor_wi_p1[2] # df
    resultli[[k]][11, 3] <- cor_wi_p2[4] # correlation coefficient
    resultli[[k]][11, 4] <- cor_wi_p2[3] # p-value
    resultli[[k]][12, 3] <- cor_wi_p3[4] # correlation coefficient
    resultli[[k]][12, 4] <- cor_wi_p3[3] # p-value

    resultli[[k]][c(7:12), 2] <- c("pearson", "spearman", "kendall")

    # 봄
    cor_sp_t1 <- cor.test(mam_t$t2m, mam_t$AVG_TEMP, method = "pearson")
    cor_sp_t2 <- cor.test(mam_t$t2m, mam_t$AVG_TEMP, method = "spearman")
    cor_sp_t3 <- cor.test(mam_t$t2m, mam_t$AVG_TEMP, method = "kendall", exact = FALSE)

    cor_sp_p1 <- cor.test(mam_p$tp, mam_p$MON_PRE, method = "pearson")
    cor_sp_p2 <- cor.test(mam_p$tp, mam_p$MON_PRE, method = "spearman")
    cor_sp_p3 <- cor.test(mam_p$tp, mam_p$MON_PRE, method = "kendall", exact = FALSE)

    resultli[[k]][c(13:15), 1] <- paste0("avg_temp.", season[2])
    resultli[[k]][13, 3] <- cor_sp_t1[4] # correlation coefficient
    resultli[[k]][13, 4] <- cor_sp_t1[3] # p-value
    resultli[[k]][13, 5] <- cor_sp_t1[2] # df
    resultli[[k]][14, 3] <- cor_sp_t2[4] # correlation coefficient
    resultli[[k]][14, 4] <- cor_sp_t2[3] # p-value
    resultli[[k]][15, 3] <- cor_sp_t3[4] # correlation coefficient
    resultli[[k]][15, 4] <- cor_sp_t3[3] # p-value

    resultli[[k]][c(16:18), 1] <- paste0("pre.", season[2])
    resultli[[k]][16, 3] <- cor_sp_p1[4] # correlation coefficient
    resultli[[k]][16, 4] <- cor_sp_p1[3] # p-value
    resultli[[k]][16, 5] <- cor_sp_p1[2] # df
    resultli[[k]][17, 3] <- cor_sp_p2[4] # correlation coefficient
    resultli[[k]][17, 4] <- cor_sp_p2[3] # p-value
    resultli[[k]][18, 3] <- cor_sp_p3[4] # correlation coefficient
    resultli[[k]][18, 4] <- cor_sp_p3[3] # p-value

    resultli[[k]][c(13:18), 2] <- c("pearson", "spearman", "kendall")

    # 여름
    cor_su_t1 <- cor.test(jja_t$t2m, jja_t$AVG_TEMP, method = "pearson")
    cor_su_t2 <- cor.test(jja_t$t2m, jja_t$AVG_TEMP, method = "spearman")
    cor_su_t3 <- cor.test(jja_t$t2m, jja_t$AVG_TEMP, method = "kendall", exact = FALSE)

    cor_su_p1 <- cor.test(jja_p$tp, jja_p$MON_PRE, method = "pearson")
    cor_su_p2 <- cor.test(jja_p$tp, jja_p$MON_PRE, method = "spearman")
    cor_su_p3 <- cor.test(jja_p$tp, jja_p$MON_PRE, method = "kendall", exact = FALSE)

    resultli[[k]][c(19:21), 1] <- paste0("avg_temp.", season[3])
    resultli[[k]][19, 3] <- cor_su_t1[4] # correlation coefficient
    resultli[[k]][19, 4] <- cor_su_t1[3] # p-value
    resultli[[k]][19, 5] <- cor_su_t1[2] # df
    resultli[[k]][20, 3] <- cor_su_t2[4] # correlation coefficient
    resultli[[k]][20, 4] <- cor_su_t2[3] # p-value
    resultli[[k]][21, 3] <- cor_su_t3[4] # correlation coefficient
    resultli[[k]][21, 4] <- cor_su_t3[3] # p-value

    resultli[[k]][c(22:24), 1] <- paste0("pre.", season[3])
    resultli[[k]][22, 3] <- cor_su_p1[4] # correlation coefficient
    resultli[[k]][22, 4] <- cor_su_p1[3] # p-value
    resultli[[k]][22, 5] <- cor_su_p1[2] # df
    resultli[[k]][23, 3] <- cor_su_p2[4] # correlation coefficient
    resultli[[k]][23, 4] <- cor_su_p2[3] # p-value
    resultli[[k]][24, 3] <- cor_su_p3[4] # correlation coefficient
    resultli[[k]][24, 4] <- cor_su_p3[3] # p-value

    resultli[[k]][c(19:24), 2] <- c("pearson", "spearman", "kendall")

    # 가을
    cor_au_t1 <- cor.test(son_t$t2m, son_t$AVG_TEMP, method = "pearson")
    cor_au_t2 <- cor.test(son_t$t2m, son_t$AVG_TEMP, method = "spearman")
    cor_au_t3 <- cor.test(son_t$t2m, son_t$AVG_TEMP, method = "kendall", exact = FALSE)

    cor_au_p1 <- cor.test(son_p$tp, son_p$MON_PRE, method = "pearson")
    cor_au_p2 <- cor.test(son_p$tp, son_p$MON_PRE, method = "spearman")
    cor_au_p3 <- cor.test(son_p$tp, son_p$MON_PRE, method = "kendall", exact = FALSE)

    resultli[[k]][c(25:27), 1] <- paste0("avg_temp.", season[4])
    resultli[[k]][25, 3] <- cor_au_t1[4] # correlation coefficient
    resultli[[k]][25, 4] <- cor_au_t1[3] # p-value
    resultli[[k]][25, 5] <- cor_au_t1[2] # df
    resultli[[k]][26, 3] <- cor_au_t2[4] # correlation coefficient
    resultli[[k]][26, 4] <- cor_au_t2[3] # p-value
    resultli[[k]][27, 3] <- cor_au_t3[4] # correlation coefficient
    resultli[[k]][27, 4] <- cor_au_t3[3] # p-value

    resultli[[k]][c(28:30), 1] <- paste0("pre.", season[4])
    resultli[[k]][28, 3] <- cor_au_p1[4] # correlation coefficient
    resultli[[k]][28, 4] <- cor_au_p1[3] # p-value
    resultli[[k]][28, 5] <- cor_au_p1[2] # df
    resultli[[k]][29, 3] <- cor_au_p2[4] # correlation coefficient
    resultli[[k]][29, 4] <- cor_au_p2[3] # p-value
    resultli[[k]][30, 3] <- cor_au_p3[4] # correlation coefficient
    resultli[[k]][30, 4] <- cor_au_p3[3] # p-value

    resultli[[k]][c(25:30), 2] <- c("pearson", "spearman", "kendall")
  }

  # 지점별 전년(ann)

  for (k in 1:l) {
    temp <- datalist[[k]]
    temp <- subset(temp, year >= S & year <= E)

    # 연도별 계절 평균값 구하기, 만약 한 달이라도 na가 있으면 그 연도는 계산 안됨.
    for (j in 0:3) {
      for (i in 1:nrow(temp)) {
        if (is.na(temp[i, 4 + j])) {
          temp[which(temp[, 2] == temp[i, 2]), 4 + j] <- NA
        }
      }
    }

    ann_t1 <- aggregate(t2m ~ ID + year, data = temp, mean) # 평균기온
    ann_t2 <- aggregate(AVG_TEMP ~ ID + year, data = temp, mean)
    ann_t <- merge(ann_t1, ann_t2, by = c("ID", "year"), all = TRUE)

    ann_p1 <- aggregate(tp ~ ID + year, data = temp, mean) # 강수량
    ann_p2 <- aggregate(MON_PRE ~ ID + year, data = temp, mean)
    ann_p <- merge(ann_p1, ann_p2, by = c("ID", "year"), all = TRUE)

    # 결측값 제거
    ann_t <- ann_t[complete.cases(ann_t[, c("AVG_TEMP")]), ]
    ann_p <- ann_p[complete.cases(ann_p[, c("MON_PRE")]), ]

    cor_t1 <- cor.test(ann_t$t2m, ann_t$AVG_TEMP, method = "pearson")
    cor_t2 <- cor.test(ann_t$t2m, ann_t$AVG_TEMP, method = "spearman")
    cor_t3 <- cor.test(ann_t$t2m, ann_t$AVG_TEMP, method = "kendall", exact = FALSE)

    cor_p1 <- cor.test(ann_p$tp, ann_p$MON_PRE, method = "pearson")
    cor_p2 <- cor.test(ann_p$tp, ann_p$MON_PRE, method = "spearman")
    cor_p3 <- cor.test(ann_p$tp, ann_p$MON_PRE, method = "kendall", exact = FALSE)

    resultli[[k]][c(31:33), 1] <- "avg_temp.ANN"
    resultli[[k]][31, 3] <- cor_t1[4] # correlation coefficient
    resultli[[k]][31, 4] <- cor_t1[3] # p-value
    resultli[[k]][31, 5] <- cor_t1[2] # df
    resultli[[k]][32, 3] <- cor_t2[4] # correlation coefficient
    resultli[[k]][32, 4] <- cor_t2[3] # p-value
    resultli[[k]][33, 3] <- cor_t3[4] # correlation coefficient
    resultli[[k]][33, 4] <- cor_t3[3] # p-value
    resultli[[k]][c(31:33), 2] <- c("pearson", "spearman", "kendall")

    resultli[[k]][c(34:36), 1] <- "pre.ANN"
    resultli[[k]][34, 3] <- cor_p1[4] # correlation coefficient
    resultli[[k]][34, 4] <- cor_p1[3] # p-value
    resultli[[k]][34, 5] <- cor_p1[2] # df
    resultli[[k]][35, 3] <- cor_p2[4] # correlation coefficient
    resultli[[k]][35, 4] <- cor_p2[3] # p-value
    resultli[[k]][36, 3] <- cor_p3[4] # correlation coefficient
    resultli[[k]][36, 4] <- cor_p3[3] # p-value

    resultli[[k]][c(34:36), 2] <- c("pearson", "spearman", "kendall")
  }

  t <- data.frame()

  for (i in 1:l) {
    temp <- as.data.frame(resultli[[i]])
    temp <- transform(temp, station = datalist[[i]][1, 1])
    temp <- temp[, c(6, 1:5)]
    t <- rbind(t, temp)
  }
  return(t)
}

# required package:
# input:
# output:

temp_anl <- function(name, datalist, S, E) { # <- 두개로 나누기(평균값 구하기, 상관분석 하기)

  # 상관분석(평균기온): 지점별, 월별. 켄달, 피어슨, 스피어만 상관분석

  resultli <- list()

  l <- length(name)

  # 지점별 시계열 분석
  for (i in 1:l) {
    mon <- datalist[[i]]
    mon <- subset(mon, year >= S & year <= E) # 분석기간 설정

    mon <- mon[complete.cases(mon[, c("AVG_TEMP", "MON_PRE")]), ]

    cor_t1 <- cor.test(mon$t2m, mon$AVG_TEMP, method = "pearson")
    cor_t2 <- cor.test(mon$t2m, mon$AVG_TEMP, method = "spearman")
    cor_t3 <- cor.test(mon$t2m, mon$AVG_TEMP, method = "kendall")

    cor_p1 <- cor.test(mon$tp, mon$MON_PRE, method = "pearson")
    cor_p2 <- cor.test(mon$tp, mon$MON_PRE, method = "spearman")
    cor_p3 <- cor.test(mon$tp, mon$MON_PRE, method = "kendall")

    bias_t <- bias(mon$AVG_TEMP, mon$t2m)
    rmse_t <- rmse(mon$AVG_TEMP, mon$t2m)

    # wet_mon <- mon[which(mon$MON_PRE != 0),]
    bias_p <- bias(mon$MON_PRE, mon$tp)
    rmse_p <- rmse(mon$MON_PRE, mon$tp)

    result <- data.frame("variable", "cor.test", "coefficient", "p-value", "df")

    n <- 1

    result[c((-4 + (5 * n)):(5 * n)), 1] <- "avg_temp.MON"
    result[-4 + (5 * n), 3] <- cor_t1[4] # correlation coefficient(pearson)
    result[-4 + (5 * n), 4] <- cor_t1[3] # p-value
    result[-3 + (5 * n), 3] <- cor_t2[4] # correlation coefficient(spearman)
    result[-3 + (5 * n), 4] <- cor_t2[3] # p-value
    result[-2 + (5 * n), 3] <- cor_t3[4] # correlation coefficient(kendall)
    result[-2 + (5 * n), 4] <- cor_t3[3] # p-value
    result[-4 + (5 * n), 5] <- cor_t1[2] # df
    result[-1 + (5 * n), 3] <- bias_t # Bias
    result[(5 * n), 3] <- rmse_t # RMSE

    result[c((1 + (5 * n)):(5 + 5 * n)), 1] <- "pre.MON"
    result[1 + (5 * n), 3] <- cor_p1[4] # correlation coefficient(pearson)
    result[1 + (5 * n), 4] <- cor_p1[3] # p-value
    result[2 + (5 * n), 3] <- cor_p2[4] # correlation coefficient(spearman)
    result[2 + (5 * n), 4] <- cor_p2[3] # p-value
    result[3 + (5 * n), 3] <- cor_p3[4] # correlation coefficient(kendall)
    result[3 + (5 * n), 4] <- cor_p3[3] # p-value
    result[1 + (5 * n), 5] <- cor_p1[2] # df
    result[4 + (5 * n), 3] <- bias_p # Bias
    result[5 + (5 * n), 3] <- rmse_p # RMSE

    result[, 2] <- c("pearson", "spearman", "kendall", "bias", "rmse")
    resultli[[i]] <- result
  }

  #-------------------------------------------------------------------------------
  # 지점별 계절별

  for (k in 1:l) {
    # 겨울의 경우, 전년 12월이 들어가므로, 계산의 편리를 위해 전년을 다음해로 바꾸기
    seasonal <- datalist[[k]]

    seasonal <- subset(seasonal, year >= S & year <= E)

    # 전년도 12월, 1월, 2월이 한쌍이 되도록 해야함... => 12월은 연도를 +1 로 하자
    seasonal[which(seasonal$month == 12), 2] <- seasonal[which(seasonal$month == 12), 2] + 1

    wi <- subset(seasonal, month >= 12 | month <= 2)
    sp <- subset(seasonal, month >= 3 & month <= 5)
    su <- subset(seasonal, month >= 6 & month <= 8)
    au <- subset(seasonal, month >= 9 & month <= 11)

    wssa <- list(wi, sp, su, au)

    # 연도별 계절 평균값 구하기, 만약 한 달이라도 na가 있으면 그 연도는 계산 안됨.

    for (n in 1:4) {
      ssn <- data.frame(wssa[n])
      for (j in 0:3) {
        for (i in 1:nrow(ssn)) {
          if (is.na(ssn[i, 4 + j])) {
            ssn[which(ssn[, 2] == ssn[i, 2]), 4 + j] <- NA
          }
        }
      }
      # 평균기온(계절)
      ssn_t1 <- aggregate(t2m ~ ID + year, data = ssn, mean)
      ssn_t2 <- aggregate(AVG_TEMP ~ ID + year, data = ssn, mean)
      ssn_t <- merge(ssn_t1, ssn_t2, by = c("ID", "year"), all = TRUE)

      ssn_t1 <- aggregate(t2m ~ ID + year, data = ssn, mean)
      ssn_t2 <- aggregate(AVG_TEMP ~ ID + year, data = ssn, mean)

      # 강수량(계절)
      ssn_p1 <- aggregate(tp ~ ID + year, data = ssn, mean)
      ssn_p2 <- aggregate(MON_PRE ~ ID + year, data = ssn, mean)
      ssn_p <- merge(ssn_p1, ssn_p2, by = c("ID", "year"), all = TRUE)

      ssn_t <- ssn_t[complete.cases(ssn_t[, c("AVG_TEMP")]), ]
      ssn_p <- ssn_p[complete.cases(ssn_p[, c("MON_PRE")]), ]

      # 계절
      season <- c("DJF", "MAM", "JJA", "SON")

      # 기온
      cor_ssn_t1 <- cor.test(ssn_t$t2m, ssn_t$AVG_TEMP, method = "pearson")
      cor_ssn_t2 <- cor.test(ssn_t$t2m, ssn_t$AVG_TEMP, method = "spearman")
      cor_ssn_t3 <- cor.test(ssn_t$t2m, ssn_t$AVG_TEMP,
        method = "kendall",
        exact = FALSE
      )

      bias_ssn_t <- bias(ssn_t$AVG_TEMP, ssn_t$t2m)
      rmse_ssn_t <- rmse(ssn_t$AVG_TEMP, ssn_t$t2m)

      # 강수
      cor_ssn_p1 <- cor.test(ssn_p$tp, ssn_p$MON_PRE, method = "pearson")
      cor_ssn_p2 <- cor.test(ssn_p$tp, ssn_p$MON_PRE, method = "spearman")
      cor_ssn_p3 <- cor.test(ssn_p$tp, ssn_p$MON_PRE,
        method = "kendall",
        exact = FALSE
      )

      # wet_ssn_p <- ssn_p[which(ssn_p$MON_PRE != 0),]

      bias_ssn_p <- bias(ssn_p$MON_PRE, ssn_p$tp)
      rmse_ssn_p <- rmse(ssn_p$MON_PRE, ssn_p$tp)

      resultli[[k]][c(1 + (5 * (2 * n)):(5 + (5 * (2 * n)))), 1] <- paste0("avg_temp.", season[n])
      resultli[[k]][1 + (5 * (2 * n)), 3] <- cor_ssn_t1[4] # correlation coefficient
      resultli[[k]][1 + (5 * (2 * n)), 4] <- cor_ssn_t1[3] # p-value
      resultli[[k]][1 + (5 * (2 * n)), 5] <- cor_ssn_t1[2] # df
      resultli[[k]][2 + (5 * (2 * n)), 3] <- cor_ssn_t2[4] # correlation coefficient
      resultli[[k]][2 + (5 * (2 * n)), 4] <- cor_ssn_t2[3] # p-value
      resultli[[k]][3 + (5 * (2 * n)), 3] <- cor_ssn_t3[4] # correlation coefficient
      resultli[[k]][3 + (5 * (2 * n)), 4] <- cor_ssn_t3[3] # p-value
      resultli[[k]][4 + (5 * (2 * n)), 3] <- bias_ssn_t
      resultli[[k]][5 + (5 * (2 * n)), 3] <- rmse_ssn_t

      resultli[[k]][c(1 + (5 * (2 * n + 1)):(5 + (5 * (2 * n + 1)))), 1] <- paste0("pre.", season[n])
      resultli[[k]][1 + (5 * (2 * n + 1)), 3] <- cor_ssn_p1[4] # correlation coefficient
      resultli[[k]][1 + (5 * (2 * n + 1)), 4] <- cor_ssn_p1[3] # p-value
      resultli[[k]][1 + (5 * (2 * n + 1)), 5] <- cor_ssn_p1[2] # df
      resultli[[k]][2 + (5 * (2 * n + 1)), 3] <- cor_ssn_p2[4] # correlation coefficient
      resultli[[k]][2 + (5 * (2 * n + 1)), 4] <- cor_ssn_p2[3] # p-value
      resultli[[k]][3 + (5 * (2 * n + 1)), 3] <- cor_ssn_p3[4] # correlation coefficient
      resultli[[k]][3 + (5 * (2 * n + 1)), 4] <- cor_ssn_p3[3] # p-value
      resultli[[k]][4 + (5 * (2 * n + 1)), 3] <- bias_ssn_p
      resultli[[k]][5 + (5 * (2 * n + 1)), 3] <- rmse_ssn_p

      resultli[[k]][c((1 + (5 * (2 * n))):(5 + (5 * (2 * n + 1)))), 2] <- c(
        "pearson", "spearman",
        "kendall", "bias",
        "rmse"
      )
    }
  }

  # 지점별 전년(ann)

  for (k in 1:l) {
    annual <- datalist[[k]]
    annual <- subset(annual, year >= S & year <= E)

    # 연도별 계절 평균값 구하기, 만약 한 달이라도 na가 있으면 그 연도는 계산 안됨.
    for (j in 0:3) {
      for (i in 1:nrow(annual)) {
        if (is.na(annual[i, 4 + j])) {
          annual[which(annual[, 2] == annual[i, 2]), 4 + j] <- NA
        }
      }
    }

    ann_t1 <- aggregate(t2m ~ ID + year, data = annual, mean) # 평균기온
    ann_t2 <- aggregate(AVG_TEMP ~ ID + year, data = annual, mean)
    ann_t <- merge(ann_t1, ann_t2, by = c("ID", "year"), all = TRUE)

    ann_p1 <- aggregate(tp ~ ID + year, data = annual, mean) # 강수량
    ann_p2 <- aggregate(MON_PRE ~ ID + year, data = annual, mean)
    ann_p <- merge(ann_p1, ann_p2, by = c("ID", "year"), all = TRUE)

    # 결측값 제거
    ann_t <- ann_t[complete.cases(ann_t[, c("AVG_TEMP")]), ]
    ann_p <- ann_p[complete.cases(ann_p[, c("MON_PRE")]), ]

    cor_ann_t1 <- cor.test(ann_t$t2m, ann_t$AVG_TEMP, method = "pearson")
    cor_ann_t2 <- cor.test(ann_t$t2m, ann_t$AVG_TEMP, method = "spearman")
    cor_ann_t3 <- cor.test(ann_t$t2m, ann_t$AVG_TEMP,
      method = "kendall",
      exact = FALSE
    )
    bias_ann_t <- bias(ann_t$AVG_TEMP, ann_t$t2m)
    rmse_ann_t <- rmse(ann_t$AVG_TEMP, ann_t$t2m)
    cor_ann_p1 <- cor.test(ann_p$tp, ann_p$MON_PRE, method = "pearson")
    cor_ann_p2 <- cor.test(ann_p$tp, ann_p$MON_PRE, method = "spearman")
    cor_ann_p3 <- cor.test(ann_p$tp, ann_p$MON_PRE,
      method = "kendall",
      exact = FALSE
    )

    # wet_ann_p <- ann_p[which(ann_p$MON_PRE != 0),]
    bias_ann_p <- bias(ann_p$MON_PRE, ann_p$tp)
    rmse_ann_p <- rmse(ann_p$MON_PRE, ann_p$tp)

    n <- 5

    resultli[[k]][c((1 + (5 * (2 * n))):(5 + (5 * (2 * n)))), 1] <- "avg_temp.ANN"
    resultli[[k]][1 + (5 * (2 * n)), 3] <- cor_ann_t1[4] # correlation coefficient
    resultli[[k]][1 + (5 * (2 * n)), 4] <- cor_ann_t1[3] # p-value
    resultli[[k]][1 + (5 * (2 * n)), 5] <- cor_ann_t1[2] # df
    resultli[[k]][2 + (5 * (2 * n)), 3] <- cor_ann_t2[4] # correlation coefficient
    resultli[[k]][2 + (5 * (2 * n)), 4] <- cor_ann_t2[3] # p-value
    resultli[[k]][3 + (5 * (2 * n)), 3] <- cor_ann_t3[4] # correlation coefficient
    resultli[[k]][3 + (5 * (2 * n)), 4] <- cor_ann_t3[3] # p-value
    resultli[[k]][4 + (5 * (2 * n)), 3] <- bias_ann_t
    resultli[[k]][5 + (5 * (2 * n)), 3] <- rmse_ann_t

    resultli[[k]][c((1 + (5 * (2 * n + 1))):(5 + (5 * (2 * n + 1)))), 1] <- "pre.ANN"
    resultli[[k]][1 + (5 * (2 * n + 1)), 3] <- cor_ann_p1[4] # correlation coefficient
    resultli[[k]][1 + (5 * (2 * n + 1)), 4] <- cor_ann_p1[3] # p-value
    resultli[[k]][1 + (5 * (2 * n + 1)), 5] <- cor_ann_p1[2] # df
    resultli[[k]][2 + (5 * (2 * n + 1)), 3] <- cor_ann_p2[4] # correlation coefficient
    resultli[[k]][2 + (5 * (2 * n + 1)), 4] <- cor_ann_p2[3] # p-value
    resultli[[k]][3 + (5 * (2 * n + 1)), 3] <- cor_ann_p3[4] # correlation coefficient
    resultli[[k]][3 + (5 * (2 * n + 1)), 4] <- cor_ann_p3[3] # p-value
    resultli[[k]][4 + (5 * (2 * n + 1)), 3] <- bias_ann_p
    resultli[[k]][5 + (5 * (2 * n + 1)), 3] <- rmse_ann_p

    resultli[[k]][c((1 + (5 * (2 * n))):(5 + (5 * (2 * n + 1)))), 2] <- c(
      "pearson", "spearman",
      "kendall", "bias",
      "rmse"
    )
  }

  t <- data.frame()

  for (i in 1:l) {
    temp <- as.data.frame(resultli[[i]])
    temp <- transform(temp, station = name[i])
    temp <- temp[, c(6, 1:5)]
    t <- rbind(t, temp)
  }
  return(t)
}

# figure 2

mean_map <- function(data) {
  map <- ggplot() +
    borders(fill = "grey99") +
    geom_sf(data = data, aes(color = AvgValue), size = 1.5, shape = 19) +
    coord_sf(xlim = c(109.5, 146), ylim = c(19.5, 55), expand = FALSE)

  # geom_point

  if (data$var[1] == "temp") {
    maxtmp <- max(data[which(data$var == "temp"), ]$AvgValue)
    mintmp <- min(data[which(data$var == "temp"), ]$AvgValue)

    map <- map +
      scale_color_gradientn(
        colors = rev(temp_div1), oob = squish,
        name = "Temperature\n(°C)", limits = c(-25, 30),
        na.value = "grey50", aesthetics = "color"
      )
  } else if (data$var[1] == "pre") {
    maxpre <- max(data[which(data$var == "pre"), ]$AvgValue)

    map <- map +
      scale_color_gradientn(
        name = "Precipitation\n(mm/day)", oob = squish,
        colors = prec_seq9, aesthetics = "color",
        na.value = "grey50", limits = c(0, 10),
        breaks = seq(0, 10, 2)
      )
  }

  map <- map +
    borders(fill = NA) +
    theme_bw() +
    guides(color = guide_colourbar(
      barheight = 15, nbin = 30, raster = FALSE,
      frame.colour = "black",
      ticks.colour = "black",
      ticks.linewidth = 0.25,
    )) +
    scale_x_continuous(breaks = seq(110, 140, 10)) +
    scale_y_continuous(breaks = seq(20, 50, 10)) +
    labs(x = "Longitude (°E)", y = "Latituide (°N)") +
    facet_wrap(~season, ncol = 5) +
    scalebar(data,
      dist = 500, dist_unit = "km", model = "WGS84", transform = TRUE,
      st.size = 2.5,
      st.bottom = TRUE,
      border.size = 0.5,
      # location = "bottomright",
      anchor = c(x = 143, y = 21),
      facet.var = c("season"), facet.lev = c("ANN")
    ) +
    ggtitle(ifelse(data$var[1] == "pre", "(b) Precipitation", "(a) Temperature")) +
    theme(
      legend.direction = "vertical",
      legend.title = element_text(size = 10),
      legend.position = "right",
      # legend.text = element_text(size = 14 , face = 'bold'),
      plot.title = element_text(
        size = 18, color = "black",
        face = "bold"
      ),
      panel.background = element_rect(fill = "azure2"),
      strip.background.x = element_blank(), strip.text.x = element_blank(),
      axis.title = element_text(),
      axis.text = element_text()
    ) +
    geom_text(
      data = data.frame(
        season = unique(data$season),
        label = unique(data$season)
      ),
      inherit.aes = FALSE, size = 5,
      aes(x = 114, y = 53.5, label = label, fontface = "bold")
    )

  return(map)
}


# figure 3, 4, 5
# inset plot ---------------------------------------------------------------
get_subplot2 <- function(data) {

  # r <- max(data.frame(table(cut(data$coefficient, seq(0,1,0.2))))$Freq) * 1.2
  # r <- round(r)
  # t <- 650
  # p <- 500

  scl <- ifelse(data$var[1] == "avg_temp", 537 * 1.15, 537 * 0.85)
  hvar <- ifelse(data$var[1] == "pre", 0.45, 0.5)

  brks <- c(-0.2, seq(0, 1, 0.2))

  Subplot <- ggplot(data = data, aes(x = coefficient)) +
    geom_histogram(
      bins = 7, breaks = brks,
      color = "#000000", fill = "grey90"
    ) +
    geom_text(aes(label = paste0(round((..count..) / 537 * 100, 1), "%")),
      vjust = -0.5, hjust = hvar,
      color = "black", size = 1.7, stat = "bin", breaks = brks
    ) +
    scale_x_continuous(
      limits = c(-0.2, 1), breaks = brks, oob = squish,
      label = c(-0.5, seq(0, 1, 0.2))
    ) +
    scale_y_continuous(
      limits = c(0, scl), breaks = seq(0, 537, 537 / 4),
      label = paste0(seq(0, 100, 25), "%")
    ) +
    theme_pubclean() +
    theme( # panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.title = element_blank(),
      axis.text.y = element_text(size = 5, colour = "black"),
      axis.text.x = element_text(size = 5, colour = "black", vjust = 5),
      axis.ticks = element_blank()
    ) +
    facet_grid(season ~ data) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      # axis.text.x =  element_blank(),
      # axis.text.y = element_text(size = 10, colour = 'black'),
      axis.ticks.x = element_blank(),
      plot.margin = unit(c(0, 0, -0.1, 0), "cm")
    )

  return(Subplot)
}

# annotation function ------------------------------------------------------
annotation_custom2 <- function(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
  layer(
    data = data, stat = StatIdentity, position = PositionIdentity,
    geom = ggplot2:::GeomCustomAnn,
    inherit.aes = TRUE, params = list(
      grob = grob,
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax
    )
  )
}

# analysis function --------------------------------------------------------
anl_map <- function(data, subplot = 0) {
  map <- ggplot(data = data) +
    borders(fill = "grey99")

  if (data$method[[1]] %in% c("pearson", "spearman", "kendall")) {
    map <- map + geom_sf(
      data = data, aes(
        fill = coefficient,
        color = significance
      ),
      size = 2.25, stroke = 0.5, shape = 21
    ) +
      scale_color_manual(
        values = c("transparent", "grey50"),
        aesthetics = "color"
      ) +
      guides(color = "none")
    if (data$method[[1]] == "pearson") {
      map <- map +
        scale_fill_gradient(
          name = "R-values", aesthetics = "fill", space = "Lab",
          limits = c(0, 1), oob = squish,
          high = "red", low = "white",
          label = function(x) sprintf("%.1f", x),
          guide = guide_colourbar(
          barheight = 16, nbin = 30, raster = FALSE,
          frame.colour = "black", ticks.colour = "black"
        ))
    } else if (data$method[[1]] == "spearman") {
      map <- map +
        scale_fill_gradient(
          name = "Spearman \n coefficient \n (ρ)",
          aesthetics = "fill", space = "Lab",
          limits = c(0, 1), oob = squish,
          high = "red", low = "white",
          label = function(x) sprintf("%.1f", x),
          guide = guide_colourbar(
          barheight = 16, nbin = 30, raster = FALSE,
          frame.colour = "black", ticks.colour = "black"
        )
        )
    } else if (data$method[[1]] == "kendall") {
      map <- map +
        scale_fill_gradient(
          name = "Kendall \n coefficient \n (τ)",
          aesthetics = "fill", space = "Lab",
          limits = c(0, 1), oob = squish,
          high = "red", low = "white",
          label = function(x) sprintf("%.1f", x),
          guide = guide_colourbar(
          barheight = 16, nbin = 30, raster = FALSE,
          frame.colour = "black", ticks.colour = "black"
        )
        )
    }
  } else if (data$method[[1]] == "bias" & data$var[[i]] == "pre") {
    lim_pre <- ifelse(data$season[1] == "JJA", 3.0, 1.5)

    map <- map + geom_sf(
      data = data, aes(color = coefficient),
      size = 2.5, shape = 19
    ) +
      scale_color_gradientn(
        name = "Bias\n(mm/day)", guide = "colorbar",
        na.value = "#660000", aesthetics = "color",
        limits = c(-1, 1) * lim_pre,
        label = function(x) sprintf("%.1f", x),
        colors = prec_div, oob = squish
      )
  } else if (data$method[[1]] == "bias" & data$var[[i]] == "avg_temp") {
    map <- map + geom_sf(
      data = data, aes(color = coefficient),
      size = 2.5, shape = 19
    ) +
      scale_color_gradientn(
        name = "Bias\n(°C)", guide = "colorbar",
        na.value = "#660000", aesthetics = "color",
        limits = c(-1, 1) * 5,
        colors = rev(temp_div1), oob = squish
      )
  } else if (data$method[[1]] == "rmse") {
    lim_pre <- ifelse(data$season[1] == "JJA", 4.0, 2)

    map <- map + geom_sf(
      data = data, aes(color = coefficient),
      size = 2.5, shape = 19
    ) +
      scale_color_gradient(
        name = ifelse(data$var[1] == "pre",
          "RMSE\n(mm/day)", "RMSE\n(°C)"
        ),
        guide = "colorbar", aesthetics = "color",
        limits = c(0, ifelse(data$var[1] == "pre",
          lim_pre, 4
        )),
        high = "red", low = "white", oob = squish,
        label = function(x) sprintf("%.1f", x)
      )
  }

  map <- map +
    borders(fill = NA) +
    coord_sf(xlim = c(109.5, 146), ylim = c(19.5, 55), expand = FALSE) +
    scale_x_continuous(breaks = seq(110, 140, 10)) +
    scale_y_continuous(breaks = seq(20, 50, 10)) +
    theme_bw()
  #+ labs(x="Longitude (°E)", y="Latitude (°N)")

  map <- map + facet_grid(season ~ data) # nrow = 2,

  if (data$method[[1]] %in% c("pearson", "spearman", "kendall")) {
    if (subplot == TRUE) {
      inset1 <- data %>%
        split(f = .$data) %>%
        purrr::map(~ annotation_custom2(
          grob = ggplotGrob(get_subplot2(.)),
          data = data.frame(data = unique(.$data)),
          xmin = 128, xmax = 145.7,
          ymin = 48, ymax = 54.7
        ))
      map <- map + inset1
    }
  } else {
    map <- map +
      guides(color = guide_colourbar(
        barheight = 16, nbin = 30, raster = FALSE,
        frame.colour = "black",
        ticks.colour = "black"
      )) +
      theme(
        legend.direction = "vertical",
        legend.title = element_text(size = 10),
        legend.position = "right"
      )
  }

  map <- map +
    scalebar(data,
      dist = 500, dist_unit = "km", model = "WGS84", transform = TRUE,
      st.size = 2.5,
      st.bottom = TRUE,
      border.size = 0.5,
      # location = "bottomright",
      anchor = c(x = 143.5, y = 21),
      facet.var = c("data"),
      facet.lev = c("ERA5")
    )

  map <- map +
    # ggtitle(paste0(data$season)) +
    theme( # plot.title = element_text(size = 15, color = "black", hjust = 0),
      panel.background = element_rect(fill = "azure2"),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_text(face = "bold", size = 14),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    )
  # geom_text(data = data.frame(data = unique(data$data),
  #                            label = unique(data$data)),
  #          aes(x = 114, y = 53.5, label = label, fontface = 'bold'))


  if (data$season[[1]] == "ANN") {
    map <- map +
      theme(strip.text.x = element_text(hjust = 0.5, face = "bold", size = 14))
  }

  test <- (data$season[[1]] == "SON" &
    data$method[[1]] %in% c("pearson", "spearman", "kendall"))
  test2 <- (((data$season[[1]] == "JJA" & data$var[[1]] == "pre") |
    (data$season[[1]] == "SON" & data$var[[1]] == "avg_temp")) &
    data$method[[1]] %in% c("rmse", "bias"))

  if (test | test2) {
    map <- map +
      theme(axis.ticks.x = element_line(), axis.text.x = element_text())
  }

  return(map)
}

# figure 6
# inset plot ---------------------------------------------------------------
get_subplot <- function(data) {
  t <- 450
  p <- 450

  Subplot <- ggplot(data = data, aes(data)) +
    stat_count(color = "#000000", fill = "grey90") +
    geom_text(aes(label = paste0(round((..count..) / 537 * 100, 1), "%")),
      stat = "count", vjust = -0.6, color = "black", size = 2.0
    ) +
    scale_y_continuous(
      limits = c(0, 450), breaks = seq(0, 537, 537 / 4),
      label = paste0(seq(0, 100, 25), "%")
    ) +
    theme_pubclean() +
    theme( # panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.title = element_blank(),
      axis.text.x = element_text(vjust = 5, size = 5, color = "black"),
      axis.text.y = element_text(size = 5, color = "black"),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0, 0, -0.1, 0), "cm")
    ) +
    facet_grid(method ~ season) +
    theme(
      strip.text = element_blank(),
      # strip.text.y = element_text(face = "bold", size = 14),
      # axis.text.x =  element_blank(),
      # axis.text.y = element_text(size = 5, color = 'black'),
      strip.background = element_blank()
    )

  return(Subplot)
}
# annotation function -----------------------------------------------------

annotation_custom2 <- function(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
  layer(
    data = data, stat = StatIdentity, position = PositionIdentity,
    geom = ggplot2:::GeomCustomAnn,
    inherit.aes = TRUE, params = list(
      grob = grob,
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax
    )
  )
}

# best stations 
best_map <- function(data, subplot = NULL) {
  map <- ggplot() +
    borders(fill = "grey99") +
    geom_sf(data = data, aes(fill = data), size = 2, stroke = 0.1, shape = 21) +
    borders(fill = NA) +
    coord_sf(xlim = c(109.5, 146), ylim = c(19.5, 55), expand = FALSE) +
    theme_bw()

  map <- map + scale_fill_manual(
    values = colorset, aesthetics = "fill",
    name = "",
    guide = guide_legend(
      nrow = 1, reverse = FALSE,
      override.aes = list(size = 4)
    )
  )

  map <- map +
    # labs(x="Longitude.(°E)", y="Latitude.(°N)") +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = seq(110, 140, 10)) +
    scale_y_continuous(breaks = seq(20, 50, 10)) +
    theme( # legend.background = element_rect(fill="grey90", size=0.1,
      #                                 linetype="blank", colour ="black"),
      legend.justification = "bottom", legend.position = "bottom",
      legend.key.size = unit(1.5, "cm"),
      legend.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 14, face = "bold")
    )

  map <- map + facet_grid(method ~ season)

  if (subplot == TRUE) {
    inset1 <- data %>%
      split(f = .$season) %>%
      purrr::map(~ annotation_custom2(
        grob = ggplotGrob(get_subplot(.)),
        data = data.frame(season = unique(.$season)),
        xmin = 130, xmax = 145.9,
        ymin = 48, ymax = 54.9
      ))
    map <- map + inset1
  }


  map <- map +
    scalebar(data,
      dist = 500, dist_unit = "km", model = "WGS84", st.size = 2.5,
      transform = TRUE, st.bottom = TRUE,
      border.size = 0.5,
      # location = "bottomright",
      anchor = c(x = 143.4, y = 21),
      facet.var = c("season"), facet.lev = c("ANN")
    )

  map <- map +
    # ggtitle(paste0(ifelse(data$method[i] == 'rmse', 'RMSE',
    #                      stringr::str_to_title(data$method[1])))) +
    theme(
      plot.title = element_text(size = 15, color = "black"),
      panel.background = element_rect(fill = "azure2"),
      strip.background = element_blank(),
      strip.text.y = element_text(hjust = 0.5, face = "bold", size = 14),
      strip.text.x = element_blank(), axis.title = element_blank(),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    ) #+
  # geom_text(data = data.frame(season = unique(data$season),
  #                            label = unique(data$season)),
  # inherit.aes = FALSE,
  #          aes(x = 113, y = 53.5, label = label, fontface = 'bold'))
}

# figure 7
# plot function -----------------------------------------------------------

get_multi <- function(allavg) {
  d.list <- factor(c("ERA5", "JRA55", "NCEP2", "CRU"),
    levels = c("ERA5", "JRA55", "NCEP2", "CRU")
  )

  X <- allavg$avg_TEMP #
  Y <- allavg$avg_t2m #

  ifelse(floor(min(X)) > floor(min(Y)),
    range_min <- floor(min(Y)), range_min <- floor(min(X))
  )
  ifelse(ceiling(max(X)) > ceiling(max(Y)),
    range_max <- ceiling(max(X)), range_max <- ceiling(max(Y))
  )
  range <- range_max - range_min + 2

  data_type <- split(allavg, allavg$data)
  bias <- vector()
  rmse <- vector()

  for (i in 1:4) {
    bias[i] <- round(bias(data_type[[i]]$avg_TEMP, data_type[[i]]$avg_t2m), 3)
    rmse[i] <- round(rmse(data_type[[i]]$avg_TEMP, data_type[[i]]$avg_t2m), 3)
  }
  annot <- data.frame(
    x = range_max - range / 10, y = range_min + range / 4,
    label1 = bias, label2 = rmse, data = d.list,
    sort = FALSE
  )

  allavg$density <- fields::interp.surface(
    MASS::kde2d(allavg$avg_t2m, allavg$avg_TEMP),
    allavg[, c("avg_t2m", "avg_TEMP")]
  )

  plot <- ggplot(data = allavg, aes(x = avg_t2m, y = avg_TEMP)) +
    labs(x = "Reanalysis (°C)", y = "Observation (°C)") +
    geom_point(aes(alpha = 1 / density, shape = country, color = country), size = 2) +
    scale_color_brewer(palette = "Dark2", labels = c(
      "Japan", "South Korea",
      "China", "North Korea"
    )) +
    scale_shape_manual(
      values = c(0, 1, 3, 2),
      labels = c(
        "Japan", "South Korea",
        "China", "North Korea"
      )
    ) +
    scale_alpha(range = c(.4, 1), guide = "none") +
    stat_smooth(method = lm, se = FALSE, linetype = "dashed", color = "black") +
    geom_text(
      data = annot, inherit.aes = FALSE,
      aes(x, y, label = paste("RMSE =", label2, "\n Bias =", label1))
    ) +
    stat_poly_eq(
      parse = TRUE, aes(label = ..eq.label..), formula = y ~ x,
      label.x = "right", label.y = "bottom", coef.digits = 3
    ) +
    coord_cartesian(
      xlim = c(range_min - 0.1, range_max + 1),
      ylim = c(range_min - 0.1, range_max + 1)
    ) +
    geom_abline(intercept = 0, slope = 1) +
    stat_cor(
      method = "pearson", label.x = range_min - 0.1,
      label.y = range_max + 0.2 - (range / 100) * 3,
      cor.coef.name = "R", r.accuracy = 0.001, p.accuracy = 0.001
    ) +
    stat_cor(
      method = "spearman", label.x = range_min - 0.1,
      label.y = range_max + 0.2 - (range / 100) * 8,
      cor.coef.name = "rho", r.accuracy = 0.001, p.accuracy = 0.001
    ) +
    stat_cor(
      method = "kendall", label.x = range_min - 0.1,
      label.y = range_max + 0.2 - (range / 100) * 13,
      cor.coef.name = "tau", r.accuracy = 0.001, p.accuracy = 0.001
    ) +
    # scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    # scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    theme_bw() +
    facet_grid(season ~ data) +
    theme(
      axis.title = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_text(face = "bold", size = 14),
      legend.direction = "horizontal", legend.key.size = unit(1.5, "cm"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold")
    )

  return(plot)
}

get_multi2 <- function(allavg) {
  d.list <- factor(c("ERA5", "JRA55", "NCEP2", "CRU"),
    levels = c("ERA5", "JRA55", "NCEP2", "CRU")
  )

  X <- allavg$avg_PRE
  Y <- allavg$avg_tp

  ifelse(floor(min(X)) > floor(min(Y)),
    range_min <- floor(min(Y)), range_min <- floor(min(X))
  )
  ifelse(ceiling(max(X)) > ceiling(max(Y)),
    range_max <- ceiling(max(X)), range_max <- ceiling(max(Y))
  )
  range <- range_max - range_min + 2

  data_type <- split(allavg, allavg$data)

  bias <- vector()
  rmse <- vector()
  for (i in 1:4) {
    bias[i] <- round(bias(data_type[[i]]$avg_PRE, data_type[[i]]$avg_tp), 3)
    rmse[i] <- round(rmse(data_type[[i]]$avg_PRE, data_type[[i]]$avg_tp), 3)
  }
  annot <- data.frame(
    x = range_max - range / 10, y = range_min + range / 4,
    label1 = bias, label2 = rmse, data = d.list,
    sort = FALSE
  )

  allavg$density <- fields::interp.surface(
    MASS::kde2d(allavg$avg_tp, allavg$avg_PRE),
    allavg[, c("avg_tp", "avg_PRE")]
  )

  plot <- ggplot(data = allavg, aes(x = avg_tp, y = avg_PRE)) +
    labs(x = "Reanalysis (mm/day)", y = "Observation (mm/day)") +
    geom_point(aes(alpha = 1 / density, shape = country, color = country), size = 2) +
    scale_color_brewer(palette = "Dark2", labels = c(
      "Japan", "South Korea",
      "China", "North Korea"
    )) +
    scale_shape_manual(
      values = c(0, 1, 3, 2),
      labels = c(
        "Japan", "South Korea",
        "China", "North Korea"
      )
    ) +
    scale_alpha(range = c(.4, 1), guide = "none") +
    stat_smooth(method = lm, se = FALSE, linetype = "dashed", color = "black") +
    geom_text(
      data = annot, aes(x, y, label = paste("RMSE =", label2, "\n Bias =", label1)),
      inherit.aes = FALSE
    ) +
    stat_poly_eq(
      parse = TRUE, aes(label = ..eq.label..), formula = y ~ x,
      label.x = "right", label.y = "bottom", coef.digits = 3
    ) +
    coord_cartesian(
      xlim = c(range_min - 0.1, range_max + 0.1),
      ylim = c(range_min - 0.1, range_max + 0.1)
    ) +
    geom_abline(intercept = 0, slope = 1) +
    stat_cor(
      method = "pearson", label.x = range_min - 0.1,
      label.y = range_max + 0.2 - (range / 100) * 3,
      cor.coef.name = "R", r.accuracy = 0.001, p.accuracy = 0.001
    ) +
    stat_cor(
      method = "spearman", label.x = range_min - 0.1,
      label.y = range_max + 0.2 - (range / 100) * 8,
      cor.coef.name = "rho", r.accuracy = 0.001, p.accuracy = 0.001
    ) +
    stat_cor(
      method = "kendall", label.x = range_min - 0.1,
      label.y = range_max + 0.2 - (range / 100) * 13,
      cor.coef.name = "tau", r.accuracy = 0.001, p.accuracy = 0.001
    ) +
    # scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    # scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    theme_bw() +
    facet_grid(season ~ data) +
    theme(
      axis.title = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_text(face = "bold", size = 14),
      legend.direction = "horizontal", legend.key.size = unit(1.5, "cm"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold")
    )

  return(plot)
}

# inset plot -------------------------------------------------------------------

get_hist <- function(data) {
  Subplot <- ggplot(data = data, aes(data)) +
    stat_count(color = "#000000", fill = "grey90") +
    geom_text(
      stat = "count",
      aes(label = ..count..), vjust = -0.6, color = "black", size = 2.0
    ) +
    theme_pubclean() +
    theme( # panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
    facet_grid(data ~ season) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black")
    )

  return(Subplot)
}

# bar plot -------------------------------------------------------------------

get_ctr_avg <- function(data) {
  boxplot <- ggplot(data = data, aes(data)) +
    geom_boxplot(aes(x = data, y = coefficient)) +
    stat_summary(aes(x = data, y = coefficient),
      geom = "point", fun = mean,
      shape = 4, color = "red"
    ) +
    facet_grid(method ~ country) +
    theme_pubclean() +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_text(face = "bold", size = 14),
      axis.text.x = element_blank()
    )

  if (data$method[1] %in% c("pearson", "spearman", "kendall")) {
    boxplot <- boxplot +
      scale_y_continuous(limits = c(0.5, 1), oob = squish)
  }


  return(boxplot)
}

# ipcc color guide  -----------------------------------------------------------

temp_div1 <- c(
  "#67001F", "#B2182B", "#D6604D", "#F4A582", "#F4A582", "#F7F7F7",
  "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"
)
temp_div2 <- c(
  "#7f3b08", "#b35806", "#e08214", "#fdb863", "#fee0b6", "#f7f7f7",
  "#d8daeb", "#b2abd2", "#8073ac", "#542788", "#2d004b"
)
prec_div <- c(
  "#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5",
  "#C7EAE5", "#80CDC1", "#35975E", "#01665E", "#003C30"
)

prec_seq5 <- c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c")
prec_seq9 <- c(
  "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6",
  "#2171b5", "#08519c", "#08306b"
)

# country label  --------------------------------------------------------------

c.label <- data.frame(
  country = c("S. Korea", "N. Korea", "China", "Japan"),
  long = c(123.5, 133, 125, 143), lat = c(35, 40, 29, 35)
)

taylor.diagram2 <-
  function(ref, model, add = FALSE, col = "red", pch = 19, pos.cor = TRUE,
           xlab = "Standard deviation", ylab = "", main = "Taylor Diagram",
           show.gamma = TRUE, ngamma = 3, gamma.col = 8, sd.arcs = 0,
           ref.sd = FALSE, sd.method = "sample",
           grad.corr.lines = c(0.2, 0.4, 0.6, 0.8, 0.9),
           pcex = 1, cex.axis = 1, normalize = FALSE,
           mar = c(4, 3, 4, 3), lwd = 1, text, ...) {
    grad.corr.full <- c(
      0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99,
      1
    )
    R <- cor(ref, model, use = "pairwise")
    if (is.list(ref)) {
      ref <- unlist(ref)
    }
    if (is.list(model)) {
      ref <- unlist(model)
    }
    SD <- function(x, subn) {
      meanx <- mean(x, na.rm = TRUE)
      devx <- x - meanx
      ssd <- sqrt(sum(devx * devx, na.rm = TRUE) / (length(x[!is.na(x)]) -
        subn))
      return(ssd)
    }
    subn <- sd.method != "sample"
    sd.r <- SD(ref, subn)
    sd.f <- SD(model, subn)
    if (normalize) {
      sd.f <- sd.f / sd.r
      sd.r <- 1
    }
    maxsd <- 1.5 * max(sd.f, sd.r)
    oldpar <- par("mar", "xpd", "xaxs", "yaxs")
    if (!add) {
      par(mar = mar)
      if (pos.cor) {
        if (nchar(ylab) == 0) {
          ylab <- "Standard deviation"
        }
        plot(0,
          xlim = c(0, maxsd * 1.1), ylim = c(0, maxsd *
            1.1), xaxs = "i", yaxs = "i", axes = FALSE,
          main = main, xlab = "", ylab = ylab, type = "n",
          cex = cex.axis, ...
        )
        mtext(xlab, side = 1, line = 2.3)
        if (grad.corr.lines[1]) {
          for (gcl in grad.corr.lines) {
            lines(c(0, maxsd *
              gcl), c(0, maxsd * sqrt(1 - gcl^2)), lty = 3)
          }
        }
        segments(c(0, 0), c(0, 0), c(0, maxsd), c(
          maxsd,
          0
        ))
        axis.ticks <- pretty(c(0, maxsd))
        axis.ticks <- axis.ticks[axis.ticks <= maxsd]
        axis(1, at = axis.ticks, cex.axis = cex.axis)
        axis(2, at = axis.ticks, cex.axis = cex.axis)
        if (sd.arcs[1]) {
          if (length(sd.arcs) == 1) {
            sd.arcs <- axis.ticks
          }
          for (sdarc in sd.arcs) {
            xcurve <- cos(seq(0, pi / 2, by = 0.03)) * sdarc
            ycurve <- sin(seq(0, pi / 2, by = 0.03)) * sdarc
            lines(xcurve, ycurve, col = "blue", lty = 3)
          }
        }
        if (show.gamma[1]) {
          if (length(show.gamma) > 1) {
            gamma <- show.gamma
          } else {
            gamma <- pretty(c(0, maxsd), n = ngamma)[-1]
          }
          if (gamma[length(gamma)] > maxsd) {
            gamma <- gamma[-length(gamma)]
          }
          labelpos <- seq(45, 70, length.out = length(gamma))
          for (gindex in 1:length(gamma)) {
            xcurve <- cos(seq(0, pi, by = 0.03)) * gamma[gindex] +
              sd.r
            endcurve <- which(xcurve < 0)
            endcurve <- ifelse(length(endcurve), min(endcurve) -
              1, 105)
            ycurve <- sin(seq(0, pi, by = 0.03)) * gamma[gindex]
            maxcurve <- xcurve * xcurve + ycurve * ycurve
            startcurve <- which(maxcurve > maxsd * maxsd)
            startcurve <- ifelse(length(startcurve), max(startcurve) +
              1, 0)
            lines(xcurve[startcurve:endcurve], ycurve[startcurve:endcurve],
              col = gamma.col
            )
            if (xcurve[labelpos[gindex]] > 0) {
              boxed.labels(xcurve[labelpos[gindex]], ycurve[labelpos[gindex]],
                gamma[gindex],
                border = FALSE
              )
            }
          }
        }
        xcurve <- cos(seq(0, pi / 2, by = 0.01)) * maxsd
        ycurve <- sin(seq(0, pi / 2, by = 0.01)) * maxsd
        lines(xcurve, ycurve)
        bigtickangles <- acos(seq(0.1, 0.9, by = 0.1))
        medtickangles <- acos(seq(0.05, 0.95, by = 0.1))
        smltickangles <- acos(seq(0.91, 0.99, by = 0.01))
        segments(cos(bigtickangles) * maxsd, sin(bigtickangles) *
          maxsd, cos(bigtickangles) * 0.97 * maxsd, sin(bigtickangles) *
          0.97 * maxsd)
        par(xpd = TRUE)
        if (ref.sd) {
          xcurve <- cos(seq(0, pi / 2, by = 0.01)) * sd.r
          ycurve <- sin(seq(0, pi / 2, by = 0.01)) * sd.r
          lines(xcurve, ycurve)
        }
        points(sd.r, 0, cex = pcex)
        text(cos(c(bigtickangles, acos(c(0.95, 0.99)))) *
          1.05 * maxsd, sin(c(bigtickangles, acos(c(
          0.95,
          0.99
        )))) * 1.05 * maxsd, c(
          seq(0.1, 0.9, by = 0.1),
          0.95, 0.99
        ), cex = cex.axis)
        text(maxsd * 0.8, maxsd * 0.8, "Correlation",
          srt = 315,
          cex = cex.axis
        )
        segments(cos(medtickangles) * maxsd, sin(medtickangles) *
          maxsd, cos(medtickangles) * 0.98 * maxsd, sin(medtickangles) *
          0.98 * maxsd)
        segments(cos(smltickangles) * maxsd, sin(smltickangles) *
          maxsd, cos(smltickangles) * 0.99 * maxsd, sin(smltickangles) *
          0.99 * maxsd)
      } else {
        x <- ref
        y <- model
        R <- cor(x, y, use = "pairwise.complete.obs")
        E <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
        xprime <- x - mean(x, na.rm = TRUE)
        yprime <- y - mean(y, na.rm = TRUE)
        sumofsquares <- (xprime - yprime)^2
        Eprime <- sqrt(sum(sumofsquares) / length(complete.cases(x)))
        E2 <- E^2 + Eprime^2
        if (add == FALSE) {
          maxray <- 1.5 * max(sd.f, sd.r)
          plot(c(-maxray, maxray), c(0, maxray),
            type = "n",
            asp = 1, bty = "n", xaxt = "n", yaxt = "n",
            xlim = c(-1.1 * maxray, 1.1 * maxray), xlab = xlab,
            ylab = ylab, main = main, cex = cex.axis
          )
          discrete <- seq(180, 0, by = -1)
          listepoints <- NULL
          for (i in discrete) {
            listepoints <- cbind(listepoints, maxray *
              cos(i * pi / 180), maxray * sin(i * pi / 180))
          }
          listepoints <- matrix(listepoints, 2, length(listepoints) / 2)
          listepoints <- t(listepoints)
          lines(listepoints[, 1], listepoints[, 2])
          lines(c(-maxray, maxray), c(0, 0))
          lines(c(0, 0), c(0, maxray))
          for (i in grad.corr.lines) {
            lines(c(0, maxray * i), c(0, maxray * sqrt(1 -
              i^2)), lty = 3)
            lines(c(0, -maxray * i), c(0, maxray * sqrt(1 -
              i^2)), lty = 3)
          }
          for (i in grad.corr.full) {
            text(1.05 * maxray * i, 1.05 * maxray * sqrt(1 -
              i^2), i, cex = cex.axis, adj = cos(i) / 2)
            text(-1.05 * maxray * i, 1.05 * maxray * sqrt(1 -
              i^2), -i, cex = cex.axis, adj = 1 - cos(i) / 2)
          }
          seq.sd <- seq.int(0, 2 * maxray, by = (maxray / 10))[-1]
          for (i in seq.sd) {
            xcircle <- sd.r + (cos(discrete * pi / 180) *
              i)
            ycircle <- sin(discrete * pi / 180) * i
            for (j in 1:length(xcircle)) {
              if ((xcircle[j]^2 + ycircle[j]^2) < (maxray^2)) {
                points(xcircle[j], ycircle[j],
                  col = "darkgreen",
                  pch = "."
                )
                if (j == 10) {
                  text(xcircle[j], ycircle[j], signif(
                    i,
                    2
                  ),
                  cex = cex.axis, col = "darkgreen",
                  srt = 90
                  )
                }
              }
            }
          }
          seq.sd <- seq.int(0, maxray, length.out = 5)
          for (i in seq.sd) {
            xcircle <- cos(discrete * pi / 180) * i
            ycircle <- sin(discrete * pi / 180) * i
            if (i) {
              lines(xcircle, ycircle, lty = 3, col = "blue")
            }
            text(min(xcircle), -0.06 * maxray, signif(
              i,
              2
            ), cex = cex.axis, col = "blue")
            text(max(xcircle), -0.06 * maxray, signif(
              i,
              2
            ), cex = cex.axis, col = "blue")
          }
          text(0, -0.14 * maxray, "Standard Deviation",
            cex = cex.axis, col = "blue"
          )
          text(0, -0.22 * maxray, "Centered RMS Difference",
            cex = cex.axis, col = "darkgreen"
          )
          points(sd.r, 0,
            pch = 22, bg = "darkgreen",
            cex = pcex
          )
          text(0, 1.2 * maxray, "Correlation Coefficient",
            cex = cex.axis
          )
        }
        S <- (2 * (1 + R)) / (sd.f + (1 / sd.f))^2
      }
    }
    points(sd.f * R, sd.f * sin(acos(R)),
      pch = pch, col = col,
      cex = pcex, lwd = lwd
    )
    text(sd.f * R, sd.f * sin(acos(R)),
      col = col,
      labels = text, cex = pcex - 0.2, pos = 3
    )
    invisible(oldpar)
  }
