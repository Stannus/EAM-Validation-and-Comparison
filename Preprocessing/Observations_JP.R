library("lubridate")
library(dplyr)
library(readr)
options(stringsAsFactors = FALSE)

#station data (JMA-AmeDAS) preprocessing

#meta revise
meta <- read.csv("C:/Users/user/Desktop/JMA/meta_e.csv", 
                 fileEncoding = 'utf-8')

#Degree, Minute, Second -> Degree
latm <- meta$latitude_minutes
lonm <- meta$longtitude_minutes

latd <- meta$latitude_degrees
lond <- meta$longtitude_degrees

lat <- latd + latm/60
lon <- lond + lonm/60

meta <- transform(meta, 'lat' = lat)
meta <- transform(meta, 'lon' = lon)

meta <- meta[,c('id', 'lat', 'lon', 'sea_level', 'rem1')]

n_occur <- data.frame(table(meta$id))
dup <- meta[meta$id %in% n_occur$Var1[n_occur$Freq >1],]

dup_list <- levels(factor(dup$id))

unq <- meta[meta$id %in% n_occur$Var1[n_occur$Freq == 1],]

for (i in 1:length(dup_list)){
  temp <- dup[which(dup$id == dup_list[i])[2],]
  unq <- rbind(unq, temp)
}

n_occur <- data.frame(table(unq$id)) #check

write.csv(unq, "C:/Users/user/Desktop/JMA/meta_lonlatel.csv")

#station revise 
station <- read.csv("C:/Users/user/Desktop/JMA/JMA_Station_Data.csv", 
                    fileEncoding = 'utf-8')

date <- strsplit(station[,'date'], "/")

year <- vector()
month <- vector()

for (i in 1:length(date)){
  year[i] <- date[[i]][1]
  month[i] <- date[[i]][2]
}

station <- transform(station, 'year' = as.numeric(year))
station <- transform(station, 'month' = as.numeric(month))


# 격자데이터의 월 강수량 단위(mm/day)와 관측데이터의 월 강수량 단위(mm/month)가 불일치

# -> 관측데이터의 강수량을 해당 월의 일수로 나누어 표준화.
station <- transform(station, length = days_in_month(station$month))

# 월합 강수량을 월평균 강수량으로 변환
station <- transform(station, avg_pre = station$total_pre/station$length)

write.csv(station, "\\\\climatelab\\Minseok\\Comparison\\Station\\all_month_amedas.csv")

#--------------------------------------------------------------------------------------------------------------------
year_avg1 <- aggregate(avg_temp ~ year, data = station, mean)
year_avg2 <- aggregate(avg_pre ~ year, data = station, mean) #월평균 강수량, 월합계 강수량을 구하려면 sum 

year_avg <- merge(year_avg1, year_avg2, by="year", all = TRUE)
colnames(year_avg) <- c("year", "avg_temp", "avg_pre")

write.csv(year_avg, "\\\\climatelab\\Minseok\\Comparison\\Station\\all_year_amedas.csv")
#--------------------------------------------------------------------------------------------------------------------
season_avg1 <- aggregate(avg_temp ~ season, data = station, mean)
season_avg2 <- aggregate(avg_pre ~ season, data = station, mean) #월평균 강수량, 월합계 강수량을 구하려면 sum 

season_avg <- merge(season_avg1, season_avg2, by="season", all = TRUE)
colnames(season_avg) <- c("season", "avg_temp", "mean_pre")

write.csv(season_avg, "\\\\climatelab\\Minseok\\Comparison\\Station\\all_season_amedas.csv")
#--------------------------------------------------------------------------------------------------------------------
