# 제목: South Korea observation preprocessing
# 목적: 재분석 자료의 비교 검증을 위한 ASOS 데이터 가공 및 전처리
# 작성자: 김민석
# 사용 패키지:
library("lubridate")

# 사용 디렉토리:
setwd("\\\\climatelab\\data\\ASOS\\")

# 1981-01-01 이후 기록이 시작된 관측소 제거--------------------------------------------------------

meta <- read.csv("\\\\climatelab\\data\\ASOS\\META_OBSERVATION_STATION_INFO.csv", header = TRUE)
meta <- meta[, c(1:9)]

head(meta)
colnames(meta) <- c("station", "s_date", "e_date", "name", "address", "md", "lat", "lon", "height")
head(meta)

str(meta)

# 공백인("") 데이터들을 na로 변환하여 lubridate 패키지를 사용하여 날짜를 계산할 수 있도록 변환
meta[] <- lapply(meta, function(x) ifelse(nchar(trimws(x)) == 0, NA, x))

# 관측 시작-종료일의 데이터를 날짜형으로 변환
meta$e_date <- as.Date(meta$e_date)
meta$s_date <- as.Date(meta$s_date)

# 1981년 이후 관측 기록이 있는 관측소들중 위치가 변경되어 시작일이 기록된 것들을 제거

# 시작일이 1981년 이후, 종료일이 없는 지점
meta_r <- subset(meta, is.na(e_date))
meta_r <- subset(meta_r, s_date > "1981-01-01")

# 시작일이 1981년 이전, 종료일이 있는 지점
meta_o <- subset(meta, is.na(e_date) == FALSE)
meta_o <- subset(meta_o, s_date < "1981-01-01")

# 시작일이 1981년 이후인 지점들 중 관측지가 이동하여 시작일이 기록된 지점들을 제외
meta_u <- setdiff(meta_r$name, meta_o$name)

# 1981년 이전에 관측기록이 없는 지점들만 추출하여 전체에서 차집합
station_name <- setdiff(meta$name, meta_u)

# 추출된 지점들을 저장
meta_c <- subset(meta, meta$name %in% station_name)
meta_c <- subset(meta_c, is.na(e_date))
meta_c <- subset(meta_c, name != "안동") # 안동의 데이터 기간의 공백(1979~1983)으로 제외함. (시작일과 별도)

write.csv(meta_c, "meta_revised.csv")

#---------------------------------------------------------------------------------------------------

# OBS_ASOS_MNH_1981_2020.csv를 관측 지점 별로 파일 분리---------------------------------------------

df <- read.csv("OBS_ASOS_MNH_1981_2020.csv")

# 앞에서 추출한 관측소의 데이터 추출

meta <- read.csv("meta_revised.csv")

df <- subset(df, df$station %in% meta$station)

colnames(df[4]) <- "avg_temp"
colnames(df[5]) <- "monthly_prep"

# 결측치 확인
sum(is.na(df$avg_temp)) # :0
sum(is.na(df$monthly_prep)) # :20

# 날짜형으로 계산하기 용이하게 변환
df$date <- as.Date(df$date)

# 날짜를 연, 월로 분리
df <- transform(df, year = as.integer(substr(date, 1, 4)))
df <- transform(df, month = as.integer(substr(date, 6, 7)))

# 월자료에서 계절 정보를 추출

season <- vector()

# 현년의 DJF -> 전년의 D + 현년의 JF
for (i in 1:length(df$month)) {
  if (df$month[i] < 3) {
    season[i] <- paste0(df$year[i], "-DJF")
  } else if (df$month[i] > 2 & df$month[i] < 6) {
    season[i] <- paste0(df$year[i], "-MAM")
  } else if (df$month[i] > 5 & df$month[i] < 9) {
    season[i] <- paste0(df$year[i], "-JJA")
  } else if (df$month[i] > 8 & df$month[i] < 12) {
    season[i] <- paste0(df$year[i], "-SON")
  } else if (df$year[i] == 2020) { # 2020년 12월의 경우 다음 관측점의 연-계절로 추출되어 지정해줌.
    season[i] <- "2021-DJF"
  } else if (df$month[i] == 12) {
    season[i] <- paste0(df$year[i + 1], "-DJF")
  }
}

df <- transform(df, season = season)

head(df)

# 격자데이터의 월 강수량 단위(mm/day)와 관측데이터의 월 강수량 단위(mm/month)가 불일치

# -> 관측데이터의 강수량을 해당 월의 일수로 나누어 표준화.
df <- transform(df, length = days_in_month(df$date))

# 월합 강수량을 월평균 강수량으로 변환
df <- transform(df, avg_prep = df$monthly_prep / df$length)

head(df)
str(df)

# 지점별로 데이터 분리
data <- split(df[, -2], df[, 1])

# 지점별로 데이터 저장
station <- meta$station

for (i in 1:length(station)) {
  write.csv(data[[i]], file = paste0("./monthly_data/", station[i], ".csv"))
}

#-------------------------------------------------------------------------------------------------------------------

# 월별 데이터를 년별 데이터로 변환(산술평균)하여 저장

for (k in 1:length(station)) {
  data <- read.csv(paste0("\\\\climatelab\\data\\ASOS\\monthly_data\\", station[k], ".csv"))

  colnames(data) <- c("X", "station", "date", "avg_temp", "monthly_prep", "year", "month", "season", "length", "avg_pre")

  year_avg1 <- aggregate(avg_temp ~ year, data = data, mean)
  year_avg2 <- aggregate(avg_prep ~ year, data = data, mean) # 월평균 강수량, 월합계 강수량을 구하려면 sum

  year_avg <- merge(year_avg1, year_avg2, by = "year", all = TRUE)
  colnames(year_avg) <- c("year", "avg_temp", "avg_prep ")

  write.csv(year_avg, paste0("\\\\climatelab\\data\\ASOS\\yearly_data\\", station[k], ".csv"))
}

head(data)

# 월별 데이터를 계절별 데이터로 변환(산술평균)하여 저장

for (k in 1:length(station)) {
  data <- read.csv(paste0("\\\\climatelab\\data\\ASOS\\monthly_data\\", station[k], ".csv"), header = TRUE)

  colnames(data) <- c("X", "station", "date", "avg_temp", "monthly_prep", "year", "month", "season", "length", "avg_pre")

  season_avg1 <- aggregate(avg_temp ~ season, data = data, mean)
  season_avg2 <- aggregate(avg_prep ~ season, data = data, mean) # 월평균 강수량, 월합계 강수량을 구하려면 sum

  season_avg <- merge(season_avg1, season_avg2, by = "season", all = TRUE)
  colnames(season_avg) <- c("season", "avg_temp", "mean_prep ")

  write.csv(season_avg, paste0("\\\\climatelab\\data\\ASOS\\seasonly_data\\", station[k], ".csv"))
}

head(data)