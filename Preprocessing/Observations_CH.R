
# cma_pre-processing ------------------------------------------------------

meta <- read.csv("\\\\climatelab\\data\\CMA\\meta_revised.csv")

t_dir <- c("\\\\climatelab\\data\\CMA\\Temperature")
p_dir <- c("\\\\climatelab\\data\\CMA\\Precipitation")

colt <- c("ID", "year", "month", "day", "MAX_TEMP", "MIN_TEMP")
colp <- c("ID", "year", "month", "day", "PREC", "TEMP")

data <- data.frame()

# 분석기간을 포함하지 않는 관측지점 확인(1981~2010)  ---------------------------------------

for (i in 1:nrow(meta)){
  temp <- read.table(paste0(t_dir, '\\', meta$province[i], '\\', 
                            meta$province[i], meta$id[i], "temps"), 
                     header = FALSE, col.names = colt)
  
  prec <- read.table(paste0(p_dir, '\\', meta$province[i], '\\', 
                            meta$province[i], meta$id[i]),
                     header = FALSE, col.names = colp)
  
  all <- merge(temp, prec, sort = FALSE, all = TRUE)
  
  l <- nrow(all)
  
  meta[i,6] <- paste(all$year[1], all$month[1], all$day[1], sep = "-")
  meta[i,7] <- paste(all$year[l], all$month[l], all$day[l], sep = "-")
  
  print(i)
}

names(meta)[6:7] <- c("s_date", "e_date")
meta[,6] <- as.Date(meta[,6])
meta[,7] <- as.Date(meta[,7])

rf_meta <- data.frame(subset(rf_meta, (rf_meta$s_date <= as.Date("1981-01-01") & 
                           rf_meta$e_date > as.Date("2011-01-01"))))

# nrow(rf_meta) : [1] 517 ->아래 subset 결과와 동일


# day2month ---------------------------------------------------------------

cnt <- 0

for (i in 1:nrow(meta)){
  #기온
  temp <- read.table(paste0(t_dir, '\\', meta$province[i], '\\', 
                            meta$province[i], meta$id[i], "temps"), 
                     header = FALSE, col.names = colt)
  
  for (k in 1:nrow(temp)){
    if (temp$MAX_TEMP[k] == 32766){
      temp$MAX_TEMP[k] <- NA
    }
    if (temp$MIN_TEMP[k] == 32766){
      temp$MIN_TEMP[k] <- NA
    }
  }
  
  temp[,5:6] <- temp[,5:6]/10 #원단위: 0.1'C
  
  #강수
  prec <- read.table(paste0(p_dir, '\\', meta$province[i], '\\', 
                            meta$province[i], meta$id[i]),
                     header = FALSE, col.names = colp)
  
  prec[,5:6] <- prec[,5:6]/10 #원단위: 0.1 mm
  
  #for(k in 1:nrow(prec)){
  #  if(prec[k,5]==0){
  #    if(is.na(temp[k,5]) | is.na(temp[k,6])){
  #      prec[k,5] <- NA
  #      print(paste(meta$ID[i], temp$year[k], temp$month[k], temp$day[k],
  #                  ": PRE 0 value has replaced to NA"))}}
  #}
  
  tp <- merge(temp, prec, sort = FALSE)
  tp <- subset(tp, (tp[,2] > 1980) & (tp[,2] < 2011), sort = FALSE)
  
  #if (length(levels(factor(tp[,2]))) != 30){
  #  print(meta$id[i])
  #  next}
  
  max <- tp[,c(1:4,5)]
  max <- na.omit(max)
  a <- as.data.frame(table(max$year, max$month), stringsAsFactors = FALSE)
  
  min <- tp[,c(1:4,6)]
  min <- na.omit(min)
  b <- as.data.frame(table(min$year, min$month), stringsAsFactors = FALSE)
  
  sum <- tp[,c(1:4,7)]
  sum <- na.omit(sum)
  c <- as.data.frame(table(sum$year, sum$month), stringsAsFactors = FALSE)
  
  # 관측치 개수
  names(a)[3] <- "max_freq"
  names(b)[3] <- "min_freq"
  names(c)[3] <- "sum_freq"
  
  freq <- merge(a, b)
  freq <- merge(freq, c)
  
  freq$Var1 <- as.numeric(freq$Var1)
  freq$Var2 <- as.numeric(freq$Var2)
  
  sm <- c(2, 4, 6, 9, 11)
  bm <- c(1, 3, 5, 7, 8, 10, 12)
  
  for(k in 1:3){
    freq[,5+k] <- 0
    
    for(l in 1:nrow(freq)){
      if(freq[l,2] %in% bm){
        if(freq[l,2+k]>=25){
          freq[l,5+k] <- 1
        }
      } else {
        if(freq[l,2] %in% sm){
          if(freq[l,2+k]>=24){
            freq[l,5+k] <- 1
          }
        }
      } #윤년만 따로 추출해서 검사 
      for(l in 1:nrow(freq)){
        if(freq[l,1]%%4>0){   #0으로 떨어지면 윤년, 아니면 평년 
          if(freq[l,2]==2){
            if(freq[l,2+k]>=23){
              freq[l,5+k] <- 1
            }}}}}}
  
  names(freq)[c(1:2,6:8)] <- c("year", "month", "max_con", "min_con", "sum_con")
  
  tp$AVG_TEMP <- (tp$MAX_TEMP + tp$MIN_TEMP)/2
  
  month_avg1 <- aggregate(AVG_TEMP ~ year + month, data = tp, mean)  # 사용자 정의함수를 쓰려면 mean 대신에 fun
  month_avg2 <- aggregate(MAX_TEMP ~ year + month, data = tp, mean)
  month_avg3 <- aggregate(MIN_TEMP ~ year + month, data = tp, mean)
  month_avg4 <- aggregate(PREC ~ year + month, data = tp, mean) #월평균 강수량, 월합계 강수량을 구하려면 sum 
  
  month_avg <- merge(month_avg1, month_avg2, by=c("year", "month"), all = TRUE)
  month_avg <- merge(month_avg, month_avg3, by=c("year", "month"), all = TRUE)
  month_avg <- merge(month_avg, month_avg4, by=c("year", "month"), all = TRUE)
  
  colnames(month_avg)[6] <- "AVG_PRE"
  
  # 관측빈도와 합병 
  month_avg <- merge(month_avg, freq, by=c("year", "month"), all = TRUE)
  
  # 관측개수(V7~V10)가 0인 것(관측값개수 80%이하)은 NA가 되도록 함  
  for(k in 1:3){
    for(l in 1:nrow(month_avg)){
      if(month_avg[l,9+k] == 0){
        month_avg[l,2+k] <- NA
        print(paste("ID :", meta$id[i], l, "th data has replaced as NA"))
      }
    }()
  }
  
  write.csv(month_avg, 
            paste0("\\\\climatelab\\data\\CMA\\monthly\\", tp$ID[i] ,".csv"))
  
  month_avg$ID <- tp$ID[i]
  
  data <- rbind(data, month_avg) #nrow = 517 * 360
  
  cnt <- cnt + 1
  print(paste(cnt, "th station file saved."))
}

write.csv(data, "\\\\climatelab\\data\\CMA\\all_month_cma.csv")

# errors ------------------------------------------------------------------

errors <- station[which(station$AVG_TEMP >= 45),-1]

colnames(meta)[2] <- "ID"
rawlist <- list()

for (i in 1:nrow(errors)){
  meta_error <- meta[which(meta$ID == errors[i,'ID']),] 
  temp <- read.table(paste0(t_dir, '\\', meta_error$province, '\\', 
                                    meta_error$province, meta_error$ID, "temps"), 
                             header = FALSE, col.names = colt)
  rawlist[[i]] <- subset(temp, (temp$year == errors[i,'year'] & temp$month == errors[i,'month']))
}

raw <- data.frame
for (i in 1:nrow(errors)){
  raw <- rbind(raw, data.frame(rawlist[[i]]))
}

write.csv(raw, "\\\\climatelab\\data\\CMA\\errors.csv")