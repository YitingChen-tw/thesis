# ==============================================================================
# 環境初始化
# ==============================================================================
graphics.off()     # 關閉舊圖表
gc()               # 釋放記憶體

#========================================
# Drinking water exposure Monte Carlo
#========================================
library(truncnorm)
library(ggplot2)

rm(list = ls())
set.seed(123)

n_sim <- 10000

# 平均餘命(year)
LT_m <- 76.94 * 365
LT_f <- 83.74 * 365
LT_a <- 80.77 * 365

# ng → mg
CF <- 1e-6

# 暴露頻率
EF <- 365

# 暴露年數
ed_vec <- c(11, 10, 10, 11)

#========================================
# IR (mL/day)
#========================================
ir_a_means <- c(1343.8, 1438.8, 1508.4, 1442.2)
ir_a_sds   <- c(781.1, 770.0, 1010.6, 906.5)

ir_m_means <- c(1429.6, 1524.9, 1727.6, 1576.1)
ir_m_sds   <- c(812.3, 807.5, 1060.4, 967.7)

ir_f_means <- c(1235.5, 1352.1, 1307.2, 1304.4)
ir_f_sds   <- c(725.6, 710.4, 875.2, 801.3)

#========================================
# BW (kg)
#========================================
bw_a_means <- c(67.10, 67.10, 65.10, 65.10)
bw_a_sds   <- c(20.32, 20.32, 16.53, 16.53)

bw_m_means <- c(75.40, 75.40, 71.80, 71.80)
bw_m_sds   <- c(21.58, 21.58, 18.88, 18.88)

bw_f_means <- c(58.70, 58.70, 59.00, 59.00)
bw_f_sds   <- c(14.66, 14.66, 12.68, 12.68)

#========================================
# Chemical concentration (ng/L)
#========================================
chem <- c("BPA","BP-1","BP-3","BP-8","MeP","DEET")

C_mean <- c(7.68, 0.77, 10.78, 1.41, 5.00, 99.83)
C_sd   <- c(3.39, 0.24, 12.40, 0.00, 7.27, 44.22)

#========================================
# Monte Carlo function
#========================================
calc_ladd <- function(C_mean, C_sd,
                      ir_means, ir_sds,
                      bw_means, bw_sds,
                      LT){
  
  # concentration
  C_sim <- rtruncnorm(
    n_sim,
    a = 0,
    b = Inf,
    mean = C_mean,
    sd = ifelse(C_sd == 0, 1e-6, C_sd)
  )
  
  total_daily_dose <- rep(0, n_sim)
  
  # 儲存敏感度分析參數
  IR_store <- matrix(0, n_sim, length(ed_vec))
  BW_store <- matrix(0, n_sim, length(ed_vec))
  
  for(i in 1:length(ed_vec)){
    
    ir_tmp <- rtruncnorm(
      n_sim,
      a = 200,
      b = 5000,
      mean = ir_means[i],
      sd = ir_sds[i]
    )
    
    bw_tmp <- rtruncnorm(
      n_sim,
      a = 50,
      b = 100,
      mean = bw_means[i],
      sd = bw_sds[i]
    )
    
    # 儲存抽樣值
    IR_store[,i] <- ir_tmp
    BW_store[,i] <- bw_tmp
    
    # mL/day → L/day
    ir_tmp <- ir_tmp / 1000
    
    dose_i <- C_sim * ir_tmp *
      EF * ed_vec[i] * CF / bw_tmp
    
    total_daily_dose <- total_daily_dose + dose_i
  }
  
  ladd <- total_daily_dose / LT
  
  # 平均暴露參數
  IR_mean <- rowMeans(IR_store)
  BW_mean <- rowMeans(BW_store)
  
  return(list(
    ladd = ladd,
    C = C_sim,
    IR = IR_mean,
    BW = BW_mean
  ))
}

#========================================
# Run simulation
#========================================
water_result_t  <- list()
water_result_m  <- list()
water_result_f  <- list()

for(j in 1:length(chem)){
  
  water_result_t[[j]] <- calc_ladd(
    C_mean[j], C_sd[j],
    ir_a_means, ir_a_sds,
    bw_a_means, bw_a_sds,
    LT_a
  )
  
  water_result_m[[j]] <- calc_ladd(
    C_mean[j], C_sd[j],
    ir_m_means, ir_m_sds,
    bw_m_means, bw_m_sds,
    LT_m
  )
  
  water_result_f[[j]] <- calc_ladd(
    C_mean[j], C_sd[j],
    ir_f_means, ir_f_sds,
    bw_f_means, bw_f_sds,
    LT_f
  )
}

names(water_result_t)   <- chem
names(water_result_m)   <- chem
names(water_result_f)   <- chem
graphics.off()


save(water_result_t, water_result_m, water_result_f, file = "Water_LADD.RData")
#========================================
# Plot function (ggplot2 version)
#========================================
plot_box <- function(result_list, group_name){
  
  # 將 list 轉成 data.frame
  df <- do.call(rbind, lapply(seq_along(result_list), function(i){
    data.frame(
      Chemical = names(result_list)[i],
      LADD = result_list[[i]]$ladd
    )
  }))
  
  # 正式標題
  plot_title <- paste("Drinking Water LADD Distribution –", group_name)
  
  # 箱型圖
  p <- ggplot(df, aes(x = Chemical, y = LADD)) +
    geom_boxplot(fill = "lightblue", color = "black") +  # 黑邊框
    
    # 平均值
    stat_summary(fun = mean, geom = "point", 
                 shape = 19, size = 3, color = "red") +
    
    # 真正最大值
    stat_summary(fun = max, geom = "errorbar",
                 width = 0.4, linewidth = 1.8, color = "darkgreen") +
  
    coord_flip() +
    
    
    labs(title = plot_title, x = "", y = "LADD (mg/kg/day)") +
    theme_minimal(base_size = 13) +
    theme(
      axis.title.y = element_text(margin = margin(r = 20)),
      axis.text.y  = element_text(size = 14),
      axis.text.x  = element_text(size = 14),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5) # 圖表外框
    ) +
    
    # 在右上角加紅點圖示 + 標註
    annotate("point", 
             x = length(unique(df$Chemical)) + 0.2, 
             y = max(df$LADD) * 0.92, 
             color = "red", size = 3) +
    annotate("text", 
             x = length(unique(df$Chemical)) + 0.2, 
             y = max(df$LADD) * 0.94, 
             label = "Mean", hjust = 0, size = 3.5)
  
  print(p)
}

# 測試
plot_box(water_result_t, "Total Population")
plot_box(water_result_m, "Male")
plot_box(water_result_f, "Female")

#========================================
# Batch Sensitivity Analysis
#========================================
make_tornado <- function(result_obj,
                         chem_name,
                         group_name){
  
  # 建立 sensitivity dataframe
  sens <- data.frame(
    C    = result_obj[[chem_name]]$C,
    IR   = result_obj[[chem_name]]$IR,
    BW   = result_obj[[chem_name]]$BW,
    LADD = result_obj[[chem_name]]$ladd
  )
  
  # Spearman correlation(每個輸入變數 vs 最終 LADD 的相關性)
  tornado_df <- data.frame(
    Chemical    = chem_name,
    Group       = group_name,
    Variable    = c("Concentration", "Ingestion Rate", "Body Weight"),
    Correlation = c(
      cor(sens$C, sens$LADD, method="spearman"),
      cor(sens$IR, sens$LADD, method="spearman"),
      cor(sens$BW, sens$LADD, method="spearman")
    )
  )
  
  print(tornado_df)   # Console 顯示
  
}

all_corr <- list()

for(chem_name in chem){
  all_corr[[chem_name]] <- rbind(
    make_tornado(water_result_t, chem_name, "Total"),
    make_tornado(water_result_m, chem_name, "Male"),
    make_tornado(water_result_f, chem_name, "Female")
  )
}

# 合併成一個大表
all_corr_df <- do.call(rbind, all_corr)
head(all_corr_df)

write.csv(all_corr_df, "Water_LADD_Sensitivity.csv", row.names = FALSE)

#========================================
# Summary table
#========================================
make_summary <- function(x){
  
  data.frame(
    Chemical = names(x),
    Mean   = sapply(x, function(z) mean(z$ladd)),
    Median = sapply(x, function(z) median(z$ladd)),
    P95    = sapply(x, function(z)
      quantile(z$ladd, 0.95))
  )
}

summary_t <- make_summary(water_result_t)
summary_m <- make_summary(water_result_m)
summary_f <- make_summary(water_result_f)

print(summary_t)
print(summary_m)
print(summary_f)

write.csv(summary_t,"Summary_t.csv",row.names=FALSE)
write.csv(summary_m,"Summary_Male.csv",row.names=FALSE)
write.csv(summary_f,"Summary_Female.csv",row.names=FALSE)