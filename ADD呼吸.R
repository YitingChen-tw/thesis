# ==============================================================================
# 環境初始化
# ==============================================================================
graphics.off()     # 關閉舊圖表
gc()               # 釋放記憶體


#========================
# 1. Monte Carlo setup
#========================
library(truncnorm)

n_sim <- 10000
set.seed(123)
AT_m <- 47 * 365  # 男性平均餘命
AT_f <- 47 * 365  # 女性平均餘命  # 終身致癌評估分母 (年)
CF <- 1e-6     # 單位換算: ng 轉 mg (因為圖中濃度是 ng/m3)

#========================
# 2. 定義各年齡層參數向量
#========================
# 暴露年數 ED
ed_vec <- c(6, 5, 5, 5, 5, 5, 5, 5, 5) 

# --- 男性參數-----
ve_m_means <- c(15.37, 16.35, 16.35, 14.11, 14.11, 14.11, 14.11, 14.06, 14.06)
ve_m_sds   <- c(1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01) 

bw_m_means <- c(75.40, 75.40, 75.40, 75.40, 75.40, 71.80, 71.80, 71.80, 71.80)
bw_m_sds   <- c(21.58, 21.58, 21.58, 21.58, 21.58, 18.88, 18.88, 18.88, 18.88)

# --- 女性參數 (依此類推) ---
ve_f_means <- c(10.13, 10.26, 10.26, 11.50, 11.50, 11.50, 11.50, 10.07, 10.07)
ve_f_sds   <- c(0.71, 0.71, 0.71, 0.71, 0.71, 0.71, 0.71, 0.71, 0.71) 

bw_f_means <- c(58.70, 58.70, 58.70, 58.70, 58.70, 59.00, 59.00, 59.00, 59.00)
bw_f_sds   <- c(14.66, 14.66, 14.66, 14.66, 14.66, 12.68, 12.68, 12.68, 12.68)

# 共同參數
C_sim <- rtruncnorm(n_sim, a = 0.000007, b = 100, mean = 33.32, sd = 112.08)
EF <- 365

#========================
# 3. 執行加權計算 (向量化迴圈)
#========================
library(truncnorm)

# 初始化存儲向量
total_daily_dose_m <- rep(0, n_sim)
total_daily_dose_f <- rep(0, n_sim)

# 儲存敏感度分析用參數
VE_m_store <- matrix(0, n_sim, length(ed_vec))
BW_m_store <- matrix(0, n_sim, length(ed_vec))

VE_f_store <- matrix(0, n_sim, length(ed_vec))
BW_f_store <- matrix(0, n_sim, length(ed_vec))


for(i in 1:length(ed_vec)) {
  
  # 男性該階段抽樣
  ve_m_tmp <- rtruncnorm(n_sim, a = 0.12, b = 23.71, mean = ve_m_means[i], sd = ve_m_sds[i])
  bw_m_tmp <- rtruncnorm(n_sim, a = 10.67, b = 128.44, mean = bw_m_means[i], sd = bw_m_sds[i])
  
  # 女性該階段抽樣
  ve_f_tmp <- rtruncnorm(n_sim, a = 0.11, b = 16.80, mean = ve_f_means[i], sd = ve_f_sds[i])
  bw_f_tmp <- rtruncnorm(n_sim, a = 14.72, b = 102.68, mean = bw_f_means[i], sd = bw_f_sds[i])
  
  # 儲存敏感度分析抽樣值
  VE_m_store[,i] <- ve_m_tmp
  BW_m_store[,i] <- bw_m_tmp
  
  VE_f_store[,i] <- ve_f_tmp
  BW_f_store[,i] <- bw_f_tmp
  
  
  # 計算該階段對一生的貢獻分量： (C * VE * EF * ED * CF) / BW
  total_daily_dose_m <- total_daily_dose_m + (C_sim * ve_m_tmp * EF * ed_vec[i] * CF / bw_m_tmp)
  total_daily_dose_f <- total_daily_dose_f + (C_sim * ve_f_tmp * EF * ed_vec[i] * CF / bw_f_tmp)
}

# 最後分別除以各自的 AT
# 注意：分子已經包含了各階段的 (ED * 365)，所以分母直接除以 AT (年) 即可
inh_BPA_ADD_m <- pmax(total_daily_dose_m / AT_m, 0)
inh_BPA_ADD_f <- pmax(total_daily_dose_f / AT_f, 0)


#========================
# 4. 敏感度分析（Spearman）
#========================
#男性
# 取平均暴露參數
VE_m_mean <- rowMeans(VE_m_store)
BW_m_mean <- rowMeans(BW_m_store)

# 建立資料框
sens_m <- data.frame(C  = C_sim, VE = VE_m_mean, BW = BW_m_mean, ADD = inh_BPA_ADD_m)

# Spearman correlation
cor_m <- cor(sens_m, method = "spearman")

round(cor_m, 3)


#女性
VE_f_mean <- rowMeans(VE_f_store)
BW_f_mean <- rowMeans(BW_f_store)

sens_f <- data.frame(C = C_sim, VE = VE_f_mean, BW = BW_f_mean, ADD = inh_BPA_ADD_f)

cor_f <- cor(sens_f, method = "spearman")

round(cor_f, 3)


#========================================
# Unified inhalation registry
#========================================

inh_VE_m <- VE_m_mean

inh_VE_f <- VE_f_mean

inh_VE_t <- (VE_m_mean + VE_f_mean) / 2

inh_BPA_ADD_t <- (inh_BPA_ADD_m + inh_BPA_ADD_f) / 2

inh_C_BPA <- C_sim


BW_m <- BW_m_mean

BW_f <- BW_f_mean

BW_t <- (BW_m_mean + BW_f_mean) / 2



save(
  inh_BPA_ADD_m, inh_BPA_ADD_f, inh_BPA_ADD_t,
  inh_C_BPA,
  inh_VE_m, inh_VE_f, inh_VE_t,
  BW_m, BW_f, BW_t,
  file = "Inhalation_BPA.RData"
)


#========================
# 6. 繪圖
#========================

m <- inh_BPA_ADD_m
f <- inh_BPA_ADD_f

options(scipen = -10)

#========================
# 男女共用X軸與breaks
#========================
xmax_all <- max(c(m, f))

common_breaks <- seq(
  0,
  xmax_all,
  length.out = 71   # 70 bins
)

#========================
# 繪圖函數
#========================
plot_dist <- function(data, title,
                      breaks_use,
                      xmax_use){
  
  # histogram（只計算，不畫）
  tmp <- graphics::hist(
    data,
    breaks = breaks_use,
    plot = FALSE
  )
  
  dens <- density(data, adjust = 1.2)
  
  ymax <- max(tmp$density, dens$y) * 1.2
  
  mean_v   <- mean(data)
  median_v <- median(data)
  p95_v    <- quantile(data, 0.95)
  
  par(
    bg = "white",
    mar = c(5, 8, 3, 2),
    mgp = c(3, 1, 0), 
    xaxs = "i",
    yaxs = "i",
    cex.main = 1.4,
    cex.lab  = 1.3,
    cex.axis = 1.2,
    las = 1
  )
  
  graphics::hist(
    data,
    breaks = breaks_use,
    probability = TRUE,
    col = "royalblue1",
    border = "white",
    xlim = c(0, xmax_use),
    ylim = c(0, ymax),
    main = title,
    xlab = "ADD (mg/kg/day)",
    ylab = "",
    las = 1,
    xaxt = "n",   # 關掉預設 x 軸
    yaxt = "n"    # 關掉預設 y 軸
  )
  
  # 自訂 X 軸
  x_ticks <- pretty(c(0, xmax_use))
  
  axis(
    side = 1,
    at = x_ticks,
    labels = sprintf("%.2e", x_ticks)
  )
  
  # 自訂 Y 軸
  y_ticks <- pretty(c(0, ymax))
  
  axis(
    side = 2,
    at = y_ticks,
    labels = sprintf("%.2e", y_ticks),
    las = 1
  )
  
  title(ylab = "Probability density", line = 6)
  
  lines(dens, lwd = 1.4, col = "blue4")
  
  abline(v = mean_v, col = "red", lty = 2)
  abline(v = median_v, col = "red", lty = 2)
  abline(v = p95_v, col = "red", lty = 2)
  
  legend("topright",
         legend = c(
           paste0("Mean = ", sprintf("%.2e", mean_v)),
           paste0("Median = ", sprintf("%.2e", median_v)),
           paste0("P95 = ", sprintf("%.2e", p95_v))
         ),
         bty = "n",
         cex = 1.2)
}



# 測試
plot_dist(m, "Male BPA Inhalation ADD Distribution",
          breaks_use = common_breaks, xmax_use = xmax_all)
plot_dist(f, "Female BPA Inhalation ADD Distribution",
          breaks_use = common_breaks, xmax_use = xmax_all)