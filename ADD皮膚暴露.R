# ==============================================================================
# 環境初始化
# ==============================================================================
graphics.off()     # 關閉舊圖表
gc()               # 釋放記憶體


#========================================
# Dermal Exposure Monte Carlo with Sensitivity Analysis
#========================================
library(ggplot2)
library(truncnorm)

rm(list = ls())
set.seed(123)
n_sim <- 10000

#========================
# Lifetime (days) -> Averaging time
#========================
AT_m <- 47 * 365
AT_f <- 47 * 365
AT_t <- 47 * 365

#========================
# Exposure parameters
#========================
EF <- 365
ED <- 46
CF <- 0.001

#========================
# Body weight (kg)
#========================
BW_m <- rtruncnorm(n_sim, a=10.67, b=140.13, mean=73.45, sd=28.67)
BW_f <- rtruncnorm(n_sim, a=14.72, b=102.68, mean=58.86, sd=19.38)
BW_t <- rtruncnorm(n_sim, a=10.67, b=140.13, mean=66.02, sd=26.19)

#========================
# PCP parameters
#========================
AQ <- list(shampoo = 8, perfume = 0.75, bodywash = 5, bodylotion = 8) # g/time
EV <- list(shampoo = 1, perfume = 1, bodywash = 1.07, bodylotion = 0.71) # times/day
F <- list(shampoo = 0.01, perfume = 1, bodywash = 0.01, bodylotion = 1) # retention factor
A_BPA <- 0.036
A_MeP <- 0.609

#========================
# Chemical concentrations in PCP (µg/g)
#========================
rnorm_pos <- function(n, mean, sd) pmax(rnorm(n, mean, sd), 0)
C_BPA_shampoo <- rnorm_pos(n_sim, 3.79, 1.5)
C_MeP_shampoo <- rnorm_pos(n_sim, 2.8, 0.7)
C_MeP_perfume <- rnorm_pos(n_sim, 796.93, 199.23)
C_MeP_bodywash <- rnorm_pos(n_sim, 676.17, 234.73)
C_MeP_bodylotion <- rnorm_pos(n_sim, 747.32, 454.08)

#========================
# ADD function (mg/kg/day)
#========================
calc_ADD <- function(Cpcp, AQ, EV, EF, ED, A, F, CF, BW, AT){
  (Cpcp * AQ * EV * EF * ED * A * F * CF) / (BW * AT)
}

#========================
# BPA
#========================
ADD_BPA_m <- calc_ADD(C_BPA_shampoo, AQ$shampoo, EV$shampoo, EF, ED, A_BPA, F$shampoo, CF, BW_m, AT_m)
ADD_BPA_f <- calc_ADD(C_BPA_shampoo, AQ$shampoo, EV$shampoo, EF, ED, A_BPA, F$shampoo, CF, BW_f, AT_f)
ADD_BPA_t <- calc_ADD(C_BPA_shampoo, AQ$shampoo, EV$shampoo, EF, ED, A_BPA, F$shampoo, CF, BW_t, AT_t)

#========================
# MeP
#========================
MeP <- function(BW, AT){
  list(
    shampoo = calc_ADD(C_MeP_shampoo, AQ$shampoo, EV$shampoo, EF, ED, A_MeP, F$shampoo, CF, BW, AT),
    perfume = calc_ADD(C_MeP_perfume, AQ$perfume, EV$perfume, EF, ED, A_MeP, F$perfume, CF, BW, AT),
    bodywash = calc_ADD(C_MeP_bodywash, AQ$bodywash, EV$bodywash, EF, ED, A_MeP, F$bodywash, CF, BW, AT),
    bodylotion = calc_ADD(C_MeP_bodylotion, AQ$bodylotion, EV$bodylotion, EF, ED, A_MeP, F$bodylotion, CF, BW, AT)
  )
}

MeP_m <- MeP(BW_m, AT_m)
MeP_f <- MeP(BW_f, AT_f)
MeP_t <- MeP(BW_t, AT_t)

#========================
# Combine into lists for plotting
#========================
dermal_result_t <- list(BPA_shampoo=ADD_BPA_t, MeP_shampoo=MeP_t$shampoo, MeP_perfume=MeP_t$perfume,
                 MeP_bodywash=MeP_t$bodywash, MeP_bodylotion=MeP_t$bodylotion)
dermal_result_m <- list(BPA_shampoo=ADD_BPA_m, MeP_shampoo=MeP_m$shampoo, MeP_perfume=MeP_m$perfume,
                 MeP_bodywash=MeP_m$bodywash, MeP_bodylotion=MeP_m$bodylotion)
dermal_result_f <- list(BPA_shampoo=ADD_BPA_f, MeP_shampoo=MeP_f$shampoo, MeP_perfume=MeP_f$perfume,
                 MeP_bodywash=MeP_f$bodywash, MeP_bodylotion=MeP_f$bodylotion)




#========================
# Sensitivity Analysis
#========================
make_tornado <- function(result_obj, chem_name, group_name, BW){
  
  if(chem_name == "BPA_shampoo"){
    C_val <- C_BPA_shampoo
    AQ_val <- AQ$shampoo
    EV_val <- EV$shampoo
    F_val <- F$shampoo
    A_val <- A_BPA
  } else {
    pcp_name <- gsub("MeP_","",chem_name)
    C_val <- switch(pcp_name,
                    shampoo = C_MeP_shampoo,
                    perfume = C_MeP_perfume,
                    bodywash = C_MeP_bodywash,
                    bodylotion = C_MeP_bodylotion)
    AQ_val <- AQ[[pcp_name]]
    EV_val <- EV[[pcp_name]]
    F_val <- F[[pcp_name]]
    A_val <- A_MeP
  }
  
  sens <- data.frame(
    Cpcp = C_val,
    AQ   = rep(AQ_val, length.out = n_sim),
    EV   = rep(EV_val, length.out = n_sim),
    EF   = rep(EF, length.out = n_sim),
    ED   = rep(ED, length.out = n_sim),
    A    = rep(A_val, length.out = n_sim),
    F    = rep(F_val, length.out = n_sim),
    BW   = rep(BW, length.out = n_sim),
    ADD = result_obj[[chem_name]]
  )
  
  vars <- c("Cpcp","AQ","EV","EF","ED","A","F","BW")
  Correlation <- sapply(vars, function(v){
    if(length(unique(sens[[v]])) <= 1){
      return(NA)
    } else {
      return(cor(sens[[v]], sens$ADD, method="spearman"))
    }
  })
  
  tornado_df <- data.frame(
    Chemical = chem_name,
    Group = group_name,
    Variable = vars,
    Correlation = Correlation
  )
  
  return(tornado_df)
}

chem_names <- names(dermal_result_t)
all_corr <- list()
for(chem_name in chem_names){
  all_corr[[chem_name]] <- rbind(
    make_tornado(dermal_result_t, chem_name, "Total", BW_t),
    make_tornado(dermal_result_m, chem_name, "Male", BW_m),
    make_tornado(dermal_result_f, chem_name, "Female", BW_f)
  )
}
all_corr_df <- do.call(rbind, all_corr)
write.csv(all_corr_df, "Dermal_ADD_Sensitivity.csv", row.names=FALSE)


save(
  dermal_result_t,
  dermal_result_m,
  dermal_result_f,
  C_BPA_shampoo,
  BW_m,
  BW_f,
  file = "Dermal_ADD.RData"
)


#========================
# Dermal data visualized with Drinking Water style
#========================
library(ggplot2)

plot_dermal_as_water <- function(result_list, group_name){
  
  # 將 list 轉成 data.frame
  df <- do.call(rbind, lapply(seq_along(result_list), function(i){
    data.frame(
      Chemical = names(result_list)[i],
      ADD = result_list[[i]]
    )
  }))
  
  plot_title <- paste("Lifetime Average Daily Dose (ADD) via Dermal Exposure –", group_name)
  
  # 計算最大值，自動調整軸範圍
  max_ADD <- max(df$ADD, na.rm = TRUE)
  y_upper <- max_ADD * 1.1  # 留 10% 緩衝
  
  # 計算平均值與 P95，用於右上角圖示標註
  overall_mean <- mean(df$ADD, na.rm = TRUE)
  overall_p95  <- quantile(df$ADD, 0.95, na.rm = TRUE)
  
  # 使用箱型圖 + 平均值 + 離群值
  p <- ggplot(df, aes(x = Chemical, y = ADD, fill = Chemical)) +
    
    geom_boxplot(width = 0.6, fill = "#C2B280", color = "black", outlier.shape = 19,
                 outlier.size = 1.8, alpha = 0.9) +

    
    stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
    
    coord_flip() +
    
    scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
    
    labs(title = plot_title, x = "", y = "ADD (mg/kg/day)") +

    
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y  = element_text(size = 14),
      axis.text.x  = element_text(size = 14),
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )+
    
    annotate("point", 
             x = length(unique(df$Chemical)) + 0.2, 
             y = max(df$ADD) * 0.92, 
             color = "red", size = 3) +
    annotate("text", 
             x = length(unique(df$Chemical)) + 0.2, 
             y = max(df$ADD) * 0.94, 
             label = "Mean", hjust = 0, size = 5)    
  
  print(p)
}

#========================
# 測試範例
plot_dermal_as_water(dermal_result_t, "Total Population")
plot_dermal_as_water(dermal_result_m, "Male")
plot_dermal_as_water(dermal_result_f, "Female")
#========================
# Summary table
#========================
make_summary <- function(x){
  data.frame(
    Chemical = names(x),
    Mean = sapply(x, mean),
    Median = sapply(x, median),
    P95 = sapply(x, function(z) quantile(z,0.95))
  )
}

summary_t <- make_summary(dermal_result_t)
summary_m <- make_summary(dermal_result_m)
summary_f <- make_summary(dermal_result_f)

write.csv(summary_t,"Dermal_Summary_Total.csv",row.names=FALSE)
write.csv(summary_m,"Dermal_Summary_Male.csv",row.names=FALSE)
write.csv(summary_f,"Dermal_Summary_Female.csv",row.names=FALSE)