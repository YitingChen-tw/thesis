install.packages("truncnorm")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("purrr")
install.packages("tidyr")
install.packages("scales")
install.packages("rlang")
install.packages("vctrs")
install.packages("tidyselect")
install.packages("pillar")
install.packages("tibble")
install.packages("forcats")
install.packages("scales")
install.packages("writexl")

# ==============================================================================
# 環境初始化
# ==============================================================================
graphics.off()     # 關閉舊圖表
gc()               # 釋放記憶體

#========================================
# Dietary ingestion rate Monte Carlo
#========================================
library(truncnorm)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(forcats)
library(scales)
library(writexl)

rm(list = ls())

set.seed(123)
n_sim <- 10000

EF <- 365
ED <- 47
CF <- 1e-6
AT <- list(
  m = 47 * 365,
  f = 47 * 365,
  t = 47 * 365
)

#========================================
# BW (kg)
#========================================
BW_m <- rtruncnorm(n_sim, a=10.67, b=140.13, mean=73.45, sd=28.67)
BW_f <- rtruncnorm(n_sim, a=14.72, b=102.68, mean=58.86, sd=19.38)
BW_t <- rtruncnorm(n_sim, a=6.14, b=128.06, mean=66.02, sd=26.19)


#========================================
# Ingestion rate parameter table
#========================================
IR_para <- list(
  
  grain = list(
    m = c(mean=431.64, sd=254.69, min=7.49, max=1927.16),
    f = c(mean=280.23, sd=194.59, min=0.06, max=1862),
    t = c(mean=355.13, sd=238.63, min=0.06, max=1927.16)
  ),
  
  vegetable = list(
    m = c(mean=273.94, sd=227.46, min=0.01, max=2021.52),
    f = c(mean=279.59, sd=220.01, min=0.0001, max=1990.77),
    t = c(mean=276.82, sd=223.69, min=0.0001, max=2021.52)
  ),
  
  fruit = list(
    m = c(mean=266.92, sd=249.72, min=0.04, max=1877.35),
    f = c(mean=235.97, sd=219.40, min=0.38, max=2754.44),
    t = c(mean=249.55, sd=233.64, min=0.04, max=2754.44)
  ),
  
  dairy = list(
    m = c(mean=182.72, sd=209.29, min=0.4, max=1693.08),
    f = c(mean=162.9, sd=162.91, min=0.54, max=1291.23),
    t = c(mean=171.54, sd=184.74, min=0.4, max=1693.08)
  ),
  
  egg = list(
    m = c(mean=64.34, sd=48.12, min=0.05, max=309.55),
    f = c(mean=56.91, sd=41.68, min=0.01, max=486.16),
    t = c(mean=60.62, sd=45.15, min=0.01, max=486.16)
  ),
  
  chicken = list(
    m = c(mean=120.1, sd=106.94, min=4.41, max=773.43),
    f = c(mean=82.36, sd=69.72, min=2.76, max=468.12),
    t = c(mean=101.54, sd=92.48, min=2.76, max=773.43)
  ),
  
  fish = list(
    m = c(mean=76.94, sd=74.04, min=0.03, max=819.40),
    f = c(mean=63.33, sd=62.37, min=0.04, max=563.17),
    t = c(mean=70.09, sd=68.74, min=0.03, max=819.40)
  ),
  
  seafood = list(
    m = c(mean=48.86, sd=57.13, min=0.16, max=444.91),
    f = c(mean=39.21, sd=40.17, min=0.48, max=255.88),
    t = c(mean=43.98, sd=49.51, min=0.16, max=444.91)
  ),
  
  pork = list(
    m = c(mean=111.55, sd=94.66, min=0.19, max=880.91),
    f = c(mean=78.25, sd=67.71, min=0.19, max=556.37),
    t = c(mean=95.57, sd=84.46, min=0.19, max=880.91)
  ),
  
  beef = list(
    m = c(mean=100.53, sd=96.46, min=0.11, max=677.29),
    f = c(mean=72.91, sd=73.52, min=3.77, max=633.44),
    t = c(mean=88.27, sd=88.02, min=0.11, max=677.29)
  ),
  
  canned_food = list(
    m = c(mean=483.35, sd=852.26, min=0.04, max=8646.88),
    f = c(mean=282.02, sd=531.11, min=0.15, max=5578.63),
    t = c(mean=394.00, sd=718.45, min=0.04, max=8646.88)
  ),
  canned_drink = list(
    m = c(mean=755.11, sd=1088.30, min=0.12, max=8646.88),
    f = c(mean=488.06, sd=716.09, min=0.20, max=5578.63),
    t = c(mean=624.59, sd=925.13, min=0.12, max=8646.88)
  ),
  canned_fruit = list(
    m = c(mean=149.53, sd=289.71, min=0.04, max=2448.51),
    f = c(mean=83.92, sd=209.25, min=0.15, max=2448.51),
    t = c(mean=112.16, sd=246.92, min=0.04, max=2448.51)
  ),
  canned_dessert = list(
    m = c(mean=139.27, sd=190.64, min=2.36, max=465.72),
    f = c(mean=107.12, sd=174.65, min=2.36, max=417.67),
    t = c(mean=123.97, sd=181.73, min=2.36, max=465.72)
  ),
  canned_vegetable = list(
    m = c(mean=22.59, sd=30.91, min=0.05, max=220.70),
    f = c(mean=19.97, sd=30.29, min=0.82, max=154.93),
    t = c(mean=22.04, sd=31.20, min=0.05, max=220.70)
  ),
  canned_peanut = list(
    m = c(mean=40.37, sd=89.12, min=0.96, max=680.00),
    f = c(mean=23.06, sd=46.91, min=0.22, max=340.00),
    t = c(mean=31.79, sd=71.76, min=0.22, max=680.00)
  ),
  canned_seafood = list(
    m = c(mean=25.15, sd=47.86, min=0.21, max=240.00),
    f = c(mean=21.11, sd=31.76, min=0.21, max=137.07),
    t = c(mean=28.63, sd=46.51, min=0.21, max=240.00)
  ),
  canned_meat = list(
    m = c(mean=70.65, sd=50.93, min=6.32, max=284.93),
    f = c(mean=53.38, sd=35.58, min=5.37, max=199.34),
    t = c(mean=62.80, sd=45.25, min=5.37, max=284.93)
  )
  
)

make_ir <- function(x){
    rtruncnorm(n_sim, a = x["min"], b = x["max"],mean = x["mean"], sd = x["sd"])
  }


IR_sim <- list()

for(food in names(IR_para)){
    IR_sim[[food]] <- list(
    m = make_ir(IR_para[[food]]$m),
    f = make_ir(IR_para[[food]]$f),
    t = make_ir(IR_para[[food]]$t)
  )
  
}

#========================================
# Example output
#========================================

head(IR_sim$grain$m)

summary(IR_sim$canned_food$t)


#========================================
# Monte Carlo simulation
#========================================

Cfood_para <- list(
    grain = list(
      BPA = c(mean = 4.80, sd = 3.54, min = 0.83, max = 9.28),
      MeP = c(mean = 20.32, sd = 81.89, min = 0.12, max = 597.45),
      'BP-1' = c(mean = 0.69, sd = 0.83, min = 0.00, max = 3.70),
      'BP-3' = c(mean = 3.65, sd = 15.30, min = 0.00, max = 117.20),
      'BP-8' = c(mean = 0.17, sd = 0.37, min = 0.00, max = 1.40)
    ),
    vegetable = list(
      BPA = c(mean = 0.37, sd = 0.19, min = 0.13, max = 0.92),
      MeP = c(mean = 5.72, sd = 19.77, min = 0.04, max = 105.57)
    ),
    fruit = list(
      BPA = c(mean = 0.71, sd = 0.32, min = 0.25, max = 0.97),
      MeP = c(mean = 8.72, sd = 10.15, min = 0.03, max = 31.51)
    ),
    dairy = list(
      BPA = c(mean = 0.85, sd = 1.50, min = 0.09, max = 7.86),
      MeP = c(mean = 4.95, sd = 9.14, min = 0.12, max = 38.53)
    ),
    egg = list(
      BPA = c(mean = 1.53, sd = 1.02, min = 0.77, max = 5.34)
    ),
    chicken = list(
      BPA = c(mean = 4.88, sd = 5.60, min = 0.75, max = 14.13)
    ),
    fish = list(
      BPA = c(mean = 13.36, sd = 18.38, min = 0.21, max = 59.21),
      MeP = c(mean = 3.57, sd = 7.23, min = 0.04, max = 32.50),
      'BP-1' = c(mean = 0.05, sd = 0.05, min = 0.00, max = 0.15),
      'BP-3' = c(mean = 0.13, sd = 0.18, min = 0.00, max = 0.90),
      'BP-8' = c(mean = 0.01, sd = 0.01, min = 0.00, max = 0.01)
    ),
    seafood = list(
      BPA = c(mean = 4.56, sd = 3.34, min = 0.29, max = 9.19),
      'BP-1' = c(mean = 0.37, sd = 0.30, min = 0.01, max = 0.70),
      'BP-3' = c(mean = 0.06, sd = 1.22, min = 0.01, max = 0.56),
      'BP-8' = c(mean = 0.57, sd = 0.45, min = 0.01, max = 0.90)
    ),
    pork = list(
      BPA = c(mean = 1.38, sd = 1.10, min = 0.27, max = 3.06)
    ),
    beef = list(
      BPA = c(mean = 5.77, sd = 4.85, min = 0.69, max = 11.92)
    ),
    canned_food = list(
      BPA = c(mean = 28.75, sd = 73.67, min = 0.09, max = 166.34)
    ),
    canned_drink = list(
      BPA = c(mean = 7.86, sd = 6.99, min = 0.09, max = 14.75)
    ),
    canned_fruit = list(
      BPA = c(mean = 2.44, sd = 11.80, min = 0.26, max = 22.10)
    ),
    canned_dessert = list(
      BPA = c(mean = 8.35, sd = 3.71, min = 5.81, max = 14.60)
    ),
    canned_vegetable = list(
      BPA = c(mean = 10.13, sd = 11.99, min = 0.87, max = 33.90)
    ),
    canned_peanut = list(
      BPA = c(mean = 94.23, sd = 13.76, min = 73.77, max = 103.66)
    ),
    canned_seafood = list(
      BPA = c(mean = 13.32, sd = 13.95, min = 1.13, max = 39.90)
    ),
    canned_meat = list(
      BPA = c(mean = 86.73, sd = 71.19, min = 1.22, max = 166.34)
    )
)

make_cfood <- function(x){
  rtruncnorm(n_sim, a = x["min"], b = x["max"], mean = x["mean"], sd = x["sd"])
}


Cfood_sim <- list()

for(food in names(Cfood_para)){
  
  Cfood_sim[[food]] <- list()
  
  for(chem in names(Cfood_para[[food]])){
    
    Cfood_sim[[food]][[chem]] <-
      make_cfood(Cfood_para[[food]][[chem]])
    
  }
  
}

#========================================
# Monte Carlo function
#========================================
calc_ADD <- function(C, IR, BW, AT){
  (C * IR * EF * ED * CF) / (BW * AT)
}
ADD_sim <- list()
for(food in names(Cfood_para)){
  
  ADD_sim[[food]] <- list()
  
  for(chem in names(Cfood_para[[food]])){
    
    # 濃度
    C <- Cfood_sim[[food]][[chem]]
    
    # ADD calculation
    ADD_sim[[food]][[chem]] <- list(
      
      m = calc_ADD(
        C = C,
        IR = IR_sim[[food]]$m,
        BW = BW_m,
        AT = AT$m
      ),
      
      f = calc_ADD(
        C = C,
        IR = IR_sim[[food]]$f,
        BW = BW_f,
        AT = AT$f
      ),
      
      t = calc_ADD(
        C = C,
        IR = IR_sim[[food]]$t,
        BW = BW_t,
        AT = AT$t
      )
      
    )
    
  }
}

get_summary <- function(x){
  data.frame(
    Mean = mean(x),
    SD = sd(x),
    Median = median(x),
    P50 = median(x),
    P95 = quantile(x, 0.95),
    P99 = quantile(x, 0.99)
  )
}

get_food_summary <- function(food, chem){
  
  data.frame(
    Male = get_summary(ADD_sim[[food]][[chem]]$m),
    Female = get_summary(ADD_sim[[food]][[chem]]$f),
    Total = get_summary(ADD_sim[[food]][[chem]]$t)
  )
}



save(
  ADD_sim,
  Cfood_sim,
  IR_sim,
  BW_m,
  BW_f,
  BW_t,
  file = "Diet_ADD.RData"
)

#========================================
# 整理成long table(畫圖用)
#========================================

all_summary <- list()

for(food in names(ADD_sim)){
  for(chem in names(ADD_sim[[food]])){
    
    tmp <- data.frame(
      Food = food,
      Chemical = chem,
      Group = c("Male","Female","Total"),
      Mean = c(
        mean(ADD_sim[[food]][[chem]]$m),
        mean(ADD_sim[[food]][[chem]]$f),
        mean(ADD_sim[[food]][[chem]]$t)
      ),
      P95 = c(
        quantile(ADD_sim[[food]][[chem]]$m,0.95),
        quantile(ADD_sim[[food]][[chem]]$f,0.95),
        quantile(ADD_sim[[food]][[chem]]$t,0.95)
      )
    )
    
    all_summary[[paste(food,chem)]] <- tmp
  }
}

all_summary_df <- do.call(rbind, all_summary)

write_xlsx(all_summary_df, "Diet_ADD_Summary.xlsx")


#========================================
# Food category grouping
#========================================

general_foods <- c(
  "grain",
  "vegetable",
  "fruit",
  "dairy",
  "egg",
  "chicken",
  "fish",
  "seafood",
  "pork",
  "beef",
  "canned_food"
)

canned_detail_foods <- c(
  "canned_drink",
  "canned_fruit",
  "canned_dessert",
  "canned_vegetable",
  "canned_peanut",
  "canned_seafood",
  "canned_meat"
)

#========================================
# 整理成 long format
#========================================

plot_df <- all_summary_df %>%
  tidyr::pivot_longer(
    cols = c(Mean, P95),
    names_to = "Metric",
    values_to = "Value"
  )

general_food_df <- subset(
  plot_df,
  Food %in% general_foods
)

canned_detail_df <- subset(
  plot_df,
  Food %in% canned_detail_foods
)



#========================================
# Plot 1 : General foods (BOXPLOT)
#========================================

# 建立 chemical 分組
chem_group_1 <- c("BP-1", "BP-3", "BP-8")
chem_group_2 <- c("BPA", "MeP")


general_dist_plot <- data.frame()

for(food in general_foods){
  
  for(chem in names(ADD_sim[[food]])){
    
    tmp <- rbind(
      data.frame(
        Food = food,
        Chemical = chem,
        Group = "Male",
        ADD = ADD_sim[[food]][[chem]]$m
      ),
      data.frame(
        Food = food,
        Chemical = chem,
        Group = "Female",
        ADD = ADD_sim[[food]][[chem]]$f
      ),
      data.frame(
        Food = food,
        Chemical = chem,
        Group = "Total",
        ADD = ADD_sim[[food]][[chem]]$t
      )
    )
    
    general_dist_plot <- rbind(
      general_dist_plot,
      tmp
    )
  }
}


dist_plot_g1 <- subset(
  general_dist_plot,
  Chemical %in% chem_group_1
)

dist_plot_g2 <- subset(
  general_dist_plot,
  Chemical %in% chem_group_2
)

# 畫圖函數
  plot_box <- function(df, title_text){
    
    ggplot(
      df,
      aes(x = Food, y = ADD, fill = Group)
      ) +

    geom_boxplot(
      position = position_dodge(width = 0.8),
      width = 0.5,
      alpha = 0.7,
      outlier.size = 1
      ) +
      
    # 平均值紅點（關鍵）
    stat_summary(
      fun = mean,
      geom = "point",
      position = position_dodge(width = 0.8),
      color = "red",
      size = 2
      ) +
      
    facet_wrap(~Chemical) +
      
    scale_y_log10(labels = scales::scientific) +

    scale_fill_manual(
      values = c(
        Male = "#A8E6A3",    # 淡綠
        Female = "#F7B6D2",  # 淡粉
        Total = "#A7C7E7"    # 淡藍
        )
      )+
                  
    theme_bw(base_size = 19) +
      theme(
        
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = 16        # ✔️ X軸字變大
        ),
        
        axis.text.y = element_text(
          size = 16        # ✔️ Y軸字變大
        ),
        
        axis.title = element_text(
          size = 18,
          face = "bold"
        ),
        
        strip.text = element_text(
          size = 18,
          face = "bold"
        ),
        
        legend.position = "bottom",
        legend.box = "horizontal",        
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        
        panel.grid.minor = element_blank(),
        
        plot.title = element_text(
          size = 20,
          face = "bold",
          hjust = 0.5   # ✔️ 標題置中
        )
      ) +
      
      labs(
        title = title_text,
        x = "Food Category",
        y = "ADD (log scale)",
        fill = "Population Group"
      )
  }

# 畫兩張圖
plot_box(dist_plot_g1,"")
plot_box(dist_plot_g2,"")
  
  
  
#========================================
# Plot 2 : Canned food details
#========================================

#-----------------------------
# 建立 canned boxplot 原始資料
#-----------------------------
canned_dist_plot <- data.frame()

for(food in canned_detail_foods){
  
  for(chem in names(ADD_sim[[food]])){
    
    tmp <- rbind(
      data.frame(
        Food = food,
        Chemical = chem,
        Group = "Male",
        ADD = ADD_sim[[food]][[chem]]$m
      ),
      data.frame(
        Food = food,
        Chemical = chem,
        Group = "Female",
        ADD = ADD_sim[[food]][[chem]]$f
      ),
      data.frame(
        Food = food,
        Chemical = chem,
        Group = "Total",
        ADD = ADD_sim[[food]][[chem]]$t
      )
    )
    
    canned_dist_plot <- rbind(canned_dist_plot, tmp)
  }
}

#-----------------------------
# Boxplot function（同你 general style）
#-----------------------------
plot_canned_box <- function(df, title_text){
  
  ggplot(df, aes(x = Food, y = ADD, fill = Group)) +
    
    geom_boxplot(
      position = position_dodge(width = 0.8),
      width = 0.5,
      alpha = 0.7,
      outlier.size = 0.8
    ) +
    
    stat_summary(
      fun = mean,
      geom = "point",
      position = position_dodge(width = 0.8),
      color = "red",
      size = 2
    ) +
    
    facet_wrap(~Chemical) +
    
    scale_y_log10(labels = scales::scientific) +
    
    scale_fill_manual(
      values = c(
        Male   = "#A8E6A3",   # 淡綠
        Female = "#F7B6D2",   # 淡粉
        Total  = "#A7C7E7"    # 淡藍
      )
    ) +
    
    theme_bw(base_size = 18) +
    
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 18, face = "bold"),
      strip.text = element_text(size = 16, face = "bold"),
      
      legend.position = "bottom",
      legend.box = "horizontal",
      
      panel.grid.minor = element_blank(),
      
      plot.title = element_text(
        hjust = 0.5,
        size = 18,
        face = "bold"
      )
    ) +
    
    labs(
      title = title_text,
      x = "Canned Food Subcategories",
      y = "ADD (log scale)",
      fill = "Population Group"
    )
}

#-----------------------------
# 畫圖
#-----------------------------
plot_canned_box(canned_dist_plot, "Canned Food Exposure Distribution")