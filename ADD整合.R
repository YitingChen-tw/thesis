install.packages("viridis")
install.packages("forcats")
install.packages("crayon")

library(scales)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(viridis)
library(forcats)

#========================================================
# PUBLICATION-GRADE CUMULATIVE EXPOSURE ENGINE
#========================================================

#========================================================
# Packages
#========================================================

required_packages <- c(
  "ppcor",
  "dplyr",
  "purrr",
  "tidyr"
)

invisible(
  lapply(
    required_packages,
    function(pkg){
      if(!require(pkg, character.only = TRUE)){
        install.packages(pkg)
        library(pkg, character.only = TRUE)
      }
    }
  )
)

#========================================================
# Global setup
#========================================================

graphics.off()
gc()

set.seed(123)

n_sim <- 10000

#========================================================
# Load pathway outputs
#========================================================

load("Inhalation_BPA.RData")
load("Water_ADD.RData")
load("Diet_ADD.RData")
load("Dermal_ADD.RData")

#========================================================
# Registry object
#========================================================

EXPOSURE_REGISTRY <- list()

#========================================================
# Universal exposure registrar
#========================================================

register_exposure <- function(
    registry,
    group,
    pathway,
    chemical,
    source,
    ADD,
    inputs = list()
){
  
  if(is.null(registry[[group]])){
    registry[[group]] <- list()
  }
  
  if(is.null(registry[[group]][[chemical]])){
    registry[[group]][[chemical]] <- list()
  }
  
  if(is.null(
    registry[[group]][[chemical]][[pathway]]
  )){
    
    registry[[group]][[chemical]][[pathway]] <- list()
  }
  
  registry[[group]][[chemical]][[pathway]][[source]] <- list(
    
    ADD = ADD,
    
    inputs = inputs
  )
  
  return(registry)
}

#========================================================
# Register inhalation
#========================================================

inh_data <- list(
  
  male = list(
    ADD = inh_BPA_ADD_m,
    VE = inh_VE_m,
    BW = BW_m,
    AT = 47 * 365
  ),
  
  female = list(
    ADD = inh_BPA_ADD_f,
    VE = inh_VE_f,
    BW = BW_f,
    AT = 47 * 365
  )
)

for(group in names(inh_data)){
  
  EXPOSURE_REGISTRY <- register_exposure(
      registry = EXPOSURE_REGISTRY,
      group = group,
      pathway = "inhalation",
      chemical = "BPA",
      source = "indoor_air",
      ADD = inh_data[[group]]$ADD,
      inputs = list(
        C = inh_C_BPA,
        VE = inh_data[[group]]$VE,
        BW = inh_data[[group]]$BW
      )
    )
}

#========================================================
# Register water
#========================================================

water_map <- list(
  
  male   = water_result_m,
  
  female = water_result_f,
  
  total  = water_result_t
)

for(group in names(water_map)){
  
  result_obj <- water_map[[group]]
  
  for(chem in names(result_obj)){
    
    EXPOSURE_REGISTRY <- register_exposure(
        registry = EXPOSURE_REGISTRY,
        group = group,
        pathway = "water",
        chemical = chem,
        source = "drinking_water",
        ADD = result_obj[[chem]]$ADD,
        inputs = list(
          C = result_obj[[chem]]$C,
          IR = result_obj[[chem]]$IR,
          BW = result_obj[[chem]]$BW
        )
      )
  }
}

#========================================================
# Register diet
#========================================================

group_map <- c(
  
  male   = "m",
  
  female = "f",
  
  total  = "t"
)

for(food in names(ADD_sim)){
  
  for(chem in names(ADD_sim[[food]])){
    
    for(group in names(group_map)){
      
      g <- group_map[group]
      
      EXPOSURE_REGISTRY <-
        register_exposure(
          registry = EXPOSURE_REGISTRY,
          group = group,
          pathway = "diet",
          chemical = chem,
          source = food,
          ADD = ADD_sim[[food]][[chem]][[g]],
          inputs = list(
            C = Cfood_sim[[food]][[chem]],
            IR = IR_sim[[food]][[g]]
          )
        )
    }
  }
}

#========================================================
# Register dermal
#========================================================

dermal_map <- list(
  
  male   = dermal_result_m,
  
  female = dermal_result_f,
  
  total  = dermal_result_t
)

for(group in names(dermal_map)){
  
  obj <- dermal_map[[group]]
  
  for(name in names(obj)){
    
    split_name <- strsplit(name, "_")[[1]]
    
    chem <- split_name[1]
    
    source <- paste(
      split_name[-1],
      collapse = "_"
    )
    
    EXPOSURE_REGISTRY <-
      register_exposure(
        registry = EXPOSURE_REGISTRY,
        group = group,
        pathway = "dermal",
        chemical = chem,
        source = source,
        ADD = obj[[name]],
        inputs = list()
      )
  }
}

#========================================================
# Automatic chemical discovery
#========================================================

ALL_CHEMICALS <- unique(
  
  unlist(
    
    lapply(
      EXPOSURE_REGISTRY,
      names
    )
  )
)

#========================================================
# Total cumulative ADD
#========================================================

calc_total_ADD <- function(
    registry,
    group,
    chemical
){
  
  chem_obj <- registry[[group]][[chemical]]
  
  if(is.null(chem_obj)){
    return(NULL)
  }
  
  ADDs <- unlist(
    
    lapply(
      chem_obj,
      function(pathway){
        
        lapply(
          pathway,
          `[[`,
          "ADD"
        )
      }
    ),
    
    recursive = FALSE
  )
  
  Reduce("+", ADDs)
}

#========================================================
# Pathway contribution
#========================================================

calc_pathway_contribution <- function(
    registry,
    group,
    chemical
){
  
  chem_obj <- registry[[group]][[chemical]]
  
  if(is.null(chem_obj)){
    return(NULL)
  }
  
  total <- calc_total_ADD(
    registry,
    group,
    chemical
  )
  
  out <- data.frame()
  
  for(pathway in names(chem_obj)){
    
    pathway_sum <- Reduce(
      "+",
      
      lapply(
        chem_obj[[pathway]],
        `[[`,
        "ADD"
      )
    )
    
    tmp <- data.frame(
      
      Pathway = pathway,
      
      Contribution =
        mean(pathway_sum / total) * 100
    )
    
    out <- rbind(out, tmp)
  }
  
  out <- out[
    order(
      out$Contribution,
      decreasing = TRUE
    ),
  ]
  
  rownames(out) <- NULL
  
  return(out)
}

#========================================================
# Extract all sensitivity inputs
#========================================================

extract_inputs <- function(
    registry,
    group,
    chemical
){
  
  chem_obj <- registry[[group]][[chemical]]
  
  if(is.null(chem_obj)){
    return(NULL)
  }
  
  out <- list()
  
  for(pathway in names(chem_obj)){
    
    for(source in names(
      chem_obj[[pathway]]
    )){
      
      tmp <- chem_obj[[pathway]][[source]]
      
      for(v in names(tmp$inputs)){
        
        out[[paste(
          pathway,
          source,
          v,
          sep = "_"
        )]] <- tmp$inputs[[v]]
      }
    }
  }
  
  as.data.frame(out)
}

#========================================================
# PRCC sensitivity analysis
#========================================================

calc_prcc <- function(
    registry,
    group,
    chemical
){
  
  X <- extract_inputs(
    registry,
    group,
    chemical
  )
  
  if(is.null(X)){
    return(NULL)
  }
  
  keep <- sapply(
    X,
    function(x)
      length(unique(x)) > 1
  )
  
  X <- X[, keep, drop = FALSE]
  
  if(ncol(X) == 0){
    return(NULL)
  }
  
  Y <- calc_total_ADD(
    registry,
    group,
    chemical
  )
  
  rank_X <- as.data.frame(
    apply(X, 2, rank)
  )
  
  rank_df <- data.frame(
    rank_X,
    Y = rank(Y)
  )
  
  prcc <- pcor(
    rank_df,
    method = "pearson"
  )
  
  out <- data.frame(
    
    Variable = colnames(rank_X),
    
    PRCC =
      prcc$estimate[
        colnames(rank_X),
        "Y"
      ]
  )
  
  out <- out[
    order(
      abs(out$PRCC),
      decreasing = TRUE
    ),
  ]
  
  rownames(out) <- NULL
  
  return(out)
}

#========================================================
# Automatic cumulative ADD
#========================================================

TOTAL_ADD <- list()

for(group in names(EXPOSURE_REGISTRY)){
  
  TOTAL_ADD[[group]] <- list()
  
  for(chem in ALL_CHEMICALS){
    
    TOTAL_ADD[[group]][[chem]] <-
      tryCatch(
        
        calc_total_ADD(
          EXPOSURE_REGISTRY,
          group,
          chem
        ),
        
        error = function(e) NULL
      )
  }
}

#========================================================
# Automatic contribution analysis
#========================================================

TOTAL_CONTRIBUTION <- list()

for(group in names(EXPOSURE_REGISTRY)){
  
  TOTAL_CONTRIBUTION[[group]] <- list()
  
  for(chem in ALL_CHEMICALS){
    
    TOTAL_CONTRIBUTION[[group]][[chem]] <-
      tryCatch(
        
        calc_pathway_contribution(
          EXPOSURE_REGISTRY,
          group,
          chem
        ),
        
        error = function(e) NULL
      )
  }
}

#========================================================
# Automatic PRCC analysis
#========================================================

TOTAL_PRCC <- list()

for(group in names(EXPOSURE_REGISTRY)){
  
  TOTAL_PRCC[[group]] <- list()
  
  for(chem in ALL_CHEMICALS){
    
    cat(
      "Running PRCC:",
      group,
      chem,
      "\n"
    )
    
    TOTAL_PRCC[[group]][[chem]] <-
      tryCatch(
        
        calc_prcc(
          EXPOSURE_REGISTRY,
          group,
          chem
        ),
        
        error = function(e){
          
          message(
            "PRCC failed:",
            group,
            chem
          )
          
          NULL
        }
      )
  }
}

for(group in names(EXPOSURE_REGISTRY)){
  for(chem in names(EXPOSURE_REGISTRY[[group]])){
    
    X <- extract_inputs(EXPOSURE_REGISTRY, group, chem)
    
    cat("\nCHEM:", chem,
        "\nN cols:", ncol(X),
        "\nUnique check:\n")
    
    print(sapply(X, function(x) length(unique(x))))
  }
}

#========================================================
# Example outputs
#========================================================

TOTAL_ADD$male$BPA

TOTAL_CONTRIBUTION$male$BPA

TOTAL_PRCC$male$BPA

#========================================================
# Save final cumulative engine
#========================================================

save(
  
  EXPOSURE_REGISTRY,
  
  TOTAL_ADD,
  
  TOTAL_CONTRIBUTION,
  
  TOTAL_PRCC,
  
  file = "Cumulative_Exposure_Engine.RData"
)


#========================================
# Rebuild all_summary_df (MUST exist before plotting)
#========================================

all_summary <- list()

for(food in names(ADD_sim)){
  for(chem in names(ADD_sim[[food]])){
    
    tmp <- data.frame(
      Food = food,
      Chemical = chem,
      Group = c("Male", "Female", "Total"),
      Mean = c(
        mean(ADD_sim[[food]][[chem]]$m),
        mean(ADD_sim[[food]][[chem]]$f),
        mean(ADD_sim[[food]][[chem]]$t)
      ),
      P95 = c(
        unname(quantile(ADD_sim[[food]][[chem]]$m, 0.95)),
        unname(quantile(ADD_sim[[food]][[chem]]$f, 0.95)),
        unname(quantile(ADD_sim[[food]][[chem]]$t, 0.95))
      )
      )
    
    all_summary[[paste(food, chem)]] <- tmp
  }
}

all_summary_df <- dplyr::bind_rows(all_summary)

#========================================
# Rebuild plotting dataset 
#========================================

plot_df <- all_summary_df %>%
  tidyr::pivot_longer(
    cols = c(Mean, P95),
    names_to = "Metric",
    values_to = "Value"
  )

plot_df <- plot_df %>%
  mutate(
    Food = as.character(Food),
    Chemical = as.character(Chemical),
    Group = as.character(Group),
    Metric = as.character(Metric),
    Value = as.numeric(Value)
  )


#========================================
# Food category definitions
#========================================

# General food categories
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

# Canned food subcategories
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
# Unified theme (centered title)
#========================================

theme_center <- theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )


#========================================
# Toxicological reference values
# Unit: mg/kg/day
#========================================

TOX <- c(
  
  BPA  = 2e-7,
  `BP-1` = 3e-5,
  `BP-3` = 2e-3,
  `BP-8` = 3e-5,
  DEET = 1e-4,
  MEP  = 1e-2
)



#========================================
# Hazard Quotient (HQ)
#========================================

HQ <- list()

for(sex in names(TOTAL_ADD)){
  
  HQ[[sex]] <- list()
  
  for(chem in names(TOTAL_ADD[[sex]])){
    
    if(!chem %in% names(TOX)) next
    
    ADD <- TOTAL_ADD[[sex]][[chem]]
    
    if(is.null(ADD) || length(ADD) == 0) next
    
    HQ[[sex]][[chem]] <- ADD / TOX[[chem]]
  }
}


#========================================
# Hazard Index (HI)
#========================================
HI <- list()

for(sex in names(HQ)){
  
  valid_hq <- HQ[[sex]][
    sapply(HQ[[sex]], function(x) !is.null(x))
  ]
  
  if(length(valid_hq) == 0){
    HI[[sex]] <- NULL
    next
  }
  
  HI[[sex]] <- Reduce(`+`, valid_hq)
}


#========================================
# HQ exceedance probability
#========================================

HI_EXCEED <- data.frame(
  
  Sex = names(HI),
  
  Exceed_Prob = sapply(HI, function(x){
    
    if(is.null(x) || length(x) == 0) return(NA)
    
    mean(x > 1, na.rm = TRUE) * 100
  })
)


#========================================
# HI exceedance probability
#========================================

HI_EXCEED <- data.frame(
  
  Sex = names(HI),
  
  Exceed_Prob = sapply(
    HI,
    function(x){
      
      mean(x > 1) * 100
    }
  )
)

HI_EXCEED


#========================================
# Risk summary function
#========================================

make_risk_summary <- function(x){
  
  data.frame(
    
    Mean = mean(x),
    
    Median = median(x),
    
    P95 = quantile(x, 0.95),
    
    P99 = quantile(x, 0.99),
    
    Max = max(x),
    
    Exceed_Prob = mean(x > 1) * 100
  )
}


#========================================
# HQ summary
#========================================

HQ_SUMMARY <- list()

for(sex in names(HQ)){
  
  HQ_SUMMARY[[sex]] <- do.call(
    
    rbind,
    
    lapply(
      
      names(HQ[[sex]]),
      
      function(chem){
        
        tmp <- make_risk_summary(
          HQ[[sex]][[chem]]
        )
        
        tmp$Chemical <- chem
        
        tmp
      }
    )
  )
}

HQ_SUMMARY$m



#========================================
# HI summary
#========================================

HI_SUMMARY <- do.call(
  
  rbind,
  
  lapply(
    
    names(HI),
    
    function(sex){
      
      tmp <- make_risk_summary(HI[[sex]])
      
      tmp$Sex <- sex
      
      tmp
    }
  )
)

HI_SUMMARY


library(ggplot2)

#========================================
# Prepare HQ plotting dataframe
#========================================

HQ_PLOT <- list()

for(sex in names(HQ)){
  
  for(chem in names(HQ[[sex]])){
    
    HQ_PLOT[[paste(sex, chem)]] <- data.frame(
      
      Sex = sex,
      
      Chemical = chem,
      
      HQ = HQ[[sex]][[chem]]
    )
  }
}

HQ_PLOT_DF <- do.call(
  rbind,
  HQ_PLOT
)

#========================================
# Figure 1: HQ boxplot
#========================================

ggplot(HQ_PLOT_DF, aes(x = Chemical, y = HQ, fill = Sex)) +
  
  geom_boxplot(width = 0.5, alpha = 0.5, color = "black", # 箱型圖

               position = position_dodge(width = 0.7),  # 將不同性別分開

               outlier.size = 1 # 界外值(outlier)黑點縮小
              ) +  
  # 平均值
  stat_summary(fun = mean, geom = "point", shape = 19, size = 1.5, color = "red",
               position = position_dodge(width = 0.7) # 跟著箱型圖位置
              ) + 
  
  # HQ = 1 閾值線
  geom_hline(yintercept = 1, linetype = 2, color = "red", linewidth = 1) +
  
  # log scale
  scale_y_log10() +
  
  # 轉向
  coord_flip() +
  
  theme_bw(base_size = 13) +

  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5), # 主標題
    axis.title = element_text(size = 18), # 座標軸標題
    axis.text = element_text(size = 16), # 座標軸刻度文字
    legend.title = element_text(size = 16), # legend標題
    legend.text = element_text(size = 14) # legend內容
  ) +
  
  labs(
    title = "Hazard Quotient Distribution",
    y = "HQ (log scale)",
    x = ""
  )

#========================================
# HI plotting dataframe
#========================================

HI_DF <- do.call(
  
  rbind,
  
  lapply(
    
    names(HI),
    
    function(sex){
      
      data.frame(
        
        Sex = sex,
        
        HI = HI[[sex]]
      )
    }
  )
)

#========================================
# Figure 2: HI distribution
#========================================
ggplot(
  HI_DF,
  aes(x = HI, fill = Sex)
) +
  
  geom_density(alpha = 0.4) +
  
  geom_vline(xintercept = 1, linetype = 2, color = "red", linewidth = 1) +
  
  scale_x_log10() +
  
  theme_bw(base_size = 14) +
  
  labs(
    title = "Hazard Index Distribution",
    x = "HI (log scale)",
    y = "Density"
  ) +
  
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5), # 主標題
    axis.title = element_text(size = 18), # 座標軸標題
    axis.text = element_text(size = 16), # 座標軸刻度文字
    legend.title = element_text(size = 16), # legend標題
    legend.text = element_text(size = 14) # legend內容
  ) 

#========================================================
# PRCC Tornado Plot
#========================================================

prcc_df <- list()

for(group in names(TOTAL_PRCC)){
  
  for(chem in names(TOTAL_PRCC[[group]])){
    
    tmp <- TOTAL_PRCC[[group]][[chem]]
    
    if(is.null(tmp)) next
    
    tmp$Group <- group
    
    tmp$Chemical <- chem
    
    prcc_df[[paste(group, chem)]] <- tmp
  }
}

prcc_df <- bind_rows(prcc_df)

#========================================================
# Keep top variables only
#========================================================

prcc_top <- prcc_df %>%
  
  group_by(
    Group,
    Chemical
  ) %>%
  
  slice_max(
    order_by = abs(PRCC),
    n = 8
  )

#========================================================
# Figure 4: PRCC Sensitivity Analysis
#========================================================
p4 <- ggplot(prcc_top, aes(
  x = reorder(Variable, abs(PRCC)),
  y = PRCC,
  fill = PRCC
)) +
  geom_col() +
  coord_flip() +
  facet_grid(Group ~ Chemical, scales = "free_y") +
  
  # 顏色 scale：低值深藍，中間黃，高值深紅
  scale_fill_gradient2(
    low = "#2166AC",      # 深藍
    mid = "#FFFF00",      # 黃色
    high = "#B2182B",     # 深紅
    midpoint = 0,
    limits = c(-0.9, 0.9),
    oob = scales::squish
  ) +
  
  # Y scale: 自動依照資料生成刻度
  scale_y_continuous(
    limits = c(-0.9, 0.9),
    breaks = pretty(c(prcc_top$PRCC, -0.9, 0.9))
  ) +
  
  labs(
    title = "PRCC Sensitivity Analysis",
    x = "",
    y = "PRCC"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p4)


