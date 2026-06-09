
library(tidyverse)
library(ggpubr)
library(rstatix)
# 主函数：分析临床变量与表达量的差异（修正版）
#' Function Title
#'
#' Function Description
#' @export
analyze_clinical_expression <- function(data, clinical_vars, expression_cols) {
  # data: 数据框，包含临床变量和表达量
  # clinical_vars: 临床变量名称向量（字符型）
  # expression_cols: 表达量列名称向量（字符型）
  
  results <- list()
  plots <- list()
  
  for (clin_var in clinical_vars) {
    for (expr_col in expression_cols) {
      # 准备数据：去除缺失值，并确保临床变量为因子
      sub_data <- data[, c(clin_var, expr_col)] %>%
        drop_na() %>%
        mutate(!!sym(clin_var) := as.factor(!!sym(clin_var)))
      
      # 获取临床变量的唯一值个数
      n_levels <- nlevels(sub_data[[clin_var]])
      
      # 如果只有一个水平，跳过
      if (n_levels < 2) {
        next
      }
      
      # 根据水平数选择统计方法
      if (n_levels == 2) {
        # Wilcoxon秩和检验
        test_result <- wilcox.test(sub_data[[expr_col]] ~ sub_data[[clin_var]], data = sub_data)
        
        # 计算效应量（r值）
        z_value <- qnorm(test_result$p.value/2)
        r_effect <- abs(z_value) / sqrt(nrow(sub_data))
        
        # 计算中位数和四分位数
        stats <- sub_data %>%
          group_by(!!sym(clin_var)) %>%
          summarise(
            median = median(!!sym(expr_col), na.rm = TRUE),
            q1 = quantile(!!sym(expr_col), 0.25, na.rm = TRUE),
            q3 = quantile(!!sym(expr_col), 0.75, na.rm = TRUE),
            n = n()
          )
        
        # 构建结果行
        stats_text <- paste0(
          levels(sub_data[[clin_var]])[1], ": median=", round(stats$median[1], 2), 
          " (Q1-Q3: ", round(stats$q1[1], 2), "-", round(stats$q3[1], 2), "), n=", stats$n[1], "\n",
          levels(sub_data[[clin_var]])[2], ": median=", round(stats$median[2], 2), 
          " (Q1-Q3: ", round(stats$q1[2], 2), "-", round(stats$q3[2], 2), "), n=", stats$n[2]
        )
        
        result_row <- data.frame(
          Clinical_Variable = clin_var,
          Expression = expr_col,
          Test_Method = "Wilcoxon",
          Levels_Compared = paste(levels(sub_data[[clin_var]]), collapse = " vs "),
          Statistic = test_result$statistic,
          P_value = test_result$p.value,
          Effect_Size = r_effect,
          Significance = case_when(
            test_result$p.value < 0.001 ~ "***",
            test_result$p.value < 0.01 ~ "**",
            test_result$p.value < 0.05 ~ "*",
            TRUE ~ "ns"
          ),
          Statistics_Summary = stats_text,
          Tukey_HSD = NA,
          stringsAsFactors = FALSE
        )
        
      } else if (n_levels >= 3) {
        # ANOVA分析 - 修正点：创建正确的公式
        formula_obj <- as.formula(paste(expr_col, "~", clin_var))
        aov_result <- aov(formula_obj, data = sub_data)
        anova_summary <- summary(aov_result)
        p_value <- anova_summary[[1]][["Pr(>F)"]][1]
        
        # 计算效应量（eta squared）
        ss_effect <- anova_summary[[1]][["Sum Sq"]][1]
        ss_total <- sum(anova_summary[[1]][["Sum Sq"]])
        eta_squared <- ss_effect / ss_total
        
        # 两两比较（Tukey HSD）
        tukey_result <- TukeyHSD(aov_result)
        pairwise_df <- as.data.frame(tukey_result[[1]])
        pairwise_df$Comparison <- rownames(pairwise_df)
        pairwise_df$P_adj <- pairwise_df$`p adj`
        pairwise_df <- pairwise_df[, c("Comparison", "diff", "lwr", "upr", "P_adj")]
        
        # 计算各组统计量
        stats <- sub_data %>%
          group_by(!!sym(clin_var)) %>%
          summarise(
            mean = mean(!!sym(expr_col), na.rm = TRUE),
            sd = sd(!!sym(expr_col), na.rm = TRUE),
            median = median(!!sym(expr_col), na.rm = TRUE),
            q1 = quantile(!!sym(expr_col), 0.25, na.rm = TRUE),
            q3 = quantile(!!sym(expr_col), 0.75, na.rm = TRUE),
            n = n()
          )
        
        # 构建Tukey HSD结果的文本摘要
        tukey_summary <- pairwise_df %>%
          mutate(Significant = ifelse(P_adj < 0.05, "Yes", "No")) %>%
          select(Comparison, diff, P_adj, Significant)
        
        result_row <- data.frame(
          Clinical_Variable = clin_var,
          Expression = expr_col,
          Test_Method = "ANOVA",
          Levels_Compared = paste(levels(sub_data[[clin_var]]), collapse = ", "),
          Statistic = anova_summary[[1]][["F value"]][1],
          P_value = p_value,
          Effect_Size = eta_squared,
          Significance = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            TRUE ~ "ns"
          ),
          Statistics_Summary = paste0(capture.output(print(stats)), collapse = "\n"),
          Tukey_HSD = paste(capture.output(print(tukey_summary)), collapse = "\n"),
          stringsAsFactors = FALSE
        )
      }
      
      # 存储结果
      results[[paste(clin_var, expr_col, sep = "_")]] <- result_row
      
      # 生成箱线图
      # 准备显著性标记数据
      if (n_levels == 2) {
        p_val_label <- paste("Wilcoxon p =", format(test_result$p.value, scientific = TRUE, digits = 3))
        if (test_result$p.value < 0.05) {
          p_val_label <- paste0(p_val_label, " ", 
                                ifelse(test_result$p.value < 0.001, "***",
                                       ifelse(test_result$p.value < 0.01, "**", "*")))
        }
      } else {
        p_val_label <- paste("ANOVA p =", format(p_value, scientific = TRUE, digits = 3))
        if (p_value < 0.05) {
          p_val_label <- paste0(p_val_label, " ",
                                ifelse(p_value < 0.001, "***",
                                       ifelse(p_value < 0.01, "**", "*")))
        }
      }
      
      p <- ggplot(sub_data, aes(x = !!sym(clin_var), y = !!sym(expr_col))) +
        geom_boxplot(aes(fill = !!sym(clin_var)), alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
        stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
        labs(
          title = paste("Expression of", expr_col, "by", clin_var),
          x = clin_var,
          y = paste(expr_col, "Expression Level"),
          caption = p_val_label
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 11),
          plot.caption = element_text(hjust = 0, face = "italic", size = 9),
          legend.position = "none"
        ) +
        scale_fill_brewer(palette = "Set2")
      
      # 对于ANOVA且显著的情况，添加两两比较的显著性标记
      if (n_levels >= 3 && p_value < 0.05) {
        sig_pairs <- pairwise_df[pairwise_df$P_adj < 0.05, ]
        if (nrow(sig_pairs) > 0) {
          # 计算y轴范围用于放置文本
          y_range <- range(sub_data[[expr_col]], na.rm = TRUE)
          y_max <- y_range[2]
          y_min <- y_range[1]
          y_step <- (y_max - y_min) * 0.05
          
          # 构建显著性标记文本
          sig_text <- paste(sig_pairs$Comparison, 
                            sprintf("(p=%.3f)", sig_pairs$P_adj), 
                            collapse = "\n")
          
          p <- p + annotate("text", 
                            x = length(levels(sub_data[[clin_var]])) / 2, 
                            y = y_max + y_step,
                            label = sig_text,
                            size = 3, 
                            color = "red",
                            hjust = 0.5)
        }
      }
      
      plots[[paste(clin_var, expr_col, sep = "_")]] <- p
    }
  }
  
  # 合并所有统计结果
  if (length(results) > 0) {
    result_df <- bind_rows(results)
    
    # 调整p值（FDR校正）
    result_df$P_adjusted <- p.adjust(result_df$P_value, method = "fdr")
    result_df$Significance_Adj <- case_when(
      result_df$P_adjusted < 0.001 ~ "***",
      result_df$P_adjusted < 0.01 ~ "**",
      result_df$P_adjusted < 0.05 ~ "*",
      TRUE ~ "ns"
    )
    
    # 重新排列列的顺序
    result_df <- result_df %>%
      select(Clinical_Variable, Expression, Test_Method, Levels_Compared,
             Statistic, P_value, P_adjusted, Effect_Size, Significance, Significance_Adj,
             everything())
  } else {
    result_df <- data.frame()
    warning("No valid comparisons were performed!")
  }
  
  # 返回结果
  return(list(
    statistics = result_df,
    plots = plots
  ))
}
