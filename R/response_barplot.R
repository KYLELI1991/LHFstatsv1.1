#' Plot ORR and DCR Rates by Group
#'
#' @export
plot_response_rates <- function(
    df,
    group_col = "group",
    response_col = "response",
    group_labels = NULL,
    colors = NULL,
    p_value_color = "black",
    p_value_shape = "bracket",
    title = "ORR and DCR Rates by Group",
    subtitle = "with 95% Confidence Intervals",
    x_label = NULL,
    y_label = "Response Rate",
    legend_title = "Group",
    font_size = 12,
    show_ci_labels = TRUE,
    show_n_labels = TRUE,
    ci_label_position = "above",
    ci_label_size = 3,
    p_adjust_method = "none",
    response_levels = c("CR", "PR", "SD", "PD"),
    orr_levels = c("CR", "PR"),
    dcr_levels = c("CR", "PR", "SD"),

    # 新增控制参数
    show_comparisons = TRUE,
    max_comparisons_show = NULL,
    comparison_start_y = 0.85,
    comparison_spacing = 0.08,
    comparison_text_size = 3,
    comparison_line_size = 0.8,

    rotate_x_labels = FALSE,
    x_label_angle = 45,

    bar_width = 0.25,
    bar_alpha = 0.8,

    show_orr_only = FALSE,
    show_dcr_only = FALSE,

    # ===== 新增参数 =====
    show_cr = FALSE,        # 是否在ORR柱子中叠加显示CR部分（深色）
    show_dcr = TRUE,        # 是否显示DCR柱子

    y_max_manual = NULL,
    y_breaks = seq(0, 1, by = 0.2),

    show_legend = TRUE,
    legend_position = "top",

    output_stats_only = FALSE,

    show_p_value_table = FALSE,
    p_value_table_position = "bottom"
) {

  # 加载必要的包
  if (!require(dplyr)) stop("需要安装dplyr包")
  if (!require(tidyr)) stop("需要安装tidyr包")
  if (!require(ggplot2)) stop("需要安装ggplot2包")
  if (!require(rlang)) stop("需要安装rlang包")
  if (show_p_value_table && !require(gridExtra)) {
    warning("显示P值表格需要安装gridExtra包，将不显示表格")
    show_p_value_table <- FALSE
  }
  if (show_p_value_table && !require(grid)) {
    warning("显示P值表格需要安装grid包，将不显示表格")
    show_p_value_table <- FALSE
  }

  # 检查必要的列
  if (!group_col %in% colnames(df)) stop(paste("组别列", group_col, "不存在于数据框中"))
  if (!response_col %in% colnames(df)) stop(paste("响应列", response_col, "不存在于数据框中"))

  # 重命名
  df <- df %>% rename(group = all_of(group_col), response = all_of(response_col))

  if (!is.character(df$response) && !is.factor(df$response)) stop("响应列必须是字符型或因子型")
  if (!is.factor(df$response)) df$response <- factor(df$response, levels = response_levels)

  groups <- unique(df$group)
  n_groups <- length(groups)
  if (n_groups < 2) stop("至少需要2个组进行比较")

  # ===== 确定显示哪些类别 =====
  # show_dcr_only 优先级最高；否则 show_dcr=FALSE 等价于只显示ORR
  if (show_dcr_only) {
    display_orr <- FALSE
    display_dcr <- TRUE
  } else if (show_orr_only || !show_dcr) {
    display_orr <- TRUE
    display_dcr <- FALSE
  } else {
    display_orr <- TRUE
    display_dcr <- TRUE
  }

  if (show_orr_only && show_dcr_only) {
    warning("show_orr_only和show_dcr_only不能同时为TRUE，将显示两者")
    display_orr <- TRUE
    display_dcr <- TRUE
  }

  # 组标签
  if (is.null(group_labels)) {
    group_labels <- setNames(as.character(groups), groups)
  } else if (length(group_labels) != n_groups) {
    warning("组标签数量与组数不匹配，将使用组名作为标签")
    group_labels <- setNames(as.character(groups), groups)
  }

  # 颜色
  if (is.null(colors)) {
    colors <- setNames(scales::hue_pal()(n_groups), groups)
  } else if (length(colors) < n_groups) {
    warning("颜色数量不足，将使用默认颜色")
    colors <- setNames(scales::hue_pal()(n_groups), groups)
  }

  # 1. 计算 ORR / DCR / CR
  df <- df %>%
    mutate(
      ORR = ifelse(response %in% orr_levels, 1, 0),
      DCR = ifelse(response %in% dcr_levels, 1, 0),
      CR  = ifelse(response %in% "CR", 1, 0)
    )

  cat("数据摘要:\n")
  cat("组别分布:\n"); print(table(df$group))
  cat("\n响应分布:\n"); print(table(df$response))
  cat("\nORR定义:", paste(orr_levels, collapse = ", "), "\n")
  cat("DCR定义:", paste(dcr_levels, collapse = ", "), "\n")
  if (show_cr) cat("CR将在ORR柱子中以深色叠加显示\n")

  # 2. 按组汇总
  summary_stats <- df %>%
    group_by(group) %>%
    summarise(
      n_total = n(),
      n_CR    = sum(CR),
      n_ORR   = sum(ORR),
      n_DCR   = sum(DCR),
      CR_rate  = mean(CR),
      ORR_rate = mean(ORR),
      DCR_rate = mean(DCR),
      ORR_se = sqrt(mean(ORR) * (1 - mean(ORR)) / n()),
      DCR_se = sqrt(mean(DCR) * (1 - mean(DCR)) / n()),
      .groups = "drop"
    )

  # 3. 置信区间
  calculate_ci <- function(success, total) {
    if (total == 0) return(c(0, 0))
    if (success == 0) return(c(0, 0))
    if (success == total) return(c(1, 1))
    test <- prop.test(success, total, conf.level = 0.95)
    c(test$conf.int[1], test$conf.int[2])
  }

  cis <- df %>%
    group_by(group) %>%
    summarise(
      ORR_ci_lower = calculate_ci(sum(ORR), n())[1],
      ORR_ci_upper = calculate_ci(sum(ORR), n())[2],
      DCR_ci_lower = calculate_ci(sum(DCR), n())[1],
      DCR_ci_upper = calculate_ci(sum(DCR), n())[2],
      .groups = "drop"
    )

  final_stats <- left_join(summary_stats, cis, by = "group") %>%
    mutate(
      PR_part_rate = ORR_rate - CR_rate   # ORR中PR部分的比例
    )

  # 4. P值（两两比较）
  group_pairs <- combn(groups, 2, simplify = FALSE)
  n_comparisons <- length(group_pairs)
  p_values_raw <- list()

  for (i in seq_along(group_pairs)) {
    pair <- group_pairs[[i]]
    pair_name <- paste(pair, collapse = "_vs_")
    stats_pair <- final_stats %>% filter(group %in% pair)

    orr_test <- tryCatch(
      prop.test(x = stats_pair$n_ORR, n = stats_pair$n_total),
      error = function(e) { warning(paste("ORR比较", pair_name, "失败:", e$message)); list(p.value = NA) }
    )
    dcr_test <- tryCatch(
      prop.test(x = stats_pair$n_DCR, n = stats_pair$n_total),
      error = function(e) { warning(paste("DCR比较", pair_name, "失败:", e$message)); list(p.value = NA) }
    )

    p_values_raw[[pair_name]] <- list(
      pair = pair,
      ORR_raw = ifelse(is.null(orr_test$p.value), NA, orr_test$p.value),
      DCR_raw = ifelse(is.null(dcr_test$p.value), NA, dcr_test$p.value),
      ORR_test = orr_test,
      DCR_test = dcr_test
    )
  }

  p_values_adjusted <- list()
  if (p_adjust_method != "none" && n_comparisons > 1) {
    all_orr_p <- sapply(p_values_raw, function(x) x$ORR_raw)
    all_dcr_p <- sapply(p_values_raw, function(x) x$DCR_raw)
    valid_orr <- !is.na(all_orr_p)
    valid_dcr <- !is.na(all_dcr_p)
    adjusted_orr_p <- rep(NA, length(all_orr_p))
    adjusted_dcr_p <- rep(NA, length(all_dcr_p))
    if (sum(valid_orr) > 0) adjusted_orr_p[valid_orr] <- p.adjust(all_orr_p[valid_orr], method = p_adjust_method)
    if (sum(valid_dcr) > 0) adjusted_dcr_p[valid_dcr] <- p.adjust(all_dcr_p[valid_dcr], method = p_adjust_method)
    for (i in seq_along(p_values_raw)) {
      pair_name <- names(p_values_raw)[i]
      p_values_adjusted[[pair_name]] <- list(
        pair = p_values_raw[[i]]$pair,
        ORR = adjusted_orr_p[i], DCR = adjusted_dcr_p[i],
        ORR_raw = p_values_raw[[i]]$ORR_raw, DCR_raw = p_values_raw[[i]]$DCR_raw,
        ORR_test = p_values_raw[[i]]$ORR_test, DCR_test = p_values_raw[[i]]$DCR_test
      )
    }
  } else {
    for (i in seq_along(p_values_raw)) {
      pair_name <- names(p_values_raw)[i]
      p_values_adjusted[[pair_name]] <- list(
        pair = p_values_raw[[i]]$pair,
        ORR = p_values_raw[[i]]$ORR_raw, DCR = p_values_raw[[i]]$DCR_raw,
        ORR_raw = p_values_raw[[i]]$ORR_raw, DCR_raw = p_values_raw[[i]]$DCR_raw,
        ORR_test = p_values_raw[[i]]$ORR_test, DCR_test = p_values_raw[[i]]$DCR_test
      )
    }
  }

  p_value_table_data <- NULL
  if (show_p_value_table) {
    p_value_table_data <- data.frame(
      比较 = names(p_values_adjusted),
      ORR_P值 = sapply(p_values_adjusted, function(x) {
        if (is.na(x$ORR)) return("NA")
        if (x$ORR < 0.001) return("<0.001")
        sprintf("%.3f", x$ORR)
      }),
      DCR_P值 = sapply(p_values_adjusted, function(x) {
        if (is.na(x$DCR)) return("NA")
        if (x$DCR < 0.001) return("<0.001")
        sprintf("%.3f", x$DCR)
      }),
      stringsAsFactors = FALSE
    )
    p_value_table_data$比较 <- sapply(p_values_adjusted, function(x) {
      paste(group_labels[x$pair[1]], "vs", group_labels[x$pair[2]])
    })
  }

  if (output_stats_only) {
    return(list(
      stats = final_stats, p_values = p_values_adjusted,
      comparisons = group_pairs, n_comparisons = n_comparisons,
      p_adjust_method = p_adjust_method, p_value_table = p_value_table_data
    ))
  }

  # 5. 准备绘图数据
  categories_to_show <- c()
  category_labels <- c()
  x_base_map <- c()

  if (display_orr) {
    categories_to_show <- c(categories_to_show, "ORR_rate")
    category_labels <- c(category_labels, "ORR")
    x_base_map["ORR_rate"] <- length(x_base_map) + 1
  }
  if (display_dcr) {
    categories_to_show <- c(categories_to_show, "DCR_rate")
    category_labels <- c(category_labels, "DCR")
    x_base_map["DCR_rate"] <- length(x_base_map) + 1
  }

  plot_data <- final_stats %>%
    pivot_longer(cols = all_of(categories_to_show), names_to = "category", values_to = "rate") %>%
    mutate(
      ci_lower = case_when(
        category == "ORR_rate" ~ ORR_ci_lower,
        category == "DCR_rate" ~ DCR_ci_lower
      ),
      ci_upper = case_when(
        category == "ORR_rate" ~ ORR_ci_upper,
        category == "DCR_rate" ~ DCR_ci_upper
      ),
      category_label = factor(
        ifelse(category == "ORR_rate", "ORR", "DCR"),
        levels = category_labels, ordered = TRUE
      ),
      n_events = ifelse(category == "ORR_rate", n_ORR, n_DCR),
      ci_text = sprintf("%.1f-%.1f%%", ci_lower * 100, ci_upper * 100),
      ci_y_position = case_when(
        ci_label_position == "middle" ~ rate * 0.5,
        ci_label_position == "above" ~ ci_upper + 0.02,
        ci_label_position == "below" ~ ci_lower - 0.02,
        ci_label_position == "none" ~ NA_real_,
        TRUE ~ ci_upper + 0.02
      ),
      x_base = x_base_map[category]
    )

  # 组内偏移
  group_positions <- data.frame(
    group = groups,
    pos_offset = seq(-bar_width, bar_width, length.out = n_groups)
  )

  plot_data <- plot_data %>%
    left_join(group_positions, by = "group") %>%
    mutate(x_position = x_base + pos_offset)

  # ===== CR 叠加数据 =====
  cr_overlay_data <- NULL
  if (show_cr && display_orr) {
    cr_overlay_data <- plot_data %>%
      filter(category == "ORR_rate") %>%
      mutate(
        cr_rate = CR_rate,
        pr_rate = PR_part_rate
      )
  }

  # y_max
  if (!is.null(y_max_manual)) {
    y_max <- y_max_manual
  } else {
    y_max <- max(plot_data$ci_upper, na.rm = TRUE) * 1.3
    if (show_comparisons && n_comparisons > 1 && p_value_shape != "none") {
      n_show <- ifelse(is.null(max_comparisons_show), n_comparisons,
                       min(n_comparisons, max_comparisons_show))
      y_max <- max(y_max, max(plot_data$ci_upper, na.rm = TRUE) * (1 + n_show * comparison_spacing * 0.8))
    }
  }

  # 6. 基础图形
  p <- ggplot() +
    geom_bar(
      data = plot_data,
      aes(x = x_position, y = rate, fill = group),
      stat = "identity",
      width = bar_width * 0.8,
      alpha = bar_alpha,
      position = position_identity()
    )

  # ===== 叠加 CR 柱子 =====
  if (!is.null(cr_overlay_data) && nrow(cr_overlay_data) > 0) {
    p <- p +
      geom_bar(
        data = cr_overlay_data,
        aes(x = x_position, y = cr_rate, fill = group),
        stat = "identity",
        width = bar_width * 0.8,
        alpha = 1,
        position = position_identity()
      )
  }

  # 误差线
  p <- p +
    geom_errorbar(
      data = plot_data,
      aes(x = x_position, ymin = ci_lower, ymax = ci_upper),
      width = bar_width * 0.3, size = 0.8, color = "black"
    )

  # 总百分比标签（柱顶）
  p <- p +
    geom_text(
      data = plot_data,
      aes(x = x_position, y = rate, label = sprintf("%.1f%%", rate * 100)),
      vjust = -1.2, size = font_size * 0.3, fontface = "bold"
    )

  # ===== CR / PR 内部标签 =====
  if (!is.null(cr_overlay_data) && nrow(cr_overlay_data) > 0) {
    # CR 标签（仅在CR比例足够大时显示）
    cr_label_data <- cr_overlay_data %>% filter(cr_rate >= 0.05)
    if (nrow(cr_label_data) > 0) {
      p <- p +
        geom_text(
          data = cr_label_data,
          aes(x = x_position, y = cr_rate / 2,
              label = sprintf("CR: %.1f%%", cr_rate * 100)),
          size = font_size * 0.22, color = "white", fontface = "bold"
        )
    }
    # PR 标签
    pr_label_data <- cr_overlay_data %>% filter(pr_rate >= 0.05)
    if (nrow(pr_label_data) > 0) {
      p <- p +
        geom_text(
          data = pr_label_data,
          aes(x = x_position, y = cr_rate + pr_rate / 2,
              label = sprintf("PR: %.1f%%", pr_rate * 100)),
          size = font_size * 0.22, color = "white", fontface = "bold"
        )
    }
  }

  # CI 标签
  if (show_ci_labels && ci_label_position != "none") {
    ci_label_data <- plot_data %>% filter(!is.na(ci_y_position))
    if (nrow(ci_label_data) > 0) {
      p <- p +
        geom_text(
          data = ci_label_data,
          aes(x = x_position, y = ci_y_position, label = ci_text),
          size = ci_label_size, color = "darkblue", fontface = "bold",
          vjust = ifelse(ci_label_position == "above", -0.2,
                         ifelse(ci_label_position == "below", 1.2, 0.5))
        )
    }
  }

  # 样本量标签
  if (show_n_labels) {
    # 如果显示CR，样本量放到柱子顶部（避免与CR/PR标签重叠）
    n_y_pos <- if (show_cr && display_orr) {
      plot_data$rate * 0.95
    } else {
      plot_data$rate / 2
    }
    p <- p +
      geom_text(
        data = plot_data,
        aes(x = x_position, y = n_y_pos, label = sprintf("%d/%d", n_events, n_total)),
        size = font_size * 0.25, color = "white", fontface = "bold"
      )
  }

  # 7. P值标注
  if (show_comparisons && p_value_shape != "none" && n_comparisons > 0) {
    n_show <- ifelse(is.null(max_comparisons_show), n_comparisons,
                     min(n_comparisons, max_comparisons_show))
    if (n_show > 0) {
      p_values_list <- p_values_adjusted[1:n_show]
      p_values_list <- p_values_list[order(
        sapply(p_values_list, function(x) min(x$ORR, x$DCR, na.rm = TRUE))
      )]
      base_y <- y_max * comparison_start_y

      for (i in seq_along(p_values_list)) {
        comparison <- p_values_list[[i]]
        pair <- comparison$pair
        pos1 <- group_positions$pos_offset[group_positions$group == pair[1]]
        pos2 <- group_positions$pos_offset[group_positions$group == pair[2]]
        y_line <- base_y - (i - 1) * comparison_spacing * y_max

        x_positions <- list()
        if (display_orr) x_positions$orr <- x_base_map["ORR_rate"]
        if (display_dcr) x_positions$dcr <- x_base_map["DCR_rate"]

        for (cat_name in names(x_positions)) {
          cat_x_base <- x_positions[[cat_name]]
          p_value <- ifelse(cat_name == "orr", comparison$ORR, comparison$DCR)
          if (is.na(p_value)) next

          x_start <- cat_x_base + pos1
          x_end <- cat_x_base + pos2
          x_mid <- (x_start + x_end) / 2

          if (p_value_shape %in% c("bracket", "line", "dashed", "arrow")) {
            line_type <- ifelse(p_value_shape == "dashed", "dashed", "solid")
            p <- p +
              annotate("segment", x = x_start, xend = x_end, y = y_line, yend = y_line,
                       size = comparison_line_size, color = p_value_color, linetype = line_type)
            if (p_value_shape == "bracket") {
              p <- p +
                annotate("segment", x = x_start, xend = x_start,
                         y = y_line - 0.01, yend = y_line + 0.01,
                         size = comparison_line_size, color = p_value_color) +
                annotate("segment", x = x_end, xend = x_end,
                         y = y_line - 0.01, yend = y_line + 0.01,
                         size = comparison_line_size, color = p_value_color)
            }
            if (p_value_shape == "arrow") {
              p <- p +
                annotate("segment", x = x_start, xend = x_end,
                         y = y_line, yend = y_line,
                         size = comparison_line_size, color = p_value_color,
                         arrow = arrow(ends = "both", length = unit(0.1, "inches")))
            }
            p <- p +
              annotate("text", x = x_mid, y = y_line + 0.02,
                       label = ifelse(p_value < 0.001, "P < 0.001",
                                      paste0("P = ", sprintf("%.3f", p_value))),
                       size = comparison_text_size, fontface = "bold", color = p_value_color)
          }
        }
      }
    }
  }

  # 8. 多重比较校正说明
  if (p_adjust_method != "none" && n_comparisons > 1 && show_comparisons) {
    method_names <- c("bonferroni" = "Bonferroni", "holm" = "Holm",
                      "hochberg" = "Hochberg", "hommel" = "Hommel",
                      "fdr" = "FDR (Benjamini-Hochberg)", "BH" = "Benjamini-Hochberg",
                      "BY" = "Benjamini-Yekutieli")
    method_label <- ifelse(p_adjust_method %in% names(method_names),
                           method_names[p_adjust_method], p_adjust_method)
    p <- p + labs(caption = paste0("P值已使用", method_label, "方法进行多重比较校正"))
  }

  # 9. X轴
  x_breaks <- sort(unique(plot_data$x_base))
  x_labels <- category_labels[order(x_base_map)]
  p <- p +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = c(min(x_breaks) - 0.5, max(x_breaks) + 0.5),
      expand = expansion(mult = 0.1)
    )

  # 10. Y轴
  p <- p +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, y_max),
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0.05))
    )

  # 11. 颜色
  p <- p +
    scale_fill_manual(values = colors, labels = group_labels, name = legend_title)

  # 12. 标题（如果显示CR，在副标题中自动加说明）
  subtitle_final <- subtitle
  if (show_cr && display_orr) {
    subtitle_final <- paste0(subtitle, "  |  ORR bar: dark = CR, light = PR")
  }
  p <- p + labs(title = title, subtitle = subtitle_final, x = x_label, y = y_label)

  # 13. 主题
  theme_args <- list(
    plot.title = element_text(hjust = 0.5, face = "bold", size = font_size + 2),
    plot.subtitle = element_text(hjust = 0.5, size = font_size - 1, color = "gray50"),
    plot.caption = element_text(hjust = 0.5, size = font_size - 2, color = "gray50", face = "italic"),
    axis.text.x = element_text(size = font_size, face = "bold", margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    axis.line.x = element_line(color = "black", size = 0.5)
  )
  if (show_legend) {
    theme_args$legend.position <- legend_position
    theme_args$legend.title <- element_text(face = "bold", size = font_size - 1)
    theme_args$legend.text <- element_text(size = font_size - 2)
  } else {
    theme_args$legend.position <- "none"
  }
  if (rotate_x_labels) {
    theme_args$axis.text.x <- element_text(
      size = font_size, face = "bold", angle = x_label_angle,
      hjust = 1, margin = margin(t = 5)
    )
  }

  p <- p + do.call(theme_minimal, list(base_size = font_size)) + do.call(theme, theme_args)

  # 14. P值表格
  if (show_p_value_table && !is.null(p_value_table_data)) {
    table_grob <- tableGrob(
      p_value_table_data, rows = NULL,
      theme = ttheme_minimal(
        core = list(fg_params = list(hjust = 0, x = 0.1), bg_params = list(fill = "white")),
        colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.1))
      )
    )
    if (p_value_table_position == "bottom") {
      p <- gridExtra::grid.arrange(p, table_grob, ncol = 1, heights = c(3, 1))
    } else if (p_value_table_position == "top") {
      p <- gridExtra::grid.arrange(table_grob, p, ncol = 1, heights = c(1, 3))
    } else if (p_value_table_position == "right") {
      p <- gridExtra::grid.arrange(p, table_grob, ncol = 2, widths = c(3, 1.5))
    }
  }

  return(list(
    plot = p, stats = final_stats, p_values = p_values_adjusted,
    comparisons = group_pairs, n_comparisons = n_comparisons,
    p_adjust_method = p_adjust_method, group_labels = group_labels,
    colors = colors, plot_data = plot_data, p_value_table = p_value_table_data
  ))
}


# 使用例子
if(F){
    
# 显示 ORR（含 CR 深色叠加）+ DCR
df <- data_plot_wf
plot_response_rates(df, group_col = "cohort", response_col = "response",
                    show_cr = TRUE, show_dcr = TRUE)

# 只显示 ORR（含 CR 叠加）
plot_response_rates(df, group_col = "cohort", response_col = "response",
                    show_cr = TRUE, show_dcr = FALSE)

# 原行为：ORR + DCR 单层柱
plot_response_rates(df, group_col = "cohort", response_col = "response",
                    show_cr = FALSE, show_dcr = TRUE)

result_custom <- plot_response_rates(
  df = data_plot_wf,
  group_col = "cohort",
  # 自定义组标签
  group_labels = c("Cohort1: SHR-A1811+Adebrelimab", 
                   "Cohort2: SHR-A1811+Adebrelimab", 
                   "Cohort3: SHR-A1811"),
  # 自定义颜色
  colors = c("#4B32C9", "#32B0C9", "#B0C932", "#C94B32"),
  title = "Best of response by cohort",
  subtitle = paste0("By LHF; ",date),
  response_col = "response",
  response_levels = c("CR", "PR", "SD", "PD", "NE"),  # 包含NE
  orr_levels = c("CR", "PR"),                          # ORR定义
  dcr_levels = c("CR", "PR", "SD"),                    # DCR定义
  p_adjust_method = "none",
  show_comparisons = F,
  show_ci_labels = F,
  show_cr = T,show_dcr = F
)

print(result_custom$plot)
    }
