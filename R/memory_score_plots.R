# Internal display helpers for the experimental-barcode memory-gene benchmark.
# Statistics and classification thresholds are supplied by CompareMemoryScores().
.iscoral_plot_memory_score_comparison <- function(
    plot_data,
    roc_curve,
    pr_curve,
    summary_table,
    plot_parameters
) {
  method_color <- "#C74B46"
  reference_color <- "#A0A0A0"
  class_colors <- c(
    TP = "#D37648", TN = "#E0E0E0", FP = "#91B35C", FN = "#4C7AA5"
  )
  class_shapes <- c(TP = 16, TN = 1, FP = 17, FN = 15)
  class_order <- c("TP", "TN", "FP", "FN")
  base_size <- plot_parameters$base_size
  base_family <- plot_parameters$base_family

  benchmark_theme <- ggplot2::theme_bw(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = base_size + 2, hjust = 0.5,
        margin = ggplot2::margin(b = 12)
      ),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "#333333", fill = NA, linewidth = 1
      ),
      panel.grid.major = ggplot2::element_line(
        color = "#E5E5E5", linewidth = 0.4, linetype = "dashed"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black", size = max(base_size - 2, 8)),
      axis.title = ggplot2::element_text(color = "black", size = base_size),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = max(base_size - 3, 8)),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )

  classified <- is.finite(summary_table$best_threshold)
  # Draw neutral genes first; retain all compared genes without subsampling.
  scatter_data <- plot_data[
    order(match(plot_data$classification, c("TN", "FN", "FP", "TP"))),
    , drop = FALSE
  ]
  scatter <- ggplot2::ggplot(
    scatter_data,
    ggplot2::aes(
      x = .data$ground_truth_memory_score,
      y = .data$inferred_memory_score
    )
  ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, color = reference_color,
      linetype = "dashed", linewidth = 0.7
    )

  if (length(unique(scatter_data$ground_truth_memory_score)) > 1L) {
    scatter <- scatter + ggplot2::geom_smooth(
      method = "lm", formula = y ~ x, se = FALSE,
      color = method_color, linewidth = 1.2
    )
  }
  if (classified) {
    scatter <- scatter +
      ggplot2::geom_vline(
        xintercept = plot_parameters$benchmark_truth_cutoff,
        color = "#333333", linetype = "dashed", linewidth = 0.7
      ) +
      ggplot2::geom_hline(
        yintercept = summary_table$best_threshold,
        color = "#333333", linetype = "dashed", linewidth = 0.7
      ) +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$classification, shape = .data$classification),
        size = 1.8, alpha = 0.8, stroke = 0.4
      ) +
      ggplot2::scale_color_manual(
        name = "Classification", values = class_colors,
        limits = class_order, drop = FALSE
      ) +
      ggplot2::scale_shape_manual(
        name = "Classification", values = class_shapes,
        limits = class_order, drop = FALSE
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))
      )
  } else {
    # A single-class comparison cannot support TP/TN/FP/FN assignments.
    scatter <- scatter + ggplot2::geom_point(
      color = "#8E8E8E", size = 1.8, alpha = 0.8
    )
  }

  if (classified && plot_parameters$label_top_n > 0L) {
    label_data <- plot_data[plot_data$classification == "TP", , drop = FALSE]
    label_data <- label_data[
      order(-label_data$inferred_memory_score, label_data$name), , drop = FALSE
    ]
    label_data <- utils::head(label_data, plot_parameters$label_top_n)
    if (nrow(label_data) > 0L) {
      scatter <- scatter + ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(label = .data$name),
        color = class_colors[["TP"]], family = base_family,
        size = base_size / 4, seed = 42,
        segment.color = "#333333", segment.linetype = "dashed",
        segment.size = 0.4, box.padding = 0.6, point.padding = 0.3,
        max.overlaps = Inf, show.legend = FALSE
      )
    }
  }
  scatter <- scatter +
    ggplot2::labs(
      title = "Ground truth vs. isCORAL",
      x = expression("Ground-truth memory score (" * omega^2 * ")"),
      y = expression("Inferred memory score (" * omega^2 * ")")
    ) +
    benchmark_theme

  # Do not render legacy placeholder coordinates as valid classification curves.
  if (!classified) {
    return(list(scatter = scatter, roc = NULL, pr = NULL))
  }

  roc <- ggplot2::ggplot(
    roc_curve,
    ggplot2::aes(x = .data$false_positive_rate, y = .data$true_positive_rate)
  ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, color = reference_color,
      linetype = "dashed", linewidth = 0.7
    ) +
    ggplot2::geom_path(
      ggplot2::aes(color = "isCORAL"), linewidth = 1.5
    ) +
    ggplot2::scale_color_manual(values = c(isCORAL = method_color)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = "Memory-gene ROC",
      x = "False positive rate", y = "True positive rate"
    ) +
    benchmark_theme

  best_f1_point <- pr_curve[which.max(pr_curve$f1), , drop = FALSE]
  pr <- ggplot2::ggplot(
    pr_curve,
    ggplot2::aes(x = .data$recall, y = .data$precision)
  ) +
    ggplot2::geom_hline(
      yintercept = summary_table$random_auprc, color = reference_color,
      linetype = "dashed", linewidth = 0.7
    ) +
    ggplot2::geom_path(
      ggplot2::aes(color = "isCORAL"), linewidth = 1.5
    ) +
    ggplot2::geom_point(
      data = best_f1_point, ggplot2::aes(shape = "Max F1"),
      fill = method_color, color = "#333333", size = 3.5, stroke = 0.6
    ) +
    ggplot2::scale_color_manual(values = c(isCORAL = method_color)) +
    ggplot2::scale_shape_manual(values = c("Max F1" = 23)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = "Precision-recall curve",
      x = "Recall", y = "Precision"
    ) +
    benchmark_theme

  list(scatter = scatter, roc = roc, pr = pr)
}
