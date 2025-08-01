
library(ggplot2)
library(dplyr)

# 读数据
deg_data <- read.table("input.res.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 过滤无效数据
deg_data <- deg_data[deg_data$KW_Pvalue != "-", ]
deg_data$KW_Pvalue <- as.numeric(deg_data$KW_Pvalue)
deg_data$LDA <- as.numeric(deg_data$LDA)

# 赋正负LDA，构造绘图数据,LDA阈值>=2,p值<0.05
plot_data <- deg_data %>%
  mutate(
    LDA_signed = case_when(
      Enriched_Groups == "ID" ~ LDA,
      Enriched_Groups == "HC" ~ -LDA,
      TRUE ~ 0
    ),
    negLog10P = -log10(KW_Pvalue),
    Regulation = case_when(
      LDA_signed >= 2 & KW_Pvalue <= 0.05 ~ "Up",
      LDA_signed <= -2 & KW_Pvalue <= 0.05 ~ "Down",
      TRUE ~ "Normal"
    )
  )

# 绘图
p <- ggplot(plot_data, aes(x = LDA_signed, y = negLog10P, color = Regulation)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = c("Down" = "cyan4", "Normal" = "grey50", "Up" = "orangered")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = -4:4) +  # 设置横坐标刻度为 -4 到 4
  labs(
    x = "LDA",
    y = "-Log10 p-value",
    color = "Regulation"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank()
  )+
  # 标记部分显著物种，阈值可调整
  geom_text(
    data = plot_data %>% filter(negLog10P > 1.3 & (LDA_signed > 2 | LDA_signed < -2)),
    aes(label = Biomarker_Names),
    size = 2,
    vjust = -0.5,
    color = "black"
  )

print(p)

