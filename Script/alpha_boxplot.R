
# 加载必要的包
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

# 读取数据
alpha <- read.delim("bracken.S.alpha", header = TRUE, row.names = 1, quote = "", check.names = FALSE)
meta <- read.delim("metadata.txt", header = TRUE, stringsAsFactors = FALSE)

# 合并数据
alpha$SampleID <- rownames(alpha)
df <- merge(alpha, meta, by = "SampleID")

# 转换为因子
df$Group <- factor(df$Group, levels = c("HC", "ID"))

# 多样性指标
alpha_indices <- c("richness", "chao1", "ACE", "shannon", "simpson", "invsimpson", "goods_coverage", "pielou_e")

# 转换为长格式
df_long <- melt(df, id.vars = c("SampleID", "Group"), measure.vars = alpha_indices,
                variable.name = "Index", value.name = "Value")

# 绘图
p <- ggplot(df_long, aes(x = Group, y = Value, color = Group)) +
  geom_boxplot(fill = NA, linewidth = 0.7, outlier.shape = NA) +
  geom_jitter(
    aes(fill = Group),
    shape = 21,
    color = "black",     # 散点边框为黑色
    width = 0.4,
    size = 2,
    alpha = 0.9,
    stroke = 0.3
  ) +
  facet_wrap(~ Index, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y.npc = "top") +
  labs(x = "Group", y = "Diversity Index Value") +
  scale_color_manual(values = c("HC" = "#990005", "ID" = "#285F7D")) +
  scale_fill_manual(values = c("HC" = "#990005", "ID" = "#285F7D")) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = "none"
  )

# 显示图形
print(p)

