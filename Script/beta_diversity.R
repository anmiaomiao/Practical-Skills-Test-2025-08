
###计算beta多样性
# 加载所需包
library(vegan)

# 读取 OTU 表格
otu_table <- read.delim("bracken.S.norm", header = TRUE, row.names = 1, quote = "", check.names = FALSE)

# 数据预处理：对OTU进行log10转化（加1防止log(0)），并过滤掉稀疏物种（存在于样本数低于10%的物种）
otu_table_log <- log10(otu_table + 1)
otu_table_filtered <- otu_table_log[, colSums(otu_table_log > 0) > (0.1 * nrow(otu_table_log))]

# 转置数据矩阵：行为样本，列为物种
otu_transposed <- t(otu_table_filtered)

# 计算 Bray-Curtis 距离
bray_dist <- vegdist(otu_transposed, method = "bray")
# 保存 Bray-Curtis 距离矩阵到文件
write.table(as.matrix(bray_dist), file = "bray_dist.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# 计算 Jaccard 距离（二值化）
jaccard_dist <- vegdist(otu_transposed, method = "jaccard", binary = TRUE)
# 保存 Jaccard 距离矩阵到文件
write.table(as.matrix(jaccard_dist), file = "jaccard_dist.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


###PCoA可视化####

# 加载必要的包
library(ggplot2)
library(rlang)
library(vegan)

# 1. 读取距离矩阵
distance_matrix <- as.matrix(read.delim("jaccard_dist.txt", row.names = 1, sep = "\t", check.names = FALSE))

# 2. 处理NA
distance_matrix[is.na(distance_matrix)] <- 0

# 3. 读取分组信息（metadata）
group <- read.delim("metadata.txt", row.names = 1, sep = "\t", check.names = FALSE)
group_col <- "Group"  # 元数据中表示分组的列名

# 4. 匹配样本名（确保分组和距离矩阵一致）
common_samples <- intersect(rownames(distance_matrix), rownames(group))
distance_matrix <- distance_matrix[common_samples, common_samples]
group <- group[common_samples, , drop = FALSE]

# 5. 计算 PCoA
pcoa <- cmdscale(as.dist(distance_matrix), k = (nrow(distance_matrix) - 1), eig = TRUE)
plot_data <- data.frame(pcoa$points)[, 1:2]
colnames(plot_data) <- c("PCoA1", "PCoA2")
plot_data$Sample <- rownames(plot_data)

# 6. 计算解释的变异百分比
eig <- pcoa$eig
PCOA1 <- format(100 * eig[1] / sum(eig), digits = 4)
PCOA2 <- format(100 * eig[2] / sum(eig), digits = 4)

# 7. 合并数据
plot_data <- merge(plot_data, group, by.x = "Sample", by.y = "row.names", all = TRUE)
plot_data <- plot_data[!is.na(plot_data[[group_col]]) & plot_data[[group_col]] != "", ]

# 8. PERMANOVA 检验
permanova_result <- adonis2(as.dist(distance_matrix) ~ Group, data = group, permutations = 999, method = "bray")
print(permanova_result)

# 提取 p 值和 R2
pval <- format(permanova_result$`Pr(>F)`[1], digits = 4)
R2 <- format(permanova_result$R2[1], digits = 3)
Fval <- format(permanova_result$F[1], digits = 4)  # 提取 F 值

# 9. 设置颜色
default_colors <- c("ID" = "#415686", "HC" = "#990005")      # 散点颜色
ellipse_fills <- c("ID" = "#9FACC0", "HC" = "#EDD0C6")       # 椭圆底色

# 10. 绘图
p <- ggplot(plot_data, aes(x = PCoA1, y = PCoA2)) +
  geom_point(
    aes(fill = !!sym(group_col)),
    shape = 21,
    size = 3.5,
    color = "black",
    stroke = 0.7
  ) +
  stat_ellipse(
    aes(color = !!sym(group_col), fill = !!sym(group_col)),
    geom = "polygon",
    linetype = "solid",
    level = 0.98,
    linewidth = 0.85,
    alpha = 0.05
  ) +
  scale_fill_manual(values = default_colors) +
  scale_color_manual(values = default_colors) +
  labs(
    x = paste0("PCoA1 (", PCOA1, "%)"),
    y = paste0("PCoA2 (", PCOA2, "%)"),
    title = paste0("PCoA (PERMANOVA: R² = ", R2,
                   ", F = ", Fval,
                   ", p = ", pval, ")")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.background = element_rect(fill = "transparent"),
    legend.position = "right"
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5)

# 11. 显示图形
print(p)
