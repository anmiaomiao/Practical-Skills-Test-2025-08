###门水平####
# 安装和加载必要的包
packages <- c("ggplot2", "tidyverse", "reshape2", "scales", "RColorBrewer")
lapply(packages, library, character.only = TRUE)

# ==== 设置参数 ====
input_file <- "bracken.P.txt"
design_file <- "metadata.txt"
group_column <- "Group"
top_n <- 10

# ==== 读取数据 ====
metadata <- read.table(design_file, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
taxonomy <- read.table(input_file, header = TRUE, row.names = 1, sep = "\t", comment.char = "", quote = "")

# ==== 计算全体物种归一化后的相对丰度 ====
taxonomy_rel <- sweep(taxonomy, 2, colSums(taxonomy), FUN = "/") * 100

# ==== 筛选前 N 丰度物种 ====
taxonomy$sum_abundance <- rowSums(taxonomy)
top_taxa <- rownames(taxonomy[order(taxonomy$sum_abundance, decreasing = TRUE)[1:top_n], ])

# ==== 合并其余为 Others ====
taxonomy_rel$Taxon <- rownames(taxonomy_rel)
taxonomy_long <- melt(taxonomy_rel, id.vars = "Taxon", variable.name = "Sample", value.name = "Abundance")
taxonomy_long$Taxon <- ifelse(taxonomy_long$Taxon %in% top_taxa, taxonomy_long$Taxon, "Others")

# ==== 添加分组信息 ====
metadata$Sample <- rownames(metadata)
taxonomy_merged <- merge(taxonomy_long, metadata, by.x = "Sample", by.y = "Sample")

# ==== 计算各组在 Phylum 水平的平均相对丰度 ====
group_avg <- taxonomy_merged %>%
  group_by(!!sym(group_column), Taxon) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")

# ==== 按组内相对丰度排序 Taxon ====
group_avg <- group_avg %>%
  group_by(!!sym(group_column)) %>%
  arrange(mean_abundance, .by_group = TRUE) %>%
  mutate(Taxon = factor(Taxon, levels = unique(Taxon)))

my_colors <- c("#BC80BD", "skyblue1", "#CCEBC5","#FCCDE5","#FDB462",
               "#104E8B","#BEBADA","#B2DFEE","#B3DE69","#DAD386",
               "#80B1D3")

# ==== 绘图 ====
p <- ggplot(group_avg, aes(x = !!sym(group_column), y = mean_abundance, fill = Taxon)) +
  geom_bar(stat = "identity",color = "black", linewidth = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  xlab(group_column) +
  ylab("Average Relative Abundance (%)") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = my_colors)

print(p)

# ==== 统计检验部分：Kruskal-Wallis检验 ====

# 计算每个门（Taxon）在两组样本中的相对丰度数据，用于检验
kruskal_results <- taxonomy_merged %>%
  group_by(Taxon) %>%
  summarise(p_value = kruskal.test(Abundance ~ !!sym(group_column))$p.value)

# 计算两组的平均相对丰度（方便导出）
avg_abundance_wide <- group_avg %>%
  pivot_wider(names_from = !!sym(group_column), values_from = mean_abundance)

# 合并平均丰度和检验结果
final_results <- merge(avg_abundance_wide, kruskal_results, by = "Taxon")

# 添加显著性星号
final_results$significance <- cut(final_results$p_value,
                                  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                  labels = c("***", "**", "*", "ns"))

# 重命名列名，方便理解（假设组名为 "HC" 和 "ID"）
colnames(final_results)[2:3] <- c("Mean_Abundance_HC", "Mean_Abundance_ID")

# 导出结果到csv
write.csv(final_results, file = "Phylum_abundance_kruskal_results.csv", row.names = FALSE)

# 打印前几行结果查看
print(head(final_results))


###种水平####
# ==== 加载必要包 ====
packages <- c("ggplot2", "tidyverse", "reshape2", "scales", "RColorBrewer")
lapply(packages, library, character.only = TRUE)

# ==== 参数设置 ====
input_file <- "bracken.S.txt"       # 种水平丰度表
design_file <- "metadata.txt"       # 分组文件
group_column <- "Group"             # metadata 中分组列名
top_n <- 15                        # 前 N 个物种展示

# ==== 读取数据 ====
metadata <- read.table(design_file, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
taxonomy <- read.table(input_file, header = TRUE, row.names = 1, sep = "\t", comment.char = "", quote = "")

# ==== 计算每个样本的相对丰度（列归一化）====
taxonomy_rel <- sweep(taxonomy, 2, colSums(taxonomy), FUN = "/") * 100

# ==== 筛选全局丰度前 top_n 的物种 ====
taxonomy$sum_abundance <- rowSums(taxonomy)
top_taxa <- rownames(taxonomy[order(taxonomy$sum_abundance, decreasing = TRUE)[1:top_n], ])

# ==== 转为长格式 + 合并为 "Others" ====
taxonomy_rel$Taxon <- rownames(taxonomy_rel)
taxonomy_long <- melt(taxonomy_rel, id.vars = "Taxon", variable.name = "Sample", value.name = "Abundance")
taxonomy_long$Taxon <- ifelse(taxonomy_long$Taxon %in% top_taxa, taxonomy_long$Taxon, "Others")

# ==== 合并 others 的 abundance ====
taxonomy_grouped <- taxonomy_long %>%
  group_by(Sample, Taxon) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# ==== 添加分组信息 ====
metadata$Sample <- rownames(metadata)
taxonomy_merged <- merge(taxonomy_grouped, metadata, by.x = "Sample", by.y = "Sample")

# ==== 计算每组的平均相对丰度 ====
group_avg <- taxonomy_merged %>%
  group_by(!!sym(group_column), Taxon) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")

# ==== 全局排序 Taxon（用于控制颜色顺序） ====
taxon_levels <- group_avg %>%
  group_by(Taxon) %>%
  summarise(overall_mean = mean(mean_abundance)) %>%
  arrange(overall_mean) %>%
  pull(Taxon)

group_avg$Taxon <- factor(group_avg$Taxon, levels = taxon_levels)

# ==== 自定义颜色（需 ≥ Taxon 类别数） ====
my_colors <- unique(c( "#be1420","#FB8072", "#E9D7CB","#FCCDE5","#DAA5BB",
                       "#BEBADA","#B893BE", "#FFFFB3","#FDB462","#B3DE69",
                       "#CCEBC5","#7ABB8B","#B2DFEE",
                       "#80B1D3","#748ebb","#104E8B"))

# ==== 绘图 ====
p <- ggplot(group_avg, aes(x = !!sym(group_column), y = mean_abundance, fill = Taxon)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +   # 添加边框
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  xlab(group_column) +
  ylab("Average Relative Abundance (%)") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = my_colors)

print(p)

# ==== 计算两组每个物种的平均相对丰度（方便导出） ====
avg_abundance_wide <- group_avg %>%
  pivot_wider(names_from = !!sym(group_column), values_from = mean_abundance)

# ==== Kruskal-Wallis 检验（两组物种丰度） ====
kruskal_results <- taxonomy_merged %>%
  group_by(Taxon) %>%
  summarise(p_value = kruskal.test(Abundance ~ !!sym(group_column))$p.value)

# ==== 合并平均丰度和检验结果 ====
final_results <- merge(avg_abundance_wide, kruskal_results, by = "Taxon")

# ==== 添加显著性星号 ====
final_results$significance <- cut(final_results$p_value,
                                  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                  labels = c("***", "**", "*", "ns"))

# ==== 只保留需要的列，并改列名方便理解 ====
colnames(final_results)[2:(1 + length(unique(metadata[[group_column]])))] <- paste0("Mean_Abundance_", colnames(final_results)[2:(1 + length(unique(metadata[[group_column]])))])
final_results <- final_results %>%
  select(Taxon, starts_with("Mean_Abundance"), p_value, significance)

# ==== 导出结果到csv ====
write.csv(final_results, file = "species_abundance_kruskal_results.csv", row.names = FALSE)

# ==== 打印前几行结果查看 ====
print(head(final_results))






