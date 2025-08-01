
# 使用 BiocManager 安装 file2meco
if (!requireNamespace("file2meco", quietly = TRUE)) {
  BiocManager::install("file2meco")
}


library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(file2meco)

# 读入MetaCyc pathway mapping
data(MetaCyc_pathway_map)
head(MetaCyc_pathway_map)

# 读入数据
df <- read.delim("associate.txt", sep = "\t", check.names = FALSE)

# 查看列名，确认是否包含括号等特殊符号
colnames(df)

# 为了避免特殊字符问题，先重命名列
colnames(df)[2] <- "Level_means"
colnames(df)[1] <- "Feature"

# 查看列名，确认是否包含括号等特殊符号
colnames(df)

# 分离 HC 和 ID 丰度，并计算 log2FC
df <- df %>%
  separate(Level_means, into = c("HC_abundance", "ID_abundance"), sep = "\\|") %>%
  mutate(
    HC_abundance = as.numeric(sub("HC:", "", HC_abundance)),
    ID_abundance = as.numeric(sub("ID:", "", ID_abundance)),
    log2FC = log2((ID_abundance + 1e-9) / (HC_abundance + 1e-9)),
    abs_log2FC = abs(log2FC)
  )

# 查看前几行
head(df)

# 提取Feature冒号后面部分命名Pathway
df <- df %>%
  mutate(
    Pathway = str_trim(str_extract(Feature, "(?<=: ).*"))
  )

# 提取Feature冒号前面的部分命名Pathway_ID
df <- df %>%
  mutate(
    Pathway_ID = str_trim(str_extract(Feature, "^[^:]+"))
  )

# 将associate.txt和MetaCyc_pathway_map进行合并，给Feature添加分类信息
df <- df %>%
  left_join(MetaCyc_pathway_map %>%
              select(pathway, Superclass1, Superclass2),
            by = c("Pathway" = "pathway"))

# 筛选显著差异通路（Q值<0.05）, |log2FC| 排名前20的行
sig_df <- df %>% 
  filter(`Q-value` < 0.05) %>%
  arrange(desc(abs_log2FC)) %>%    # 按绝对log2FC降序排列
  slice_head(n = 20)               # 取前20条记录

# 画气泡图，只显示排名前20的路径
ggplot(sig_df, aes(x = log2FC, y = reorder(Feature, log2FC))) +
  geom_point(aes(size = abs_log2FC, fill = `Q-value`), 
             color = "black",    
             shape = 21,         
             stroke = 0.5,       
             alpha = 0.8) +
  scale_fill_gradientn(
    name = "Q-value",
    colors = c("#B73D30","#F08C57","#F4F4E4","#acd2e5","#6090c1","#0D5875"),   # 三段颜色，可以改成你喜欢的
    limits = c(0, 0.05)                    # 限制范围
  )+
  scale_size_continuous(name = "|log2FC|", range = c(1.5, 4.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Top 20 Differentially Abundant MetaCyc Pathways by |log2FC|",
    x = "log2 Fold Change (ID / HC)",
    y = "MetaCyc Pathway"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # 用linewidth
  )

# 查看前几行数据
head(df)
head(sig_df)

# 保存完整 df
write.csv(df, "df_all_results.csv", row.names = FALSE)

# 保存显著差异 sig_df
write.csv(sig_df, "sig_df_top20.csv", row.names = FALSE)


