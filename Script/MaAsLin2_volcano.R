
# 简介 Introduction #----

# 程序功能：运用MaAsLin2方法计算组间差异
# Functions: Difference analysis using MaAsLin2


# 设置工作目录
setwd("C:/Users/11251/Desktop/RAI/202508")
opts <- list(
  input = "MaAsLin2_overall_difference.csv",
  trend = "MaAsLin2_enriched_depleted.csv",
  output = "C:/Users/11251/Desktop/RAI/202508/"
)

print(opts)


# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","magrittr","dplyr")) 
}
# load related packages
suppressWarnings(suppressMessages(library("ggplot2")))
suppressWarnings(suppressMessages(library("dplyr")))
suppressWarnings(suppressMessages(library("magrittr")))

data_species01 <- read.csv(opts$trend, row.names = 1)
data_species01$species <- rownames(data_species01)
data_species01 <- data_species01[, -c(1:4)]

data_MWAS <- read.csv(opts$input, row.names = 1)
data_MWAS$species <- rownames(data_MWAS)

data_species02 <- merge(data_species01, data_MWAS, by = "species")

#devtools::install_github("BioSenior/ggvolcano", force = TRUE)
library(ggvolcano)
data_vol <- data_species02
data_vol = as.data.frame(data_vol)

data_vol2 <- data_vol
data_vol2$padj2 <- -log10(data_vol2$FDR)

#P.Value = 0.05
library(ggplot2)

# 设置颜色分组：FDR >= 0.05 或 |Beta| <= 0.3 都标记为 NotSig
data_vol2$color_group <- ifelse(abs(data_vol2$Beta) <= 0.3 | data_vol2$FDR >= 0.05, "NotSig", data_vol2$NPC.association)


# 重新绘图（增加 x = -0.3 和 x = 0.3 的虚线）
p_volcano1 <- ggplot(data = data_vol2, aes(x = -Beta, y = padj2)) +
  geom_point(alpha = 0.7, size = 1.0, aes(color = color_group)) + 
  ylab("-log10(Pvalue)") +
  scale_color_manual(values = c("NotSig"="grey50",
                                "enriched"="orangered",
                                "depleted"="cyan4")) + 
  geom_vline(xintercept = 0, lty = 4, col = "black", lwd = 0.4) + 
  # 新增 Beta 阈值虚线
  geom_vline(xintercept = -0.3, lty = 2, col = "darkgrey", lwd = 0.5) +
  geom_vline(xintercept =  0.3, lty = 2, col = "darkgrey", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.4) + 
  labs(x = "Coef. (by MaAsLin2)", y = bquote(atop(-Log[10]~italic(FDR)))) +
  theme_bw()


# add labels
library(dplyr)

# using geom_text_repel() to add labels
library(ggrepel)

# 选取 |Beta| > 0.3 且 FDR < 0.05 的点
beta_label_data <- filter(data_vol2, abs(Beta) > 0.3 & FDR < 0.05)

# 新增：对 |Beta| > 0.3 且 FDR < 0.05 的点加标注
p_volcano2 <- p_volcano1 +  
  geom_text_repel(
    data = beta_label_data, 
    aes(x = -Beta, y = padj2, label = Feature),
    size = 2,
    color = "black",
    max.overlaps = 100
  ) +
  theme(
    legend.position.inside = c(0.84, 0.85),
    panel.grid = element_blank()
  )

p_volcano2

