#!/usr/bin/env Rscript

# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"

# 设置用户自己的R库路径，避免写入系统目录
.libPaths("/data/gaoyunyun/Rlibs")

# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}

# 命令行参数设置
option_list <- list(
  make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
              help="Input reads count file; such as OTU table, kraken2 taxonomy counts table [default %default]"),
  make_option(c("-d", "--depth"), type="numeric", default=0,
              help="Rarefaction depth; default 0 to auto using minimum [default %default]"),
  make_option(c("-s", "--seed"), type="numeric", default=1,
              help="Random sample seed, can set any integer [default %default]"),
  make_option(c("-n", "--normalize"), type="character", default="result/otutab_rare.txt",
              help="Normalized filename [default %default]"),
  make_option(c("-o", "--output"), type="character", default="result/alpha/vegan.txt",
              help="Output alpha diversity filename [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

print(paste("The input feature table is ", opts$input,  sep = ""))

# 加载依赖包
package_list <- c("vegan")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 创建输出目录，避免写入错误
suppressWarnings(dir.create(dirname(opts$normalize), recursive = TRUE))
suppressWarnings(dir.create(dirname(opts$output), recursive = TRUE))

# 读取物种表
species = read.table(opts$input, header=T, sep="\t", quote = "", row.names=1, comment.char="") 

# 计算样本最小数目作为 rarefaction depth
print(paste0("Samples size are:"))
colSums(species)
min = min(colSums(species))
if (opts$depth==0){
  opts$depth=min
}
print(paste("Rarefaction depth ", opts$depth, ". If depth set 0 will using sample minimum size ", min, sep = ""))

print(paste("Random sample number: ", opts$seed,  sep = ""))
set.seed(opts$seed)
otu = vegan::rrarefy(t(species), opts$depth)

idx = rowSums(otu) >= opts$depth
print(paste0("The discard samples: ",rownames(otu[!idx,])))
suppressWarnings(write.table(rownames(otu[!idx,]), file=paste(opts$normalize,".discard",sep=""), append = F, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F))

otu = otu[idx,]

# 计算 alpha 多样性指数
estimateR = t(estimateR(otu))[,c(1,2,4)]
colnames(estimateR) = c("richness", "chao1", "ACE")
shannon = diversity(otu, index = "shannon")
simpson = diversity(otu, index = "simpson")
invsimpson = diversity(otu, index = "invsimpson")
goods_coverage = apply(otu, 1, function(x) 1 - sum(x == 1)/sum(x))
pielou_e = shannon / log(rowSums(otu > 0))

alpha_div = cbind(estimateR, shannon, simpson, invsimpson, goods_coverage, pielou_e)
colnames(alpha_div) = c("richness", "chao1", "ACE", "shannon", "simpson", "invsimpson", "goods_coverage", "pielou_e")

print(paste0("Calculate 8 alpha diversities: richness, chao1, ACE, shannon, simpson, invsimpson, goods_coverage, pielou_e"))
head(alpha_div, n=1)

# 输出 normalized table
write.table("#OTUID\t", file=paste(opts$normalize,sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(t(otu), file=paste(opts$normalize,sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

# 输出 alpha diversity table
write.table("SampleID\t", file=paste(opts$output,sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(alpha_div, file=paste(opts$output,sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

print(paste("Name of rarefaction file ", opts$normalize,  sep = ""))
print(paste("Output alpha diversity filename ", opts$output, sep = ""))
print(paste("The discard samples file ",  opts$normalize, ".discard", sep = ""))
