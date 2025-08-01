#!/bin/bash

# 设置包含 SRA 访问号的文件路径
input_file="SRR_Acc_List.txt"

# 设置输出目录
output_dir="/data/gaoyunyun/project/anmiaomiao/download/SRA_seq"
sra_cache_dir="/data/gaoyunyun/project/anmiaomiao/download/SRA_cache"

# 创建输出目录（如果不存在）
mkdir -p "$output_dir"
mkdir -p "$sra_cache_dir"

# 设置 SRA 工具的缓存目录
export NCBI_SETTINGS="$sra_cache_dir"

# 检查必要的工具是否存在
command -v prefetch >/dev/null 2>&1 || { echo "错误: prefetch 未安装或不在 PATH 中" >&2; exit 1; }
command -v fastq-dump >/dev/null 2>&1 || { echo "错误: fastq-dump 未安装或不在 PATH 中" >&2; exit 1; }

# 检查输入文件是否存在
if [ ! -f "$input_file" ]; then
    echo "错误: 输入文件 $input_file 不存在"
    exit 1
fi

echo "开始下载和转换 SRA 数据..."
echo "输入文件: $input_file"
echo "输出目录: $output_dir"
echo "SRA 缓存目录: $sra_cache_dir"
echo "----------------------------------------"

# 统计总数
total_count=$(wc -l < "$input_file")
current_count=0

# 循环读取每个 SRA 访问号
while IFS= read -r sra_id || [ -n "$sra_id" ]
do
    # 跳过空行和注释行
    if [[ -z "$sra_id" || "$sra_id" =~ ^#.* ]]; then
        continue
    fi
    
    current_count=$((current_count + 1))
    echo "[$current_count/$total_count] 处理 $sra_id..."
    
    # 步骤1: 使用 prefetch 下载 SRA 文件
    echo "  正在下载 SRA 文件..."
    if prefetch --output-directory "$sra_cache_dir" "$sra_id"; then
        echo "  SRA 文件下载成功"
    else
        echo "  警告: SRA 文件下载失败，跳过 $sra_id"
        continue
    fi
    
    # 步骤2: 使用 fastq-dump 转换为 FASTQ 格式
    echo "  正在转换为 FASTQ 格式..."
    sra_file="$sra_cache_dir/$sra_id/$sra_id.sra"
    
    if [ -f "$sra_file" ]; then
        if fastq-dump --gzip --split-files --outdir "$output_dir" "$sra_file"; then
            echo "  FASTQ 转换成功"
            
            # 可选: 删除 SRA 文件以节省空间（取消注释下面的行来启用）
            # rm -rf "$sra_cache_dir/$sra_id"
            # echo "  已删除临时 SRA 文件"
        else
            echo "  错误: FASTQ 转换失败 $sra_id"
        fi
    else
        echo "  错误: 找不到 SRA 文件 $sra_file"
    fi
    
    echo "  完成 $sra_id"
    echo "----------------------------------------"
done < "$input_file"

echo "所有任务完成！"
echo "FASTQ 文件保存在: $output_dir"
echo "SRA 缓存文件保存在: $sra_cache_dir"
echo "如果不再需要 SRA 文件，可以删除缓存目录以节省空间"

# 显示输出目录内容
echo ""
echo "输出文件列表:"
ls -lh "$output_dir"
