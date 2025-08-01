import os
import json

# 目录和阈值设置
qc_dir = "temp/qc"
output_path = os.path.join(qc_dir, "qc_summary.tsv")

# 质控要求
q20_threshold = 0.95
q30_threshold = 0.85
min_reads = 10_000_000
min_retention = 0.80
max_n_content = 0.02      # N碱基占比 ≤ 2%
min_read_length = 50      # 最短read长度 ≥ 50
min_adapter_trimmed = 1   # 去接头读数 > 0

results = []

for file in sorted(os.listdir(qc_dir)):
    if file.endswith("_fastp.json"):
        sample = file.replace("_fastp.json", "")
        file_path = os.path.join(qc_dir, file)

        try:
            with open(file_path) as f:
                data = json.load(f)

            # 获取过滤前后数据
            before = data['summary']['before_filtering']
            after = data['summary']['after_filtering']

            raw_reads = before['total_reads']
            clean_reads = after['total_reads']
            retention = clean_reads / raw_reads if raw_reads > 0 else 0
            q20 = after['q20_rate']
            q30 = after['q30_rate']

            # 从read1和read2后的content_curves中取N含量曲线
            n_curve_read1 = data.get('read1_after_filtering', {}).get('content_curves', {}).get('N', [])
            n_curve_read2 = data.get('read2_after_filtering', {}).get('content_curves', {}).get('N', [])

            max_n_read1 = max(n_curve_read1) if n_curve_read1 else 0.0
            max_n_read2 = max(n_curve_read2) if n_curve_read2 else 0.0

            # 取最大值作为N含量
            n_content = max(max_n_read1, max_n_read2)

            # 获取平均 read 长度
            read_length = after.get('read1_mean_length', 0)

            # 获取去接头数
            adapter_trimmed = data['adapter_cutting']['adapter_trimmed_reads']

            # 判断是否合格
            passed = (
                q20 >= q20_threshold and
                q30 >= q30_threshold and
                clean_reads >= min_reads and
                retention >= min_retention and
                n_content <= max_n_content and
                read_length >= min_read_length and
                adapter_trimmed >= min_adapter_trimmed
            )

            results.append([
                sample,
                q20,
                q30,
                clean_reads,
                retention,
                n_content,    # 保持原始浮点数，不做四舍五入
                int(read_length),
                adapter_trimmed,
                "PASS" if passed else "FAIL"
            ])

        except Exception as e:
            print(f"[ERROR] Failed to process {file}: {e}")

# 输出结果
with open(output_path, "w") as out:
    header = "Sample\tQ20\tQ30\tCleanReads\tRetention\tN_Content\tAvgReadLen\tAdapterTrimmed\tPass"
    out.write(header + "\n")
    print(header)

    for row in results:
        line = "\t".join(map(str, row))
        out.write(line + "\n")
        print(line)

print(f"\n✅ QC summary saved to: {output_path}")