#!/bin/bash
source config_1.ini
#echo "[DEBUG] 工作目录结构:"
#tree -L 3 "$WORKDIR" || ls -lR "$WORKDIR" | grep -E "alignment/|bam$"

mkdir -p "$WORKDIR"/{clean_data,alignment/{raw,unique},read_position,specific_tags,vcf/{raw,filtered},snp,logs}

function run_fastp {
    echo "[$(date)] 数据质控..."
    for sample in "$MALE_SAMPLE" "$FEMALE_SAMPLE"; do
        fastp -i "$RAW_DATA_DIR/${sample}_R1.fq.gz" \
              -I "$RAW_DATA_DIR/${sample}_R2.fq.gz" \
              -o "$WORKDIR/clean_data/${sample}_clean_R1.fq.gz" \
              -O "$WORKDIR/clean_data/${sample}_clean_R2.fq.gz" \
              --detect_adapter_for_pe \
              --cut_front --cut_tail \
              --average_qual 20 \
              --length_required "$MIN_LENGTH" \
              --thread "$THREADS" \
              --html "$WORKDIR/clean_data/${sample}_fastp_report.html" \
              --json "$WORKDIR/clean_data/${sample}_fastp_report.json" \
              2> "$WORKDIR/logs/${sample}_fastp.log"
        # 验证输出
        echo "${sample} reads数:"
        gzcat "$WORKDIR/clean_data/${sample}_clean_R1.fq.gz" | echo $((`wc -l`/4))
    done
}

function run_alignment {
    echo "[$(date)] 比对参考基因组..."
    for sample in "$MALE_SAMPLE" "$FEMALE_SAMPLE"; do
        bam_file="$WORKDIR/alignment/raw/${sample}.bam"
        # 原始比对
        bwa mem -t "$THREADS" -v 3 \
            -T 30 -B 4 -O 6,1 \
            -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
            "$REF_GENOME" \
            "$WORKDIR/clean_data/${sample}_clean_R1.fq.gz" \
            "$WORKDIR/clean_data/${sample}_clean_R2.fq.gz" \
            2> "$WORKDIR/logs/${sample}_bwa.log" \
            | samtools sort -@ "$THREADS" -O BAM \
            -o "$WORKDIR/alignment/raw/${sample}.bam"
        # 验证比对结果
        echo "=== ${sample}比对统计 ==="
        samtools flagstat "$WORKDIR/alignment/raw/${sample}.bam"
        # 生成原始BAM的索引
        samtools index -c "$bam_file" && \
        chmod 644 "$bam_file" "$bam_file.csi" || {
            echo "错误: 无法为原始BAM生成索引"
            exit 1
        }
        echo "染色体比对分布:"
        samtools idxstats "$WORKDIR/alignment/raw/${sample}.bam" | head
    done
}

function bam_filter {
    echo "[$(date)] 过滤唯一比对..."
    for sample in "$MALE_SAMPLE" "$FEMALE_SAMPLE"; do
        # 定义清晰的文件路径变量
        input_bam="$WORKDIR/alignment/raw/${sample}.bam"
        output_bam="$WORKDIR/alignment/unique/${sample}.unique.bam"
        
        # 执行过滤
        bamtools filter \
            -in "$input_bam" \
            -out "$output_bam" \
            -mapQuality ">=${MAPQ_THRESHOLD}" \
            -tag "NM:<=${MAX_MISMATCH}" || {
                echo "错误: bamtools过滤失败"
                exit 1
            }
        
        # 为过滤后的BAM生成索引
        samtools index -c "$output_bam" && \
        chmod 644 "$output_bam" "$output_bam.csi" || {
            echo "错误: 无法生成/设置索引权限"
            exit 1
        }
        
        # 验证文件
        samtools quickcheck "$output_bam" || {
            echo "错误: 输出的BAM文件无效"
            exit 1
        }
    done
    
    echo "[$(date)] 过滤完成，已生成.csi索引"
}

function extract_read_positions {
    echo "[$(date)] 提取reads位置信息..."
    for sample in "$MALE_SAMPLE" "$FEMALE_SAMPLE"; do
        samtools view "$WORKDIR/alignment/unique/${sample}.unique.bam" \
            | awk 'BEGIN{OFS="\t"} {
                print $1, $3, $4, $4+length($10)-1
            }' > "$WORKDIR/read_position/${sample}.position.tsv"
    done
}

# 添加Python脚本内容作为函数
function run_tag_analysis {
    echo "[$(date)] 分析tag位置信息..."
    
    # 定义路径（全部在read_position目录）
    INPUT_DIR="$WORKDIR/read_position"
    CLUSTER_FILE="$INPUT_DIR/tag_clusters.tsv"
    MALE_TAG_FILE="$INPUT_DIR/male_specific_tags.tsv"
    FEMALE_TAG_FILE="$INPUT_DIR/female_specific_tags.tsv"
    SHARED_TAG_FILE="$INPUT_DIR/shared_tags.tsv"
    
    # 嵌入式Python聚类分析
    python3 << EOF | tee "$CLUSTER_FILE"
import csv
from collections import defaultdict

def load_data(male_file, female_file):
    data = defaultdict(list)
    for sample, file in [('male', male_file), ('female', female_file)]:
        with open(file) as f:
            for line in csv.reader(f, delimiter='\t'):
                if len(line) == 4:
                    chrom, start, end = line[1], int(line[2]), int(line[3])
                    data[chrom].append((start, end, sample))
    return data

def cluster_tags(data, max_distance):
    results = []
    for chrom in data:
        tags = sorted(data[chrom])
        current_cluster = []
        
        for tag in tags:
            if not current_cluster:
                current_cluster.append(tag)
            elif tag[0] - current_cluster[-1][1] <= max_distance:
                current_cluster.append(tag)
            else:
                results.append((chrom, process_cluster(current_cluster)))
                current_cluster = [tag]
        
        if current_cluster:
            results.append((chrom, process_cluster(current_cluster)))
    
    return results

def process_cluster(cluster):
    counts = defaultdict(int)
    for _, _, sample in cluster:
        counts[sample] += 1
    
    category = "shared" if len(counts) == 2 else f"{list(counts.keys())[0]}_specific"
    
    return {
        'start': min(t[0] for t in cluster),
        'end': max(t[1] for t in cluster),
        'male_depth': counts.get('male', 0),
        'female_depth': counts.get('female', 0),
        'category': category
    }

# 主程序
male_file = "$INPUT_DIR/${MALE_SAMPLE}.position.tsv"
female_file = "$INPUT_DIR/${FEMALE_SAMPLE}.position.tsv"
max_distance = $CLUSTER_DISTANCE

print("Chromosome\tStart\tEnd\tMale_Depth\tFemale_Depth\tCategory")
for chrom, cluster in cluster_tags(load_data(male_file, female_file), max_distance):
    print(f"{chrom}\t{cluster['start']}\t{cluster['end']}\t"
          f"{cluster['male_depth']}\t{cluster['female_depth']}\t"
          f"{cluster['category']}")
EOF

    # 分类输出结果文件
    echo "分类输出tag结果..."
    awk -F'\t' -v male="$MALE_TAG_FILE" -v female="$FEMALE_TAG_FILE" -v shared="$SHARED_TAG_FILE" '
    BEGIN {
        # 写入表头
        print "Chromosome\tStart\tEnd\tMale_Depth\tFemale_Depth\tCategory" > male
        print "Chromosome\tStart\tEnd\tMale_Depth\tFemale_Depth\tCategory" > female
        print "Chromosome\tStart\tEnd\tMale_Depth\tFemale_Depth\tCategory" > shared
    }
    NR > 1 {  # 跳过Python输出的表头
        if ($6 == "male_specific") print >> male
        else if ($6 == "female_specific") print >> female
        else if ($6 == "shared") print >> shared
    }
    ' "$CLUSTER_FILE"
    
    # 统计结果
    echo "Tag分析结果统计:"
    echo " - 总cluster数: $(tail -n +2 "$CLUSTER_FILE" | wc -l)"
    echo " - 雄性特异tag数: $(tail -n +2 "$MALE_TAG_FILE" | wc -l)"
    echo " - 雌性特异tag数: $(tail -n +2 "$FEMALE_TAG_FILE" | wc -l)"
    echo " - 共享tag数: $(tail -n +2 "$SHARED_TAG_FILE" | wc -l)"
    
    # 验证输出
    echo "=== 分类文件前2行 ==="
    for file in "$MALE_TAG_FILE" "$FEMALE_TAG_FILE" "$SHARED_TAG_FILE"; do
        echo "[$(basename "$file")]"
        head -n 2 "$file"
    done
}

# 筛选包含转座子末端序列的雌雄特有tag
function filter_specific_tags {
    echo "[$(date)] 高效筛选特异序列..."
    
    source config_1.ini
    
    # 分析三类tag：male, female和shared
    for category in male female shared; do
        echo "处理${category} tag..."
        
        # 定义输出文件路径
        READS_FILE="$WORKDIR/specific_tags/${category}_with_sequence.tsv"
        STATS_FILE="$WORKDIR/specific_tags/${category}_tag_stats.tsv"
        
        # 初始化输出文件
        echo -e "Chromosome\tStart\tEnd\tReadID\tSequence\tMatchPos" > "$READS_FILE"
        echo -e "Chromosome\tStart\tEnd\tMatchedReads" > "$STATS_FILE"
        
        # 确定输入文件
        if [ "$category" == "shared" ]; then
            TAG_FILE="$WORKDIR/read_position/shared_tags.tsv"
            BAM_FILE="$WORKDIR/alignment/unique/male_merged.unique.bam"  # 使用任意一个样本的BAM
        else
            TAG_FILE="$WORKDIR/read_position/${category}_specific_tags.tsv"
            BAM_FILE="$WORKDIR/alignment/unique/${category}_merged.unique.bam"
        fi
        
        # 步骤1：提取候选reads
        samtools view -F 4 "$BAM_FILE" | \
            awk -v seq="$SPECIFIC_SEQUENCE" '
            BEGIN {OFS="\t"}
            $10 ~ seq {
                pos = index($10, seq);
                if (pos > 0) print $3, $4, $4+length($10)-1, $1, $10, pos;
            }' > tmp_candidates.tsv
        
        # 步骤2：生成reads明细文件
        awk -F'\t' '
        NR==FNR {tags[$1,$2,$3]=1; next}
        ($1,$2,$3) in tags {print $1, $2, $3, $4, $5, $6}
        ' "$TAG_FILE" tmp_candidates.tsv >> "$READS_FILE"
        
        # 步骤3：生成tag统计文件
        awk -F'\t' '
        NR==FNR {
            if (FNR > 1) {  # 跳过表头
                region = $1 "\t" $2 "\t" $3
                tag_regions[region] = 1
            }
            next
        }
        {
            region = $1 "\t" $2 "\t" $3
            if (region in tag_regions) {
                counts[region]++
            }
        }
        END {
            for (region in tag_regions) {
                print region "\t" (counts[region]+0)
            }
        }
        ' "$TAG_FILE" tmp_candidates.tsv | \
        sort -k1,1 -k2,2n >> "$STATS_FILE"
        
        # 清理临时文件
        rm tmp_candidates.tsv
        
        # 输出统计信息
        total_tags=$(tail -n +2 "$TAG_FILE" | wc -l)
        matched_tags=$(awk -F'\t' '$4 > 0' "$STATS_FILE" | wc -l)
        echo " - ${category} tag总数: $total_tags"
        echo " - 包含序列的tag数: $matched_tags"
        echo " - 匹配reads数: $(tail -n +2 "$READS_FILE" | wc -l)"
    done
}

# 在SNP calling前添加完整性检查
function check_bam_files {
    echo "[$(date)] 检查BAM文件完整性..."
    for bam in "$WORKDIR/alignment/unique/"*.bam; do
        # 检查主文件
        if [ ! -r "$bam" ]; then
            echo "错误: 无法读取文件 $bam"
            exit 1
        fi
        
        # 检查索引文件
        if [ ! -f "${bam}.csi" ] && [ ! -f "${bam}.bai" ]; then
            echo "警告: 缺少索引文件，尝试重新生成..."
            samtools index -M "$bam" || {
                echo "错误: 无法生成索引"
                exit 1
            }
        fi
        
        # 验证文件完整性
        samtools quickcheck "$bam" || {
            echo "错误: BAM文件损坏 $bam"
            exit 1
        }
    done
}

function run_snp_calling {
    echo "[$(date)] SNP calling..."
    for sample in "$MALE_SAMPLE" "$FEMALE_SAMPLE"; do
        # 1. 生成未过滤的VCF用于调试
        # bcftools mpileup -f "$REF_GENOME" \
            "$WORKDIR/alignment/unique/${sample}.unique.bam" \
            > "$WORKDIR/vcf/raw/${sample}.debug.vcf"
        
        # 2. 正式调用
        bcftools mpileup -Ou \
            -f "$REF_GENOME" \
            -a "AD,DP" \
            -q "$SNP_BASEQ" -Q "$SNP_MAPQ" \
            --max-depth "$MAX_DEPTH" \
            "$WORKDIR/alignment/unique/${sample}.unique.bam" \
            | bcftools call -mv \
            -o "$WORKDIR/vcf/raw/${sample}.raw.vcf.gz" \
            --ploidy "$PLOIDY" \
            2> "$WORKDIR/logs/${sample}_mpileup.log"
        
        # 3. 检查是否生成有效变异
        if [ $(bcftools view -H "$WORKDIR/vcf/raw/${sample}.raw.vcf.gz" | wc -l) -eq 0 ]; then
            echo "警告: ${sample} 未检测到变异，尝试降低过滤阈值"
            bcftools filter -i "QUAL >= 10 && FORMAT/DP >= 2" \
                -Oz -o "$WORKDIR/vcf/filtered/${sample}.filtered.vcf.gz" \
                "$WORKDIR/vcf/raw/${sample}.raw.vcf.gz"
        else
            bcftools filter -i "QUAL >= $SNP_QUALITY && FORMAT/DP >= $MIN_DP" \
                -Oz -o "$WORKDIR/vcf/filtered/${sample}.filtered.vcf.gz" \
                "$WORKDIR/vcf/raw/${sample}.raw.vcf.gz"
        fi
        
        bcftools index -f "$WORKDIR/vcf/filtered/${sample}.filtered.vcf.gz"
        bcftools stats "$WORKDIR/vcf/filtered/${sample}.filtered.vcf.gz" \
            > "$WORKDIR/snp/${sample}.stats.txt"

    done
}

function classify_snps {
    echo "[$(date)] SNP分类..."
    bcftools isec -n =2 -w1 \
        "$WORKDIR/vcf/filtered/${MALE_SAMPLE}.filtered.vcf.gz" \
        "$WORKDIR/vcf/filtered/${FEMALE_SAMPLE}.filtered.vcf.gz" \
        -p "$WORKDIR/snp/shared"
    
    for sample in "$MALE_SAMPLE" "$FEMALE_SAMPLE"; do
        bcftools view -i "MAF >= $MAF_THRESHOLD" \
            "$WORKDIR/vcf/filtered/${sample}.filtered.vcf.gz" \
            | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%DP\t[%AD]\n" \
            > "$WORKDIR/snp/${sample}_specific.tsv"
    done
}

function analyze_specific_snps {
    echo "[$(date)] 分析特异性SNP位点..."
    
    # 使用修正后的文件名
    male_file="$WORKDIR/snp/$MALE_SPECIFIC_FILE"
    female_file="$WORKDIR/snp/$FEMALE_SPECIFIC_FILE"
    
    # 验证输入文件
    if [[ ! -f "$male_file" || ! -f "$female_file" ]]; then
        echo "错误: 输入文件缺失！"
        echo " - 雄性文件: $male_file " $(test -f "$male_file" && echo "存在" || echo "缺失")
        echo " - 雌性文件: $female_file " $(test -f "$female_file" && echo "存在" || echo "缺失")
        exit 1
    fi

    # 筛选共有多态性位点（增加等位基因频率计算）
    echo "筛选雌雄共有多态性位点..."
    awk -v OFS='\t' '
    BEGIN {
        # 加载雌性数据到内存
        while (getline < "'"$female_file"'" > 0) {
            if (NF == 6) {
                female[$1,$2] = $0
                # 预计算雌性MAF
                split($6, f_ad, ",");
                f_total = f_ad[1] + f_ad[2];
                f_maf = (f_total > 0) ? f_ad[2]/f_total : 0;
                female_maf[$1,$2] = f_maf
            }
        }
    }
    NF == 6 && ($1,$2) in female {
        split($6, m_ad, ",");
        m_total = m_ad[1] + m_ad[2];
        m_maf = (m_total > 0) ? m_ad[2]/m_total : 0;
        
        split(female[$1,$2], f, "\t");
        split(f[6], f_ad, ",");
        f_total = f_ad[1] + f_ad[2];
        
        print $1, $2, $3, $4, $5, f[5], 
              $6, f[6], 
              m_maf, female_maf[$1,$2],
              (m_maf + female_maf[$1,$2])/2  # 平均MAF
    }' "$male_file" > tmp_both.tsv
    
    # 添加表头（包含MAF信息）
    { echo -e "Chrom\tPos\tRef\tAlt\tMale_DP\tFemale_DP\tMale_AD\tFemale_AD\tMale_MAF\tFemale_MAF\tAvg_MAF";
      cat tmp_both.tsv; } > "$WORKDIR/snp/both_polymorphic.tsv"

    # 筛选性别特异位点（保持不变）
    for gender in male female; do
        echo "筛选仅在${gender}中多态的位点..."
        current_file="$WORKDIR/snp/$([ "$gender" == "male" ] && echo "$MALE_SPECIFIC_FILE" || echo "$FEMALE_SPECIFIC_FILE")"
        opposite_file="$WORKDIR/snp/$([ "$gender" == "male" ] && echo "$FEMALE_SPECIFIC_FILE" || echo "$MALE_SPECIFIC_FILE")"
        
        awk -v OFS='\t' '
        BEGIN {
            # 加载相反性别数据
            while (getline < "'"$opposite_file"'" > 0) {
                if (NF == 6) opposite[$1,$2] = 1
            }
        }
        NF == 6 && !(($1,$2) in opposite) {
            split($6, ad, ",");
            total = ad[1] + ad[2];
            maf = (total > 0) ? ad[2]/total : 0;
            print $1, $2, $3, $4, $5, $6, maf
        }' "$current_file" > "tmp_${gender}.tsv"
        
        # 添加表头
        { echo -e "Chrom\tPos\tRef\tAlt\tDP\tAD\tMAF"; cat "tmp_${gender}.tsv"; } > "$WORKDIR/snp/${gender}_only_polymorphic.tsv"
    done
    
    # 清理临时文件
    rm tmp_*.tsv 2>/dev/null
    
    # 结果统计（排除表头）
    echo "SNP分类结果:"
    for file in both male_only female_only; do
        count=$(( $(wc -l < "$WORKDIR/snp/${file}_polymorphic.tsv") - 1 ))
        echo " - ${file}_polymorphic.tsv: $count 个位点"
    done
    
    # 输出MAF统计摘要
    echo "共有多态性位点MAF统计:"
    awk -F'\t' '
    NR>1 {
        male_sum += $9; 
        female_sum += $10;
        count++
    } 
    END {
        printf " - 雄性平均MAF: %.4f\n", male_sum/count;
        printf " - 雌性平均MAF: %.4f\n", female_sum/count;
    }' "$WORKDIR/snp/both_polymorphic.tsv"
}



# 主流程
run_fastp
run_alignment
bam_filter
extract_read_positions
run_tag_analysis
filter_specific_tags
check_bam_files
run_snp_calling
classify_snps
analyze_specific_snps

echo "分析完成! 结果目录: $WORKDIR"
