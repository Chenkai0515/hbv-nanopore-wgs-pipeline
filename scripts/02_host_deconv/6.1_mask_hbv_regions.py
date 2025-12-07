#!/usr/bin/env python3
import sys, gzip, argparse
from collections import defaultdict

def read_paf_intervals(paf_path, min_mapq=20, merge_distance=10, min_match_len=0):
    """
    解析PAF，返回 {read_id: [(qs, qe), ...]} ，按读段(query)坐标。
    仅保留 mapq >= min_mapq 的命中，可设置 merge_distance 合并间隔。
    """
    intervals = defaultdict(list)
    # 读PAF（可为 .gz 或纯文本）
    _open = gzip.open if paf_path.endswith(".gz") else open
    with _open(paf_path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:  # PAF基础字段不足
                continue
            qname = parts[0]
            try:
                qstart = int(parts[2]); qend = int(parts[3])
                mapq  = int(parts[11])
            except Exception:
                continue

            if mapq < min_mapq:
                continue
            if min_match_len > 0:
                # 第10列为 #matches, 第11列为 aln block length；此处用 qend-qstart 作为近似过滤
                if (qend - qstart) < min_match_len:
                    continue
            if qend > qstart:
                intervals[qname].append((qstart, qend))

    # 合并重叠/近邻区间
    merged = {}
    for rid, ivs in intervals.items():
        ivs.sort()
        out = []
        for s, e in ivs:
            if not out:
                out.append([s, e])
            else:
                ps, pe = out[-1]
                if s <= pe + merge_distance:  # 近邻也并入
                    out[-1][1] = max(pe, e)
                else:
                    out.append([s, e])
        merged[rid] = [(s, e) for s, e in out]
    return merged

def mask_fastq(fq_in, fq_out, iv_dict):
    """
    按 iv_dict 中给定的区间，把FASTQ中对应read的序列置为'N'。
    质量值保持不变。
    """
    _in = gzip.open(fq_in, "rt") if fq_in.endswith(".gz") else open(fq_in, "rt")
    _out = gzip.open(fq_out, "wt") if fq_out.endswith(".gz") else open(fq_out, "wt")

    masked_reads = 0
    masked_bases = 0
    total_reads  = 0

    with _in as fin, _out as fout:
        while True:
            name = fin.readline()
            if not name:
                break
            seq  = fin.readline().rstrip("\n")
            plus = fin.readline()
            qual = fin.readline().rstrip("\n")
            total_reads += 1
            rid = name.strip().split()[0]
            rid = rid[1:] if rid.startswith("@") else rid  # 去掉开头的'@'

            if rid in iv_dict and iv_dict[rid]:
                s_list = list(seq)
                # 统计被mask的碱基数量
                for s, e in iv_dict[rid]:
                    # 防御性截断
                    s0 = max(0, min(len(s_list), s))
                    e0 = max(0, min(len(s_list), e))
                    if e0 > s0:
                        for i in range(s0, e0):
                            s_list[i] = 'N'
                        masked_bases += (e0 - s0)
                seq = "".join(s_list)
                masked_reads += 1

            fout.write(name)
            fout.write(seq + "\n")
            fout.write(plus)
            fout.write(qual + "\n")

    sys.stderr.write(f"[mask_hbv_regions] total_reads={total_reads}, masked_reads={masked_reads}, masked_bases={masked_bases}\n")

def main():
    ap = argparse.ArgumentParser(description="Mask HBV-matching regions in FASTQ using PAF alignments.")
    ap.add_argument("--paf", required=True, help="PAF file from minimap2 (reads vs HBV panel).")
    ap.add_argument("--fastq", required=True, help="Input FASTQ(.gz).")
    ap.add_argument("--out", required=True, help="Output masked FASTQ(.gz).")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--merge-distance", type=int, default=10)
    ap.add_argument("--min-match-len", type=int, default=0, help="Minimum (qend-qstart) to keep an alignment block.")
    args = ap.parse_args()

    ivs = read_paf_intervals(args.paf, args.min_mapq, args.merge_distance, args.min_match_len)
    mask_fastq(args.fastq, args.out, ivs)

if __name__ == "__main__":
    main()
