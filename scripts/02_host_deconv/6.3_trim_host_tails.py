#!/usr/bin/env python3
# 6.3_trim_host_tails.py
import sys, argparse, gzip
from typing import List, Tuple, Dict

def gzopen(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def parse_paf(paf_path: str, min_mapq: int, min_match_len: int, min_frac: float, min_pid: float) -> Dict[str, List[Tuple[int,int,int]]]:
    hits = {}
    with gzopen(paf_path, "rt") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            qlen  = int(parts[1])
            qs    = int(parts[2])
            qe    = int(parts[3])
            mapq  = int(parts[11])
            try:
                block_len = int(parts[10])
            except Exception:
                block_len = max(0, qe - qs)
            pid = None
            try:
                matches = int(parts[9])
                if block_len > 0:
                    pid = matches / block_len
            except Exception:
                pid = None
            for tag in parts[12:]:
                if tag.startswith("NM:i:"):
                    try:
                        nm = int(tag[5:])
                        if block_len > 0:
                            pid_nm = max(0.0, 1.0 - nm / block_len)
                            if pid is None:
                                pid = pid_nm
                    except Exception:
                        pass
            if pid is None:
                pid = 0.0
            if mapq < min_mapq:
                continue
            aln_q = max(0, qe - qs)
            if aln_q < min_match_len:
                continue
            if qlen > 0 and (aln_q / qlen) < min_frac:
                continue
            if pid < min_pid:
                continue
            hits.setdefault(qname, []).append((qs, qe, qlen))
    return hits

def merge_intervals(ints: List[Tuple[int,int]], merge_dist: int) -> List[Tuple[int,int]]:
    if not ints:
        return []
    ints = sorted(ints, key=lambda x: (x[0], x[1]))
    out = []
    cs, ce = ints[0]
    for s, e in ints[1:]:
        if s <= ce + merge_dist:
            if e > ce: ce = e
        else:
            out.append((cs, ce))
            cs, ce = s, e
    out.append((cs, ce))
    return out

def classify_and_clip(qlen: int, hit_qints: List[Tuple[int,int]], tail_near_end: int, tail_max_len: int, tail_max_frac: float, merge_dist: int, min_keep_len: int, discard_internal: bool):
    if not hit_qints:
        return True, 0, 0, "no_trusted_human_hits"
    merged = merge_intervals(hit_qints, merge_dist=merge_dist)
    left_ints, right_ints, internal = [], [], []
    for s, e in merged:
        left = (s <= tail_near_end)
        right = ((qlen - e) <= tail_near_end)
        if left and not right:
            left_ints.append((s, e))
        elif right and not left:
            right_ints.append((s, e))
        elif left and right:
            internal.append((s, e))
        else:
            internal.append((s, e))
    if discard_internal and internal:
        return False, 0, 0, "internal_human_segments"
    trim_left = 0
    trim_right = 0
    if left_ints:
        le = max(e for _, e in left_ints)
        trim_left = le
    if right_ints:
        rs = min(s for s, _ in right_ints)
        trim_right = qlen - rs
    if trim_left + trim_right >= qlen:
        return False, 0, 0, "clipped_to_empty"
    max_allowed = min(tail_max_len, int(tail_max_frac * qlen))
    if (trim_left > 0 and trim_left > max_allowed) or (trim_right > 0 and trim_right > max_allowed):
        return False, 0, 0, "tail_too_long"
    kept_len = qlen - trim_left - trim_right
    if kept_len < min_keep_len:
        return False, 0, 0, "kept_len_below_min"
    return True, trim_left, trim_right, ("tails_clipped" if (trim_left or trim_right) else "kept_unchanged")

def iter_fastq(fh):
    while True:
        name = fh.readline()
        if not name:
            return
        seq = fh.readline()
        plus = fh.readline()
        qual = fh.readline()
        if not qual:
            return
        yield name.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")

def main():
    ap = argparse.ArgumentParser(description="Trim human tails from Nanopore reads using PAF (human).")
    ap.add_argument("--paf", required=True, help="PAF (human). .gz ok")
    ap.add_argument("--fastq", required=True, help="Original UNMASKED FASTQ(.gz)")
    ap.add_argument("--out-viral", required=True, help="Output viral_enriched.unmasked.fastq.gz")
    ap.add_argument("--out-host", required=True, help="Output host_only.unmasked.fastq.gz")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-match-len", type=int, default=100)
    ap.add_argument("--min-frac", type=float, default=0.05)
    ap.add_argument("--min-pid", type=float, default=0.80)
    ap.add_argument("--tail-near-end", type=int, default=80)
    ap.add_argument("--tail-max-len", type=int, default=300)
    ap.add_argument("--tail-max-frac", type=float, default=0.20)
    ap.add_argument("--merge-distance", type=int, default=10)
    ap.add_argument("--min-keep-len", type=int, default=800)
    ap.add_argument("--discard-internal", type=int, default=1)
    args = ap.parse_args()

    hits = parse_paf(args.paf, args.min_mapq, args.min_match_len, args.min_frac, args.min_pid)
    total = kept = trimmed = host_only = trimmed_bp = 0
    reasons = {}

    with gzopen(args.fastq, "rt") as fin, gzip.open(args.out_viral, "wt") as fv, gzip.open(args.out_host, "wt") as fh:
        for name, seq, plus, qual in iter_fastq(fin):
            total += 1
            rid = name[1:].split()[0]
            qlen = len(seq)
            qints = [(s, e) for (s, e, L) in hits.get(rid, [])]
            keep, tl, tr, why = classify_and_clip(
                qlen, qints,
                args.tail_near_end, args.tail_max_len, args.tail_max_frac,
                args.merge_distance, args.min_keep_len, bool(args.discard_internal)
            )
            reasons[why] = reasons.get(why, 0) + 1
            if keep:
                if tl or tr:
                    new_seq = seq[tl: qlen - tr]
                    new_qual = qual[tl: qlen - tr]
                    new_name = f"{name.rstrip()} clip=L{tl},R{tr}\n"
                    fv.write(new_name)
                    fv.write(new_seq + "\n")
                    fv.write(plus + "\n")
                    fv.write(new_qual + "\n")
                    trimmed += 1
                    trimmed_bp += (tl + tr)
                else:
                    fv.write(name if name.endswith("\n") else name + "\n")
                    fv.write(seq if seq.endswith("\n") else seq + "\n")
                    fv.write(plus if plus.endswith("\n") else plus + "\n")
                    fv.write(qual if qual.endswith("\n") else qual + "\n")
                kept += 1
            else:
                fh.write(name if name.endswith("\n") else name + "\n")
                fh.write(seq if seq.endswith("\n") else seq + "\n")
                fh.write(plus if plus.endswith("\n") else plus + "\n")
                fh.write(qual if qual.endswith("\n") else qual + "\n")
                host_only += 1

    def pct(n): 
        return 0.0 if total == 0 else (100.0 * n / total)
    sys.stderr.write(f"[trim_host_tails] total_reads={total}\tkept={kept}\ttrimmed_reads={trimmed}\ttrimmed_bp={trimmed_bp}\thost_only={host_only}\n")
    for k in sorted(reasons):
        sys.stderr.write(f"[trim_host_tails] reason_count\t{k}={reasons[k]}\n")

if __name__ == "__main__":
    main()