#!/usr/bin/env python3
# 6.4_paf_partition_ids.py
import sys, argparse, gzip

def gzopen(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def main():
    ap = argparse.ArgumentParser(description="Emit host read IDs based on PAF thresholds (reads vs human).")
    ap.add_argument("--paf", required=True, help="PAF (reads vs human). .gz ok")
    ap.add_argument("--out-host-ids", required=True, help="Output file of host read IDs")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-match-len", type=int, default=100)
    ap.add_argument("--min-frac", type=float, default=0.05)
    ap.add_argument("--min-pid", type=float, default=0.80)
    ap.add_argument("--min-total-aln", type=int, default=200)
    ap.add_argument("--min-total-frac", type=float, default=0.05)
    args = ap.parse_args()

    totals = {}   # qname -> (qlen, total_aln_bases)
    with gzopen(args.paf, "rt") as f:
        for line in f:
            if not line or line[0] == "#": 
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            qname = p[0]
            qlen  = int(p[1])
            qs    = int(p[2]); qe = int(p[3])
            mapq  = int(p[11])
            try:
                block_len = int(p[10])
            except Exception:
                block_len = max(0, qe - qs)
            pid = None
            try:
                matches = int(p[9])
                if block_len > 0:
                    pid = matches / block_len
            except Exception:
                pid = None
            for tag in p[12:]:
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
            aln_q = max(0, qe - qs)
            if mapq < args.min_mapq: 
                continue
            if aln_q < args.min_match_len:
                continue
            if qlen > 0 and (aln_q / qlen) < args.min_frac:
                continue
            if pid < args.min_pid:
                continue
            if qname not in totals:
                totals[qname] = [qlen, 0]
            totals[qname][1] += aln_q

    out = open(args.out_host_ids, "wt")
    host = 0
    for qname, (qlen, tot) in totals.items():
        if tot >= args.min_total_aln and (qlen == 0 or (tot / qlen) >= args.min_total_frac):
            out.write(qname + "\n")
            host += 1
    out.close()
    sys.stderr.write(f"[paf_partition_ids] host_ids={host}\n")

if __name__ == "__main__":
    main()