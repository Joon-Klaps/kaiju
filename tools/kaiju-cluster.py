#!/usr/bin/env python3
"""
kaiju-cluster: single-command tool to
  - run Kaiju (-v) and cluster by Jaccard over top accessions, or
  - cluster an existing Kaiju verbose TSV.

Usage examples:
  # Run Kaiju and cluster
  kaiju-cluster run \
    --nodes nodes.dmp --fmi kaiju_db.fmi -i contigs.fasta \
    --out-clusters membership.tsv [kaiju options...] [cluster options...]

  # Cluster an existing Kaiju -v TSV
  kaiju-cluster cluster -i kaiju.out -o membership.tsv [cluster options...]
"""
import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from collections import defaultdict


def parse_accessions_field(field: str):
    if not field:
        return []
    parts = [p for p in field.strip().split(',') if p]
    return parts


def load_evidence(tsv, min_metric=None):
    acc_to_reads = defaultdict(set)
    read_to_accs = defaultdict(set)
    read_metric = {}
    for line in tsv:
        if not line.strip():
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 7:
            continue
        status, read_id, _taxon, metric, _ids, accs, _frags = cols[:7]
        if status not in ('C', 'U'):
            continue
        try:
            metric_val = float(metric)
        except ValueError:
            metric_val = None
        if min_metric is not None and metric_val is not None and metric_val < min_metric:
            continue
        acc_list = parse_accessions_field(accs)
        if not acc_list:
            continue
        s = read_to_accs[read_id]
        for a in acc_list:
            s.add(a)
            acc_to_reads[a].add(read_id)
        if metric_val is not None:
            read_metric[read_id] = max(metric_val, read_metric.get(read_id, metric_val))
    return acc_to_reads, read_to_accs, read_metric


class DSU:
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            return x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return False
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1
        return True


def jaccard(a, b):
    if not a or not b:
        return 0.0
    inter = len(a & b)
    if inter == 0:
        return 0.0
    union = len(a | b)
    return inter / union


def do_cluster(acc_to_reads, read_to_accs, min_jaccard=0.5, min_shared=None, max_degree=None):
    dsu = DSU()
    for reads in acc_to_reads.values():
        if max_degree is not None and len(reads) > max_degree:
            continue
        rl = list(reads)
        if len(rl) < 2:
            continue
        anchor = rl[0]
        for r in rl[1:]:
            if min_shared is not None:
                shared = len(read_to_accs[anchor] & read_to_accs[r])
                if shared < min_shared:
                    continue
            jac = jaccard(read_to_accs[anchor], read_to_accs[r])
            if jac >= min_jaccard:
                dsu.union(anchor, r)
    clusters = defaultdict(list)
    for r in read_to_accs.keys():
        root = dsu.find(r)
        clusters[root].append(r)
    return clusters


def write_membership(clusters, out):
    for _idx, (_root, members) in enumerate(clusters.items(), 1):
        centroid = min(members)
        for m in members:
            out.write(f"{centroid}\t{m}\n")


def default_kaiju_bin(repo_root: Path) -> Path:
    return repo_root / "bin" / "kaiju"


def build_kaiju_cmd(args, kaiju_bin: Path):
    cmd = [str(kaiju_bin), "-t", args.nodes, "-f", args.fmi, "-i", args.input]
    if args.input2:
        cmd += ["-j", args.input2]
    if args.mode:
        cmd += ["-a", args.mode]
    if args.protein:
        cmd += ["-p"]
    if args.disable_seg:
        cmd += ["-X"]
    elif args.enable_seg:
        cmd += ["-x"]
    if args.threads is not None:
        cmd += ["-z", str(args.threads)]
    if args.min_match_len is not None:
        cmd += ["-m", str(args.min_match_len)]
    if args.min_score is not None:
        cmd += ["-s", str(args.min_score)]
    if args.mismatches is not None:
        cmd += ["-e", str(args.mismatches)]
    if args.min_evalue is not None:
        cmd += ["-E", str(args.min_evalue)]
    # Exposed caps
    if args.max_matches_si is not None:
        cmd += ["-S", str(args.max_matches_si)]
    if args.max_match_ids is not None:
        cmd += ["-I", str(args.max_match_ids)]
    if args.max_match_acc is not None:
        cmd += ["-A", str(args.max_match_acc)]
    # Always verbose
    cmd.append("-v")
    return cmd


def main():
    repo_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description="Run Kaiju and/or cluster by Jaccard in one command")
    sub = ap.add_subparsers(dest="cmd", required=True)

    # run: run Kaiju then cluster
    apr = sub.add_parser("run", help="Run Kaiju (-v) and cluster the results")
    apr.add_argument("--nodes", "-t", required=True)
    apr.add_argument("--fmi", "-f", required=True)
    apr.add_argument("--input", "-i", required=True)
    apr.add_argument("--input2", "-j", default=None)
    apr.add_argument("--kaiju-bin", default=str(default_kaiju_bin(repo_root)))
    apr.add_argument("--mode", "-a", choices=["mem", "greedy"], default=None)
    apr.add_argument("--threads", "-z", type=int, default=None)
    apr.add_argument("--protein", "-p", action="store_true")
    apr.add_argument("--enable-seg", "-x", action="store_true")
    apr.add_argument("--disable-seg", "-X", action="store_true")
    apr.add_argument("--min-match-len", "-m", type=int, default=None)
    apr.add_argument("--min-score", "-s", type=int, default=None)
    apr.add_argument("--mismatches", "-e", type=int, default=None)
    apr.add_argument("--min-evalue", "-E", type=float, default=None)
    apr.add_argument("--max-matches-si", "-S", type=int, default=None)
    apr.add_argument("--max-match-ids", "-I", type=int, default=None)
    apr.add_argument("--max-match-acc", "-A", type=int, default=None)
    # cluster outputs
    apr.add_argument("--out-clusters", "-o", required=True)
    # cluster options
    apr.add_argument("--min-jaccard", type=float, default=0.5)
    apr.add_argument("--min-metric", type=float, default=None)
    apr.add_argument("--min-shared", type=int, default=None)
    apr.add_argument("--max-degree", type=int, default=None)

    # cluster: only cluster an existing Kaiju TSV
    apc = sub.add_parser("cluster", help="Cluster an existing Kaiju verbose TSV (-v)")
    apc.add_argument("--input", "-i", required=True, help="Kaiju -v TSV")
    apc.add_argument("--output", "-o", required=True, help="Output centroid\\tmember TSV")
    apc.add_argument("--min-jaccard", type=float, default=0.5)
    apc.add_argument("--min-metric", type=float, default=None)
    apc.add_argument("--min-shared", type=int, default=None)
    apc.add_argument("--max-degree", type=int, default=None)

    args = ap.parse_args()

    if args.cmd == "run":
        kaiju_bin = Path(args.kaiju_bin)
        if not kaiju_bin.exists():
            ap.error(f"kaiju executable not found at {kaiju_bin}")
        with tempfile.NamedTemporaryFile(prefix="kaiju_out_", suffix=".tsv", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            cmd = build_kaiju_cmd(args, kaiju_bin)
            with open(tmp_path, "w", encoding="utf-8") as outfh:
                subprocess.run(cmd, check=True, stdout=outfh)
            with open(tmp_path, "r", encoding="utf-8") as infh:
                acc_to_reads, read_to_accs, _rm = load_evidence(infh, args.min_metric)
            clusters = do_cluster(acc_to_reads, read_to_accs,
                                  min_jaccard=args.min_jaccard,
                                  min_shared=args.min_shared,
                                  max_degree=args.max_degree)
            os.makedirs(Path(args.out_clusters).parent, exist_ok=True)
            with open(args.out_clusters, "w", encoding="utf-8") as out:
                write_membership(clusters, out)
        finally:
            try:
                os.remove(tmp_path)
            except FileNotFoundError:
                pass
    elif args.cmd == "cluster":
        with open(args.input, "r", encoding="utf-8") as fh:
            acc_to_reads, read_to_accs, _rm = load_evidence(fh, args.min_metric)
        clusters = do_cluster(acc_to_reads, read_to_accs,
                              min_jaccard=args.min_jaccard,
                              min_shared=args.min_shared,
                              max_degree=args.max_degree)
        os.makedirs(Path(args.output).parent, exist_ok=True)
        with open(args.output, "w", encoding="utf-8") as out:
            write_membership(clusters, out)


if __name__ == "__main__":
    main()
