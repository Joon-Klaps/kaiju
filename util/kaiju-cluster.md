# kaiju-cluster

Single-command tool to run Kaiju (with `-v`) and cluster reads/contigs by Jaccard similarity over top accessions, or to cluster an existing Kaiju verbose TSV.

Two modes:

- Run Kaiju and cluster:
  - `utils/kaiju-cluster.py run -t nodes.dmp -f kaiju_db.fmi -i contigs.fasta -o membership.tsv [kaiju opts] [cluster opts]`
- Cluster an existing Kaiju TSV:
  - `utils/kaiju-cluster.py cluster -i kaiju.out -o membership.tsv [cluster opts]`

Input (Kaiju `-v` TSV) columns:

1) C/U
2) read_id
3) taxon_id
4) metric (score or length)
5) ids
6) accessions
7) fragments

Output (membership.tsv):
  centroid\tmember
  acc1\tacc1
  acc1\tacc2
  acc3\tacc3
  acc3\tacc4

Tips:

- Use `--max-degree` to ignore ubiquitous proteins.
- Use `--min-metric` to filter weak matches (e.g., 65 in Greedy mode).
