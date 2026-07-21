
process SATIVASCORE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    // All N per-sequence jplace files for one dataset are staged under placements/
    tuple val(meta), path(jplace_files, stageAs: "placements/*"), path(taxonomy)

    output:
    tuple val(meta), path("*.mislabels.tsv"), emit: mislabels
    tuple val(meta), path("*.summary.txt"),   emit: summary
    path "versions.yml",                      emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ''
    """
    python3 - "${taxonomy}" "${prefix}.mislabels.tsv" "${prefix}.summary.txt" ${args} << 'PYEOF'
import sys
import glob
import gzip
import json
import argparse
from collections import Counter

def parse_jplace_tree(tree_str):
    # Extended Newick used by EPA-ng: each node's branch length is followed by
    # "{edge_num}". Returns {edge_num: set(leaf names in the subtree below it)}.
    tree_str = tree_str.strip().rstrip(';')
    pos = 0
    edge_leaves = {}

    def parse_clade():
        nonlocal pos
        leaves = set()
        if tree_str[pos] == '(':
            pos += 1
            while True:
                leaves |= parse_clade()
                if tree_str[pos] == ',':
                    pos += 1
                    continue
                if tree_str[pos] == ')':
                    pos += 1
                    break
        else:
            start = pos
            while tree_str[pos] not in ':,(){};':
                pos += 1
            leaves = {tree_str[start:pos]}
        if pos < len(tree_str) and tree_str[pos] == ':':
            pos += 1
            while pos < len(tree_str) and tree_str[pos] not in ',(){};':
                pos += 1
        if pos < len(tree_str) and tree_str[pos] == '{':
            pos += 1
            start = pos
            while tree_str[pos] != '}':
                pos += 1
            edge_leaves[tree_str[start:pos]] = leaves
            pos += 1
        return leaves

    parse_clade()
    return edge_leaves


def load_taxonomy(path):
    tax = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\\n')
            if not line:
                continue
            name, _, ranks = line.partition('\\t')
            tax[name] = [r.strip() for r in ranks.split(';')]
    return tax


def majority_taxonomy(names, tax):
    # Per-rank majority vote among the leaves neighbouring the placement edge.
    rows = [tax[n] for n in names if n in tax]
    if not rows:
        return []
    n_ranks = min(len(r) for r in rows)
    return [Counter(r[i] for r in rows).most_common(1)[0][0] for i in range(n_ranks)]


def first_mismatch_rank(predicted, original):
    for i, (p, o) in enumerate(zip(predicted, original), start=1):
        if p != o:
            return i
    return None


parser = argparse.ArgumentParser()
parser.add_argument('taxonomy')
parser.add_argument('mislabels_tsv')
parser.add_argument('summary_txt')
parser.add_argument('--min-lwr', type=float, default=0.5,
                     help='Minimum like_weight_ratio of the best placement required to trust it enough to flag a mismatch.')
opts = parser.parse_args()

tax = load_taxonomy(opts.taxonomy)

rows = []
for jplace_file in sorted(glob.glob('placements/*.jplace.gz')):
    with gzip.open(jplace_file, 'rt') as fh:
        data = json.load(fh)

    edge_leaves = parse_jplace_tree(data['tree'])
    fields = data['fields']
    idx_edge = fields.index('edge_num')
    idx_lwr = fields.index('like_weight_ratio')

    for placement in data['placements']:
        query_name = placement['n'][0] if 'n' in placement else placement['nm'][0][0]

        best = max(placement['p'], key=lambda p: p[idx_lwr])
        edge_num = str(best[idx_edge])
        lwr = best[idx_lwr]

        neighbours = edge_leaves.get(edge_num, set())
        predicted = majority_taxonomy(neighbours, tax)
        original = tax.get(query_name, [])

        mismatch_rank = first_mismatch_rank(predicted, original)
        is_mislabel = mismatch_rank is not None and lwr >= opts.min_lwr

        rows.append({
            'seq_name': query_name,
            'original_label': ';'.join(original),
            'predicted_label': ';'.join(predicted),
            'lwr': lwr,
            'mismatch_rank': mismatch_rank,
            'is_mislabel': is_mislabel,
        })

with open(opts.mislabels_tsv, 'w') as fh:
    print('seq_name\\toriginal_label\\tpredicted_label\\tlwr\\tmismatch_rank', file=fh)
    for row in rows:
        if row['is_mislabel']:
            print(f"{row['seq_name']}\\t{row['original_label']}\\t{row['predicted_label']}\\t{row['lwr']:.6f}\\t{row['mismatch_rank']}", file=fh)

with open(opts.summary_txt, 'w') as fh:
    print(f"sequences scored: {len(rows)}", file=fh)
    print(f"min_lwr threshold: {opts.min_lwr}", file=fh)
    print(f"putative mislabels: {sum(1 for r in rows if r['is_mislabel'])}", file=fh)
    no_neighbours = sum(1 for r in rows if not r['predicted_label'])
    if no_neighbours:
        print(f"queries with no neighbouring reference leaves (best edge not found in tree): {no_neighbours}", file=fh)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mislabels.tsv ${prefix}.summary.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
