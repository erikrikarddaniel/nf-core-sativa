
process SATIVASCORE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // quay.io/biocontainers/python:3.11 has no build-hash-suffixed tag to pin to (unlike
    // real bioconda-recipe images); pin by digest instead so the underlying image can't
    // silently drift and shift floating-point tie-breaks in SATIVASCORE's rank voting.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python@sha256:b322907f8e52b2055ccad4e46848d28a4a5631b403116cc80ddf61ec8601e05e' }"

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
from collections import Counter, defaultdict

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
    # Per-rank majority vote among the leaves neighbouring one placement edge.
    # `names` is a set (from parse_jplace_tree), so iteration order depends on
    # PYTHONHASHSEED; sort it so Counter.most_common()'s tie-breaking (first
    # encountered wins) is deterministic across runs instead of picking a
    # different equally-weighted rank value each time.
    rows = [tax[n] for n in sorted(names) if n in tax]
    if not rows:
        return []
    n_ranks = min(len(r) for r in rows)
    return [Counter(r[i] for r in rows).most_common(1)[0][0] for i in range(n_ranks)]


def weighted_rank_votes(placement, edge_leaves, tax, fields, n_ranks):
    # A query's placement is a probability distribution (like_weight_ratio) over
    # candidate edges, not a single best edge. At each rank, sum the LWR of every
    # candidate whose neighbouring-leaf majority taxonomy agrees on a value; the
    # winner is the value with the most total weight, and that weight is the
    # confidence at that rank. Summing across the full candidate set (not just
    # the top-1 edge) matters: a query can be >95% confidently outside its
    # declared phylum while no single candidate edge holds more than ~10% of the
    # weight, because that ~95% is spread thinly across many edges within the
    # correct (different) clade.
    idx_edge = fields.index('edge_num')
    idx_lwr = fields.index('like_weight_ratio')

    rank_weights = [defaultdict(float) for _ in range(n_ranks)]
    for candidate in placement['p']:
        edge_num = str(candidate[idx_edge])
        lwr = candidate[idx_lwr]
        neighbours = edge_leaves.get(edge_num, set())
        local_taxonomy = majority_taxonomy(neighbours, tax)
        for i in range(min(n_ranks, len(local_taxonomy))):
            rank_weights[i][local_taxonomy[i]] += lwr

    predicted, confidence = [], []
    for weights in rank_weights:
        if not weights:
            predicted.append(None)
            confidence.append(0.0)
            continue
        best_value, best_weight = max(weights.items(), key=lambda kv: kv[1])
        predicted.append(best_value)
        confidence.append(best_weight)
    return predicted, confidence


def first_mismatch(predicted, original):
    for i, (p, o) in enumerate(zip(predicted, original)):
        if p != o:
            return i
    return None


parser = argparse.ArgumentParser()
parser.add_argument('taxonomy')
parser.add_argument('mislabels_tsv')
parser.add_argument('summary_txt')
parser.add_argument('--min-lwr', type=float, default=0.5,
                     help='Minimum aggregated placement-weight confidence at the mismatching rank required to trust it enough to flag a mislabel.')
opts = parser.parse_args()

tax = load_taxonomy(opts.taxonomy)

rows = []
for jplace_file in sorted(glob.glob('placements/*.jplace.gz')):
    with gzip.open(jplace_file, 'rt') as fh:
        data = json.load(fh)

    edge_leaves = parse_jplace_tree(data['tree'])
    fields = data['fields']

    for placement in data['placements']:
        query_name = placement['n'][0] if 'n' in placement else placement['nm'][0][0]
        original = tax.get(query_name, [])

        predicted, confidence = weighted_rank_votes(placement, edge_leaves, tax, fields, len(original))

        mismatch_rank = first_mismatch(predicted, original)
        rank_confidence = confidence[mismatch_rank] if mismatch_rank is not None else None
        is_mislabel = mismatch_rank is not None and rank_confidence >= opts.min_lwr

        rows.append({
            'seq_name': query_name,
            'original_label': ';'.join(original),
            'predicted_label': ';'.join(str(p) for p in predicted),
            'lwr': rank_confidence,
            'mismatch_rank': mismatch_rank + 1 if mismatch_rank is not None else None,
            'is_mislabel': is_mislabel,
        })

with open(opts.mislabels_tsv, 'w') as fh:
    print('seq_name\\toriginal_label\\tpredicted_label\\tlwr\\tmismatch_rank\\tmethod', file=fh)
    for row in rows:
        if row['is_mislabel']:
            print(f"{row['seq_name']}\\t{row['original_label']}\\t{row['predicted_label']}\\t{row['lwr']:.6f}\\t{row['mismatch_rank']}\\tsativa", file=fh)

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
