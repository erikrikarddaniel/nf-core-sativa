
process SATIVALOOSPLIT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(alignment), path(tree)

    output:
    tuple val(meta), path("queryaln/*.fasta"), path("referencealn/*.fasta"), path("referencetree/*.nwk"), emit: loo
    path "versions.yml",                                                                                  emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 - "${alignment}" "${tree}" << 'PYEOF'
import sys
import os
import Bio
from Bio import SeqIO, Phylo

alignment_file, tree_file = sys.argv[1], sys.argv[2]

os.makedirs('queryaln', exist_ok=True)
os.makedirs('referencealn', exist_ok=True)
os.makedirs('referencetree', exist_ok=True)

# The subworkflow normalises whatever format the pipeline was given (FASTA,
# Clustal, PHYLIP) to FASTA via EMBOSS_SEQRET before calling this module, so
# there is exactly one format to parse here. Not PHYLIP: EMBOSS's phylip writer
# truncates sequence names to 10 characters, silently colliding (and corrupting
# the alignment) for anything with longer real-world identifiers.
records = list(SeqIO.parse(alignment_file, 'fasta'))
if not records:
    sys.exit('No sequences found in: ' + alignment_file)

for idx, record in enumerate(records, start=1):
    tag = f"{idx:04d}_{record.id}"

    # Suffixes must differ across the three output types: EPANG_PLACE stages
    # queryaln/referencealn/referencetree side by side in the same task work
    # directory, and same-named files (e.g. two "<tag>.fasta") collide there.
    SeqIO.write([record], f"queryaln/{tag}.query.fasta", "fasta")
    SeqIO.write([r for r in records if r.id != record.id], f"referencealn/{tag}.reference.fasta", "fasta")

    pruned = Phylo.read(tree_file, "newick")
    pruned.prune(record.id)
    Phylo.write(pruned, f"referencetree/{tag}.tree.nwk", "newick")

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
    print('    biopython: ' + Bio.__version__, file=fh)
PYEOF
    """

    stub:
    """
    mkdir -p queryaln referencealn referencetree
    touch queryaln/0001_stub.query.fasta referencealn/0001_stub.reference.fasta referencetree/0001_stub.tree.nwk
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
