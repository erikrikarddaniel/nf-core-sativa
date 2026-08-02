process RESOLVETAXONOMY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'quay.io/biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(taxonomy), path(sequences)

    output:
    tuple val(meta), path("*.resolved.tax"),                  emit: taxonomy
    tuple val(meta), path("*.resolved.${sequences.extension}"), emit: sequences
    tuple val(meta), path("*.warnings.txt"),                  emit: warnings
    path "versions.yml",                                      emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Nextflow stages an absent optional path(taxonomy) as an empty list -- falsy in
    // Groovy -- rather than as a file, so this correctly distinguishes "no --taxonomy
    // file given" from a real one, without ever interpolating the literal text "[]"
    // into the command line.
    def taxonomy_in = taxonomy ? "${taxonomy}" : ''
    """
    python3 - "${taxonomy_in}" "${sequences}" "${prefix}.resolved.tax" "${prefix}.resolved.${sequences.extension}" "${prefix}.warnings.txt" << 'PYEOF'
import sys
from Bio import SeqIO

taxonomy_in, sequences_in, taxonomy_out, sequences_out, warnings_out = sys.argv[1:6]

# Format sniffed from content (FASTA, Clustal or PHYLIP -- --sequences can be any of
# these; this runs before the pipeline's own EMBOSS_SEQRET normalisation to FASTA).
# Matches CHECKNAMECONSISTENCY's own sniffing logic, which runs right after this.
with open(sequences_in) as fh:
    first_line = next((l.strip() for l in fh if l.strip()), '')
if first_line.startswith('>'):
    sequences_format = 'fasta'
elif first_line.upper().startswith('CLUSTAL'):
    sequences_format = 'clustal'
else:
    sequences_format = 'phylip-relaxed'

records = list(SeqIO.parse(sequences_in, sequences_format))
warnings = []

def embedded_taxonomy(record):
    # record.description is the *whole* header line (id + any trailing text);
    # record.id is just its first token -- GTDB's own single-file convention puts
    # the taxonomy string right after the id, space-separated.
    return record.description[len(record.id):].strip()

if taxonomy_in:
    # An explicit --taxonomy file always wins. Warn (not fail) -- surfaced by the
    # caller via log.warn, not just buried in this task's own stderr -- rather than
    # silently using the wrong source if the sequences also happen to carry
    # embedded text.
    if any(embedded_taxonomy(record) for record in records):
        warnings.append(
            '--taxonomy was provided; ignoring embedded taxonomy text found in '
            '--sequences record headers.'
        )
    with open(taxonomy_in) as fh_in, open(taxonomy_out, 'w') as fh_out:
        fh_out.write(fh_in.read())
else:
    missing = [record.id for record in records if not embedded_taxonomy(record)]
    if missing:
        sys.exit(
            'No --taxonomy file was provided, and these --sequences records have no '
            'embedded taxonomy in their header either: ' + ', '.join(missing)
        )
    with open(taxonomy_out, 'w') as fh:
        for record in records:
            print(f"{record.id}\\t{embedded_taxonomy(record)}", file=fh)

# Always strip headers down to a bare id -- downstream tools (IQTREE, EPA-ng) keep
# the whole header line as the leaf name, not just the first token, so leftover
# embedded-taxonomy text risks the same kind of name-mangling bug already fixed for
# '~'/'#' characters in CHECKNAMECONSISTENCY.
for record in records:
    record.description = record.id
SeqIO.write(records, sequences_out, sequences_format)

with open(warnings_out, 'w') as fh:
    for warning in warnings:
        print(warning, file=fh)

with open('versions.yml', 'w') as fh:
    print('"${task.process}":', file=fh)
    print('    python: ' + sys.version.split()[0], file=fh)
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.resolved.tax ${prefix}.resolved.${sequences.extension} ${prefix}.warnings.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
