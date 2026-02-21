#!/bin/bash

INDEX="data/salmon_index/index"
OUTDIR="results/salmon_quant"

mkdir -p $OUTDIR

for fn in data/filtered/SRR105516{57..65}.clean.fastq.gz; do
    # Extract sample name from path
    samp=$(basename $fn .clean.fastq.gz)

    # Identify which file is processing
    echo "Processing sample ${samp}"

    # Run Salmon
    salmon quant -i $INDEX -l A \
        -r $fn \
        --validateMappings \
        --softclip \
	-p 12 \
        -o $OUTDIR/${samp}
done
