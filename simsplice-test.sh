#!/usr/bin/env bash
set -eu -o pipefail

mkdir -p /data/iarpa/analysis/PRJNA607948
cd /data/iarpa/analysis/PRJNA607948

#fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip --origfmt

@sra-toolkit bash -x <<'EOS'
esearch -db sra -query PRJNA607948 |efetch -format native |jx -xPo- -e '_.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE' -e 'a=_.RUN_SET.RUN.$accession; t=_.EXPERIMENT.TITLE.replace(/[ ]/g,"_"); o=t.replace(/_[^_]+$/,""); t.match(/SISPA/)? [] : `fastq-dump --outdir ${q(o)} --gzip --skip-technical --split-3 --clip --dumpbase --read-filter pass ${q(a)} && (cd ${q(o)} && t=${q(t)} rename -f '\''s/^/$ENV{t}./'\'' ${q(a)}*.fastq.gz)`' |parallel -q bash -xc
EOS

find * -maxdepth 0 -type d |parallel -q @unicycler bash -xc 'unicycler -1 "$0"/*_illumina.*_1.fastq.gz -2 "$0"/*_illumina.*_2.fastq.gz -l "$0"/*_ONT.*.fastq.gz -o "$0" -t "$(nproc)"'

find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 --eqx -aYx sr -t "$(nproc)" assembly.fasta <(gzip -cd *_illumina.*_1.fastq.gz) <(gzip -cd *_illumina.*_2.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o "$0"_illumina.bam && samtools index "$0"_illumina.bam'

find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 --eqx -aYx map-ont -t "$(nproc)" assembly.fasta <(gzip -cd *_ONT.*.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o "$0"_ONT.bam && samtools index "$0"_ONT.bam'

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && ~/src/simsplice/target/debug/genvcf --genome assembly.fasta --vcf modifications.vcf'

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && ~/src/simsplice/target/debug/simsplice --bam "$(basename "$0")_illumina.bam" --genome assembly.fasta --vcf modifications.vcf --outgenome modified.fasta --outreads illumina_modified_1.fastq.gz --outreads2 illumina_modified_2.fastq.gz && ~/src/simsplice/sort-fastq.sh illumina_modified_1.fastq.gz && ~/src/simsplice/sort-fastq.sh illumina_modified_2.fastq.gz'

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && ~/src/simsplice/target/debug/simsplice --bam "$(basename "$0")_ONT.bam" --genome assembly.fasta --vcf modifications.vcf --outgenome modified.fasta --outreads ONT_modified.fastq.gz && ~/src/simsplice/sort-fastq.sh ONT_modified.fastq.gz'

find -name 'ONT_modified.fastq.gz' -o -name 'illumina_modified_1.fastq.gz' -o -name 'illumina_modified_2.fastq.gz' |parallel -q bash -xc '~/src/simsplice/sort-hash-fastq.sh "${0%.fastq.gz}".{fastq.gz,hashed.fastq.gz}'

find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 --eqx -aYx sr -t "$(nproc)" assembly.fasta <(gzip -cd illumina_modified_1.fastq.gz) <(gzip -cd illumina_modified_2.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o illumina_modified.bam && samtools index illumina_modified.bam'

find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 --eqx -aYx map-ont -t "$(nproc)" assembly.fasta <(gzip -cd ONT_modified.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o ONT_modified.bam && samtools index ONT_modified.bam'