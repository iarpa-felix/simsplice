#!/usr/bin/env bash
set -eu -o pipefail

mkdir -p /data/iarpa/analysis/PRJNA607948
cd /data/iarpa/analysis/PRJNA607948

#fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip --origfmt

@sra-toolkit bash -x <<'EOS'
esearch -db sra -query PRJNA607948 |efetch -format native |jx -xPo- -e '_.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE' -e 'a=_.RUN_SET.RUN.$accession; t=_.EXPERIMENT.TITLE.replace(/[ ]/g,"_"); o=t.replace(/_[^_]+$/,""); t.match(/SISPA/)? [] : `fastq-dump --outdir ${q(o)} --gzip --skip-technical --split-3 --clip --dumpbase --read-filter pass ${q(a)} && (cd ${q(o)} && t=${q(t)} rename -f '\''s/^/$ENV{t}./'\'' ${q(a)}*.fastq.gz)`' |parallel -q bash -xc
EOS

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && ln -sfv *_pass_1.fastq.gz illumina_original_1.fastq.gz && ln -sfv *_pass_2.fastq.gz illumina_original_2.fastq.gz && ln -sfv *ONT*_pass.fastq.gz ONT_original.fastq.gz'

find * -maxdepth 0 -type d |parallel -q @unicycler bash -xc 'unicycler -1 "$0"/illumina_original_1.fastq.gz -2 "$0"/illumina_original_2.fastq.gz -l "$0"/ONT_original.fastq.gz -o "$0" -t "$(nproc)"'

find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 -aYx sr -t "$(nproc)" assembly.fasta <(gzip -cd illumina_original_1.fastq.gz) <(gzip -cd illumina_original_2.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o illumina_original.bam && samtools index illumina_original.bam'

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && bedtools bamtobed -i illumina_original.bam -bed12 >illumina_original.bed && bedToBigBed illumina_original.bed assembly.sizes illumina_original.bb'

find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 -aYx map-ont -t "$(nproc)" assembly.fasta <(gzip -cd ONT_original.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o ONT_original.bam && samtools index ONT_original.bam'

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && bedtools bamtobed -i ONT_original.bam -bed12 >ONT_original.bed && bedToBigBed ONT_original.bed assembly.sizes ONT_original.bb'

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && ~/src/simsplice/target/release/genvcf --genome assembly.fasta --vcf modifications.vcf'

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && ~/src/simsplice/target/release/simsplice --bam illumina_original.bam --genome assembly.fasta --vcf modifications.vcf --outgenome modified.fasta --outreads illumina_modified_1.fastq.gz --outreads2 illumina_modified_2.fastq.gz && ~/src/simsplice/sort-fastq.sh illumina_modified_1.fastq.gz && ~/src/simsplice/sort-fastq.sh illumina_modified_2.fastq.gz'

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && ~/src/simsplice/target/release/simsplice --bam ONT_original.bam --genome assembly.fasta --vcf modifications.vcf --outgenome modified.fasta --outreads ONT_modified.fastq.gz && ~/src/simsplice/sort-fastq.sh ONT_modified.fastq.gz'

find -name 'ONT_modified.fastq.gz' -o -name 'illumina_modified_1.fastq.gz' -o -name 'illumina_modified_2.fastq.gz' |parallel -q bash -xc '~/src/simsplice/sort-hash-fastq.sh "${0%.fastq.gz}".{fastq.gz,hashed.fastq.gz}'

find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 -aYx sr -t "$(nproc)" modified.fasta <(gzip -cd illumina_modified_1.fastq.gz) <(gzip -cd illumina_modified_2.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o illumina_modified.bam && samtools index illumina_modified.bam'

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && bedtools bamtobed -i illumina_modified.bam -bed12 >illumina_modified.bed && bedToBigBed illumina_modified.bed modified.sizes illumina_modified.bb'

find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 -aYx map-ont -t "$(nproc)" modified.fasta <(gzip -cd ONT_modified.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o ONT_modified.bam && samtools index ONT_modified.bam'

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && bedtools bamtobed -i ONT_modified.bam -bed12 >ONT_modified.bed && bedToBigBed ONT_modified.bed modified.sizes ONT_modified.bb'

find -name '*.bam' -not -name '*.collated.bam' |parallel -q bash -xc 'bam2bedgraph --bigwig --out "${0%.bam}" "$0"'

find -name '*.fasta' |parallel -q bash -xc '~/src/scripts/fasta_sizes.sh "$0" >"${0%.fasta}.sizes" && faToTwoBit "$0" "${0%.fasta}.2bit"'

find -name '*.vcf' |parallel -q bash -c '
set -eux -o pipefail
~/src/simsplice/comparative/vcf2chainbed.sh "$0" |LC_COLLATE=C sort -k1,1 -k2,2n >"${0%.vcf}.orig2mod.bed" 
bedToBigBed -type=bed3+9 -as="$HOME/src/simsplice/comparative/chain.as" "${0%.vcf}.orig2mod.bed" "$(dirname "$0")/modified.sizes" "${0%.vcf}.orig2mod.bb"
~/src/simsplice/comparative/invert-chainbed.sh "${0%.vcf}.orig2mod.bed" |LC_COLLATE=C sort -k1,1 -k2,2n >"${0%.vcf}.mod2orig.bed" 
bedToBigBed -type=bed3+9 -as="$HOME/src/simsplice/comparative/chain.as" "${0%.vcf}.mod2orig.bed" "$(dirname "$0")/assembly.sizes" "${0%.vcf}.mod2orig.bb"
~/src/simsplice/comparative/vcf2bed.sh "$0" > >(LC_COLLATE=C sort -k1,1 -k2,2n >"${0%.vcf}.qvcf.bed") 2> >(LC_COLLATE=C sort -k1,1 -k2,2n >"${0%.vcf}.rvcf.bed") 
bedToBigBed "${0%.vcf}.qvcf.bed" "$(dirname "$0")/assembly.sizes" "${0%.vcf}.qvcf.bb"
bedToBigBed "${0%.vcf}.rvcf.bed" "$(dirname "$0")/modified.sizes" "${0%.vcf}.rvcf.bb"
'

find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && ln -sfv ~/src/simsplice/comparative/*.{html,js,css} .'
