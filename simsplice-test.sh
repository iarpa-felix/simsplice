#!/usr/bin/env bash
set -eu -o pipefail

mkdir -p /data/iarpa/analysis/PRJNA607948
cd /data/iarpa/analysis/PRJNA607948

#fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip --origfmt

# download FASTQ files from SRA
@sra-toolkit bash -x <<'EOS'
esearch -db sra -query PRJNA607948 |efetch -format native |jx -xPo- -e '_.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE' -e 'a=_.RUN_SET.RUN.$accession; t=_.EXPERIMENT.TITLE.replace(/[ ]/g,"_"); o=t.replace(/_[^_]+$/,""); t.match(/SISPA/)? [] : `fastq-dump --outdir ${q(o)} --gzip --skip-technical --split-3 --clip --dumpbase --read-filter pass ${q(a)} && (cd ${q(o)} && t=${q(t)} rename -f '\''s/^/$ENV{t}./'\'' ${q(a)}*.fastq.gz)`' |parallel -q bash -xc
EOS

# make generic symlinks
find * -maxdepth 0 -type d |parallel -q @unicycler bash -xc 'ln -s *_pass_1.fastq.gz illumina_original_1.fastq.gz && ln -s *_pass_2.fastq.gz illumina_original_2.fastq.gz && ln -s *ONT*_pass.fastq.gz ONT_original.fastq.gz'

# unicycler assembly
find * -maxdepth 0 -type d |parallel -q @unicycler bash -xc 'unicycler -1 "$0"/illumina_original_1.fastq.gz -2 "$0"/illumina_original_2.fastq.gz -l "$0"/ONT_original.fastq.gz -o "$0" -t "$(nproc)"'

# run minimap with illumina reads against assembly
find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 -aYx sr -t "$(nproc)" assembly.fasta <(gzip -cd *_illumina.*_1.fastq.gz) <(gzip -cd *_illumina.*_2.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o illumina_original.bam && samtools index illumina_original.bam'

# run minimap with ONT reads against assembly
find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 -aYx map-ont -t "$(nproc)" assembly.fasta <(gzip -cd ONT_original.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o ONT_original.bam && samtools index ONT_original.bam'

# randomly generate VCF files
find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && ~/src/simsplice/target/release/genvcf --genome assembly.fasta --vcf modifications.vcf'

# run simsplice on illumina, sort reads by name
find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && ~/src/simsplice/target/release/simsplice --bam illumina_original.bam" --genome assembly.fasta --vcf modifications.vcf --outgenome modified.fasta --outreads illumina_modified_1.fastq.gz --outreads2 illumina_modified_2.fastq.gz && ~/src/simsplice/sort-fastq.sh illumina_modified_1.fastq.gz && ~/src/simsplice/sort-fastq.sh illumina_modified_2.fastq.gz'

# run simsplice on ONT, sort reads by name
find * -maxdepth 0 -type d |parallel -q bash -xc 'cd "$0" && ~/src/simsplice/target/release/simsplice --bam ONT_original.bam" --genome assembly.fasta --vcf modifications.vcf --outgenome modified.fasta --outreads ONT_modified.fastq.gz && ~/src/simsplice/sort-fastq.sh ONT_modified.fastq.gz'

# hash and sort modified FASTQ read names
find -name 'ONT_modified.fastq.gz' -o -name 'illumina_modified_1.fastq.gz' -o -name 'illumina_modified_2.fastq.gz' |parallel -q bash -xc '~/src/simsplice/sort-hash-fastq.sh "${0%.fastq.gz}".{fastq.gz,hashed.fastq.gz}'

# run minimap with illumina reads against modified assembly
find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 -aYx sr -t "$(nproc)" modified.fasta <(gzip -cd illumina_modified_1.fastq.gz) <(gzip -cd illumina_modified_2.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o illumina_modified.bam && samtools index illumina_modified.bam'

# run minimap with ONT reads against modified assembly
find * -maxdepth 0 -type d |parallel -q @minimap2 bash -xc 'cd "$0" && minimap2 -aYx map-ont -t "$(nproc)" modified.fasta <(gzip -cd ONT_modified.fastq.gz) |samtools view -Sbu - |samtools sort -@ "$(nproc)" - -o ONT_modified.bam && samtools index ONT_modified.bam'

# generate bigWig files
find -name '*.bam' -not -name '*.collated.bam' |parallel -q bash -xc 'bam2bedgraph --bigwig --out "${0%.bam}" "$0"'

# generate assembly FASTA sizes files
find -name '*.fasta' |parallel -q bash -xc '~/src/scripts/fasta_sizes.sh "$0" >"${0%.fasta}.sizes" && faToTwoBit "$0" "${0%.fasta}.2bit"'

# generate bigBed files
find -name '*.bam' -not -name '*.collated.bam' |parallel -q bash -xc 'samtools view -H "$0" |perl -F\\t -lane '\''/^\@SQ/ && print join "\t",$F[1]=~s/^SN://r,$F[2]=~s/^LN://r'\'' >"$0.sizes" && bedtools bamtobed -i "$0" -bed12 >"${0%.bam}.bed" && bedToBigBed "${0%.bam}.bed" "$0.sizes" "${0%.bam}.bb"'

# generate browser files
find -name '*.vcf' |parallel -q bash -c '
set -eux -o pipefail
#
# generate chainbed bb file
~/src/simsplice/comparative/vcf2chainbed.sh "$0" > >(LC_COLLATE=C sort -k1,1 -k2,2n >"${0%.vcf}.orig2mod.bed") 2> >(LC_COLLATE=C sort -k1,1 -k2,2n >"${0%.vcf}.mod2orig.bed") 
bedToBigBed -type=bed3+9 -as="$HOME/src/simsplice/comparative/chain.as" "${0%.vcf}.orig2mod.bed" "$(dirname "$0")/modified.sizes" "${0%.vcf}.orig2mod.bb"
bedToBigBed -type=bed3+9 -as="$HOME/src/simsplice/comparative/chain.as" "${0%.vcf}.mod2orig.bed" "$(dirname "$0")/assembly.sizes" "${0%.vcf}.mod2orig.bb"
#
# generate VCF-only chainbed bb file, for comparison track
~/src/simsplice/comparative/vcf2chainbed-vcfonly.sh "$0" > >(LC_COLLATE=C sort -k1,1 -k2,2n >"${0%.vcf}.orig2mod.vcfonly.bed") 2> >(LC_COLLATE=C sort -k1,1 -k2,2n >"${0%.vcf}.mod2orig.vcfonly.bed") 
bedToBigBed -type=bed3+9 -as="$HOME/src/simsplice/comparative/chain.as" "${0%.vcf}.orig2mod.vcfonly.bed" "$(dirname "$0")/modified.sizes" "${0%.vcf}.orig2mod.vcfonly.bb"
bedToBigBed -type=bed3+9 -as="$HOME/src/simsplice/comparative/chain.as" "${0%.vcf}.mod2orig.vcfonly.bed" "$(dirname "$0")/assembly.sizes" "${0%.vcf}.mod2orig.vcfonly.bb"
#
# generate VCF bb file, to view VCF as a feature
~/src/simsplice/comparative/vcf2bed.sh "$0" > >(LC_COLLATE=C sort -k1,1 -k2,2n >"${0%.vcf}.qvcf.bed") 2> >(LC_COLLATE=C sort -k1,1 -k2,2n >"${0%.vcf}.rvcf.bed") 
bedToBigBed "${0%.vcf}.qvcf.bed" "$(dirname "$0")/assembly.sizes" "${0%.vcf}.qvcf.bb"
bedToBigBed "${0%.vcf}.rvcf.bed" "$(dirname "$0")/modified.sizes" "${0%.vcf}.rvcf.bb"
'
