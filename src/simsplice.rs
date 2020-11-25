use structopt::StructOpt;

use std::vec::Vec;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::io::Write;
use std::io;

use rust_htslib::bam::record::Record;
use rust_htslib::{bam, bam::Read as BamRead};
use rust_htslib::{bcf, bcf::Read as BcfRead};
use rust_htslib::bam::ext::BamRecordExtensions;
use flate2::write::GzEncoder;
use flate2::Compression;
use bigtools::bigwigread::BigWigRead;

use bio::data_structures::interval_tree::IntervalTree;
use bio::alphabets::dna::revcomp;

use bio::io::fasta;
use bio::io::fastq;

use slog::{info, warn};
use slog::Drain;
use slog_term;

use rand::Rng;
use rand::seq::SliceRandom;
use std::str::{from_utf8 as utf8};

// INPUTS: reference fasta file, input vcf file, input bam file
// OUTPUTS: modified reference fasta file, fastq file with modified reads

// Case 1: insertion is smaller than original
// NOTE: If the start of the read is inside the spliced out region, throw out the read

// Case 2: insertion larger than original
// NOTE: Newly inserted to-the-right region bases get randomly assigned read start site values from -1000 .. insertion ... 1000 bp

// create histogram of start sites for filled in region
// fill in with

// 1. read gets passed through unchanged
// 2. read gets removed
// 3. read sequence gets changed:
//   a. read is moved due to previous splice(s)
//   b. read is moved to fill in start site histogram for insertion

// example data:
//   /data/iarpa/analysis/SRR11140744
//   /data/iarpa/TE/437/{assembly,illumina}


// TODO:
// - left-most is always 5prime
// - sample read start sites, not read depths
// - position moved reads by start site

use eyre::{bail, eyre};
use eyre::Result;
trait ToResult<T> {
    fn r(self) -> Result<T>;
}
impl<T> ToResult<T> for Option<T> {
    fn r(self) -> Result<T> {
        self.ok_or_else(|| eyre!("NoneError"))
    }
}

#[derive(StructOpt, Debug)]
#[structopt(name = "simsplice", about = "Apply simulated modifications to a reference genome and a set of aligned reads")]
struct Options {
    #[structopt(short="g", long="genome", help = "Input reference genome FASTA file", name="FASTAFILE")]
    reference: String,
    #[structopt(short="v", long="vcf", help = "Input VCF file with mutations", name="VCFFILE")]
    vcffile: String,
    #[structopt(short="b", long="collated-bam", help = "Input BAM file of reads, *INCLUDING* unaligned reads, collated using samtools collate", name="BAMFILE")]
    collated_bamfile: String,
    #[structopt(short="w", long="start-bigwig", help = "bigWig file of read start sites, e.g. min genomic position, NOT strand specific. Use bedtools genomecov -5 on unstranded input (all +)", name="BIGWIG")]
    start_bigwig: String,
    #[structopt(short="o", long = "outgenome", help = "Output modified reference FASTA file", name="OUTFASTAFILE")]
    outfastafile: String,
    #[structopt(short="r", long = "outreads", help = "Output modified FASTQ reads file", name="OUTFASTQFILE")]
    outfastqfile: String,
    #[structopt(short="R", long = "outreads2", help = "Output modified FASTQ reads file 2, for paired-end reads", name="OUTFASTQFILE2", default_value="")]
    outfastqfile2: String,
}

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone)]
struct Splice {
    start: i64,
    stop: i64,
    replacement: String,
}

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone)]
struct Replacement {
    origpos: i64,
    modpos: i64,
    replacement: String,
}

#[derive(Debug)]
struct ReadPair {
    collated_bamfile: String,
    collated_bam: bam::Reader,
    is_paired: bool,
    record: bam::Record,
    read_qname: Vec<u8>,
    eof: bool,
}

impl ReadPair {
    fn from_bamfile(collated_bamfile: &str) -> Result<ReadPair> {
        let mut collated_bam = bam::Reader::from_path(&collated_bamfile)?;
        let mut record = bam::Record::new();
        collated_bam.read(&mut record)?;
        let read_qname = Vec::from(record.qname());
        Ok(ReadPair {
            collated_bamfile: collated_bamfile.to_string(),
            collated_bam,
            is_paired: false,
            record,
            read_qname,
            eof: false,
        })
    }
    fn header(&mut self) -> bam::HeaderView {
        self.collated_bam.header().clone()
    }
    fn read_pair(&mut self) -> Result<[Option<bam::Record>; 2]> {
        let mut read1 = Option::<bam::Record>::None;
        let mut read2 = Option::<bam::Record>::None;
        while !self.eof {
            if !self.record.is_last_in_template() {
                read1 = match &read1 {
                    Some(_) => if !self.record.is_secondary() { Some(self.record.clone()) } else { read1 },
                    None => Some(self.record.clone()),
                }
            } else {
                read2 = match &read2 {
                    Some(_) => if !self.record.is_secondary() { Some(self.record.clone()) } else { read2 },
                    None => Some(self.record.clone()),
                };
                self.is_paired = true;
            }
            if !self.collated_bam.read(&mut self.record)? {
                self.eof = true;
                break
            }
            if self.record.qname() != self.read_qname.as_slice() {
                self.read_qname = Vec::from(self.record.qname());
                break
            }
        }
        if read1.is_some() || read2.is_some() {
            if self.is_paired && read2.is_none() {
                bail!("Read 2 not found for paired-end BAM file {}: read={:#?}", self.collated_bamfile, utf8(read1.as_ref().unwrap().qname())?)
            }
            if read1.is_none() && !read2.is_none() {
                bail!("No read 1 found for corresponding read 2 in paired-end BAM file {}: read={:#?}", self.collated_bamfile, utf8(read2.as_ref().unwrap().qname())?)
            }
            if self.is_paired && (read1.is_none() || read2.is_none()) {
                let r1name = if let Some(r)=&read1 {String::from(utf8(r.qname())?)} else {"None".to_string()};
                let r2name = if let Some(r)=&read2 {String::from(utf8(r.qname())?)} else {"None".to_string()};
                bail!("Expected paired-end reads, but only one read found: read1={:#?}, read2={:#?}",
        &r1name, &r2name);
            }
        }
        Ok([read1, read2])
    }
}

fn write_fastq_records(
    read1: &fastq::Record,
    read2: Option<&fastq::Record>,
    out: &mut fastq::Writer<Box<dyn Write>>,
    out2: &mut Option<fastq::Writer<Box<dyn Write>>>,
) -> Result<()>
{
    if let Some(read2) = read2 {
        if read1.id() != read2.id() {
            bail!("Read 1 name {} does not match read 2 name {}", read1.id(), read2.id());
        }
    }
    out.write(
        read1.id(),
        None,
        read1.seq(),
        read1.qual(),
    )?;
    if let Some(out2) = out2 {
        if let Some(read2) = read2 {
            out2.write(
                read2.id(),
                None,
                &read2.seq(),
                &read2.qual(),
            )?;
        }
        else {
            bail!("Read 2 was not given for read {}", read1.id());
        }
    }
    Ok(())
}

fn fillin_aligned_pairs(record: &bam::Record, modpos: i64, origseq: &[u8], modseq: &[u8]) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let ap = record.aligned_pairs();
    let oldseq = record.seq().as_bytes();
    let mut seq = oldseq.clone();
    for a in &ap {
        let qpos = a[0];
        let rpos = a[1];
        let newrpos = rpos-record.pos()+modpos;
        if 0 <= newrpos && newrpos < modseq.len() as i64 {
            // old match
            if oldseq[qpos as usize].to_ascii_uppercase() == origseq[rpos as usize].to_ascii_uppercase() {
                seq[qpos as usize] = modseq[newrpos as usize];
            }
            // old mismatch, new match
            else if oldseq[qpos as usize].to_ascii_uppercase() == modseq[newrpos as usize].to_ascii_uppercase() {
                let nucs = [b'A', b'C', b'G', b'T'].iter().filter(|n| **n != modseq[newrpos as usize].to_ascii_uppercase()).collect::<Vec<_>>();
                let randnuc = rng.gen_range(0, nucs.len());
                seq[qpos as usize] = *nucs[randnuc];
            }
            // match case
            seq[qpos as usize] = if oldseq[qpos as usize].is_ascii_uppercase()
                { seq[qpos as usize].to_ascii_uppercase() }
            else if oldseq[qpos as usize].is_ascii_lowercase()
                { seq[qpos as usize].to_ascii_lowercase() }
            else { seq[qpos as usize] };
        }
    }
    seq
}

fn main() -> Result<()> {
    std::env::set_var("RUST_BACKTRACE", "1");
    let options: Options = Options::from_args();

    let plain = slog_term::PlainSyncDecorator::new(std::io::stderr());
    let log = slog::Logger::root(
        slog_term::FullFormat::new(plain).build().fuse(),
        slog::o!()
    );

    let mut rng = rand::thread_rng();

    // read in the reference FASTA file
    info!(log, "Reading in FASTA assembly file {}", &options.reference);
    let mut origrefseqs = BTreeMap::<String,Vec<u8>>::new();
    let fasta = fasta::Reader::from_file(&options.reference)?;
    for r in fasta.records() {
        let r = r?;
        origrefseqs.insert(String::from(r.id()), Vec::from(r.seq()));
    }

    // get the list of splice positions from the VCF file
    let mut splices = BTreeMap::<String,BTreeSet<Splice>>::new();
    for refname in origrefseqs.keys() {
        if !splices.contains_key(refname) {
            splices.insert(String::from(refname), BTreeSet::new());
        }
    }
    info!(log, "Reading VCF file {}", &options.vcffile);
    let mut vcf = bcf::Reader::from_path(&options.vcffile)?;
    for r in vcf.records() {
        let r = r?;
        info!(log, "Processing VCF record: {:#?}", &r);
        if r.allele_count() > 2 {
            bail!("Can't handle multiple alleles in VCF file: {:#?}", r);
        }
        if r.allele_count() > 0 {
            let rid = r.rid().r()?;
            let refname = utf8(r.header().rid2name(rid)?)?;
            let pos = r.pos(); // 0-based
            let alleles = r.alleles();
            let ref_allele = alleles[0];
            for allele in &alleles[1..] {
                let splice = Splice{
                    start: pos,
                    stop: pos+ref_allele.len() as i64,
                    replacement: String::from(utf8(allele)?),
                };
                info!(log, "Storing splice: {:#?}", splice);
                splices.get_mut(refname).r()?.insert(splice);
            }
        }
    }

    {
        // Open output FASTA and FASTQ
        info!(log, "Writing output FASTA: {}", &options.outfastafile);
        let mut outfasta = bio::io::fasta::Writer::to_file(&options.outfastafile)?;

        info!(log, "Writing output FASTQ1: {}", &options.outfastqfile);
        let mut outfastq: bio::io::fastq::Writer<Box<dyn Write>> = if options.outfastqfile.ends_with(".gz") {
            bio::io::fastq::Writer::new(Box::new(GzEncoder::new(std::fs::File::create(&options.outfastqfile)?, Compression::default())))
        } else {
            bio::io::fastq::Writer::new(Box::new(io::BufWriter::new(std::fs::File::create(&options.outfastqfile)?)))
        };

        info!(log, "Writing output FASTQ2: {}", &options.outfastqfile2);
        let mut outfastq2: Option<fastq::Writer<Box<dyn Write>>> = if options.outfastqfile2.is_empty() {
            None
        } else if options.outfastqfile2.ends_with(".gz") {
            Some(bio::io::fastq::Writer::new(Box::new(GzEncoder::new(std::fs::File::create(&options.outfastqfile2)?, Compression::default()))))
        } else {
            Some(bio::io::fastq::Writer::new(Box::new(io::BufWriter::new(std::fs::File::create(&options.outfastqfile2)?))))
        };

        let sample_distance = 1000;

        let mut start_bigwig = BigWigRead::from_file_and_attach(&options.start_bigwig).
            map_err(|e| eyre!("Error Loading bigWig file {}: {:#?}", options.start_bigwig, e))?;

        let mut read_pair = ReadPair::from_bamfile(&options.collated_bamfile)?;
        let header = read_pair.header();

        // Build an interval tree using the splice positions to map between original and modified reference coordinates.
        // Also build and write the output reference fasta file.
        let mut modrefseqs = BTreeMap::<String, Vec<u8>>::new();
        let mut tree = BTreeMap::<String, IntervalTree<i64, Replacement>>::new();
        // iterate through each reference name
        for (refname, splices) in &mut splices {
            info!(log, "Building interval tree and writing fasta record for refname {}", refname);
            // add the arm to the tree
            tree.insert(String::from(refname), IntervalTree::new());
            let mut modpos = 0; // position in modified genome
            let mut origpos = 0; // position in original genome
            let mut modseq = Vec::<u8>::new(); // sequence of modified genome
            let origlen = origrefseqs[refname].len() as i64; // reference length in original genome
            // iterate through the splice records for each arm
            for splice in splices.iter() {
                info!(log, "Processing splice: {:#?}", &splice);
                // pre-replacement
                if origpos < splice.start {
                    let replacement_seq = &origrefseqs[refname][origpos as usize..splice.start as usize];
                    let replacement = Replacement {
                        origpos,
                        modpos,
                        replacement: String::from(utf8(replacement_seq)?),
                    };
                    info!(log, "Pre-Replacement at {}:{}..{}: {:#?}", refname, origpos, splice.start, &replacement);
                    tree.get_mut(refname).r()?.insert(origpos..splice.start, replacement);
                    modpos += splice.start - origpos;
                    origpos = splice.start;
                    info!(log, "modpos={}, origpos={}", modpos, origpos);
                    modseq.extend(replacement_seq);
                    info!(log, "modseq.len()={}", modseq.len());
                }
                // replacement
                let replacement = Replacement {
                    origpos,
                    modpos,
                    replacement: splice.replacement.clone(),
                };
                info!(log, "Replacement at {}:{}..{}: {:#?}", refname, splice.start, splice.stop, &replacement);
                tree.get_mut(refname).r()?.insert(splice.start..splice.stop, replacement);
                modpos += splice.replacement.len() as i64;
                origpos = splice.stop;
                info!(log, "modpos={}, origpos={}", modpos, origpos);
                modseq.extend(&mut splice.replacement.as_str().as_bytes().iter());
                info!(log, "modseq.len()={}", modseq.len());
            }
            // post-replacement
            if origpos < origlen {
                let replacement_seq = &origrefseqs[refname][origpos as usize..origlen as usize];
                let replacement = Replacement {
                    origpos,
                    modpos,
                    replacement: String::from(utf8(replacement_seq)?),
                };
                info!(log, "Post-Replacement at {}:{}..{}: {:#?}", refname, origpos, origlen, &replacement);
                tree.get_mut(refname).r()?.insert(origpos..origlen, replacement);
                modseq.extend(replacement_seq);
                info!(log, "modseq.len()={}", modseq.len());
            }
            info!(log, "Writing refname {}, {} bases", refname, modseq.len());
            outfasta.write(refname, None, &modseq)?;
            modrefseqs.insert(refname.clone(), modseq);
        }

        // fill in newly spliced in genomic regions using random samples of existing read pairs
        let mut record_buffer = Vec::<[Option<Record>; 2]>::new();
        // for each refname
        for (refname, splices) in &mut splices {
            info!(log, "Filling in regions and writing FASTQ records for refname {}", refname);
            let mut modpos = 0; // position in modified genome
            let mut origpos = 0; // position in original genome
            let origlen = origrefseqs[refname].len() as i64; // reference length in original genome

            // iterate through each splice record
            for splice in splices.iter() {
                if origpos < splice.start {
                    modpos += splice.start - origpos;
                    origpos = splice.start;
                }
                if splice.stop - splice.start < splice.replacement.len() as i64 {
                    let tid = header.tid(refname.as_bytes()).r()?;
                    info!(log, "Filling in region: {}:{}..{}", refname, splice.start, splice.stop);

                    origpos = splice.stop;
                    modpos += splice.stop-splice.start;

                    // fill up a sample histogram using -1000..1000 bp window around the area of insertion
                    // from the original genome
                    let sample_region_beg = std::cmp::max(0, origpos - sample_distance) as u64;
                    let sample_region_end = std::cmp::min(origlen, origpos + sample_distance) as u64;
                    info!(log, "Building sample region: {}:{}..{}", refname, sample_region_beg, sample_region_end);
                    let mut sample_region_histo = vec![0u32; (sample_region_end - sample_region_beg) as usize];

                    for value in start_bigwig.get_interval(
                        refname,
                        sample_region_beg as u32,
                        sample_region_end as u32)?
                    {
                        let value = value?;
                        let start = value.start as u64 - sample_region_beg;
                        let end = value.end as u64 - sample_region_beg;
                        for i in std::cmp::max(0, start)..std::cmp::min(sample_region_histo.len() as u64, end)
                        {
                            sample_region_histo[i as usize] = value.value as u32;
                        }
                    }

                    // we need to add reads to the new region
                    let fill_region_len = splice.replacement.len() as i64 - (splice.stop - splice.start);
                    let fill_region_pos = modpos + (splice.replacement.len() as i64) - fill_region_len;
                    let mut fill_region_histo = vec![0u64; fill_region_len as usize];
                    'FILL_REGION:
                    for fr_i in 0..fill_region_len {
                        // find a random depth value from the sample region
                        let sr_i = rng.gen_range(0, sample_region_histo.len());
                        let sr_depth = &sample_region_histo[sr_i];

                        // now fill the current fill_region base up to at least value_at_sr_i
                        'FILL_REGION_HISTO:
                        while fill_region_histo[fr_i as usize] < *sr_depth as u64 {
                            let reads = read_pair.read_pair()?;
                            if reads[0].is_none() && reads[1].is_none() { break 'FILL_REGION }

                            let mut read_nums = (0..reads.len()).collect::<Vec<_>>();
                            read_nums.shuffle(&mut rng);

                            let mut fastq_records: [Option<fastq::Record>; 2] = [None, None];

                            // skip read pairs not on the same refname or where after translating the positions
                            // either read would go past the edge of genome space
                            let mut read_pos = -1;
                            let mut read_tid = -1;
                            for r in read_nums {
                                let record = &reads[r];
                                if let Some(record) = record {

                                    if !record.is_unmapped() {
                                        if read_tid < 0 { read_tid = record.tid() }
                                        if read_pos < 0 { read_pos = record.pos() }

                                        // if mapped read is on the same refname as the longest block's read
                                        // it gets translated to the modified genome without
                                        // needing the tree lookup
                                        if record.tid() == read_tid {
                                            let refname = utf8(header.tid2name(record.tid() as u32))?;
                                            let origseq = origrefseqs.get(refname).r()?;
                                            let modrefname = utf8(header.tid2name(tid))?;
                                            let modseq = modrefseqs.get(modrefname).r()?;

                                            let modrecstart = fill_region_pos+fr_i+(record.pos()-read_pos);
                                            let modrecend = modrecstart+(record.reference_end()-record.pos());
                                            if modrecstart < 0 || modrecend > modseq.len() as i64 {
                                                record_buffer.push(reads);
                                                continue 'FILL_REGION_HISTO
                                            }
                                            let seq = fillin_aligned_pairs(record, modrecstart, origseq, modseq);
                                            let (seq, qual) = if record.is_reverse() {
                                                (revcomp(seq),
                                                 record.qual().iter().rev().map(|q| q+33).collect::<Vec<u8>>())
                                            } else {
                                                (seq,
                                                 record.qual().iter().map(|q| q+33).collect::<Vec<u8>>())
                                            };
                                            let fastq_record = fastq::Record::with_attrs(
                                                utf8(record.qname())?,
                                                None,
                                                &seq,
                                                &qual,
                                            );
                                            fastq_records[r] = Some(fastq_record);
                                        }
                                        // if read is on a different strand than the longest block's read,
                                        // then pass it through the interval tree
                                        else {
                                            // make sure the other read's new 5' end maps to the new
                                            let mut found_entry = false;
                                            let refname = utf8(header.tid2name(record.tid() as u32))?;
                                            for entry in tree.get(refname).r()?.find(record.pos()..record.pos() + 1) {
                                                let replacement = entry.data();
                                                let origseq = origrefseqs.get(refname).r()?;
                                                let modseq = modrefseqs.get(refname).r()?;

                                                let modpos = record.pos()-replacement.origpos+replacement.modpos;
                                                let seq = fillin_aligned_pairs(record, modpos, origseq, modseq);
                                                let (seq, qual) = if record.is_reverse() {
                                                    (revcomp(seq),
                                                     record.qual().iter().rev().map(|q| q+33).collect::<Vec<u8>>())
                                                } else {
                                                    (seq,
                                                     record.qual().iter().map(|q| q+33).collect::<Vec<u8>>())
                                                };
                                                let fastq_record = fastq::Record::with_attrs(
                                                    utf8(record.qname())?,
                                                    None,
                                                    &seq,
                                                    &qual,
                                                );
                                                fastq_records[r] = Some(fastq_record);
                                                found_entry = true;
                                                break;
                                            }
                                            if !found_entry {
                                                //info!(log, "Dropped record {}", utf8(record.qname())?);
                                                continue 'FILL_REGION_HISTO;
                                            }
                                        }
                                    }
                                    // unmapped reads get passed through unchanged
                                    else {
                                        let (seq, qual) = if record.is_reverse() {
                                            (revcomp(record.seq().as_bytes()),
                                             record.qual().iter().rev().map(|q| q+33).collect::<Vec<u8>>())
                                        } else {
                                            (record.seq().as_bytes(),
                                             record.qual().iter().map(|q| q+33).collect::<Vec<u8>>())
                                        };
                                        let fastq_record = fastq::Record::with_attrs(
                                            utf8(record.qname())?,
                                            None,
                                            &seq,
                                            &qual,
                                        );
                                        fastq_records[r] = Some(fastq_record);
                                    }
                                }
                            }

                            if let Some(read1) = &fastq_records[0] {
                                if reads[1].is_none() || fastq_records[1].is_some() {
                                    // fill in histogram
                                    for record in reads.iter() {
                                        if let Some(record) = record {
                                            if !record.is_unmapped() &&
                                                read_tid >= 0 &&
                                                record.tid() == read_tid &&
                                                read_pos >= 0
                                            {
                                                let pos = record.pos()-read_pos+fr_i;
                                                if 0 <= pos && pos < fill_region_len {
                                                    fill_region_histo[pos as usize] += 1
                                                }
                                            }
                                        }
                                    }
                                    write_fastq_records(
                                        &read1,
                                        fastq_records[1].as_ref(),
                                        &mut outfastq,
                                        &mut outfastq2,
                                    )?;
                                }
                            }
                        }
                    }
                }
                modpos += splice.replacement.len() as i64-(splice.stop-splice.start);
                origpos = splice.stop;
            }
        }

        // translate the coordinates and write out the rest of the reads
        info!(log, "Translating and writing the rest of the FASTQ reads");
        'READ_PAIR:
        loop {
            let reads = if !record_buffer.is_empty() {
                record_buffer.pop().r()?
            } else {
                read_pair.read_pair()?
            };
            if reads[0].is_none() && reads[1].is_none() && record_buffer.is_empty() {
                break 'READ_PAIR;
            }

            let mut fastq_records: [Option<fastq::Record>; 2] = [None, None];
            for (r, record) in reads.iter().enumerate() {
                if let Some(record) = record {
                    if !record.is_unmapped() {
                        let record_pos = record.pos();
                        let refname = utf8(header.target_names().get(record.tid() as usize).r()?)?;
                        let origseq = origrefseqs.get(refname).r()?;
                        let modseq = modrefseqs.get(refname).r()?;
                        let mut found_entry = false;
                        for entry in tree.get(refname).r()?.find(record_pos..record_pos + 1) {
                            let replacement = entry.data();
                            let modpos = record.pos()-replacement.origpos+replacement.modpos;
                            let seq = fillin_aligned_pairs(record, modpos, origseq, modseq);
                            let (seq, qual) = if record.is_reverse() {
                                (revcomp(seq),
                                    record.qual().iter().rev().map(|q| q+33).collect::<Vec<u8>>())
                            } else {
                                (seq,
                                    record.qual().iter().map(|q| q+33).collect::<Vec<u8>>())
                            };
                            let fastq_record = fastq::Record::with_attrs(
                                utf8(record.qname())?,
                                None,
                                &seq,
                                &qual,
                            );
                            fastq_records[r] = Some(fastq_record);
                            found_entry = true;
                            break;
                        }
                        if !found_entry {
                            //info!(log, "Dropped record {}", utf8(record.qname())?);
                            continue 'READ_PAIR;
                        }
                    }
                    // unmapped reads get passed through unchanged
                    else {
                        let (seq, qual) = if record.is_reverse() {
                            (revcomp(record.seq().as_bytes()),
                             record.qual().iter().rev().map(|q| q+33).collect::<Vec<u8>>())
                        } else {
                            (record.seq().as_bytes(),
                             record.qual().iter().map(|q| q+33).collect::<Vec<u8>>())
                        };
                        let fastq_record = fastq::Record::with_attrs(
                           utf8(record.qname())?,
                           None,
                           &seq,
                           &qual,
                        );
                       fastq_records[r] = Some(fastq_record);
                    }
                }
            }
            if let Some(read1) = &fastq_records[0] {
                write_fastq_records(
                    &read1,
                    fastq_records[1].as_ref(),
                    &mut outfastq,
                    &mut outfastq2,
                )?;
            }
            else {
                bail!("No read1 found for read!");
            }
        }
    }
    Ok(())
}
