use structopt::StructOpt;

use std::str;
use std::vec::Vec;
use std::path::Path;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::io::Write;
use std::io;

use rust_htslib::bam::record::Record;
use rust_htslib::bam::record::Cigar;
use rust_htslib::{bam, bam::Read as BamRead};
use rust_htslib::{bcf, bcf::Read as BcfRead};
use rust_htslib::bam::ext::BamRecordExtensions;
use flate2::write::GzEncoder;
use flate2::Compression;

use bio::data_structures::interval_tree::IntervalTree;

use anyhow::anyhow;

use bio::io::fasta;
use bio::io::fastq;
use duct::cmd;
use interpol::{format as f};
use num_cpus;
use scopeguard::defer;
use regex::Regex;
use lazy_static::lazy_static;
use lazy_regex::{regex as re};

use slog::info;
use sloggers::Build;
use sloggers::terminal::{TerminalLoggerBuilder, Destination};
use sloggers::types::Severity;

use rand::Rng;
use rand::distributions::Alphanumeric;
use shell_words;
use shell_words::quote;

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

pub type Result<T, E = anyhow::Error> = core::result::Result<T, E>;
trait ToResult<T> {
    fn r(self) -> Result<T>;
}
impl<T> ToResult<T> for Option<T> {
    fn r(self) -> Result<T> {
        self.ok_or_else(|| anyhow!("NoneError"))
    }
}

#[derive(StructOpt, Debug)]
#[structopt(name = "bam2bedgraph", about = "Convert bam files to bedgraph/bigWig format")]
struct Options {
    #[structopt(help = "Input reference FASTA file", name="FASTAFILE")]
    reference: String,
    #[structopt(help = "Input VCF file", name="VCFFILE")]
    vcffile: String,
    #[structopt(help = "Input BAM file", name="BAMFILE")]
    bamfile: String,
    #[structopt(long = "outfasta", help = "Output reference FASTA file", name="OUTFASTAFILE")]
    outfastafile: String,
    #[structopt(long = "outfastq", help = "Output FASTQ file", name="OUTFASTQFILE")]
    outfastqfile: String,
    #[structopt(long = "outfastq2", help = "Output FASTQ file, read 2", name="OUTFASTQFILE2", default_value="")]
    outfastqfile2: String,
}

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone)]
struct Splice {
    start: i64,
    stop: i64,
    replacement: Vec<u8>,
}

struct Replacement<'a> {
    offset: i64,
    replacement: &'a[u8],
}

//437

fn main() -> Result<()> {
    std::env::set_var("RUST_BACKTRACE", "full");
    std::env::set_var("RUST_LIB_BACKTRACE", "1");
    let options: Options = Options::from_args();

    let mut builder = TerminalLoggerBuilder::new();
    builder.level(Severity::Info);
    builder.destination(Destination::Stderr);

    let log = builder.build()?;

    // read in the reference FASTA file
    let mut reference = BTreeMap::<String,Vec<u8>>::new();
    let fasta = fasta::Reader::from_file(&options.reference)?;
    for r in fasta.records() {
        let r = r?;
        reference.insert(String::from(r.id()), Vec::from(r.seq()));
    }

    // make a bam index file if it does not exist
    let index_file = f!("{options.bamfile}.bai");
    if !Path::new(&index_file).exists() {
        let index_cmd = vec![ "samtools", "index", &options.bamfile ];
        info!(log, "Running command: {}", shell_words::join(&index_cmd));
        cmd(index_cmd[0],&index_cmd[1..]).run()?;
    }
    let mut indexed_bam = bam::IndexedReader::from_path_and_index(
        &options.bamfile,
    &index_file)?;

    // collate the bam file by name
    let collated_bamfile = f!(r#"{re!(r"\.bam$").replace_all(&options.bamfile,"")}.collated.bam"#);
    let cpus = num_cpus::get().to_string();
    if !Path::new(&collated_bamfile).exists() {
        let tmpid = rand::thread_rng().sample_iter(&Alphanumeric).take(10).collect::<String>();
        let tmpfile = f!("{collated_bamfile}.{tmpid}.bam.tmp");
        let collate_cmd = vec![
            "samtools","collate",
            "-@",&cpus,
            "-o",&tmpfile,
            &options.bamfile,
        ];
        info!(log, "Running command: {}", shell_words::join(&collate_cmd));
        cmd(collate_cmd[0],&collate_cmd[1..]).run()?;
        info!(log, "Moving {} to {}", &tmpfile, &collated_bamfile);
        std::fs::rename(&tmpfile, &collated_bamfile)?;
    }
    let mut collated_bam = bam::Reader::from_path(&collated_bamfile)?;

    // get the list of splice positions from the VCF file
    let mut splices = BTreeMap::<String,BTreeSet<Splice>>::new();
    let mut vcf = bcf::Reader::from_path(&options.vcffile)?;
    for r in vcf.records() {
        let r = r?;
        if r.allele_count() > 2 {
            Err(anyhow!("Can't handle multiple alleles in VCF file: {:?}", r))?;
        }
        if r.allele_count() > 0 {
            let rid = r.rid().r()?;
            let refname = str::from_utf8(r.header().rid2name(rid)?)?;
            // TODO: Is API 1-based or 0-based?
            let pos = r.pos();
            let alleles = r.alleles();
            let ref_allele = alleles.get(0).r()?;
            for allele in &alleles[1..] {
                if !splices.contains_key(refname) {
                    splices.insert(String::from(refname), BTreeSet::new());
                }
                splices.get_mut(refname).r()?.insert(Splice{
                    start: pos,
                    stop: pos+ref_allele.len() as i64,
                    replacement: Vec::from(*allele),
                });
            }
        }
    }

    let mut outfasta = bio::io::fasta::Writer::to_file(&options.outfastafile)?;
    let mut outfastq: bio::io::fastq::Writer<Box<dyn Write>> = if options.outfastqfile.ends_with(".gz") {
        bio::io::fastq::Writer::new(Box::new(GzEncoder::new(std::fs::File::create(&options.outfastqfile)?, Compression::default())))
    } else {
        bio::io::fastq::Writer::new(Box::new(io::BufWriter::new(std::fs::File::create(&options.outfastqfile)?)))
    };

    let mut outfastq2: Option<fastq::Writer<Box<dyn Write>>> = if options.outfastqfile2.is_empty() {
        None
    } else if options.outfastqfile2.ends_with(".gz") {
        Some(bio::io::fastq::Writer::new(Box::new(GzEncoder::new(std::fs::File::create(&options.outfastqfile2)?, Compression::default()))))
    } else {
        Some(bio::io::fastq::Writer::new(Box::new(io::BufWriter::new(std::fs::File::create(&options.outfastqfile2)?))))
    };
    let header = indexed_bam.header().clone();
    let sample_distance = 1000;

    let mut record = bam::Record::new();
    collated_bam.read(&mut record)?;
    let mut read_qname = Vec::from(record.qname());
    let mut read_pair = || -> Result<(Option<bam::Record>, Option<bam::Record>)> {
        let mut read1 = Option::<bam::Record>::None;
        let mut read2 = Option::<bam::Record>::None;
        loop {
            if !read_qname.is_empty() {
                if !record.is_last_in_template() {
                    read1 = match &read1 {
                        Some(_) => if !record.is_secondary() {Some(record.clone())} else {read1},
                        None => Some(record.clone()),
                    }
                }
                else {
                    read2 = match &read2 {
                        Some(_) => if !record.is_secondary() {Some(record.clone())} else {read2},
                        None => Some(record.clone()),
                    }
                }
            }
            if !collated_bam.read(&mut record)? || record.qname() != read_qname.as_slice() {
                read_qname = Vec::from(record.qname());
                break
            }
        }
        if !options.outfastqfile2.is_empty() && read2.is_none() {
            Err(anyhow!("Read 2 not found for paired-end BAM file {}: record={:?}", collated_bamfile, record))?
        }
        if read1.is_none() && !read2.is_none() {
            Err(anyhow!("No read 1 found for corresponding read 2 in paired-end BAM file {}: record={:?}", collated_bamfile, record))?
        }
        Ok((read1, read2))
    };

    // build an interval tree using the splice positions to map between original and modified reference coordinates
    // also write the output reference fasta file
    let mut tree = BTreeMap::<String,IntervalTree<i64,Replacement>>::new();
    for (chr, splices) in &mut splices {
        tree.insert(String::from(chr), IntervalTree::new());
        let mut pos = 0; // position in modified genome
        let mut refpos = 0; // position in original genome
        let mut sequence = Vec::<u8>::new(); // sequence of modified genome
        let reflen = reference.get(chr).r()?.len() as i64; // reference length in original genome
        for splice in splices.iter() {
            if refpos < splice.start {
                let replacement = &reference.get(chr).r()?[refpos as usize..splice.start as usize];
                tree.get_mut(chr).r()?.insert(refpos..splice.start, Replacement {
                    offset: pos-refpos,
                    replacement,
                });
                pos += splice.start-refpos;
                refpos = splice.start;
                sequence.extend(replacement);
            }
            tree.get_mut(chr).r()?.insert(splice.start..splice.stop, Replacement {
                offset: pos-refpos,
                replacement: &splice.replacement,
            });
            pos += splice.replacement.len() as i64;
            refpos = splice.stop;
            sequence.append(&mut splice.replacement.clone());
        }
        if refpos < reflen {
            let replacement = &reference.get(chr).r()?[refpos as usize..reflen as usize];
            tree.get_mut(chr).r()?.insert(refpos..reflen, Replacement {
                offset: pos-refpos,
                replacement,
            });
            sequence.extend(replacement);
        }
        outfasta.write(chr, None, &sequence)?;
    }

    for (chr, splices) in &mut splices {
        tree.insert(String::from(chr), IntervalTree::new());
        let mut pos = 0; // position in modified genome
        let mut refpos = 0; // position in original genome
        let mut sequence = Vec::<u8>::new(); // sequence of modified genome
        let reflen = reference.get(chr).r()?.len() as i64; // reference length in original genome
        for splice in splices.iter() {
            if refpos < splice.start {
                let replacement = &reference.get(chr).r()?[refpos as usize..splice.start as usize];
                tree.get_mut(chr).r()?.insert(refpos..splice.start, Replacement {
                    offset: pos-refpos,
                    replacement,
                });
                pos += splice.start-refpos;
                refpos = splice.start;
                sequence.extend(replacement);
            }
            if splice.stop-splice.start < splice.replacement.len() as i64 {
                // we need to add reads to the new region
                let fill_region_len = splice.replacement.len() as i64 - (splice.stop-splice.start);
                let fill_region_pos = pos+(splice.replacement.len() as i64)-fill_region_len;
                let fill_region_histo = vec![0u64; fill_region_len as usize];
                let tid = header.tid(chr.as_bytes()).r()?;
                let sample_region_beg = std::cmp::max(0, refpos-sample_distance) as u64;
                let sample_region_end = std::cmp::min(reflen, refpos+sample_distance) as u64;
                let mut sample_region_histo = vec![0u32; (sample_region_end-sample_region_beg) as usize];
                indexed_bam.fetch(tid, sample_region_beg, sample_region_end)?;
                for pileup in indexed_bam.pileup() {
                    let pileup = pileup?;
                    // TODO: Is API 1-based or 0-based?
                    let sr_pos = pileup.pos() - sample_region_beg as u32;
                    if 0 <= sr_pos && (sr_pos as usize) < sample_region_histo.len()
                    {
                        sample_region_histo[(pileup.pos() - sample_region_beg as u32) as usize] = pileup.depth();
                    }
                }
                'FILL_REGION:
                for fr_i in 0..fill_region_len {
                    // find a random depth value from the sample region
                    let mut rng = rand::thread_rng();
                    let sr_i = rng.gen_range(0, sample_region_histo.len());
                    let sr_depth = sample_region_histo[sr_i];
                    // now fill the current fill_region base up to at least srh_depth
                    while fill_region_histo[fr_i as usize] < sr_depth as u64 {
                        let (r1, r2) = read_pair()?;
                        if r1.is_none() && r2.is_none() { break 'FILL_REGION }
                        // find the longest block
                        let mut longest_block = None;
                        let mut longest_block_size = 0;
                        let mut longest_block_b = -1i64;
                        let mut longest_block_r = -1i64;
                        for (r, record)  in vec![r1, r2].iter().enumerate() {
                            match record {
                                Some(record) => {
                                    let blocks = record.aligned_blocks();
                                    for (b, block) in blocks.iter().enumerate() {
                                        if longest_block_size < block[1]-block[0] {
                                            longest_block_size = block[1]-block[0];
                                            longest_block_b = b as i64;
                                            longest_block_r = r as i64;
                                            longest_block = Some(block);
                                        }
                                    }
                                },
                                None => (),
                            }
                        }
                        if let Some(longest_block) = longest_block {
                            // determine the offset
                            let pos = fr_i+fill_region_pos-((longest_block_size / 2)+longest_block[0]-blocks[0][0]);
                            if 0 <= pos && pos < reflen {

                            }
                            // apply the cigar string to the offset to get the new sequence
                        }
                        // center the longest block on the current fill region index
                        // write the reads to the fastq file(s)
                        // update the fill_region_histogram
                    }
                }

            }
            tree.get_mut(chr).r()?.insert(splice.start..splice.stop, Replacement {
                offset: pos-refpos,
                replacement: &splice.replacement,
            });
            pos += splice.replacement.len() as i64;
            refpos = splice.stop;
            sequence.append(&mut splice.replacement.clone());
        }
        if refpos < reflen {
            let replacement = &reference.get(chr).r()?[refpos as usize..reflen as usize];
            tree.get_mut(chr).r()?.insert(refpos..reflen, Replacement {
                offset: pos-refpos,
                replacement,
            });
            sequence.extend(replacement);
        }
        outfasta.write(chr, None, &sequence)?;
    }

    let mut process_reads = |reads: &mut Vec<bam::Record>| -> Result<()> {
        for r in reads {
            let mut newr = r.clone();
            let tid = r.tid();
            let refname = str::from_utf8(header.target_names().get(tid as usize).r()?)?;
            let blocks = r.aligned_blocks();
            for block in &blocks {
                for replacement in tree.get(refname).r()?.find(block[0]..block[1]+1) {
                    newr.set_pos(0);
                    newr.set_mpos(0);
                    let cigar = Vec::from(r.cigar().iter().map(|c| c.clone()).collect::<Vec<Cigar>>());
                    newr.set(
                        r.qname(),
                        Some(&bam::record::CigarString(cigar)),
                        r.seq().encoded,
                        r.qual());
                }
            }
            outfastq.write(
                str::from_utf8(newr.qname())?,
                None,
                newr.seq().encoded,
                newr.qual(),
            )?;
        }
        Ok(())
    };

    let mut reads = Vec::<bam::Record>::new();
    let mut record = Record::new();
    while collated_bam.read(&mut record)? {
        if !reads.is_empty() && record.qname() != reads.get(0).r()?.qname() {
            process_reads(&mut reads)?;
            reads.clear();
        }
        reads.push(record.clone());
    }
    if !reads.is_empty() {
        process_reads(&mut reads)?;
    }
    Ok(())
}
