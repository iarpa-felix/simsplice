use structopt::StructOpt;

use std::str;
use std::vec::Vec;
use std::ops::Range;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::HashSet;
use std::path::Path;

use rust_htslib::bam::record::Record;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::{bam, bam::Read as BamRead};
use rust_htslib::{bcf, bcf::Read as BcfRead};
use rust_htslib::bam::ext::BamRecordExtensions;

use bio::data_structures::interval_tree::IntervalTree;

use anyhow::Result;
use anyhow::anyhow;

use bio::io::fasta;
use duct::cmd;
use interpol::{format as f, println, eprintln};
use num_cpus;
use scopeguard::defer;
use regex::Regex;
use lazy_static::lazy_static;
use lazy_regex::regex;

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
    #[structopt(long = "outfasta", help = "Output reference FASTA file", name="OUTFASTAFILE", default_value="")]
    outfastafile: String,
    #[structopt(long = "outfastq", help = "Output FASTQ file", name="OUTFASTQFILE", default_value="")]
    outfastqfile: String,
}

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone)]
struct Splice {
    start: i64,
    stop: i64,
    replacement: Vec<u8>,
}

struct Replacement<'a> {
    pos: i64,
    replacement: &'a[u8],
}

//437

fn main() -> Result<()> {
    std::env::set_var("RUST_BACKTRACE", "full");
    std::env::set_var("RUST_LIB_BACKTRACE", "1");
    let options = Options::from_args();

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

    // collate the bam file by name
    let collated_bam = f!(r#"{regex!(r"\.bam$").replace_all(&options.bamfile,"")}.collated.bam"#);
    let cpus = num_cpus::get().to_string();
    if !Path::new(&collated_bam).exists() {
        let tmpid = rand::thread_rng().sample_iter(&Alphanumeric).take(10).collect::<String>();
        let tmpfile = f!("{collated_bam}.{tmpid}");
        let collate_cmd = vec![
            "samtools","collate",
            "-@",&cpus,
            "-o",&tmpfile,
            &options.bamfile,
        ];
        info!(log, "Running command: {}", shell_words::join(&collate_cmd));
        cmd(collate_cmd[0],&collate_cmd[1..]).run()?;
        info!(log, "Moving {} to {}", &tmpfile, &collated_bam);
        std::fs::rename(&tmpfile, &collated_bam)?;
    }
    let mut bam = bam::Reader::from_path(&collated_bam)?;

    // build an interval tree using the splice positions to map between original and modified reference coordinates
    let mut reflen = 0i64;
    let mut refpos = 0i64;
    let mut pos = 0i64;
    let mut tree = BTreeMap::<String,IntervalTree<i64,Replacement>>::new();
    for (chr, splices) in &splices {
        if !tree.contains_key(chr) {
            tree.insert(String::from(chr), IntervalTree::new());
            pos = 0;
            refpos = 0;
            reflen = reference.get(chr).r()?.len() as i64;
        }
        for splice in splices {
            if refpos < splice.start {
                tree.get_mut(chr).r()?.insert(refpos..splice.start, Replacement {
                    pos: pos-refpos,
                    replacement: &reference.get(chr).r()?[refpos as usize..splice.start as usize],
                });
                pos += splice.start-refpos;
            }
            tree.get_mut(chr).r()?.insert(splice.start..splice.stop, Replacement {
                pos: pos-refpos,
                replacement: &splice.replacement,
            });
            pos += splice.replacement.len() as i64;
            refpos = splice.stop;
        }
        if refpos < reflen {
            tree.get_mut(chr).r()?.insert(refpos..reflen, Replacement {
                pos: pos-refpos,
                replacement: &reference.get(chr).r()?[refpos as usize..reflen as usize],
            });
        }
    }

    let mut outfasta = bio::io::fasta::Writer::to_file(&options.outfastafile)?;
    let mut outfastq = bio::io::fastq::Writer::to_file(&options.outfastqfile)?;
    let header = bam.header().clone();
    let mut reads1 = Vec::<bam::Record>::new();
    let mut reads2 = Vec::<bam::Record>::new();

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
    while bam.read(&mut record)? {
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
