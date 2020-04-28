use structopt::StructOpt;

use std::str;
use std::vec::Vec;
use std::path::Path;
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

use bio::data_structures::interval_tree::IntervalTree;

use anyhow::anyhow;

use bio::io::fasta;
use bio::io::fastq;
use duct::cmd;
use interpol::{format as f};
use num_cpus;
use regex::Regex;
use lazy_static::lazy_static;
use lazy_regex::{regex as re};

use slog::info;
use sloggers::Build;
use sloggers::terminal::{TerminalLoggerBuilder, Destination};
use sloggers::types::Severity;

use rand::Rng;
use rand::seq::SliceRandom;
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

// example data:
//   /data/iarpa/analysis/SRR11140744
//   /data/iarpa/TE/437/{assembly,illumina}


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
#[structopt(name = "simsplice", about = "Apply simulated modifications to a reference genome and a set of aligned reads")]
struct Options {
    #[structopt(short="g", long="genome", help = "Input reference genome FASTA file", name="FASTAFILE")]
    reference: String,
    #[structopt(short="v", long="vcf", help = "Input VCF file with mutations", name="VCFFILE")]
    vcffile: String,
    #[structopt(short="b", long="bam", help = "Input BAM file of reads, *INCLUDING* unaligned reads", name="BAMFILE")]
    bamfile: String,
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
    replacement: Vec<u8>,
}

struct Replacement {
    offset: i64,
    replacement: Vec<u8>,
}

fn main() -> Result<()> {
    std::env::set_var("RUST_BACKTRACE", "full");
    std::env::set_var("RUST_LIB_BACKTRACE", "1");
    let options: Options = Options::from_args();

    let mut builder = TerminalLoggerBuilder::new();
    builder.level(Severity::Info);
    builder.destination(Destination::Stderr);

    let log = builder.build()?;

    let mut rng = rand::thread_rng();
    let tmpid = rng.sample_iter(&Alphanumeric).take(10).collect::<String>();

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
        let index_cmd = [ "samtools", "index", &options.bamfile ];
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
        let tmpfile = f!("{collated_bamfile}.{tmpid}.bam.tmp");
        let collate_cmd = [
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
        info!(log, "Processing VCF record: {:?}", &options.vcffile);
        let r = r?;
        if r.allele_count() > 2 {
            Err(anyhow!("Can't handle multiple alleles in VCF file: {:?}", r))?;
        }
        if r.allele_count() > 0 {
            let rid = r.rid().r()?;
            let refname = str::from_utf8(r.header().rid2name(rid)?)?;
            let pos = r.pos(); // 0-based
            let alleles = r.alleles();
            let ref_allele = alleles.get(0).r()?;
            for allele in &alleles[1..] {
                if !splices.contains_key(refname) {
                    splices.insert(String::from(refname), BTreeSet::new());
                }
                let splice = Splice{
                    start: pos,
                    stop: pos+ref_allele.len() as i64,
                    replacement: Vec::from(*allele),
                };
                info!(log, "Storing splice: {:?}", splice);
                splices.get_mut(refname).r()?.insert(splice);
            }
        }
    }

    {
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

        let mut is_paired = false;
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
                            Some(_) => if !record.is_secondary() { Some(record.clone()) } else { read1 },
                            None => Some(record.clone()),
                        }
                    } else {
                        read2 = match &read2 {
                            Some(_) => if !record.is_secondary() { Some(record.clone()) } else { read2 },
                            None => Some(record.clone()),
                        };
                        is_paired = true;
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
            if is_paired && (read1.is_none() || read2.is_none()) {
                Err(anyhow!("Expected paired-end reads, but only one read found: read1={:?}, read2={:?}", read1, read2))?;
            }
            Ok((read1, read2))
        };

        // build an interval tree using the splice positions to map between original and modified reference coordinates
        // also build and write the output reference fasta file
        let mut newref = BTreeMap::<String, Vec<u8>>::new();
        let mut tree = BTreeMap::<String, IntervalTree<i64, Replacement>>::new();
        for (chr, splices) in &mut splices {
            info!(log, "Building interval tree and writing fasta record for chr {}", chr);
            tree.insert(String::from(chr), IntervalTree::new());
            let mut pos = 0; // position in modified genome
            let mut refpos = 0; // position in original genome
            let mut sequence = Vec::<u8>::new(); // sequence of modified genome
            let reflen = reference.get(chr).r()?.len() as i64; // reference length in original genome
            for splice in splices.iter() {
                if refpos < splice.start {
                    let replacement = &reference.get(chr).r()?[refpos as usize..splice.start as usize];
                    tree.get_mut(chr).r()?.insert(refpos..splice.start, Replacement {
                        offset: pos - refpos,
                        replacement: Vec::from(replacement),
                    });
                    pos += splice.start - refpos;
                    refpos = splice.start;
                    sequence.extend(replacement);
                }
                tree.get_mut(chr).r()?.insert(splice.start..splice.stop, Replacement {
                    offset: pos - refpos,
                    replacement: splice.replacement.clone(),
                });
                pos += splice.replacement.len() as i64;
                refpos = splice.stop;
                sequence.append(&mut splice.replacement.clone());
            }
            if refpos < reflen {
                let replacement = &reference.get(chr).r()?[refpos as usize..reflen as usize];
                tree.get_mut(chr).r()?.insert(refpos..reflen, Replacement {
                    offset: pos - refpos,
                    replacement: Vec::from(replacement),
                });
                sequence.extend(replacement);
            }
            outfasta.write(chr, None, &sequence)?;
            newref.insert(chr.clone(), sequence);
        }

        // fill in newly spliced in genomic regions using random samples of existing read pairs
        let mut record_buffer = Vec::<(Option<Record>, Option<Record>)>::new();
        for (chr, splices) in &mut splices {
            info!(log, "Filling in regions and writing FASTQ records for chr {}", chr);
            tree.insert(String::from(chr), IntervalTree::new());
            let mut pos = 0; // position in modified genome
            let mut refpos = 0; // position in original genome
            let reflen = reference.get(chr).r()?.len() as i64; // reference length in original genome
            for splice in splices.iter() {
                if refpos < splice.start {
                    pos += splice.start - refpos;
                    refpos = splice.start;
                }
                if splice.stop - splice.start < splice.replacement.len() as i64 {
                    let tid = header.tid(chr.as_bytes()).r()?;
                    info!(log, "Filling in region: {}:{}..{}", chr, splice.start, splice.stop);

                    // fill up a sample histogram using -1000..1000 bp window around the area of insertion
                    // from the original genome
                    let sample_region_beg = std::cmp::max(0, refpos - sample_distance) as u64;
                    let sample_region_end = std::cmp::min(reflen, refpos + sample_distance) as u64;
                    info!(log, "Building sample region: {}:{}..{}", chr, sample_region_beg, sample_region_end);
                    let mut sample_region_histo = vec![0u32; (sample_region_end - sample_region_beg) as usize];
                    indexed_bam.fetch(tid, sample_region_beg, sample_region_end)?;
                    for pileup in indexed_bam.pileup() {
                        let pileup = pileup?;
                        let sr_pos = pileup.pos() - sample_region_beg as u32;
                        if 0 <= sr_pos && (sr_pos as usize) < sample_region_histo.len()
                        {
                            sample_region_histo[(pileup.pos() - sample_region_beg as u32) as usize] = pileup.depth();
                        }
                    }

                    // we need to add reads to the new region
                    let fill_region_len = splice.replacement.len() as i64 - (splice.stop - splice.start);
                    let fill_region_pos = pos + (splice.replacement.len() as i64) - fill_region_len;
                    let mut fill_region_histo = vec![0u64; fill_region_len as usize];
                    let mut fr_is = (0..fill_region_len).collect::<Vec<_>>();
                    fr_is.shuffle(&mut rng);
                    'FILL_REGION:
                    for fr_i in fr_is {
                        // find a random depth value from the sample region
                        let sr_i = rng.gen_range(0, sample_region_histo.len());
                        let sr_depth = sample_region_histo[sr_i];
                        // now fill the current fill_region base up to at least srh_depth
                        'FILL_REGION_HISTO:
                        while fill_region_histo[fr_i as usize] < sr_depth as u64 {
                            let (r1, r2) = read_pair()?;
                            if r1.is_none() && r2.is_none() { break 'FILL_REGION }
                            let reads = [&r1, &r2];
                            let blocks = reads.iter().map(
                                |&r| r.as_ref().map(
                                    |rr| rr.aligned_blocks())).collect::<Vec<Option<Vec<[i64; 2]>>>>();
                            // find the longest block
                            let mut longest_block_b = -1i64;
                            let mut longest_block_r = -1i64;
                            for (r, record) in reads.iter().enumerate() {
                                if let Some(_) = record {
                                    for (b, block) in blocks[r].as_ref().r()?.iter().enumerate() {
                                        if longest_block_r >= 0 &&
                                            longest_block_b >= 0 &&
                                            blocks[longest_block_r as usize].as_ref().r()?[longest_block_b as usize][1] - blocks[longest_block_r as usize].as_ref().r()?[longest_block_b as usize][0] < block[1] - block[0]
                                        {
                                            longest_block_b = b as i64;
                                            longest_block_r = r as i64;
                                        }
                                    }
                                }
                            }
                            if longest_block_r >= 0 && longest_block_b >= 0 {
                                let longest_block = blocks[longest_block_r as usize].as_ref().r()?[longest_block_b as usize];
                                let longest_block_len = longest_block[1] - longest_block[0];
                                let fill_offset = (longest_block[0] + (longest_block_len / 2)) +
                                    (fill_region_len / 2);
                                // skip read pairs not on the same chr or where after translating the positions
                                // either read would go past the edge of genome space
                                for (r, record) in [&r1, &r2].iter().enumerate() {
                                    if let Some(record) = record {
                                        if blocks[r].as_ref().r()?[0][0] - fill_offset + fill_region_pos < 0 ||
                                            record.tid() != reads[longest_block_r as usize].as_ref().r()?.tid() ||
                                            blocks[r].as_ref().r()?.last().r()?[1] - fill_offset + fill_region_pos > newref.get(str::from_utf8(header.target_names()[record.tid() as usize])?).r()?.len() as i64
                                        {
                                            record_buffer.push((r1, r2));
                                            continue 'FILL_REGION_HISTO
                                        }
                                    }
                                }
                                for (r, record) in [&r1, &r2].iter_mut().enumerate() {
                                    if let Some(record) = record {
                                        for block in blocks[r].as_ref().r()?.iter() {
                                            let fill_start = block[0] - fill_offset;
                                            let fill_end = block[1] - fill_offset;
                                            // fill in the fillin_region_histo with the blocks
                                            for i in std::cmp::max(0, fill_start)..
                                                std::cmp::min(fill_region_len, fill_end)
                                            {
                                                fill_region_histo[i as usize] += 1
                                            }
                                        }
                                        let chr = str::from_utf8(header.target_names().get(
                                            record.tid() as usize).r()?)?;
                                        let oldrefseq = reference.get(chr).r()?;
                                        let refseq = newref.get(chr).r()?;
                                        let ap = record.aligned_pairs();
                                        let oldseq = record.seq().encoded;
                                        let mut seq = Vec::from(record.seq().encoded);
                                        for a in &ap {
                                            let qpos = a[0];
                                            let rpos = a[1];
                                            let newrpos = rpos - fill_offset + fill_region_pos;
                                            if oldseq[qpos as usize].to_ascii_uppercase() == oldrefseq[rpos as usize].to_ascii_uppercase() {
                                                seq[qpos as usize] = refseq[newrpos as usize];
                                            }
                                            else if oldseq[qpos as usize].to_ascii_uppercase() == refseq[newrpos as usize].to_ascii_uppercase() {
                                                let nucs = [b'A', b'C', b'G', b'T'].iter().filter(|n| **n != refseq[newrpos as usize].to_ascii_uppercase()).collect::<Vec<_>>();
                                                let randnuc = rng.gen_range(0, 3);
                                                seq[qpos as usize] = *nucs[randnuc];
                                            }
                                            seq[qpos as usize] = if oldseq[qpos as usize].is_ascii_uppercase()
                                                { seq[qpos as usize].to_ascii_uppercase() }
                                            else if oldseq[qpos as usize].is_ascii_lowercase()
                                                { seq[qpos as usize].to_ascii_lowercase() }
                                            else { seq[qpos as usize] };
                                        }
                                        let out = if r == 0 { &mut outfastq } else { outfastq2.as_mut().r()? };
                                        out.write(
                                            str::from_utf8(record.qname())?,
                                            None,
                                            seq.as_slice(),
                                            record.qual(),
                                        )?;
                                    }
                                }
                            }
                        }
                    }
                }
                pos += splice.replacement.len() as i64;
                refpos = splice.stop;
            }
        }

        // translate the coordinates and write out the rest of the reads
        info!(log, "Translating and writing the rest of the FASTQ reads");
        loop {
            let (r1, r2) = if !record_buffer.is_empty() {
                record_buffer.pop().r()?
            } else {
                read_pair()?
            };
            if r1.is_none() && r2.is_none() && record_buffer.is_empty() { break }
            for (r, record) in [&r1, &r2].iter().enumerate() {
                if let Some(record) = record {
                    let fprime = if record.is_reverse() { record.reference_end() } else { record.pos() };
                    let refname = str::from_utf8(header.target_names().get(record.tid() as usize).r()?)?;
                    let oldrefseq = reference.get(refname).r()?;
                    let refseq = newref.get(refname).r()?;
                    for entry in tree.get(refname).r()?.find(fprime..fprime + 1) {
                        let replacement = entry.data();
                        let ap = record.aligned_pairs();
                        let oldseq = record.seq().encoded;
                        let mut seq = Vec::from(record.seq().encoded);
                        for a in &ap {
                            let qpos = a[0];
                            let rpos = a[1];
                            let newrpos = rpos + replacement.offset;
                            if 0 <= newrpos && newrpos < refseq.len() as i64 {
                                if oldseq[qpos as usize].to_ascii_uppercase() == oldrefseq[rpos as usize].to_ascii_uppercase()
                                {
                                    seq[qpos as usize] = refseq[newrpos as usize];
                                }
                                else if oldseq[qpos as usize].to_ascii_uppercase() == refseq[newrpos as usize].to_ascii_uppercase() {
                                    let nucs = [b'A', b'C', b'G', b'T'].iter().filter(|n| **n != refseq[newrpos as usize].to_ascii_uppercase()).collect::<Vec<_>>();
                                    let randnuc = rng.gen_range(0, 3);
                                    seq[qpos as usize] = *nucs[randnuc];
                                }
                                seq[qpos as usize] = if oldseq[qpos as usize].is_ascii_uppercase()
                                { seq[qpos as usize].to_ascii_uppercase() }
                                else if oldseq[qpos as usize].is_ascii_lowercase()
                                { seq[qpos as usize].to_ascii_lowercase() }
                                else { seq[qpos as usize] };
                            }
                        }
                        let out = if r == 0 { &mut outfastq } else { outfastq2.as_mut().r()? };
                        out.write(
                            str::from_utf8(record.qname())?,
                            None,
                            seq.as_slice(),
                            record.qual(),
                        )?;
                        break;
                    }
                }
            }
        }
    }

    // sort the fastq files by read name
    info!(log, "Sorting the FASTQ files by read name");
    let tmp1 = f!("{options.outfastqfile}.{tmpid}.tmp");

    let fastq_sort = f!(r#"shopt -s nocasematch; f={quote(&options.outfastqfile)}; if [[ ${{f##*.}} = gz ]]; then gzip -cd "$f"; else cat "$f"; fi | paste - - - - | sort -k1,1 | tr "\t" "\n" |if [[ ${{f##*.}} = gz ]]; then gzip; else cat; fi >{quote(&tmp1)} && mv -v {quote(&tmp1)} {quote(&options.outfastqfile)}"#);
    let fastq_sort_cmd = ["bash","-c",&fastq_sort];
    info!(log, "Running command: {}", shell_words::join(&fastq_sort_cmd));
    cmd(fastq_sort_cmd[0],&fastq_sort_cmd[1..]).run()?;
    if !options.outfastqfile2.is_empty() {
        let tmp2 = f!("{options.outfastqfile2}.{tmpid}.tmp");
        let fastq_sort = f!(r#"shopt -s nocasematch; f={quote(&options.outfastqfile2)}; if [[ ${{f##*.}} = gz ]]; then gzip -cd "$f"; else cat "$f"; fi | paste - - - - | sort -k1,1 | tr "\t" "\n" |if [[ ${{f##*.}} = gz ]]; then gzip; else cat; fi >{quote(&tmp2)} && mv -v {quote(&tmp2)} {quote(&options.outfastqfile2)}"#);
        let fastq_sort_cmd = ["bash","-c",&fastq_sort];
        info!(log, "Running command: {}", shell_words::join(&fastq_sort_cmd));
        cmd(fastq_sort_cmd[0],&fastq_sort_cmd[1..]).run()?;
    }
    Ok(())
}