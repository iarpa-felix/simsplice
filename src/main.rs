use structopt::StructOpt;

use std::str;
use std::vec::Vec;
use std::ops::Range;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::HashSet;
use std::path::Path;

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::{bam, bam::Read as BamRead};
use rust_htslib::{bcf, bcf::Read as BcfRead};

use bio::data_structures::interval_tree::IntervalTree;

use anyhow::{Result, anyhow};

use bio::io::fasta;
use duct::cmd;
use interpol::{format as f};
use num_cpus;
use scopeguard::defer;

// Examples:

// Suppose that we have sequenced a genome of length 10 bp. Let n1,
// n2,..., n8 be the number of reads beginning at each base. Suppose
// the reads are each 3 bases long, the bases in the genome are
// b1,b2,b3,b4,b5,b6,b7,b8,b9,b10, and the insertion that is to be
// replaced in the simulation (i.e., the "sequenced insertion") is at
// bases b3,b4,b5
//
// Case 1: simulated insertion shorter than sequenced insertion
//
// Let the simulated insertion be b3',b4'; that is, we wish to
// simulate reads from the following genome: b1,b2,b3',b4',b6,b7,b8,b9,b10.
// The reads that would be simulated are as follows (frequencies in
// parentheses):
//
// (n1) b1,b2,b3'
// (n2) b2,b3',b4'
// (n3) b3',b4',b6
// (n4) b4',b6,b7
// (n6) b6,b7,b8
// (n7) b7,b8,b9
// (n8) b8,b9,b10
//
// Case 2: simulated insertion longer than sequenced insertion
//
// Let the simulated insertion be b3',b4',b5',b6',b7',b8',b9'; that
// is, we wish to simulate reads from the following genome:
// b1,b2,b3',b4',b5',b6',b7',b8',b9',b6,b7,b8,b9,b10. The reads that
// would be simulated are as follows:
//
// (n1) b1,b2,b3'
// (n2) b2,b3',b4'
// (n3) b3',b4',b5'
// (n4) b4',b5',b6'
// (n5) b5',b6',b7'
// *(n3) b6',b7',b8'
// (n4) b7',b8',b9'
// (n5) b8',b9',b6
// *(n3) b9',b6,b7
// (n6) b6,b7,b8
// (n7) b7,b8,b9
// (n8) b8,b9,b10
//
// Starting at the starred locations, a variant of the algorithm
// that reduces the periodicity of the simulated read depths would be
// to randomly choose between the order n3,n4,n5 and n5,n4,n3. However,
// fundamentally it is probably better to start with a sequenced genome
// with an insertion that is not much shorter than the one that is to
// be simulated.

trait ToResult<T> {
    fn r(self) -> Result<T>;
}
impl<T> ToResult<T> for Option<T> {
    fn r(self) -> Result<T> {
        self.ok_or_else(|| anyhow!("NoneError"))
    }
}

pub fn cigar2exons(cigar: &CigarStringView, pos: i64) -> Result<Vec<Range<i64>>> {
    let mut exons = Vec::<Range<i64>>::new();
    let mut pos = pos;
    for op in cigar {
        match op {
            &Cigar::Match(length) |
            &Cigar::Equal(length) |
            &Cigar::Diff(length) => {
                pos += length as i64;
                if length > 0 {
                    exons.push(Range{start: pos - length as i64, end: pos});
                }
            }
            &Cigar::RefSkip(length) |
            &Cigar::Del(length) => {
                pos += length as i64;
            }
            &Cigar::Ins(_) |
            &Cigar::SoftClip(_) |
            &Cigar::HardClip(_) |
            &Cigar::Pad(_) => (),
        };
    }
    Ok(exons)
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
    #[structopt(long = "outbam", help = "Output BAM file", name="OUTBAMFILE", default_value="")]
    outbamfile: String,
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

    // sort the bam file by name
    let name_sorted_bam = f!("{options.bamfile}.tmp.{std::process::id()}");
    defer! {
        if Path::new(&name_sorted_bam).exists() {
            let _ = std::fs::remove_file(&name_sorted_bam);
        }
    }
    cmd!("samtools","sort","-@",num_cpus::get().to_string(),"-n","-o",&name_sorted_bam,&options.bamfile).run()?;

    // get the list of unmapped read names from the BAM file. These will be used to construct new reads for large insertions
    let mut unmapped = HashSet::<String>::new();
    let mut bam = bam::Reader::from_path(&options.bamfile)?;
    for r in bam.records() {
        let r = r?;
        let tid = r.tid();
        if tid < 0 {
            unmapped.insert(String::from(str::from_utf8(r.qname())?));
        }
    }

    // read the BAM file again and write the output BAM and FASTQ files
    let name_sorted_outbam = f!("{options.outbamfile}.name_sorted.bam");
    defer! {
        if Path::new(&name_sorted_outbam).exists() {
            let _ = std::fs::remove_file(&name_sorted_outbam);
        }
    }
    let mut bam = bam::Reader::from_path(&name_sorted_bam)?;
    let mut outbam = bam::Writer::from_path(
        &name_sorted_outbam,
        &bam::Header::from_template(bam.header()),
        bam::Format::BAM)?;
    let mut outfastq = bio::io::fastq::Writer::to_file(&options.outfastqfile)?;
    let header = bam.header().clone();
    let mut reads = Vec::<bam::Record>::new();

    let mut process_reads = |reads: &mut Vec<bam::Record>| -> Result<()> {
        if reads.len() == 0 {
            return Err(anyhow!("reads is empty!"));
        }
        if reads.len() > 2 {
            return Err(anyhow!("More than two reads in reads!: {}", reads.len()));
        }
        for r in reads {
            let mut newr = r.clone();
            let tid = r.tid();
            let refname = str::from_utf8(header.target_names().get(tid as usize).r()?)?;
            let exons = cigar2exons(&r.cigar(), r.pos() as i64)?;
            for exon in &exons {
                for replacement in tree.get(refname).r()?.find(exon.start..exon.start+1) {
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
            outbam.write(&newr)?;
            outfastq.write(
                str::from_utf8(newr.qname())?,
                None,
                newr.seq().encoded,
                newr.qual(),
            )?;
        }
        Ok(())
    };

    for r in bam.records() {
        let r = r?;
        if !reads.is_empty() && r.qname() != reads.get(0).r()?.qname() {
            process_reads(&mut reads)?;
            reads.clear();
        }
        reads.push(r);
    }
    if !reads.is_empty() {
        process_reads(&mut reads)?;
    }
    cmd!("samtools","sort","-@",num_cpus::get().to_string(),"-o",&options.outbamfile,&name_sorted_outbam).run()?;
    Ok(())
}
