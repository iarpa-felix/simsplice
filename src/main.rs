use structopt::StructOpt;

use std::vec::Vec;
use std::ops::Range;

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::bam;
use rust_htslib::bcf;

use anyhow::{Result, anyhow};

use bio::io::fasta;

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

pub fn cigar2exons(cigar: &CigarStringView, pos: u64) -> Result<Vec<Range<u64>>> {
    let mut exons = Vec::<Range<u64>>::new();
    let mut pos = pos;
    for op in cigar {
        match op {
            &Cigar::Match(length) |
            &Cigar::Equal(length) |
            &Cigar::Diff(length) => {
                pos += length as u64;
                if length > 0 {
                    exons.push(Range{start: pos - length as u64, end: pos});
                }
            }
            &Cigar::RefSkip(length) |
            &Cigar::Del(length) => {
                pos += length as u64;
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
    #[structopt(long = "out", help = "Output file prefix", name="PREFIX", default_value="")]
    out: String,
}

fn main() -> Result<()> {
    std::env::set_var("RUST_BACKTRACE", "full");
    std::env::set_var("RUST_LIB_BACKTRACE", "1");
    let options = Options::from_args();
    let mut bam = bam::Reader::from_path(&options.bamfile)?;
    let mut vcf = bcf::Reader::from_path(&options.vcffile)?;
    let fasta = fasta::Reader::from_file(&options.reference)?;

    Ok(())
}
