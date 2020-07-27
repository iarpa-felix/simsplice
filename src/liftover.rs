use structopt::StructOpt;
use bio::data_structures::interval_tree::IntervalTree;
use bio_types::strand::Strand;

use regex::{Regex};
use lazy_static::lazy_static;
use lazy_regex::{regex as re};

use slog::info;
use slog_term;
use slog::Drain;

use eyre::{eyre, bail};
use bio::io::gff::GffType;
use std::fs::File;
use std::io::{BufReader, BufRead};
use linked_hash_map::LinkedHashMap;
use std::ops::{Range, Neg};
use bio::io::gff;

pub type Report = eyre::Report<color_eyre::Context>;
pub type Result<T, E = Report> = core::result::Result<T, E>;
trait ToResult<T> {
    fn r(self) -> Result<T>;
}
impl<T> ToResult<T> for Option<T> {
    fn r(self) -> Result<T> {
        self.ok_or_else(|| eyre!("NoneError"))
    }
}

#[derive(Clone)]
struct Liftover {
    chr1: String,
    range1: Range<u64>,
    strand1: Strand,
    chr2: String,
    range2: Range<u64>,
    strand2: Strand,
}

#[derive(StructOpt, Debug)]
#[structopt(name = "liftover", about = "Convert genomic coordinates from one system to another")]
struct Options {
    #[structopt(help = "file input (in GFF/GTF format)", name="IN")]
    input: String,
    #[structopt(help = "Liftover file of format: ", name="LIFTOVER")]
    liftover: String,
    #[structopt(help = "Output file", name="OUT")]
    out: String,
}

fn main() -> Result<()> {
    std::env::set_var("RUST_BACKTRACE", "full");
    std::env::set_var("RUST_LIB_BACKTRACE", "1");
    let options: Options = Options::from_args();

    let plain = slog_term::PlainSyncDecorator::new(std::io::stderr());
    let log = slog::Logger::root(
        slog_term::FullFormat::new(plain).build().fuse(),
        slog::o!()
    );

    if options.input.len() == 0 {
        bail!("No input files given!")
    }

    let liftover = File::open(&options.liftover)?;
    let mut liftover = BufReader::new(liftover);
    let mut line = String::new();
    let mut lineno=1;
    let mut tree = LinkedHashMap::<String,IntervalTree<u64,Liftover>>::new();
    while liftover.read_line(&mut line)? > 0 {
        let line = line.trim_end_matches("\n");
        let split = line.split("\t").collect::<Vec<_>>();
        if split.len() < 8 {
            bail!("Liftover file {}, line {} badly formatted: {}", &options.liftover, lineno, line);
        }
        let l = Liftover {
            chr1: split[0].to_string(),
            range1: split[1].parse::<u64>()?..split[2].parse::<u64>()?,
            strand1: Strand::from_char(&split[3].chars().nth(0).r()?)?,
            chr2: split[4].to_string(),
            range2: split[5].parse::<u64>()?..split[6].parse::<u64>()?,
            strand2: Strand::from_char(&split[7].chars().nth(0).r()?)?,
        };
        if l.range1.end-l.range1.start != l.range2.end-l.range2.start {
            bail!("Range1 and range2 lengths unequal: range1={:?}, range2={:?} in line {}: {}",
            &l.range1, &l.range2, lineno, line);
        }
        if !tree.contains_key(&l.chr1) {
            tree.insert(l.chr1.clone(), IntervalTree::new());
        }
        tree[&l.chr1].insert(&l.range1, l.clone());
        lineno += 1;
    }

    let input = &options.input;
    let lower = input.to_ascii_lowercase();
    if re!(r"(?i)[.]{gtf,gff,gff3,gff2}$").is_match(input) {
        let gfftype = if lower.ends_with(".gtf") { GffType::GTF2 }
        else if lower.ends_with(".gff2") { GffType::GFF2 }
        else { GffType::GFF3 };
        let mut reader = bio::io::gff::Reader::from_file(input, gfftype)?;
        let mut out = gff::Writer::to_file(&options.out, gfftype)?;

        for record in &mut reader.records() {
            let mut record = record?.clone();
            let mut chr = Option::<&str>::None;
            let mut start = Option::<u64>::None;
            let mut end = Option::<u64>::None;
            let mut strand = Option::<Strand>::None;

            if let Some(tree) = tree.get(record.seqname()) {
                for item in tree.find(*record.start()-1..*record.start()) {
                    chr = Some(&item.data().chr2);
                    if let Some(s) = record.strand() {
                        if item.data().strand1 != item.data().strand2 && item.data().strand1 != Strand::Unknown && item.data().strand2 != Strand::Unknown {
                            start = Some(item.data().range2.end-((record.start()-1)-item.data().range1.start));
                            strand = Some(s.neg());
                        }
                        else {
                            start = Some((record.start()-1)-item.data().range1.start+item.data().range2.start);
                            strand = Some(s);
                        }
                    }
                    else {
                        start = Some((record.start()-1)-item.data().range1.start+item.data().range2.start);
                    }
                    for item in tree.find(*record.end()-1..*record.end()) {
                        if let Some(s) = record.strand() {
                            if item.data().strand1 != item.data().strand2 && item.data().strand1 != Strand::Unknown && item.data().strand2 != Strand::Unknown {
                                end = Some(item.data().range2.end-((record.end()-1)-item.data().range1.start));
                                if s != Strand::Unknown && strand.is_some() && Some(s.neg()) != strand { continue }
                                strand = Some(s.neg());
                            }
                            else {
                                end = Some((record.end()-1)-item.data().range1.start+item.data().range2.start);
                                if s != Strand::Unknown && strand.is_some() && Some(s) != strand { continue }
                                strand = Some(s);
                            }
                        }
                        else {
                            end = Some((record.end()-1)-item.data().range1.start+item.data().range2.start);
                        }
                        break;
                    }
                    break;
                }
            }
            if let (Some(chr), Some(start), Some(end)) = (chr, start, end) {
                *record.seqname_mut() = chr.to_string();
                *record.start_mut() = start;
                *record.end_mut() = end;
                if let Some(s) = strand {
                    *record.strand_mut() = s.strand_symbol().to_string();
                }
                out.write(&mut record)?;
            }
        }
    }
    else {
        bail!("Unknown file extension: {}", &options.input);
    }

    Ok(())
}