use structopt::StructOpt;
use bio::data_structures::interval_tree::IntervalTree;

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
use std::ops::Range;
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
    chr2: String,
    range2: Range<u64>,
}

#[derive(StructOpt, Debug)]
#[structopt(name = "liftover", about = "Convert genomic coordinates from one system to another")]
struct Options {
    #[structopt(short="l", long="liftover", help = "Liftover file of format: ", name="LIFTOVER")]
    liftover: String,
    #[structopt(short="i", long="in", help = "file input (can be GFF/GTF/BED format", name="IN")]
    input: Vec<String>,
    #[structopt(short="o", long="output", help = "Output file", name="OUT")]
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
        if split.len() < 6 {
            bail!("Liftover file {}, line {} badly formatted: {}", &options.liftover, lineno, line);
        }
        let l = Liftover {
            chr1: split[0].to_string(),
            range1: split[1].parse::<u64>()?..split[2].parse::<u64>()?,
            chr2: split[3].to_string(),
            range2: split[4].parse::<u64>()?..split[5].parse::<u64>()?,
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

    for input in &options.input {
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

                if let Some(tree) = tree.get(record.seqname()) {
                    for item in tree.find(*record.start()-1..*record.start()) {
                        chr = Some(&item.data().chr2);
                        start = Some(item.data().range2.start+1);
                        break;
                    }
                    if let Some(chr) = chr {
                        for item in tree.find(*record.end()-1..*record.end()) {
                            if chr == &item.data().chr2 {
                                end = Some(item.data().range2.end);
                                break;
                            }
                        }
                    }
                }
                if let (Some(chr), Some(start), Some(end)) = (chr, start, end) {
                    *record.seqname_mut() = chr.to_string();
                    *record.start_mut() = start;
                    *record.end_mut() = end;
                    out.write(&mut record)?;
                }
            }
        }
        else {
            let mut reader = bio::io::bed::Reader::from_file(input)?;
            let mut out = bio::io::bed::Writer::to_file(&options.out)?;
            for record in &mut reader.records() {
                let mut record = record?.clone();
                let mut chr = Option::<&str>::None;
                let mut start = Option::<u64>::None;
                let mut end = Option::<u64>::None;

                if let Some(tree) = tree.get(record.chrom()) {
                    for item in tree.find(record.start()..record.start()+1) {
                        chr = Some(&item.data().chr2);
                        start = Some(item.data().range2.start);
                        break;
                    }
                    if let Some(chr) = chr {
                        for item in tree.find(record.end()-1..record.end()) {
                            if chr == &item.data().chr2 {
                                end = Some(item.data().range2.end);
                                break;
                            }
                        }
                    }
                }
                if let (Some(chr), Some(start), Some(end)) = (chr, start, end) {
                    record.set_chrom(chr);
                    record.set_start(start);
                    record.set_end(end);
                    out.write(&mut record)?;
                }
            }
        }
    }

    Ok(())
}