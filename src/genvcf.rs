use structopt::StructOpt;
use std::collections::{BTreeSet, HashMap};
use linked_hash_map::LinkedHashMap;
use anyhow::anyhow;
use bio::io::fasta;
use regex::{Regex, Captures};
use lazy_static::lazy_static;
use lazy_regex::{regex as re};

use slog::info;
use sloggers::Build;
use sloggers::terminal::{TerminalLoggerBuilder, Destination};
use sloggers::types::Severity;

use rand::Rng;
use rand::seq::SliceRandom;
use std::ops::Range;
use rand::distributions::Distribution;
use rust_htslib::bcf;
use rust_htslib::bcf::Format;

pub type Result<T, E = anyhow::Error> = core::result::Result<T, E>;
trait ToResult<T> {
    fn r(self) -> Result<T>;
}
impl<T> ToResult<T> for Option<T> {
    fn r(self) -> Result<T> {
        self.ok_or_else(|| anyhow!("NoneError"))
    }
}

#[derive(Debug)]
pub struct Nucleotide;

impl Distribution<u8> for Nucleotide {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> u8 {
        const RANGE: u32 = 4;
        const GEN_NUCLEOTIDE: &[u8] = b"ACGT";
        loop {
            let var = rng.next_u32() >> (32 - 2);
            if var < RANGE {
                return GEN_NUCLEOTIDE[var as usize] as u8
            }
        }
    }
}

#[derive(StructOpt, Debug)]
#[structopt(name = "simsplice", about = "Apply simulated modifications to a reference genome and a set of aligned reads")]
struct Options {
    #[structopt(long="genome", help = "Input reference genome FASTA file", name="FASTAFILE")]
    reference: String,
    #[structopt(long="vcf", help = "Output VCF file with mutations", name="VCFFILE")]
    vcffile: String,
    #[structopt(long="prob-insert", help = "Probability of an insertion event (0.0-1.0)", name="PROB", default_value="0.3")]
    prob_insert: f64,
    #[structopt(long="prob-delete", help = "Probability of a deletion event (0.0-1.0)", name="PROB", default_value="0.3")]
    prob_delete: f64,
    #[structopt(long="prob-splice", help = "Probability of a splicing event (0.0-1.0)", name="PROB", default_value="0.3")]
    prob_splice: f64,
    #[structopt(long="delete-range", help = "Range in bp for deletions", name="BP", default_value="200-1000")]
    delete_range: String,
    #[structopt(long="insert-range", help = "Range in bp for insertions", name="BP", default_value="200-1000")]
    insert_range: String,
    #[structopt(long="num-modifications", help = "The number of modification to make", name="NUM")]
    num_modifications: Option<u64>,
    #[structopt(long="modifications", help = "Specify modifications. A series of chromosome choordinate/ranges:lengths. Examples: chr1:10000, chr1:1000-2000, chr1:1000-2000:400", name="MODIFICATOINS", default_value="1")]
    modifications: Vec<String>,
}

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Clone)]
struct Splice {
    chr: String,
    start: i64,
    end: i64,
    replacement: Vec<u8>,
}

#[derive(Debug, PartialEq)]
enum ModType {
    Insert,
    Delete,
    Splice
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

    // read in the reference FASTA file
    let mut reference = LinkedHashMap::<String,Vec<u8>>::new();
    let fasta = fasta::Reader::from_file(&options.reference)?;
    for r in fasta.records() {
        let r = r?;
        reference.insert(String::from(r.id()), Vec::from(r.seq()));
    }
    let refids = reference.keys().enumerate().
        map(|(i,n)| (n.as_str(), i)).collect::<HashMap<_,_>>();

    // validate options
    let total_prob = options.prob_insert + options.prob_delete + options.prob_splice;
    if total_prob == 0.0 {
        Err(anyhow!("prob-insert, prob-delete and prob-splice sum to zero!"))?;
    }

    let delete_range: Vec<Range<i64>> = re!(r"^([0-9]+)(-|\.\.)([0-9]+)$").captures(&options.delete_range).
        iter().map(|c| Ok((c[1].parse::<i64>()?)..(c[3].parse::<i64>()?))).collect::<Result<Vec<_>>>()?;
    let delete_range = delete_range.first().ok_or_else(||
        anyhow!("Delete range could not be parsed: {}", &options.delete_range)
    )?;

    let insert_range: Vec<Range<i64>> = re!(r"^([0-9]+)(-|\.\.)([0-9]+)$").captures(&options.insert_range).
        iter().map(|c| Ok(c[1].parse::<i64>()?..c[3].parse::<i64>()?)).collect::<Result<Vec<_>>>()?;
    let insert_range = insert_range.first().ok_or_else(||
        anyhow!("Insert range could not be parsed: {}", &options.insert_range)
    )?;
    if options.num_modifications.is_some() && !options.modifications.is_empty() {
        Err(anyhow!("Specifying both num-modifications and modifications is disallowed"))?;
    }
    let modifications = (0..options.num_modifications.unwrap_or(0)).map(|_| "".to_string()).collect::<Vec<_>>();
    let mut splices = BTreeSet::<Splice>::new();
    for m in if modifications.is_empty() {&options.modifications} else {&modifications} {
        let splice = re!(r"(?xi)
            ^(?P<chr>[^:]+)?
            (:((?P<start>[0-9]+)(([.][.]|-)(?P<end>[0-9]+))?)?
            (:(?P<replace>[0-9]+|[ACGTN]+)?)?)?$
            ").captures(m).iter().map(|c: &Captures|
            {
                let chr = c.name("chr").map(|c| c.as_str());
                let start = c.name("start").map(|s| s.as_str().parse::<i64>()).transpose()?;
                let end = c.name("end").map(|e| e.as_str().parse::<i64>()).transpose()?;
                let replace = c.name("replace").map(|r| r.as_str());

                let mod_type = [
                    (&ModType::Insert, &options.prob_insert),
                    (&ModType::Delete, &options.prob_delete),
                    (&ModType::Splice, &options.prob_splice),
                ].choose_weighted(&mut rng, |t| t.1)?.0;

                let chrs = reference.keys().filter(|n| {
                    let min_len = if let (Some(start), Some(end)) = (start, end) {end-start}
                        else if let Some(end)=end {end}
                    else if let Some(start)=start {start+delete_range.start}
                    else {delete_range.start};
                    reference[n.as_str()].len() as i64 <= min_len
                }).map(|c| c.as_str()).collect::<Vec<_>>();
                if chrs.is_empty() {
                    Err(anyhow!("Could not create modification: {}: No suitable refseqs found", m))?;
                }
                let chr = match chr {
                    Some(chr) => chr,
                    None => chrs.choose_weighted(&mut rng, |c| reference[*c].len())?,
                };
                let start = match start {
                    Some(start) => start,
                    None => rng.gen_range(0, reference[chr].len() as i64 - delete_range.start),
                };
                let end = match end {
                    Some(end) => end,
                    None => if *mod_type == ModType::Insert { start }
                    else {
                        rng.gen_range(start+delete_range.start, reference[chr].len() as i64)
                    },
                };
                let replace = match replace {
                    Some(replace) => match replace.parse::<usize>() {
                        Ok(num_replace) => rng.sample_iter(&Nucleotide).take(num_replace).collect::<Vec<u8>>(),
                        Err(_) => replace.as_bytes().to_vec(),
                    }
                    None => if *mod_type == ModType::Delete { vec![] }
                    else {
                        let num_replace = rng.gen_range(insert_range.start, insert_range.end);
                        rng.sample_iter(&Nucleotide).take(num_replace as usize).collect::<Vec<u8>>()
                    },
                };
                Ok(Splice {
                    chr: chr.to_string(),
                    start,
                    end,
                    replacement: replace,
                })
            }).collect::<Result<Vec<Splice>>>()?;
        splices.extend(splice);
    }

    let mut header = bcf::header::Header::new();
    for (name, seq) in &reference {
        header.push_record(format!("##contig=<ID={},length={}>", name, seq.len()).as_bytes());
    }
    let mut vcf = bcf::Writer::from_path(&options.vcffile, &header, false, Format::VCF)?;
    let mut lastsplice: Option<&Splice> = None;
    for splice in &splices {
        if let Some(lastsplice) = &lastsplice {
            if lastsplice.chr == splice.chr && splice.start < lastsplice.end {
                Err(anyhow!("Overlapping modifications found, cannot continue: {:?}: {:?}", &splice, &lastsplice))?;
            }
        }
        lastsplice = Some(&splice);

        let mut record = vcf.empty_record();
        record.set_rid(Some(refids[splice.chr.as_str()] as u32));
        if splice.start == 0 {
            record.set_pos(splice.start);
            record.set_alleles(&[
                &reference[&splice.chr][splice.start as usize..std::cmp::min(reference[&splice.chr].len() as i64, splice.end+1) as usize],
                &[splice.replacement.to_vec(),
                    reference[&splice.chr][splice.end as usize..std::cmp::min(reference[&splice.chr].len() as i64, splice.end as i64) as usize].to_vec()].concat(),
            ])?;
        }
        else {
            record.set_pos(splice.start-1);
            record.set_alleles(&[
                &reference[&splice.chr][(splice.start-1) as usize..splice.end as usize],
                &[reference[&splice.chr][(splice.start-1) as usize..splice.start as usize].to_vec(),
                    splice.replacement.to_vec()].concat(),
            ])?;
        }
        vcf.write(&record)?;
    }
    Ok(())
}

