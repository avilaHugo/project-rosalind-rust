use rosalind::*;
use structopt::StructOpt;
use std::path::PathBuf;

#[derive(Debug, StructOpt)]
#[structopt(name = "rosalind", about = "Project Rosalind problem solver.")]
struct Opt {
    /// Problem id
    #[structopt(short, long)]
    problem_id: String,

    /// Input file
    #[structopt(short, long, parse(from_os_str))]
    input_file: PathBuf,
}


fn main() {
    let opt = Opt::from_args();
    match opt.problem_id.as_str() {
        "DNA" => counting_dna_nucleotides(opt.input_file),
	"RNA" => transcribing_dna_into_rna(opt.input_file),
        _ => panic!("Problem ID unknown"),

    }
}
