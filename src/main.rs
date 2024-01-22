use rosalind;
use std::fs;
use std::path::PathBuf;
use structopt::StructOpt;

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

fn counting_dna_nucleotides(input_file: PathBuf) {
    let seq = fs::read_to_string(input_file).expect("File not found");
    let kmer_counter = rosalind::Seq::from_str(seq).count_kmers(1);
    let result = ["A", "C", "G", "T"]
        .map(|x| kmer_counter.get(x).unwrap().to_string())
        .join(" ");

    println!("{}", result)
}

fn main() {
    let opt = Opt::from_args();
    match opt.problem_id.as_str() {
        "DNA" => counting_dna_nucleotides(opt.input_file),
        _ => panic!("Problem ID unknown"),
    }
}
