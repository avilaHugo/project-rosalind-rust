use std::fs;
use std::path::PathBuf;
use std::collections::HashMap;
use itertools::Itertools;

pub struct Seq {
    seq: String 
}

impl Seq {
    pub fn from_str(seq: String) -> Self {
	Self {
	    seq
	}
    }

    pub fn count_kmers(&self, k: usize) -> HashMap<String, u64> {
	let mut kmer_counter = HashMap::new();
	for chunk in &self.seq.chars().chunks(k) {
	    let kmer: String = chunk.collect();

	    if !kmer_counter.contains_key(&kmer) {
		kmer_counter.insert(kmer.clone(), 0);
	    }

	    *kmer_counter.get_mut(&kmer).unwrap() += 1;

	};

	kmer_counter
    }

    pub fn get_nuc_pair(nuc: char) -> String {
	match nuc {
	    'A' => String::from("T"),
	    'C' => String::from("G"),
	    'T' => String::from("A"),
	    'G' => String::from("C"),
	    _ => panic!("Unknown nuc: '{}'", nuc),
	}
    }

    pub fn reverse_complement_dna(&self) -> String {
	self.seq
	    .chars()
	    .map(|x| { Seq::get_nuc_pair(x) })
	    .join("")
	    .chars()
	    .rev()
	    .join("")
    }
}

pub fn counting_dna_nucleotides(input_file: PathBuf) {
    let seq = fs::read_to_string(input_file).expect("File not found");
    let kmer_counter = Seq::from_str(seq).count_kmers(1);
    let result = ["A", "C", "G", "T"]
        .map(|x| kmer_counter.get(x).unwrap().to_string())
        .join(" ");

    println!("{}", result)
}

pub fn transcribing_dna_into_rna(input_file: PathBuf) {
	let content = fs::read_to_string(input_file).expect("File not found !.");
	let seq = Seq::from_str(content);
	let result = seq.seq.replace("T", "U");
	println!("{}", result)
}
pub fn complementing_a_strand_of_dna(input_file: PathBuf) {
    let content = fs::read_to_string(input_file).expect("File not found !.").trim().to_string();
    let seq = Seq::from_str(content);
    println!("{}", seq.reverse_complement_dna())
}

pub fn rabbits_and_recurrence_relations(input_file: PathBuf) {
    let (n, k): (usize, usize) = fs::read_to_string(input_file)
	.expect("File not found")
	.trim()
	.split(" ")
	.map(|x| { x.parse::<usize>().unwrap()})
	.collect_tuple()
	.unwrap();

    // 1    2    3    4    5
    // 1b0a 0b1a 3b1a 3b4a 12b7a
    // 1    1    4    7    19
    let mut b = 1;
    let mut a = 0;

    for _ in 2..n+1 {
	let temp_b = b;
	b = a * k;
	a = temp_b + a;
    }

    println!("{}", a+b)
}

#[cfg(test)]
mod tests {
    use crate::Seq;

    #[test]
    fn test_count_kmers() {
	let seq = Seq::from_str("ACTG".to_string());
	let kmer_counter = seq.count_kmers(1);
	assert_eq!(1, *kmer_counter.get("A").unwrap())
    }

    #[test]
    fn test_reverse_complement(){
	let seq = Seq::from_str(String::from("AAAACCCGGT"));
	assert_eq!(String::from("ACCGGGTTTT"), seq.reverse_complement_dna())
    }

}
