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
    fn test_transcribing_dna_into_rna(){
	let seq = Seq::from_str("GATGGAACTTGACTACGTAAATT".to_string());
	let result = seq.seq.replace("T", "U");
	assert_eq!("GAUGGAACUUGACUACGUAAAUU", result)
    }
}
