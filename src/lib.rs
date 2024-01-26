use itertools::Itertools;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;


pub struct Seq {
    seq: String,
}

impl Seq {
    pub fn from_str(seq: String) -> Self {
        Self { seq }
    }

    pub fn count_kmers(&self, k: usize) -> HashMap<String, u64> {
        let mut kmer_counter = HashMap::new();
        for chunk in &self.seq.chars().chunks(k) {
            let kmer: String = chunk.collect();

            if !kmer_counter.contains_key(&kmer) {
                kmer_counter.insert(kmer.clone(), 0);
            }

            *kmer_counter.get_mut(&kmer).unwrap() += 1;
        }

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
            .map(Seq::get_nuc_pair)
            .join("")
            .chars()
            .rev()
            .join("")
    }
}

pub fn counting_dna_nucleotides(input_file: PathBuf) -> String {
    let seq = fs::read_to_string(input_file).expect("File not found");
    let kmer_counter = Seq::from_str(seq).count_kmers(1);
    let result = ["A", "C", "G", "T"]
        .map(|x| kmer_counter.get(x).unwrap().to_string())
        .join(" ");

    result
}

pub fn transcribing_dna_into_rna(input_file: PathBuf) -> String {
    let content = fs::read_to_string(input_file).expect("File not found !.");
    Seq::from_str(content)
	.seq
	.replace('T', "U")
}
pub fn complementing_a_strand_of_dna(input_file: PathBuf) -> String {
    let content = fs::read_to_string(input_file)
        .expect("File not found !.")
        .trim()
        .to_string();
    let seq = Seq::from_str(content);

    seq.reverse_complement_dna()
}

pub fn rabbits_and_recurrence_relations(input_file: PathBuf) -> String {
    let (n, k): (usize, usize) = fs::read_to_string(input_file)
        .expect("File not found")
        .trim()
        .split(' ')
        .map(|x| x.parse::<usize>().unwrap())
        .collect_tuple()
        .unwrap();

    // 1    2    3    4    5
    // 1b0a 0b1a 3b1a 3b4a 12b7a
    // 1    1    4    7    19
    let mut b = 1;
    let mut a = 0;

    for _ in 2..n + 1 {
        let temp_b = b;
        b = a * k;
        a += temp_b;
    }

    (a + b).to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_tmp_file_from_string(content: String) -> NamedTempFile {
	let mut temp_file = NamedTempFile::new().expect("Could not create new file !");
	writeln!(temp_file, "{}", content).expect("Could not write to to tmp file.");
	temp_file
    }

    #[test]
    fn test_count_kmers() {
        let seq = Seq::from_str("ACTG".to_string());
        let kmer_counter = seq.count_kmers(1);
        assert_eq!(1, *kmer_counter.get("A").unwrap())
    }

    #[test]
    fn test_reverse_complement() {
        let seq = Seq::from_str(String::from("AAAACCCGGT"));
        assert_eq!(String::from("ACCGGGTTTT"), seq.reverse_complement_dna())
    }

    #[test]
    fn test_counting_dna_nucleotides() {
	let strand = String::from("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC");
	let temp_file = create_tmp_file_from_string(strand);
	let result = counting_dna_nucleotides(temp_file.into_temp_path().to_path_buf());
	let expected = String::from("20 12 17 21");
	assert_eq!(expected, result)
    }
    
    #[test]
    fn test_transcribing_dna_into_rna() {
	let expected = String::from("GAUGGAACUUGACUACGUAAAUU");
	let temp_file = create_tmp_file_from_string(String::from("GATGGAACTTGACTACGTAAATT"));
	let result = transcribing_dna_into_rna(temp_file.into_temp_path().to_path_buf());
	assert_eq!(expected, result.trim())
    }
    
    #[test]
    fn test_complementing_a_strand_of_dna() {
	let expected = String::from("ACCGGGTTTT");
	let temp_file = create_tmp_file_from_string(String::from("AAAACCCGGT"));
	let result = complementing_a_strand_of_dna(temp_file.into_temp_path().to_path_buf());
	assert_eq!(expected, result)
    }

    #[test]
    fn test_rabbits_and_recurrence_relations() {
	let temp_file = create_tmp_file_from_string(String::from("5 3"));
	let result = rabbits_and_recurrence_relations(temp_file.into_temp_path().to_path_buf())
	    .trim()
	    .parse::<usize>()
	    .unwrap();
	let expected: usize = 19;
	assert_eq!(result, expected)
    }
}
