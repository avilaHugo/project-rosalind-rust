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

#[cfg(test)]
mod tests {
    use crate::Seq;

    #[test]
    fn test_count_kmers() {
	let seq = Seq::from_str("ACTG".to_string());
	let kmer_counter = seq.count_kmers(1);
	assert_eq!(1, *kmer_counter.get("A").unwrap())
    }
}
