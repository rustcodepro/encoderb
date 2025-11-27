use linfa::prelude::*;
use ndarray::{Array1, Array2};
use smartcore::svm::svc::SVC;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

/*
Gaurav Sablok
codeprog@icloud.com
*/

#[derive(Debug, PartialEq, PartialOrd)]
pub struct Sequences {
    pub path: String,
}

impl Sequences {
    pub fn sequences(&self) -> Result<(Vec<String>, Vec<i32>), Box<dyn Error>> {
        let fileopen = File::open(&self.path).expect("file not present");
        let fileread = BufReader::new(fileopen);
        let mut sequences: Vec<String> = Vec::new();
        let mut classlabels: Vec<i32> = Vec::new();
        for i in fileread.lines() {
            let line = i.expect("file not present");
            if line.starts_with("#") {
                let linesplit = line.split("\t").collect::<Vec<_>>()[1];
                classlabels.push(linesplit.parse::<i32>().unwrap());
            } else if !line.starts_with("#") {
                sequences.push(line.clone());
            }
        }
        Ok((sequences, classlabels))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct VariantSeq {
    pub sequences: Vec<Vec<(usize, f32)>>,
    pub labels: Vec<i32>,
    pub kmer_map: HashMap<String, usize>,
    pub feature_size: usize,
}

impl VariantSeq {
    pub fn generatesparse(pathfile: &str, labels: Vec<i32>, k: usize) -> Self {
        let nucleotides = ['A', 'C', 'G', 'T'];
        let mut kmer_map = HashMap::new();
        let mut kmer_idx = 0;
        for kmer in kmerssparse(k, &nucleotides) {
            kmer_map.insert(kmer, kmer_idx);
            kmer_idx += 1;
        }
        let pathseq = Sequences {
            path: pathfile.to_string(),
        };
        let sequencesunclone = pathseq.sequences().unwrap();
        let sequences: Vec<Vec<(usize, f32)>> = sequencesunclone
            .0
            .into_iter()
            .map(|seq| sparsefrequencies(&seq, k, &kmer_map))
            .collect();
        VariantSeq {
            sequences,
            labels,
            kmer_map,
            feature_size: 4usize.pow(k as u32),
        }
    }

    pub fn densematrix(&self, sparse: &[(usize, f32)]) -> Array1<f32> {
        let mut densematrix = Array1::zeros(self.feature_size);
        for &(idx, freq) in sparse {
            densematrix[idx] = freq;
        }
        densematrix
    }

    pub fn dataset(&self) -> (Array2<f32>, Array1<i32>) {
        let features: Array2<f32> = Array2::from_shape_vec(
            (self.sequences.len(), self.feature_size),
            self.sequences
                .iter()
                .flat_map(|sparse| self.densematrix(sparse))
                .collect(),
        )
        .unwrap();
        let targets = Array1::from_vec(self.labels.clone());
        (features, targets)
    }
}

pub fn kmerssparse(kmersize: usize, nucleotides: &[char]) -> Vec<String> {
    let mut kmers = vec!["".to_string()];
    for _ in 0..kmersize {
        let mut finalkmers = vec![];
        for kmer in kmers {
            for &nuc in nucleotides {
                let new_kmer = kmer.clone();
                finalkmers.push(nuc.to_string());
                kmers.push(new_kmer.clone());
            }
        }
        kmers = finalkmers;
    }
    kmers
}

pub fn sparsefrequencies(
    seq: &str,
    k: usize,
    kmer_map: &HashMap<String, usize>,
) -> Vec<(usize, f32)> {
    let mut counts = HashMap::new();
    let seq_len = seq.len();
    if seq_len < k {
        return vec![];
    }
    for i in 0..=(seq_len - k) {
        let kmer = &seq[i..i + k];
        if let Some(&idx) = kmer_map.get(kmer) {
            *counts.entry(idx).or_insert(0.0) += 1.0;
        }
    }
    let total = (seq_len - k + 1) as f32;
    let mut sparse = vec![];
    if total > 0.0 {
        for (idx, count) in counts {
            sparse.push((idx, count / total));
        }
    }
    sparse
}
