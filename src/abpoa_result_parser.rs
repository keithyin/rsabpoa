use crate::abpoa_sys::{abpoa_cons_t, abpoa_para_t};

#[derive(Debug)]
pub struct MsaResult {
    pub n_seq: i32,  // origin seq num
    pub n_cons: i32, // consensus num

    pub clu_n_seq: Vec<i32>,         // origin seq num in n-cluster
    pub clu_read_ids: Vec<Vec<i32>>, // read_ids in n-cluster
    pub cons_len: Vec<i32>,
    pub cons_seq: Vec<String>,
    pub cons_cov: Vec<Vec<i32>>,
    pub msa_len: i32,
    pub msa_seq: Vec<String>,
}

impl MsaResult {
    /// n_seqs, n_cons, clu_n_seq, clu_read_ids, cons_len, cons_seq, cons_cov, msa_len, msa_seq
    pub fn new(
        n_seq: i32,
        n_cons: i32,
        clu_n_seq: Vec<i32>,
        clu_read_ids: Vec<Vec<i32>>,
        cons_len: Vec<i32>,
        cons_seq: Vec<String>,
        cons_cov: Vec<Vec<i32>>,
        msa_len: i32,
        msa_seq: Vec<String>,
    ) -> Self {
        Self {
            n_seq,
            n_cons,
            clu_n_seq,
            clu_read_ids,
            cons_len,
            cons_seq,
            cons_cov,
            msa_len,
            msa_seq,
        }
    }

    pub fn n_seq(&self) -> i32 {
        self.n_seq
    }

    pub fn n_cons(&self) -> i32 {
        self.n_cons
    }

    pub fn clu_n_seq(&self) -> &Vec<i32> {
        &self.clu_n_seq
    }

    pub fn clu_read_ids(&self) -> &Vec<Vec<i32>> {
        &self.clu_read_ids
    }

    pub fn cons_len(&self) -> &Vec<i32> {
        &self.cons_len
    }

    pub fn cons_seq(&self) -> &Vec<String> {
        &self.cons_seq
    }

    pub fn cons_cov(&self) -> &Vec<Vec<i32>> {
        &self.cons_cov
    }

    pub fn msa_len(&self) -> i32 {
        self.msa_len
    }

    pub fn msa_seq(&self) -> &Vec<String> {
        &self.msa_seq
    }

    pub fn first_consensus_seq(&self) -> Option<&str> {
        if self.msa_seq.len() == self.n_seq as usize {
            None
        } else {
            Some(self.msa_seq[self.n_seq as usize].as_str())
        }
    }

    pub fn max_len_consensus_seq(&self) -> Option<&str> {
        self.cons_seq
            .iter()
            .map(|v| v.as_str())
            .max_by_key(|v| v.len())
    }

    pub fn max_num_of_seq_consensus_seq(&self) -> Option<&str> {
        let idx = self
            .clu_read_ids
            .iter()
            .enumerate()
            .map(|(idx, ids)| (idx, ids.len()))
            .max_by_key(|v| v.1);
        return if let Some(idx) = idx {
            Some(self.cons_seq[idx.0].as_str())
        } else {
            None
        };
    }

    pub fn print_msa(&self) {
        if self.msa_seq.len() == 0 {
            return;
        }

        self.msa_seq.iter().enumerate().for_each(|(idx, seq)| {
            if idx < self.n_seq as usize {
                print!("seq_{}  ", idx + 1);
            } else {
                let cons_id = if self.n_cons > 1 {
                    format!("_{}", idx - self.n_seq as usize + 1)
                } else {
                    "".to_string()
                };
                print!("consensus{}", cons_id);
            }
            println!("{}", seq);
        });
    }
}

pub fn abpoa_result_parser(
    abc: &abpoa_cons_t,
    idx2seq_ele: &Vec<char>,
    abpt: &abpoa_para_t,
) -> Option<MsaResult> {
    unsafe {
        let (
            n_cons,
            mut clu_n_seq,
            mut clu_read_ids,
            mut cons_len,
            mut cons_seq,
            mut cons_cov,
            msa_len,
            mut msa_seq,
        ) = (
            abc.n_cons,
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            abc.msa_len,
            vec![],
        );

        if n_cons == 0 {
            return None;
        }
        // there might be multiple result consensus sequences
        // clu_n_seq: number of subreads in the n-cluster
        // cons_len: length of consensus sequences

        if abpt.out_cons() == 1 {
            for i in 0..n_cons {
                clu_n_seq.push(*abc.clu_n_seq.add(i as usize));
                cons_len.push(*abc.cons_len.add(i as usize));

                let (mut clu_read_ids1, mut cons_seq1, mut cons_cov1) =
                    (vec![], String::new(), vec![]);
                for j in 0..*abc.clu_n_seq.add(i as usize) {
                    clu_read_ids1.push(*(*abc.clu_read_ids.add(i as usize)).add(j as usize));
                }

                clu_read_ids.push(clu_read_ids1);

                for j in 0..*abc.cons_len.add(i as usize) {
                    let c = *(*abc.cons_base.add(i as usize)).add(j as usize);
                    cons_seq1.push(idx2seq_ele[c as usize]);
                    cons_cov1.push(*(*abc.cons_cov.add(i as usize)).add(j as usize));
                }

                cons_seq.push(cons_seq1);
                cons_cov.push(cons_cov1);
            }
        }

        if msa_len > 0 && abpt.out_msa() == 1 {
            for i in 0..(abc.n_seq + n_cons) {
                let mut msa_seq1 = String::new();
                for j in 0..msa_len {
                    let c = *(*abc.msa_base.add(i as usize)).add(j as usize);
                    msa_seq1.push(idx2seq_ele[c as usize]);
                }

                msa_seq.push(msa_seq1);
            }
        }

        Some(MsaResult::new(
            abc.n_seq,
            n_cons,
            clu_n_seq,
            clu_read_ids,
            cons_len,
            cons_seq,
            cons_cov,
            msa_len,
            msa_seq,
        ))
    }
}
