use std::{ffi::c_int, marker::PhantomData};

use crate::{
    abpoa_sys::{
        abpoa_add_graph_alignment, abpoa_align_sequence_to_graph, abpoa_free,
        abpoa_generate_consensus, abpoa_generate_rc_msa, abpoa_init, abpoa_para_t,
        abpoa_post_set_para, abpoa_res_t, abpoa_reset,
    },
    IDX2NT, SEQ_NT4_TABLE,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignMode {
    GLOBAL,
    LOCAL,
    EXTEND,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqType {
    DNA,
    AminoAcid,
}

extern "C" {
    pub fn free(ptr: *mut u64);
}

#[derive(Debug)]
pub struct MsaResult {
    n_seq: i32,
    n_cons: i32,

    clu_n_seq: Vec<i32>,
    clu_read_ids: Vec<Vec<i32>>,
    cons_len: Vec<i32>,
    cons_seq: Vec<String>,
    cons_cov: Vec<Vec<i32>>,
    msa_len: i32,
    msa_seq: Vec<String>, 
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

#[derive(Debug)]
pub struct AbpoaParam {
    pub match_score: i32,
    pub mismatch_score: i32,

    pub gap_open1: i32,
    pub gap_open2: i32,

    pub gap_ext1: i32,
    pub gap_ext2: i32,

    pub wb: i32,
    pub wf: f32,
    pub cons_algrm: String,

    pub out_consensus: bool,
    pub out_msa: bool,
    pub max_n_cons: i32,
    pub min_freq: f64,

    align_mode: AlignMode,
    seq_type: SeqType,
    mat: Vec<i32>,
    _marker: PhantomData<*const ()>,
}

impl Default for AbpoaParam {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: 4,
            gap_open1: 4,
            gap_open2: 24,
            gap_ext1: 2,
            gap_ext2: 1,
            wb: 10,
            wf: 0.01,
            cons_algrm: "HB".to_string(),

            out_consensus: false,
            out_msa: true,
            max_n_cons: 1,
            min_freq: 0.25,

            align_mode: AlignMode::GLOBAL,
            seq_type: SeqType::DNA,
            mat: vec![0i32; 5 * 5],
            _marker: PhantomData,
        }
    }
}

impl AbpoaParam {
    pub fn set_align_mode(&mut self, align_mode: AlignMode) -> &mut Self {
        self.align_mode = align_mode;
        self
    }

    pub fn set_seq_type(&mut self, seq_type: SeqType) -> &mut Self {
        self.seq_type = seq_type;
        match self.seq_type {
            SeqType::DNA => {
                if self.mat.len() != 5 * 5 {
                    self.mat = vec![0i32; 5 * 5];
                }
            }
            SeqType::AminoAcid => {
                if self.mat.len() != 27 * 27 {
                    self.mat = vec![0i32; 27 * 27];
                }
            }
        }
        self
    }

    fn to_abpoa_para_t(&self) -> abpoa_para_t {
        let mut abpoa_para: abpoa_para_t = unsafe { std::mem::zeroed() };

        abpoa_para.match_ = self.match_score;
        abpoa_para.mismatch = self.mismatch_score;
        abpoa_para.gap_open1 = self.gap_open1;
        abpoa_para.gap_open2 = self.gap_open2;
        abpoa_para.gap_ext1 = self.gap_ext1;
        abpoa_para.gap_ext2 = self.gap_ext2;

        abpoa_para.wb = self.wb;
        abpoa_para.wf = self.wf;

        abpoa_para.zdrop = -1;
        abpoa_para.end_bonus = -1;

        abpoa_para.set_disable_seeding(1);
        abpoa_para.set_progressive_poa(0);
        abpoa_para.set_ret_cigar(1);
        abpoa_para.set_use_qv(0);

        match self.cons_algrm.to_uppercase().as_str() {
            "MF" => abpoa_para.cons_algrm = 1,
            "HB" => abpoa_para.cons_algrm = 0,
            algrm => panic!("Unknown conseneus calling mode: {}", algrm),
        }

        match self.align_mode {
            AlignMode::GLOBAL => abpoa_para.align_mode = 0,
            AlignMode::LOCAL => abpoa_para.align_mode = 1,
            AlignMode::EXTEND => abpoa_para.align_mode = 2,
        }

        match self.seq_type {
            SeqType::DNA => abpoa_para.m = 5,
            SeqType::AminoAcid => abpoa_para.m = 27,
        }

        abpoa_para.set_out_cons(if self.out_consensus {1} else {0});
        abpoa_para.set_out_msa(if self.out_msa {1} else {0});


        abpoa_para.max_n_cons = self.max_n_cons;
        abpoa_para.min_freq = self.min_freq;

        abpoa_para.incr_fn = std::ptr::null_mut();
        abpoa_para.out_pog = std::ptr::null_mut();

        abpoa_para.mat = self.mat.as_ptr() as *mut i32;

        // println!("{:?}", abpoa_para);
        unsafe {
            abpoa_post_set_para(&mut abpoa_para as *mut abpoa_para_t);
        }
        // println!("{:?}", abpoa_para);


        abpoa_para
    }
}


pub fn msa(param: &AbpoaParam, seqs: &Vec<&str>) -> Option<MsaResult> {
    let n_seqs = seqs.len();
    if n_seqs == 0 {
        return None;
    }

    let seq_ele2idx = match param.seq_type {
        SeqType::DNA => &SEQ_NT4_TABLE,
        _ => panic!("not implement yet"),
    };

    let idx2seq_ele = match param.seq_type {
        SeqType::DNA => &IDX2NT,
        _ => panic!("not implement yet"),
    };

    let mut abpoa_param = param.to_abpoa_para_t();

    let res = unsafe {
        let ab = abpoa_init();
        abpoa_reset(ab, &mut abpoa_param, seqs[0].len() as c_int);
        let abs = &mut *(*ab).abs;
        abs.n_seq = n_seqs as i32;

        let mut abpoa_res: abpoa_res_t = std::mem::zeroed();

        seqs.iter().enumerate().for_each(|(read_idx, seq)| {
            let mut seq_encoded = seq
                .as_bytes()
                .iter()
                .map(|base| seq_ele2idx[*base as usize])
                .collect::<Vec<_>>();
            // println!("seq-encoded:{:?}", seq_encoded);
            abpoa_res.n_cigar = 0;
            // println!("before:{:?}", abpoa_res);

            abpoa_align_sequence_to_graph(
                ab,
                &mut abpoa_param,
                seq_encoded.as_mut_ptr(),
                seq_encoded.len() as c_int,
                &mut abpoa_res,
            );
            // println!("after: {:?}", abpoa_res);
            abpoa_add_graph_alignment(
                ab,
                &mut abpoa_param,
                seq_encoded.as_mut_ptr(),
                std::ptr::null_mut(),
                seq_encoded.len() as c_int,
                std::ptr::null_mut(),
                abpoa_res,
                read_idx as c_int,
                n_seqs as c_int,
                1,
            );

            if abpoa_res.n_cigar > 0 {
                free(abpoa_res.graph_cigar);
            }
        });


        if abpoa_param.out_msa() == 1 {
            abpoa_generate_rc_msa(ab, &mut abpoa_param);
        } else if abpoa_param.out_cons() == 1 {
            abpoa_generate_consensus(ab, &mut abpoa_param);
        }

        let abc = &mut *(*ab).abc;
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
        for i in 0..n_cons {
            clu_n_seq.push(*abc.clu_n_seq.add(i as usize));
            cons_len.push(*abc.cons_len.add(i as usize));

            let (mut clu_read_ids1, mut cons_seq1, mut cons_cov1) = (vec![], String::new(), vec![]);
            for j in 0..*abc.clu_n_seq.add(i as usize) {
                clu_read_ids1.push(*(*abc.clu_read_ids.add(i as usize)).add(j as usize));
            }

            clu_read_ids.push(clu_read_ids1);

            for j in 0..*abc.cons_len.add(i as usize) {
                let c = *(*abc.cons_base.add(i as usize)).add(j as usize);
                cons_seq1.push(c as char);
                cons_cov1.push(*(*abc.cons_cov.add(i as usize)).add(j as usize));
            }

            cons_seq.push(cons_seq1);
            cons_cov.push(cons_cov1);
        }

        if msa_len > 0 {
            for i in 0..(abc.n_seq + n_cons) {
                let mut msa_seq1 = String::new();
                for j in 0..msa_len {
                    let c = *(*abc.msa_base.add(i as usize)).add(j as usize);
                    msa_seq1.push(idx2seq_ele[c as usize]);
                }

                msa_seq.push(msa_seq1);
            }
        }
        abpoa_free(ab);
        MsaResult::new(
            n_seqs as i32,
            n_cons,
            clu_n_seq,
            clu_read_ids,
            cons_len,
            cons_seq,
            cons_cov,
            msa_len,
            msa_seq,
        )
    };

    Some(res)
}

#[cfg(test)]
mod test {
    use crate::abpoa::{msa, AbpoaParam};

    #[test]
    fn test_poa_msa() {
        let align_param = AbpoaParam::default();
        let seqs = vec!["AAC", "AC", "C"];
        let res = msa(&align_param, &seqs).unwrap();
        res.print_msa();
        assert_eq!(res.msa_seq[0], "AAC");
        assert_eq!(res.msa_seq[1], "-AC");
        assert_eq!(res.msa_seq[2], "--C");
        
        let align_param = AbpoaParam::default();
        let seqs = vec!["AAC", "GGAC", "GC"];
        let res = msa(&align_param, &seqs).unwrap();
        res.print_msa();
        // println!("{:?}", res);
        assert_eq!(res.msa_seq()[0], "-AAC");
        assert_eq!(res.msa_seq()[1], "GGAC");
        assert_eq!(res.msa_seq()[2], "G--C");
    }
}
