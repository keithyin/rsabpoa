use std::{ffi::c_int, marker::PhantomData, ops::Deref};

use crate::{
    abpoa_result_parser::abpoa_result_parser,
    abpoa_sys::{
        abpoa_add_graph_alignment, abpoa_align_sequence_to_graph, abpoa_cons_t, abpoa_free,
        abpoa_generate_consensus, abpoa_generate_rc_msa, abpoa_init, abpoa_msa, abpoa_msa1,
        abpoa_para_t, abpoa_post_set_para, abpoa_res_t, abpoa_reset, abpoa_t,
    },
    utils::reverse_complement,
    IDX2NT, SEQ_NT4_TABLE,
};

use crate::abpoa_result_parser::MsaResult;

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
pub struct AbPoaResT {
    ab_poa_res: abpoa_res_t,
}

impl Deref for AbPoaResT {
    type Target = abpoa_res_t;
    fn deref(&self) -> &Self::Target {
        &self.ab_poa_res
    }
}

impl From<abpoa_res_t> for AbPoaResT {
    fn from(value: abpoa_res_t) -> Self {
        Self { ab_poa_res: value }
    }
}

impl Drop for AbPoaResT {
    fn drop(&mut self) {
        if self.ab_poa_res.n_cigar > 0 {
            unsafe {
                free(self.ab_poa_res.graph_cigar);
            };
        }
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
    pub fn channel_draft_param() -> Self {
        let mut align_param = Self::default();
        align_param.match_score = 2;
        align_param.mismatch_score = 5;
        align_param.gap_open1 = 2;
        align_param.gap_open2 = 24;
        align_param.gap_ext1 = 1;
        align_param.gap_ext2 = 0;
        align_param.out_msa = false;
        align_param.out_consensus = true;
        align_param.set_align_mode(AlignMode::LOCAL);
        align_param
    }

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

    // don't call abpoa_free_para, the mat is owned by rust!!!
    pub fn to_abpoa_para_t(&self) -> abpoa_para_t {
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
        abpoa_para.set_amb_strand(1);

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

        abpoa_para.set_out_cons(if self.out_consensus { 1 } else { 0 });
        abpoa_para.set_out_msa(if self.out_msa { 1 } else { 0 });

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
            abpoa_res.n_matched_bases = 0;

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
            // println!("{}", abpoa_res.n_cigar);
            println!("{}", abpoa_res.n_matched_bases);

            if abpoa_res.n_cigar > 0 {
                free(abpoa_res.graph_cigar);
            }
        });

        if abpoa_param.out_msa() == 1 {
            abpoa_generate_rc_msa(ab, &mut abpoa_param);
        } else if abpoa_param.out_cons() == 1 {
            abpoa_generate_consensus(ab, &mut abpoa_param);
        }

        let abc: &mut abpoa_cons_t = &mut *(*ab).abc;
        let res = abpoa_result_parser(abc, idx2seq_ele, &abpoa_param);
        abpoa_free(ab);
        res
    };

    res
}

/// seqs order matters
///
pub fn msa_with_adaptive_fwd_rev(param: &AbpoaParam, seqs: &Vec<&str>) -> Option<MsaResult> {
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

        seqs.iter().enumerate().for_each(|(read_idx, seq)| {
            let mut fwd_seq_encoded = seq
                .as_bytes()
                .iter()
                .map(|base| seq_ele2idx[*base as usize])
                .collect::<Vec<_>>();

            let mut rev_seq_encoded = reverse_complement(seq.as_bytes())
                .iter()
                .map(|base| seq_ele2idx[*base as usize])
                .collect::<Vec<_>>();

            let fwd_res = abpoa_align_sequence_to_graph_with_adaptive_fwd_rev(
                ab,
                &mut abpoa_param,
                fwd_seq_encoded.as_mut_ptr(),
                fwd_seq_encoded.len() as c_int,
            );
            let rev_res = abpoa_align_sequence_to_graph_with_adaptive_fwd_rev(
                ab,
                &mut abpoa_param,
                rev_seq_encoded.as_mut_ptr(),
                rev_seq_encoded.len() as c_int,
            );

            let (mut final_seq, final_res) = if fwd_res.n_matched_bases > rev_res.n_matched_bases {
                (fwd_seq_encoded, fwd_res)
            } else {
                (rev_seq_encoded, rev_res)
            };

            abpoa_add_graph_alignment(
                ab,
                &mut abpoa_param,
                final_seq.as_mut_ptr(),
                std::ptr::null_mut(),
                final_seq.len() as c_int,
                std::ptr::null_mut(),
                final_res.ab_poa_res,
                read_idx as c_int,
                n_seqs as c_int,
                1,
            );
        });

        if abpoa_param.out_msa() == 1 {
            abpoa_generate_rc_msa(ab, &mut abpoa_param);
        } else if abpoa_param.out_cons() == 1 {
            abpoa_generate_consensus(ab, &mut abpoa_param);
        }

        let abc = &mut *(*ab).abc;
        let res = abpoa_result_parser(abc, idx2seq_ele, &abpoa_param);
        abpoa_free(ab);
        res
    };

    res
}

unsafe fn abpoa_align_sequence_to_graph_with_adaptive_fwd_rev(
    ab: *mut abpoa_t,
    abpt: *mut abpoa_para_t,
    query: *mut u8,
    qlen: ::std::os::raw::c_int,
) -> AbPoaResT {
    let mut abpoa_res: abpoa_res_t = std::mem::zeroed();
    abpoa_align_sequence_to_graph(ab, abpt, query, qlen, &mut abpoa_res);
    abpoa_res.into()
}

pub fn abpoa_consensus_dna_core(
    ab: *mut abpoa_t,
    abpt: *mut abpoa_para_t,
    seqs: &mut Vec<&mut [u8]>,
) -> Option<MsaResult> {
    if seqs.len() == 0 {
        return None;
    }
    let seq_ele2idx = &SEQ_NT4_TABLE;

    let n_seqs = seqs.len();
    let mut seq_lens = seqs
        .iter()
        .map(|seq| seq.len() as c_int)
        .collect::<Vec<_>>();
    let mut seqs = seqs
        .iter_mut()
        .map(|seq| {
            seq.iter()
                .map(|&base| seq_ele2idx[base as usize])
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut seqs = seqs
        .iter_mut()
        .map(|seq_encoded| seq_encoded.as_mut_ptr())
        .collect::<Vec<_>>();

    unsafe {
        abpoa_msa(
            ab,
            abpt,
            n_seqs as c_int,
            std::ptr::null_mut(),
            seq_lens.as_mut_ptr(),
            seqs.as_mut_ptr(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        );
        let abc = &*(*ab).abc;

        let res = abpoa_result_parser(abc, &IDX2NT, &(*abpt));

        return res;
    }
}

pub fn abpoa_consensus_dna(param: &AbpoaParam, seqs: &mut Vec<&mut [u8]>) -> Option<MsaResult> {
    if seqs.len() == 0 {
        return None;
    }
    let seq_ele2idx = &SEQ_NT4_TABLE;

    let idx2seq_ele = &IDX2NT;

    let n_seqs = seqs.len();
    let mut seq_lens = seqs
        .iter()
        .map(|seq| seq.len() as c_int)
        .collect::<Vec<_>>();
    let mut seqs = seqs
        .iter_mut()
        .map(|seq| {
            seq.iter()
                .map(|&base| seq_ele2idx[base as usize])
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut seqs = seqs
        .iter_mut()
        .map(|seq_encoded| seq_encoded.as_mut_ptr())
        .collect::<Vec<_>>();

    unsafe {
        let mut abpt = param.to_abpoa_para_t();

        let ab = abpoa_init();
        abpoa_msa(
            ab,
            &mut abpt,
            n_seqs as c_int,
            std::ptr::null_mut(),
            seq_lens.as_mut_ptr(),
            seqs.as_mut_ptr(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        );
        let abc = &mut *(*ab).abc;
        let result = abpoa_result_parser(abc, idx2seq_ele, &abpt);
        abpoa_free(ab);
        return result;
    }
}

pub fn parse_consensus_seq_from_abpoa_cons_t_dna(abc: &abpoa_cons_t) -> Option<String> {
    let n_cons = abc.n_cons;
    if n_cons <= 0 {
        return None;
    }
    let mut cons_seq = String::new();
    let idx2seq_ele = &IDX2NT;

    unsafe {
        let max_cluster_idx = (0..n_cons)
            .into_iter()
            .map(|i| i as usize)
            .map(|i| (i, *abc.clu_n_seq.add(i)))
            .max_by_key(|v| v.1)
            .unwrap()
            .0;

        for j in 0..*abc.cons_len.add(max_cluster_idx) {
            let c = *(*abc.cons_base.add(max_cluster_idx)).add(j as usize);
            cons_seq.push(idx2seq_ele[c as usize]);
        }

        return Some(cons_seq);
    }
}

#[cfg(test)]
mod test {
    use crate::{
        abpoa::{abpoa_consensus_dna_core, msa, AbpoaParam},
        abpoa_sys::{abpoa_free, abpoa_init},
    };

    use super::{abpoa_consensus_dna, msa_with_adaptive_fwd_rev};

    fn get_seqs1() -> Vec<String> {
        let seqs = vec![
            "AAGAAAAAG",
            "AATGAAAAAG",
            "AAGAAAAAG",
            "AAGAAAAG",
            "TAGAAAAAAAAAAAAG",
            "CTTTTTTTTTTTTCTA",
            "AGAAAAG",
            "AAAGAAAAG",
        ];
        seqs.into_iter().map(|v| v.to_string()).collect::<Vec<_>>()
    }

    fn get_seqs2() -> Vec<String> {
        let seqs = vec![
            "AAGAAAAAG",
            "AATGAAAAAG",
            "AAGAAAAAG",
            "AAGAAAAG",
            "TAGAAAAAAAAAAAAG",
            "CTTTTTTTTTTTTCTA",
            "AGAAAAG",
            "AAAGAAAAG",
            "CTTTTTTTT",
        ];
        seqs.into_iter().map(|v| v.to_string()).collect::<Vec<_>>()
    }

    fn get_seqs3() -> Vec<String> {
        let seqs = vec![
            "AAGAAAAAG",
            "AATGAAAAAG",
            "AAGAAAAAG",
            "AAGAAAAG",
            "TAGAAAAAAAAAAAAG",
            "CTTTTTTTTTTTTCTA",
            "AGAAAAG",
        ];
        seqs.into_iter().map(|v| v.to_string()).collect::<Vec<_>>()
    }

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

    #[test]
    fn test_poa_msa2() {
        let mut align_param = AbpoaParam::default();
        align_param.mismatch_score = 6;
        align_param.gap_open1 = 2;
        align_param.gap_open2 = 24;
        align_param.gap_ext1 = 1;
        align_param.gap_ext2 = 0;

        let seqs = vec![
            "AAAAAGG",
            "AAAAAGG",
            "AAAAAGG",
            "AAAAAGG",
            "AAAAAGG",
            "AAAAAAG",
            "AAAAAAG",
            "AAAAAAG",
            "AAAAAAAAGG",
        ];
        let res = msa(&align_param, &seqs).unwrap();
        res.print_msa();
    }

    #[test]
    fn test_poa_consensus() {
        let mut align_param = AbpoaParam::default();
        align_param.mismatch_score = 6;
        align_param.gap_open1 = 2;
        align_param.gap_open2 = 24;
        align_param.gap_ext1 = 1;
        align_param.gap_ext2 = 0;
        align_param.out_consensus = true;

        let seqs = vec![
            "AAGAAAAAG",
            "AATGAAAAAG",
            "AAGAAAAAG",
            "AAGAAAAG",
            "TAGAAAAAAAAAAAAG",
            "AGAAAAG",
            "AAAGAAAAG",
        ];
        let res = msa(&align_param, &seqs).unwrap();
        res.print_msa();
    }

    #[test]
    fn test_poa_consensus2() {
        let mut align_param = AbpoaParam::default();
        align_param.match_score = 2;
        align_param.mismatch_score = 5;
        align_param.gap_open1 = 2;
        align_param.gap_open2 = 24;
        align_param.gap_ext1 = 1;
        align_param.gap_ext2 = 0;
        align_param.out_consensus = true;

        let seqs = vec![
            "AAGAAAAAG",
            "AATGAAAAAG",
            "AAGAAAAAG",
            "AAGAAAAG",
            "TAGAAAAAAAAAAAAG",
            "CTTTTTTTTTTTTCTA",
            "AGAAAAG",
            "AAAGAAAAG",
        ];
        let res = msa(&align_param, &seqs).unwrap();
        res.print_msa();
    }

    #[test]
    fn test_poa_consensus_adaptive() {
        let mut align_param = AbpoaParam::default();
        align_param.match_score = 2;
        align_param.mismatch_score = 5;
        align_param.gap_open1 = 2;
        align_param.gap_open2 = 24;
        align_param.gap_ext1 = 1;
        align_param.gap_ext2 = 0;
        align_param.out_consensus = true;

        let seqs = vec![
            "AAGAAAAAG",
            "AATGAAAAAG",
            "AAGAAAAAG",
            "AAGAAAAG",
            "TAGAAAAAAAAAAAAG",
            "CTTTTTTTTTTTTCTA",
            "AGAAAAG",
            "AAAGAAAAG",
        ];
        let res = msa_with_adaptive_fwd_rev(&align_param, &seqs).unwrap();
        res.print_msa();
    }

    #[test]
    fn test_poa_consensus_a() {
        let mut align_param = AbpoaParam::default();
        align_param.match_score = 2;
        align_param.mismatch_score = 5;
        align_param.gap_open1 = 2;
        align_param.gap_open2 = 24;
        align_param.gap_ext1 = 1;
        align_param.gap_ext2 = 0;
        align_param.out_consensus = true;

        let mut seqs = get_seqs1();
        let mut seqs = seqs
            .iter_mut()
            .map(|v| unsafe { v.as_bytes_mut() })
            .collect::<Vec<_>>();
        let res = abpoa_consensus_dna(&align_param, &mut seqs).unwrap();
        res.print_msa();
        println!("{:?}", res.cons_seq());
    }

    #[test]
    fn test_poa_consensus_b() {
        let mut align_param = AbpoaParam::default();
        align_param.match_score = 2;
        align_param.mismatch_score = 5;
        align_param.gap_open1 = 2;
        align_param.gap_open2 = 24;
        align_param.gap_ext1 = 1;
        align_param.gap_ext2 = 0;
        align_param.out_consensus = true;

        let mut seqs = get_seqs1();
        let mut seqs = seqs
            .iter_mut()
            .map(|v| unsafe { v.as_bytes_mut() })
            .collect::<Vec<_>>();
        let ab = unsafe { abpoa_init() };
        let mut abpt = align_param.to_abpoa_para_t();
        let res = abpoa_consensus_dna_core(ab, &mut abpt, &mut seqs).unwrap();
        res.print_msa();
        println!("{:?}", res.cons_seq());

        let mut seqs = get_seqs2();
        let mut seqs = seqs
            .iter_mut()
            .map(|v| unsafe { v.as_bytes_mut() })
            .collect::<Vec<_>>();
        let res = abpoa_consensus_dna_core(ab, &mut abpt, &mut seqs).unwrap();
        res.print_msa();
        println!("{:?}", res.cons_seq());

        let mut seqs = get_seqs3();
        let mut seqs = seqs
            .iter_mut()
            .map(|v| unsafe { v.as_bytes_mut() })
            .collect::<Vec<_>>();
        let res = abpoa_consensus_dna_core(ab, &mut abpt, &mut seqs).unwrap();
        res.print_msa();
        unsafe { abpoa_free(ab) };
        println!("{:?}", res.cons_seq());
    }
}
