use std::{ffi::c_int, marker::PhantomData, ops::Deref};

use crate::{
    abpoa_result_parser::abpoa_result_parser,
    abpoa_sys::{
        abpoa_add_graph_alignment, abpoa_align_sequence_to_graph, abpoa_cons_t, abpoa_free,
        abpoa_generate_consensus, abpoa_generate_rc_msa, abpoa_init, abpoa_msa, abpoa_msa1,
        abpoa_para_t, abpoa_post_set_para, abpoa_res_t, abpoa_reset, abpoa_t,
    },
    utils::reverse_complement,
    wrapper::Abpoa,
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

// pub fn msa(param: &AbpoaParam, seqs: &Vec<&str>) -> Option<MsaResult> {
//     let n_seqs = seqs.len();
//     if n_seqs == 0 {
//         return None;
//     }

//     let seq_ele2idx = match param.seq_type {
//         SeqType::DNA => &SEQ_NT4_TABLE,
//         _ => panic!("not implement yet"),
//     };

//     let idx2seq_ele = match param.seq_type {
//         SeqType::DNA => &IDX2NT,
//         _ => panic!("not implement yet"),
//     };

//     let mut abpoa_param = param.to_abpoa_para_t();

//     let res = unsafe {
//         let ab = abpoa_init();
//         abpoa_reset(ab, &mut abpoa_param, seqs[0].len() as c_int);
//         let abs = &mut *(*ab).abs;
//         abs.n_seq = n_seqs as i32;

//         let mut abpoa_res: abpoa_res_t = std::mem::zeroed();

//         seqs.iter().enumerate().for_each(|(read_idx, seq)| {
//             let mut seq_encoded = seq
//                 .as_bytes()
//                 .iter()
//                 .map(|base| seq_ele2idx[*base as usize])
//                 .collect::<Vec<_>>();
//             // println!("seq-encoded:{:?}", seq_encoded);
//             abpoa_res.n_cigar = 0;
//             abpoa_res.n_matched_bases = 0;

//             // println!("before:{:?}", abpoa_res);

//             abpoa_align_sequence_to_graph(
//                 ab,
//                 &mut abpoa_param,
//                 seq_encoded.as_mut_ptr(),
//                 seq_encoded.len() as c_int,
//                 &mut abpoa_res,
//             );
//             // println!("after: {:?}", abpoa_res);
//             abpoa_add_graph_alignment(
//                 ab,
//                 &mut abpoa_param,
//                 seq_encoded.as_mut_ptr(),
//                 std::ptr::null_mut(),
//                 seq_encoded.len() as c_int,
//                 std::ptr::null_mut(),
//                 abpoa_res,
//                 read_idx as c_int,
//                 n_seqs as c_int,
//                 1,
//             );
//             // println!("{}", abpoa_res.n_cigar);
//             println!("{}", abpoa_res.n_matched_bases);

//             if abpoa_res.n_cigar > 0 {
//                 free(abpoa_res.graph_cigar);
//             }
//         });

//         if abpoa_param.out_msa() == 1 {
//             abpoa_generate_rc_msa(ab, &mut abpoa_param);
//         } else if abpoa_param.out_cons() == 1 {
//             abpoa_generate_consensus(ab, &mut abpoa_param);
//         }

//         let abc: &mut abpoa_cons_t = &mut *(*ab).abc;
//         let res = abpoa_result_parser(abc, idx2seq_ele, &abpoa_param);
//         abpoa_free(ab);
//         res
//     };

//     res
// }

/// seqs order matters
///
// pub fn msa_with_adaptive_fwd_rev(param: &AbpoaParam, seqs: &Vec<&str>) -> Option<MsaResult> {
//     let n_seqs = seqs.len();
//     if n_seqs == 0 {
//         return None;
//     }

//     let seq_ele2idx = match param.seq_type {
//         SeqType::DNA => &SEQ_NT4_TABLE,
//         _ => panic!("not implement yet"),
//     };

//     let idx2seq_ele = match param.seq_type {
//         SeqType::DNA => &IDX2NT,
//         _ => panic!("not implement yet"),
//     };

//     let mut abpoa_param = param.to_abpoa_para_t();

//     let res = unsafe {
//         let ab = abpoa_init();
//         abpoa_reset(ab, &mut abpoa_param, seqs[0].len() as c_int);
//         let abs = &mut *(*ab).abs;
//         abs.n_seq = n_seqs as i32;

//         seqs.iter().enumerate().for_each(|(read_idx, seq)| {
//             let mut fwd_seq_encoded = seq
//                 .as_bytes()
//                 .iter()
//                 .map(|base| seq_ele2idx[*base as usize])
//                 .collect::<Vec<_>>();

//             let mut rev_seq_encoded = reverse_complement(seq.as_bytes())
//                 .iter()
//                 .map(|base| seq_ele2idx[*base as usize])
//                 .collect::<Vec<_>>();

//             let fwd_res = abpoa_align_sequence_to_graph_with_adaptive_fwd_rev(
//                 ab,
//                 &mut abpoa_param,
//                 fwd_seq_encoded.as_mut_ptr(),
//                 fwd_seq_encoded.len() as c_int,
//             );
//             let rev_res = abpoa_align_sequence_to_graph_with_adaptive_fwd_rev(
//                 ab,
//                 &mut abpoa_param,
//                 rev_seq_encoded.as_mut_ptr(),
//                 rev_seq_encoded.len() as c_int,
//             );

//             let (mut final_seq, final_res) = if fwd_res.n_matched_bases > rev_res.n_matched_bases {
//                 (fwd_seq_encoded, fwd_res)
//             } else {
//                 (rev_seq_encoded, rev_res)
//             };

//             abpoa_add_graph_alignment(
//                 ab,
//                 &mut abpoa_param,
//                 final_seq.as_mut_ptr(),
//                 std::ptr::null_mut(),
//                 final_seq.len() as c_int,
//                 std::ptr::null_mut(),
//                 final_res.ab_poa_res,
//                 read_idx as c_int,
//                 n_seqs as c_int,
//                 1,
//             );
//         });

//         if abpoa_param.out_msa() == 1 {
//             abpoa_generate_rc_msa(ab, &mut abpoa_param);
//         } else if abpoa_param.out_cons() == 1 {
//             abpoa_generate_consensus(ab, &mut abpoa_param);
//         }

//         let abc = &mut *(*ab).abc;
//         let res = abpoa_result_parser(abc, idx2seq_ele, &abpoa_param);
//         abpoa_free(ab);
//         res
//     };

//     res
// }

// unsafe fn abpoa_align_sequence_to_graph_with_adaptive_fwd_rev(
//     ab: *mut abpoa_t,
//     abpt: *mut abpoa_para_t,
//     query: *mut u8,
//     qlen: ::std::os::raw::c_int,
// ) -> AbPoaResT {
//     let mut abpoa_res: abpoa_res_t = std::mem::zeroed();
//     abpoa_align_sequence_to_graph(ab, abpt, query, qlen, &mut abpoa_res);
//     abpoa_res.into()
// }

pub fn abpoa_consensus_dna(
    abpt: *mut abpoa_para_t,
    seqs: &mut Vec<&mut [u8]>,
) -> Option<MsaResult> {
    let ab = Abpoa::new();
    abpoa_consensus_dna_core(ab.ptr(), abpt, seqs)
}

fn abpoa_consensus_dna_core(
    ab: *mut abpoa_t,
    abpt: *mut abpoa_para_t,
    seqs: &mut Vec<&mut [u8]>,
) -> Option<MsaResult> {
    if seqs.len() == 0 {
        return None;
    }
    let seq_ele2idx = &SEQ_NT4_TABLE;

    let mut seqs = seqs
        .iter_mut()
        .take(11)
        .map(|seq| {
            seq.iter()
                .map(|&base| seq_ele2idx[base as usize])
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut seq_lens = seqs
        .iter()
        .map(|seq| seq.len() as c_int)
        .collect::<Vec<_>>();
    let mut seqs = seqs
        .iter_mut()
        .map(|seq_encoded| seq_encoded.as_mut_ptr())
        .collect::<Vec<_>>();

    let n_seqs = seqs.len();

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

// pub fn abpoa_consensus_dna(param: &AbpoaParam, seqs: &mut Vec<&mut [u8]>) -> Option<MsaResult> {
//     if seqs.len() == 0 {
//         return None;
//     }
//     let seq_ele2idx = &SEQ_NT4_TABLE;

//     let idx2seq_ele = &IDX2NT;

//     let n_seqs = seqs.len();
//     let mut seq_lens = seqs
//         .iter()
//         .map(|seq| seq.len() as c_int)
//         .collect::<Vec<_>>();
//     let mut seqs = seqs
//         .iter_mut()
//         .map(|seq| {
//             seq.iter()
//                 .map(|&base| seq_ele2idx[base as usize])
//                 .collect::<Vec<_>>()
//         })
//         .collect::<Vec<_>>();
//     let mut seqs = seqs
//         .iter_mut()
//         .map(|seq_encoded| seq_encoded.as_mut_ptr())
//         .collect::<Vec<_>>();

//     unsafe {
//         let mut abpt = param.to_abpoa_para_t();

//         let ab = abpoa_init();
//         abpoa_msa(
//             ab,
//             &mut abpt,
//             n_seqs as c_int,
//             std::ptr::null_mut(),
//             seq_lens.as_mut_ptr(),
//             seqs.as_mut_ptr(),
//             std::ptr::null_mut(),
//             std::ptr::null_mut(),
//         );
//         let abc = &mut *(*ab).abc;
//         let result = abpoa_result_parser(abc, idx2seq_ele, &abpt);
//         abpoa_free(ab);
//         return result;
//     }
// }

#[cfg(test)]
mod test {
    use crate::{
        abpoa::abpoa_consensus_dna_core,
        abpoa_sys::{abpoa_free, abpoa_init},
        wrapper::AbpoaParam,
    };

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
            "AAGAAAACG",
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
    fn test_poa_consensus_b() {
        let mut align_param = AbpoaParam::new();
        align_param.modify_as_channel_consensus_param();
        align_param.post_set();
        let abpt = align_param.ptr();

        let mut seqs = get_seqs1();
        let mut seqs = seqs
            .iter_mut()
            .map(|v| unsafe { v.as_bytes_mut() })
            .collect::<Vec<_>>();
        let ab = unsafe { abpoa_init() };
        let res = abpoa_consensus_dna_core(ab, abpt, &mut seqs).unwrap();
        res.print_msa();
        println!("{:?}", res.cons_seq());

        let mut seqs = get_seqs2();
        let mut seqs = seqs
            .iter_mut()
            .map(|v| unsafe { v.as_bytes_mut() })
            .collect::<Vec<_>>();
        let res = abpoa_consensus_dna_core(ab, abpt, &mut seqs).unwrap();
        res.print_msa();
        println!("{:?}", res.cons_seq());

        let mut seqs = get_seqs3();
        let mut seqs = seqs
            .iter_mut()
            .map(|v| unsafe { v.as_bytes_mut() })
            .collect::<Vec<_>>();
        let res = abpoa_consensus_dna_core(ab, abpt, &mut seqs).unwrap();
        res.print_msa();
        println!("{:?}", res.cons_cov());
        unsafe { abpoa_free(ab) };
        println!("{:?}", res.cons_seq());
    }
}
