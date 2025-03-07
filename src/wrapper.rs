use std::marker::PhantomData;

use crate::{abpoa::AlignMode, abpoa_sys::{
    abpoa_free, abpoa_free_para, abpoa_init, abpoa_init_para, abpoa_para_t, abpoa_post_set_para,
    abpoa_t,
}};

/// !Copy, !Clone, !Send, !Sync
/// let mut ap = AbpoaParam::new();
/// ap.modify_as_channel_consensus_param();
/// ap.post_set();
/// ap.ptr();
pub struct AbpoaParam {
    abpt: *mut abpoa_para_t,
    post_set_done: bool,
    _marker: PhantomData<*const ()>,
}

impl AbpoaParam {
    pub fn new() -> Self {
        Self {
            abpt: unsafe { abpoa_init_para() },
            post_set_done: false,
            _marker: PhantomData,
        }
    }

    pub fn post_set(&mut self) {
        assert!(!self.post_set_done);
        unsafe {
            abpoa_post_set_para(self.abpt);
        }
        self.post_set_done = true;
    }

    pub fn ptr(&self) -> *mut abpoa_para_t {
        assert!(self.post_set_done);
        self.abpt
    }

    pub fn modify_as_channel_consensus_param(&self) {
        assert!(!self.post_set_done);
        unsafe {
            let abpt = &mut (*self.abpt);
            abpt.set_out_msa(1);
            abpt.set_out_cons(1);
            // abpt.set_ret_cigar(1);
            abpt.set_amb_strand(1); // adaptive strand
            // abpt.align_mode = AlignMode::LOCAL as i32;
            abpt.match_ = 2;
            abpt.mismatch = 5;
            abpt.gap_open1 = 2;
            abpt.gap_open2 = 24;
            abpt.gap_ext1 = 1;
            abpt.gap_ext2 = 0;
            // abpt.set_disable_seeding(0);
            // abpt.k = 7;
            // abpt.w = 5;
        }
    }
}

impl Drop for AbpoaParam {
    fn drop(&mut self) {
        unsafe {
            abpoa_free_para(self.abpt);
        }
    }
}

pub struct Abpoa {
    ab: *mut abpoa_t,
    _marker: PhantomData<*const ()>,
}

impl Abpoa {
    pub fn new() -> Self {
        Self {
            ab: unsafe { abpoa_init() },
            _marker: PhantomData,
        }
    }

    pub fn ptr(&self) -> *mut abpoa_t {
        self.ab
    }
}

impl Drop for Abpoa {
    fn drop(&mut self) {
        unsafe {
            abpoa_free(self.ab);
        }
    }
}

#[cfg(test)]
mod test {
    use super::AbpoaParam;

    #[test]
    fn test_abpoa_param() {
        let abpoa = AbpoaParam::new();
    }
}
