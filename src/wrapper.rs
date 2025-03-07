use std::marker::PhantomData;

use crate::abpoa_sys::{abpoa_free_para, abpoa_init_para, abpoa_para_t, abpoa_post_set_para};

/// !Copy, !Clone, !Send, !Sync
///
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
    
}

impl Drop for AbpoaParam {
    fn drop(&mut self) {
        unsafe {
            abpoa_free_para(self.abpt);
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
