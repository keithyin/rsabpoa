use std::{marker::PhantomData, ops::{Deref, DerefMut}, sync::Arc};

use crate::abpoa_sys::{abpoa_free_para, abpoa_init_para, abpoa_para_t};



/// !Copy, !Clone, !Send, !Sync
/// 
pub struct AbpoaParam {
    abpt: *mut abpoa_para_t,
    _marker: PhantomData<*const ()>,
}

impl AbpoaParam {
    pub fn new() -> Self {
        Self {
            abpt: unsafe { abpoa_init_para() },
            _marker: PhantomData,
        }
    }

    pub fn ptr(&self) -> *mut abpoa_para_t {
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