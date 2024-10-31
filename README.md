rust binding for abpoa


abpoa v1.5.3

bindgen:
```
cargo install bindgen-cli
bindgen abPOA-v1.5.3/include/abpoa.h -o src/abpoa_sys.rs --allowlist-function 'abpoa.*'
```

basic environ
```
apt-get update
apt-get install build-essential make libz-dev clang 
```



```
use raabpoa::abpoa::{msa, AbpoaParam};


fn main() {
    let align_param = AbpoaParam::default();
    let seqs = vec!["AAC", "AC", "C"];
    let res = msa(&align_param, &seqs).unwrap();
    res.print_msa();
}
```