rust binding for abpoa


abpoa v1.5.3

Basic Usage:
```rust
use rsabpoa::abpoa::{msa, AbpoaParam};


fn main() {
    let align_param = AbpoaParam::default();
    let seqs = vec!["AAC", "AC", "C"];
    let res = msa(&align_param, &seqs).unwrap();
    res.print_msa();
}
```


if you can't build rsabpoa, try install the following libs
```bash
apt-get update
apt-get install build-essential make libz-dev clang 
```


just a memo: I am using the following command to automatically generate the src/abpoa_sys.rs.
```
cargo install bindgen-cli
bindgen abPOA-v1.5.3/include/abpoa.h -o src/abpoa_sys.rs --allowlist-function 'abpoa.*'
```

