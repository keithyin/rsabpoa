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

requirements:
    * requirements: 