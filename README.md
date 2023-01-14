# kwantrace-r
KwanTrace re-implemented in Rust

## Package structure
This is a single package, which defines a library with several modules, and a binary for each scene to be rendered.

* Package
  * Library crate (in src/lib.rs)
    * vector module (in src/vector.rs)
      * test submodule 
    * matrix module (in src/matrix.rs)
  * Binary crate