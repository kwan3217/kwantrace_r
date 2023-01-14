# Phase 1

Since we are learning Rust, we will use the most basic data types to begin
with, lots of bare functions, and in general simple stuff. We will not be
implementing a full ray-tracer in phase 1. To begin with, we will however
implement the intersection between a ray and a sphere.

## `main.rs`
We aren't going to use any fancy modules, etc. We aren't going to define
our own linear algebra stuff. It looks like there are two competing linear
algebra packages which are suitable for us:

* `ndarray`, a library which attempts to fill the same space as Numpy in Python.
     This tends to deal with large multi-dimensional matrices. It doesn't have 
     any special handling for small cases like 3-or-4-element vectors, 3x3 or
     4x4 matrices, etc. It uses the same code for these as it would use for
     a 5-dimensional tensor with a million elements.
* `cgmath`, a library which is designed for *computer graphics* math. This
     is advertised as being optimized towards the use case of small matrices
     and vectors.
* The third option is to roll our own. I started to do that with `vector.rs`,
     `matrix.rs`, and `hmatrix.rs` but I don't know enough about new Rust
     concepts to feel like I could do a good job at designing the library
     using the native Rust concepts.

We will use `ndarray` first, because I found it first, and because it appears
to still be actively developed.

One notable thing about `ndarray` is that vectors are really 1-dimensional
arrays, neither row nor column vectors. It almost never really matters, and
when it does, it is almost always clear from context which is a row and which
is a column. This is handled in `ndarray` in the `dot()` function, which
calculates the dot product of two vectors, the matrix product of two compatible
matrices, the matrix product of a matrix and column vector, and the matrix product
of a row vector and matrix. In all cases, the vector is treated as having a
single index, never something like `[0,i]` or `[i,0]` where the vector is treated
as a matrix with a single row or column.

### `ray_sphere()`
The actual math is really straightforward. As with any ray-tracing math,
we are trying to find the ray parameter value where the ray hits any particular
primitive. The algebra is in the code comments for `ray_sphere()`. 

The Rust parts are where things get interesting. Just by looking through
the `ndarray` section of the Rust cookbook, it looks like we should use
ArrayView to pass the initial point and direction of the ray. We calculate
the intersection, and return an `Option<f64>` which is `None` if the ray
doesn't hit the sphere, and `Some(closest valid t)` if it does.

