//! Black-box sparse polynomial interpolation over word-sized prime fields.
#![allow(
    clippy::cast_possible_truncation,
    clippy::missing_panics_doc,
    clippy::must_use_candidate,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::module_name_repetitions
)]

pub mod modp;
pub mod poly;
pub mod univariate;
pub mod vandermonde;
pub mod zippel;

pub use modp::{Primes, Rng};
pub use poly::{Exps, ModPoly};
pub use zippel::{interpolate, BlackBox};
