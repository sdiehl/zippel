//! Black-box sparse polynomial and rational function interpolation over word-sized prime fields.
#![allow(
    clippy::cast_possible_truncation,
    clippy::missing_panics_doc,
    clippy::must_use_candidate,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::module_name_repetitions
)]

pub mod benor;
pub mod factors;
pub mod geometric;
pub mod poly;
pub mod rational;
pub mod zippel;

pub use poly::{dense, Dense, Exps, ModPoly};
pub use rational::{reconstruct, RatFunc};
pub use zippel::interpolate;
