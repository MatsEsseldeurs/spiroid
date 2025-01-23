// Precomputed factorials as f64 for 0..=4
const FACTORIALS: &[f64] = &[1.0, 1.0, 2.0, 6.0, 24.0];

pub(crate) fn factorial(x: usize) -> f64 {
    FACTORIALS[x]
}

// x == y -> 1.0
// x != y -> 0.0
pub(crate) fn kronecker_delta(x: usize, y: usize) -> f64 {
    f64!(x == y)
}

// Maps a set of 3D coordinates into a 1D index.
pub(crate) fn map_3d_to_1d(x: usize, x_max: usize, y: usize, y_max: usize, z: usize) -> usize {
    (z * x_max * y_max) + (y * x_max) + x
}

#[cfg(test)]
mod tests;
