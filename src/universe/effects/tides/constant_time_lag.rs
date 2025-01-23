



// This is a re-write of Eq. 3 and 19 from Benbakoura et al. 2019
// without the factors that are in the function semi_major_axis_13_div_2_derivative in physics.rs
// The a^-6 is here to compensate the a^6 in physics.rs
fn tidal_torque_ctl(star: &Star, tidal_quality: f64, planet: &Planet) -> f64 {
    let depth = 1E-08; // Smoothing parameter when tidal frequency is 0
    -(9. / 4.)
        * planet.mass.powi(2)
        * GRAVITATIONAL
        * planet.semi_major_axis.powi(-6)
        * tanh!(star.tidal_frequency / depth)
        * star.radius.powi(5)
        / tidal_quality
}
