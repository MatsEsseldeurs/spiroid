use super::*;
use crate::universe::effects::tides::kaula::love_number::tests::test_love_number;
use crate::universe::effects::tides::kaula::polynomials::tests::test_polynomials;
use crate::universe::particles::planet::tests::test_planet_kaula;

use pretty_assertions::assert_eq;

#[cfg(test)]
pub fn test_kaula() -> Kaula {
    let mut kaula = Kaula {
        particle_type: ParticleComposition::Solid {
            solid_file: "dummy".into(),
        },
        polynomials: test_polynomials(),
        love_number: test_love_number(),
        sum_over_m_imaginary_mfactor: 0.,
        sum_over_m_imaginary_pfactor: 0.,
        sum_over_m_real_2pq_2mp_dt: 0.,
    };
    kaula.sum_over_m_real_2pq_2mp_dt = kaula.sum_over_m_real(
        &kaula.polynomials.eccentricity_2pq_squared,
        &kaula.polynomials.inclination_2mp_squared_derivative,
    );
    kaula.sum_over_m_imaginary_mfactor = kaula.sum_over_m_imaginary_mfactor(
        &kaula.polynomials.eccentricity_2pq_squared,
        &kaula.polynomials.inclination_2mp_squared,
    );
    kaula.sum_over_m_imaginary_pfactor = kaula.sum_over_m_imaginary_pfactor(
        &kaula.polynomials.eccentricity_2pq_squared,
        &kaula.polynomials.inclination_2mp_squared,
    );
    kaula
}

#[test]
fn _summation_of_longitudinal_modes_semi_major_axis() {
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_semi_major_axis();
    let expected = 0.5630675119283106;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_spin() {
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_spin();
    let expected = 0.5641312760456983;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_eccentricity() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_eccentricity(planet.semi_minor_axis_ratio);
    let expected = 7.73439335780939e-5;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_inclination() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_inclination(&planet);
    let expected = -7423.755336647935;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_longitude_ascending_node() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_longitude_ascending_node(&planet);
    let expected = 5.076723278969285e-33;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_spin_axis_inclination() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_spin_axis_inclination(&planet);
    let expected = -0.05924449315462055;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_pericentre_eccentricity() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_pericentre_eccentricity(&planet);
    let expected = -3.2159666539736226e-25;
    assert_eq!(expected, result);
}

#[test]
fn _summation_of_longitudinal_modes_pericentre_inclination() {
    let planet = test_planet_kaula();
    let kaula = test_kaula();
    let result = kaula.summation_of_longitudinal_modes_pericentre_inclination(&planet);
    let expected = -5.258897087031897e-34;
    assert_eq!(expected, result);
}
