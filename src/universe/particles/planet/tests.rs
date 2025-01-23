use super::*;
use crate::universe::particles::star::tests::test_star;
use pretty_assertions::assert_eq;

// Test setup functions
#[cfg(test)]
fn test_planet_base() -> Planet {
    let mut planet = Planet::new();

    planet.mass = 1.8979999999999997e26;
    planet.radius = 32550000.0;
    planet.semi_major_axis = 2.0583049171152424e60_f64.powf(2. / 13.);

    planet
}

#[cfg(test)]
pub fn test_planet_magnetic() -> Planet {
    let mut planet = test_planet_base();
    planet.magnetic_field = 0.001;
    planet.magnetic_pressure = magnetic_pressure(planet.magnetic_field);
    let star = test_star();
    planet.refresh(planet.semi_major_axis, &star);

    planet
}

#[cfg(test)]
pub fn test_planet() -> Planet {
    let mut planet = test_planet_base();
    let star = test_star();
    planet.refresh(planet.semi_major_axis, &star);

    planet
}

#[cfg(test)]
pub fn test_planet_kaula() -> Planet {
    let mut planet = test_planet_base();
    planet.spin = 8.0e-7;
    planet.eccentricity = 0.005;
    planet.inclination = 0.35;
    planet.longitude_ascending_node = 1.0;
    planet.pericentre_omega = 0.0;
    planet.spin_inclination = 0.34;
    planet.radius_of_gyration_2 = 0.33;
    planet.moment_of_inertia = 5.9e37;

    let star = test_star();
    planet.refresh(planet.semi_major_axis, &star);

    planet
}

// Tests below
#[test]
fn _density_ratio() {
    let expected = 0.5587272987074543;
    let star = test_star();
    let planet = test_planet_magnetic();
    let result = planet.density_ratio(star.mass, star.radius);
    assert_eq!(expected, result);
}

#[test]
fn _mean_motion() {
    let expected = 0.0001243224991067286;
    let star = test_star();
    let planet = test_planet_magnetic();
    let result = planet.mean_motion(star.mass);
    assert_eq!(expected, result);
}

#[test]
fn _inside_alfven_radius() {
    let alfven_radius = 15078958988.616783;
    let expected = true;
    let planet = test_planet_magnetic();
    let result = planet.inside_alfven_radius(alfven_radius);
    assert_eq!(expected, result);
}

#[test]
fn _roche_limit() {
    let expected = 1600111307.3665166;
    let star = test_star();
    let planet = test_planet_magnetic();
    let result = planet.roche_limit(star.radius);
    assert_eq!(expected, result);
}

#[test]
fn _orbit_lower_limit() {
    let expected = 1600111307.3665166;
    let star = test_star();
    let planet = test_planet_magnetic();
    let result = planet.orbit_lower_limit(star.radius);
    assert_eq!(expected, result);
}

#[test]
fn _destroy_planet() {
    let expected = false;
    let planet = test_planet_magnetic();
    let result = planet.crossed_orbital_lower_limit();
    assert_eq!(expected, result);
}

#[test]
fn _magnetic_pressure() {
    let magnetic_field = 0.0;
    let expected = 0.0;
    let result = magnetic_pressure(magnetic_field);
    assert_eq!(expected, result);
}
