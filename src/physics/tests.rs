use super::*;
use crate::universe::effects::magnetism::{IsothermalWind, MagneticModel};
use crate::universe::effects::tides::TidalModel;
use crate::universe::effects::tides::kaula::tests::test_kaula;
use crate::universe::particles::planet::tests::{
    test_planet, test_planet_kaula, test_planet_magnetic,
};
use crate::universe::particles::star::tests::{
    TEST_DISK_LIFETIME, TEST_TIME, test_star, test_star_evolving,
};
use crate::universe::{Particle, ParticleType};
use pretty_assertions::assert_eq;

#[test]
fn _force_magnetic() {
    let star = test_star_evolving();
    let planet = test_planet_magnetic();

    let y = [
        star.radiative_zone_angular_momentum,
        star.convective_zone_angular_momentum,
        planet.semi_major_axis.powf(6.5),
    ];

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Wind(IsothermalWind::default()),
        },
        disk_lifetime: TEST_DISK_LIFETIME,
    };
    let mut result = y.to_vec();
    universe.update(TEST_TIME, &y).unwrap();
    let _ = force(TEST_TIME, &y, &mut result, &mut universe).unwrap();
    let expected = vec![
        -6.348994811695528e22,
        1.6351930535408648e22,
        -1.4634701453519956e43,
    ];
    assert_eq!(expected, result);
}

#[test]
fn _force_tides() {
    let star = test_star_evolving();
    let planet = test_planet();

    let y = [
        star.radiative_zone_angular_momentum,
        star.convective_zone_angular_momentum,
        planet.semi_major_axis.powf(6.5),
    ];
    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::ConstantTimeLag(1e-6),
            magnetism: MagneticModel::Disabled,
        },
        disk_lifetime: TEST_DISK_LIFETIME,
    };
    let mut result = y.to_vec();
    universe.update(TEST_TIME, &y).unwrap();
    let _ = force(TEST_TIME, &y, &mut result, &mut universe).unwrap();
    let expected = vec![
        -6.348994811695528e22,
        6.020027165936562e23,
        -1.9848639097150575e44,
    ];
    assert_eq!(expected, result);
}

#[test]
fn _force_magnetic_tides() {
    let star = test_star_evolving();
    let planet = test_planet_magnetic();

    let y = [
        star.radiative_zone_angular_momentum,
        star.convective_zone_angular_momentum,
        planet.semi_major_axis.powf(6.5),
    ];

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::ConstantTimeLag(1e-6),
            magnetism: MagneticModel::Wind(IsothermalWind::default()),
        },
        disk_lifetime: TEST_DISK_LIFETIME,
    };
    let mut result = y.to_vec();
    universe.update(TEST_TIME, &y).unwrap();
    let _ = force(TEST_TIME, &y, &mut result, &mut universe).unwrap();
    let expected = vec![
        -6.348994811695528e22,
        6.486208599049978e23,
        -2.131210924250257e44,
    ];
    assert_eq!(expected, result);
}

#[test]
fn _force_kaula() {
    let star = test_star();
    let planet = test_planet_kaula();
    let y = [
        0.0,
        0.0,
        planet.semi_major_axis.powf(6.5),
        8.062093352143078e-7,
        2.500000000179822e-5,
        0.34999207817863753,
        1.0465602799892118,
        -0.11536773671287792,
        0.31581363067032314,
    ];

    let mut universe = Universe {
        orbiting_body: Particle {
            kind: ParticleType::Planet(planet),
            tides: TidalModel::KaulaTides {
                kaula: test_kaula(),
            },
            magnetism: MagneticModel::Disabled,
        },
        central_body: Particle {
            kind: ParticleType::Star(star),
            tides: TidalModel::Disabled,
            magnetism: MagneticModel::Disabled,
        },
        disk_lifetime: TEST_DISK_LIFETIME,
    };
    let mut result = y.to_vec();
    universe.update(TEST_TIME, &y).unwrap();
    let _ = force(TEST_TIME, &y, &mut result, &mut universe).unwrap();
    let expected = vec![
        0.0,
        0.0,
        3.0436830856707734e49,
        -1.543298637727839e-9,
        5.129250585324416e-16,
        0.0007011714730129824,
        -0.0018601514573160808,
        5.876804234930107e-6,
        0.00035271712433657695,
    ];
    assert_eq!(expected, result);
}

#[test]
fn _force0() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let result = star_radiative_zone_angular_momentum_derivative(&star);
    let expected = -6.348994811695822e22;
    assert_eq!(expected, result);
}

#[test]
fn _force1() {
    let mut star = test_star();
    let planet = test_planet();
    star.refresh_tidal_frequency(&planet);
    let result =
        star_convective_zone_angular_momentum_derivative(TEST_TIME, &star, TEST_DISK_LIFETIME);
    let expected = -3.0079737689074846e22;
    assert_eq!(expected, result);
}

#[test]
fn _force2_magnetic_tides() {
    let mut star = test_star();
    let planet = test_planet_magnetic();
    star.refresh_tidal_frequency(&planet);
    let tides = TidalModel::ConstantTimeLag(1e-6);
    let mut magnetism = MagneticModel::Wind(IsothermalWind::default());
    let tidal_torque = tides.tidal_torque(&star, &planet);
    let magnetic_torque = magnetism.magnetic_torque(&planet, &star);
    star.update_torques(tidal_torque, magnetic_torque);

    let result = planet_semi_major_axis_13_div_2_derivative(&planet, &star);
    let expected = -2.1307551258578705e44;
    assert_eq!(expected, result);
}

#[test]
fn _force2_kaula() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);

    let tides = TidalModel::ConstantTimeLag(1e-6);
    let mut magnetism = MagneticModel::Wind(IsothermalWind::default());
    let tidal_torque = tides.tidal_torque(&star, &planet);
    let magnetic_torque = magnetism.magnetic_torque(&planet, &star);
    star.update_torques(tidal_torque, magnetic_torque);

    let result = planet_semi_major_axis_13_div_2_derivative(&planet, &star);
    let expected = -1.984887979983568e44;
    assert_eq!(expected, result);
}

#[test]
fn _force4() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_spin_derivative(&planet, &star, &kaula);
    let expected = 1.2501842317327892e-6;
    assert_eq!(expected, result);
}

#[test]
fn _force5() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_eccentricity_derivative(&planet, &star, &kaula);
    let expected = -2.371917363949444e-13;
    assert_eq!(expected, result);
}

#[test]
fn _force6() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_inclination_derivative(&planet, &star, &kaula);
    let expected = -0.26702870058883815;
    assert_eq!(expected, result);
}

#[test]
fn _force7() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_longitude_ascending_node_derivative(&planet, &star, &kaula);
    let expected = 0.6637879161156786;
    assert_eq!(expected, result);
}

#[test]
fn _force8() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_argument_pericentre_derivative(&planet, &star, &kaula);
    let expected = -42049166.159453586;
    assert_eq!(expected, result);
}

#[test]
fn _force9() {
    let mut star = test_star();
    let planet = test_planet_kaula();
    star.refresh_tidal_frequency(&planet);
    let kaula = test_kaula();
    let result = planet_spin_axis_inclination_derivative(&planet, &star, &kaula);
    let expected = -0.1641163464603492;
    assert_eq!(expected, result);
}
