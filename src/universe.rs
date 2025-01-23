pub(crate) mod effects;
pub(crate) mod particles;

use crate::SECONDS_IN_YEAR;
pub use effects::Kaula;
pub use particles::{Particle, ParticleType, Planet, Star, StarCsv};

use anyhow::Result;
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize, Debug)]
pub struct Universe {
    pub disk_lifetime: f64,
    pub orbiting_body: Particle,
    pub central_body: Particle,
}

impl Universe {
    #[allow(clippy::missing_errors_doc)]
    // Apply the unit conversions to initial input values.
    pub fn initialise(&mut self, time: f64) -> Result<()> {
        self.disk_lifetime *= SECONDS_IN_YEAR;

        if let ParticleType::Star(star) = &mut self.central_body.kind {
            star.initialise(time)?;
        };

        if let ParticleType::Planet(planet) = &mut self.orbiting_body.kind {
            planet.initialise();
        };

        Ok(())
    }

    // Creates a vector of initial quantities to be integrated, depending on the simulation configuration.
    pub fn integration_quantities(&self) -> Vec<f64> {
        let mut vec = vec![];
        if let ParticleType::Star(star) = &self.central_body.kind {
            vec.append(&mut vec![
                star.initial_spin * star.radiative_moment_of_inertia,
                star.initial_spin * star.convective_moment_of_inertia,
            ]);
        };

        if let ParticleType::Planet(planet) = &self.orbiting_body.kind {
            vec.append(&mut vec![
                // The actual integrated quantity is sma^6.5.
                // See comment for fn planet_semi_major_axis_13_div_2_derivative
                planet.semi_major_axis.powf(6.5),
            ]);

            if self.orbiting_body.tides.kaula_enabled() {
                vec.append(&mut vec![
                    planet.spin,
                    // The actual integrated quantity is eccentricity^2 to avoid singularities when eccentricity goes to 0.
                    planet.eccentricity.powi(2),
                    planet.inclination,
                    planet.longitude_ascending_node,
                    planet.pericentre_omega,
                    planet.spin_inclination,
                ]);
            }
        };

        vec
    }

    // Update the planet and star values from the integrator prior to the derivation step.
    pub(crate) fn update(&mut self, time: f64, y: &[f64]) -> Result<()> {
        // ***WARNING!***
        // Stateful function.planet.longitude_ascending_node
        // The order of these calculations is important.
        // Lower order calculations depend on previous values.
        // ***WARNING!***
        if let ParticleType::Star(star) = &mut self.central_body.kind {
            // Update radiative zone (y[0]) and convective zone (y[1]) angular momentum
            // and recompute independent values.
            star.refresh(time, y[0], y[1])?;

            if let ParticleType::Planet(planet) = &mut self.orbiting_body.kind {
                // Nothing to compute if the planet is already destroyed.
                if planet.is_destroyed {
                    return Ok(());
                }
                // Invert the exponent of sma^6.5 to normalise the semi major axis.
                // Recompute planet values, including those depending on star.
                planet.refresh(y[2].powf(2. / 13.), star);

                // The planet may have been destroyed in the current iteration.
                if planet.is_destroyed {
                    star.update_torques(0., 0.);
                    return Ok(());
                }

                // Compute the enabled effects (magnetism, stellar tides, planet tides)
                if time < self.disk_lifetime {
                    star.update_torques(0., 0.);
                } else {
                    // Recompute star values that depend on planet (tidal and magnetic torque).
                    star.refresh_tidal_frequency(planet);
                    let tidal_torque = self.central_body.tides.tidal_torque(star, planet);
                    let magnetic_torque = self.central_body.magnetism.magnetic_torque(planet, star);
                    star.update_torques(tidal_torque, magnetic_torque);
                }

                if self.orbiting_body.tides.kaula_enabled() {
                    //(spin, eccentricity, inclination, longitude_ascending_node, pericentre_omega, spin_inclination)
                    // Invert the exponent of e^2 to normalise the eccentricity.
                    planet.refresh_orbital_elements(y[3], sqrt!(y[4]), y[5], y[6], y[7], y[8]);
                    // Recompute the kaula tidal effects.
                    self.orbiting_body.tides.refresh_kaula(time, star, planet)?;
                }
            };
        };

        Ok(())
    }
}
