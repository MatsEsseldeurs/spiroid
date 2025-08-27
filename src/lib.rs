use anyhow::Result;
#[macro_use]
extern crate math_macros;
pub use astro_const::constants;
pub use constants::*;
pub use simulation::{Simulation, System};

mod universe;
mod utils;

use universe::physics::force;
pub use universe::{ParticleType, Planet, Star, StarCsv, Universe};

impl System<f64> for Universe {
    // This `derive` function is called by the integrator.
    // It should call the function that calculates the derivatives of the integration quantities.
    // i.e. fill `dy` with the derivatives of `y` with respect to x (`time`).
    fn derive(
        &mut self,
        time: f64,
        y: &[f64],
        dy: &mut [f64],
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        // Update the state of the universe based on the current integration values.
        Universe::update(self, time, y)?;
        // Compute the derivatives using the updated values.
        force(dy, self)?;

        Ok(())
    }

    // Update the state of the universe, prior to solution output.
    // Only called after successful integration step(s).
    fn update(
        &mut self,
        time: f64,
        y: &[f64],
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        // Update the state of the universe based on the current integration values.
        Universe::update(self, time, y)?;
        // Permanently clear destroyed particles.
        Universe::clear_destroyed_particles(self);

        Ok(())
    }
}

#[cfg(test)]
mod tests;
