use anyhow::Result;
use serde_json::{Value, json};

#[macro_use]
extern crate math_macros;
pub use astro_const::constants;
pub use constants::*;
use sci_file::OutputFile;
pub use simulation::{Simulation, System};

mod universe;
mod utils;

use universe::physics::force;
pub use universe::{ParticleType, Planet, Star, StarCsv, Universe};

#[derive(Debug)]
pub struct Spiroid {
    /// [`Universe`] contains parameters and workspace used by physics module.
    pub data: Universe,
    // Data output from the integrator is written into [`OutputFile`].
    pub output: OutputFile,
}

impl System for Spiroid {
    type Data = Universe;
    type Output = OutputFile;

    fn new(output: OutputFile, data: Universe) -> Self {
        Self { data, output }
    }

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
        self.data.update(time, y)?;
        // Compute the derivatives using the updated values.
        force(dy, &mut self.data)?;

        Ok(())
    }

    // Function that outputs each step of the solution.
    // TODO requires feature to be an interrupt to terminate the caller integration process.
    fn solout(
        &mut self,
        time: f64,
        y: &[f64],
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        let line = if self.data.orbiting_body.tides.kaula_enabled() {
            self.solout_kaula(time, y)
        } else {
            self.solout(time, y)
        };

        // Convert to JSON.
        self.output.write_json_line(&line)?;
        // TODO add option for binary output.
        Ok(())
    }
}

impl Spiroid {
    #[allow(clippy::unused_self)]
    // Called by the integrator to output the intermediate solutions.
    fn solout(&mut self, time: f64, y: &[f64]) -> Value {
        json!({
            "time": time / SECONDS_IN_YEAR,
            "radiative_zone_angular_momentum": y[0],
            "convective_zone_angular_momentum": y[1],
            "planet_semi_major_axis": y[2].powf(2. / 13.) / AU,
        })
    }

    #[allow(clippy::unused_self)]
    fn solout_kaula(&mut self, time: f64, y: &[f64]) -> Value {
        json!({
            "time": time / SECONDS_IN_YEAR,
            "radiative_zone_angular_momentum": y[0],
            "convective_zone_angular_momentum": y[1],
            "planet_semi_major_axis": y[2].powf(2. / 13.) / AU,
            "planet_spin": y[3],
            "planet_eccentricity": sqrt!(y[4]),
            "planet_inclination": y[5],
            "planet_longitude_ascending_node": y[6],
            "planet_pericentre_omega": y[7],
            "planet_spin_inclination": y[8],
        })
    }
}

#[cfg(test)]
mod tests;
