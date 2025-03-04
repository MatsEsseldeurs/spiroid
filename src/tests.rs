use super::*;
use simulation::{Integrator, System};
use sci_file::{
    OutputFile, deserialize_csv_column_vectors_from_path, deserialize_csv_rows_from_path,
    deserialize_json_from_path,
};
use simulation::InputConfig;
use std::path::PathBuf;

use pretty_assertions::assert_eq;

struct Test {
    pub data: Universe,
}

// Mock System implementation without writing any output.
impl System for Test {
    type Data = Universe;
    type Output = OutputFile;

    fn new(_: OutputFile, data: Universe) -> Self {
        Self { data }
    }

    fn derive(
        &mut self,
        time: f64,
        y: &[f64],
        dy: &mut [f64],
    ) -> Result<(), Box<(dyn std::error::Error + Send + Sync + 'static)>> {
        // Update the state of the universe based on the current integration values.
        self.data.update(time, y)?;
        // Compute the derivatives using the updated values.
        force(time, y, dy, &mut self.data)?;
        Ok(())
    }

    fn solout(
        &mut self,
        _time: f64,
        _y: &[f64],
    ) -> Result<(), Box<(dyn std::error::Error + Send + Sync + 'static)>> {
        Ok(())
    }
}

fn test_simulation(config: PathBuf) -> Vec<f64> {
    // Parse the config file.
    let mut config: InputConfig<Universe> = deserialize_json_from_path(&config).unwrap();
    config.initial_time *= SECONDS_IN_YEAR;
    config.final_time *= SECONDS_IN_YEAR;

    // Load stellar evolution data from file.
    if let ParticleType::Star(star) = &mut config.universe.central_body.kind {
        // Load stellar evolution data from file if stellar evolution is enabled.
        if let Some(star_file) = star.evolution_file() {
            let mut stellar_data = deserialize_csv_rows_from_path::<StarCsv>(star_file).unwrap();
            // Configure the stellar evolution interpolator.
            let (star_ages, star_values) = StarCsv::initialise(&mut stellar_data);
            star.initialise_evolution(&star_ages, &star_values);
        };
    }

    // Load love number data from file(s) if kaula tides are enabled.
    if let Some(kaula) = config.universe.orbiting_body.tides.kaula_get_mut() {
        if let Some(solid_file) = kaula.solid_file() {
            let love_solid = deserialize_csv_column_vectors_from_path::<f64>(solid_file).unwrap();
            kaula.initialise_love_number_solid(&love_solid);
        }
        if let Some(ocean_file) = kaula.ocean_file() {
            let love_ocean = deserialize_csv_column_vectors_from_path::<f64>(ocean_file).unwrap();
            kaula.initialise_love_number_ocean(&love_ocean);
        }
    }

    // Initialise the universe (star, planet, etc).
    config.universe.initialise(config.initial_time).unwrap();

    // Initial values for the integrator.
    let y = config.universe.integration_quantities();
    let mut system = Test {
        data: config.universe,
    };
    config
        .integrator
        .initialise(config.initial_time, config.final_time, &y)
        .unwrap();

    // Run the full integration.
    let _ = config.integrator.integrate(&mut system).unwrap();

    // Collect the final y values.
    config.integrator.y_final()
}

#[test]
fn example_no_effects() {
    let result = test_simulation("examples/no_effects.conf".into());
    let expected = vec![
        4.787118213002163e40,
        1.3794130696558446e40,
        2.8113413766640534e61,
    ];
    assert_eq!(expected, result);
}

#[test]
fn example_tides() {
    let result = test_simulation("examples/tides.conf".into());
    let expected = vec![
        4.787143163297013e40,
        1.3794202173212958e40,
        2.9069618463376594e60,
    ];
    assert_eq!(expected, result);
}

#[test]
fn example_magnetic() {
    let result = test_simulation("examples/magnetic.conf".into());
    let expected = vec![
        5.394133407982592e40,
        1.5559592497137857e40,
        7.565340545548777e59,
    ];
    assert_eq!(expected, result);
}

#[test]
fn example_magnetic_tides() {
    let result = test_simulation("examples/magnetic_tides.conf".into());
    let expected = vec![
        4.787136126036125e40,
        1.3794182013101542e40,
        7.069089651047147e60,
    ];
    assert_eq!(expected, result);
}

#[test]
fn example_kaula_solid() {
    let result = test_simulation("examples/kaula_solid.conf".into());
    let expected = vec![
        5.194432297602171e-5,
        5.194432297602171e-5,
        5.281115912876036e71,
        8.062144319508726e-7,
        2.5000000001796287e-5,
        0.3499902860956476,
        1.0464928313280928,
        -0.11519339765695413,
        0.3158644836762971,
    ];
    assert_eq!(expected, result);
}

#[test]
fn example_all_effects() {
    let result = test_simulation("examples/all_effects.conf".into());
    let expected = vec![
        1.855676942283104e42,
        4.990089513734989e41,
        5.281115912659979e71,
        8.07398442384458e-7,
        2.5000000001118837e-5,
        0.3495999235632858,
        1.0302720218046357,
        -0.07395182155411484,
        0.32784380302637944,
    ];
    assert_eq!(expected, result);
}
