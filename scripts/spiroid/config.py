"""
Generates dictionaries matching the structure of input configs into spiroid.

"""

import itertools


def make_config(planet, star, disk_lifetime, integrator, start_time, final_time):
    """Generates the simulation template required to launch the simulation in Rust."""
    return {
        "resume": False,
        "initial_time": start_time,
        "final_time": final_time,
        "integrator": integrator,
        "universe": {
            "disk_lifetime": disk_lifetime,
            "central_body": star,
            "orbiting_body": planet,
        },
    }


def make_integrator(effects):
    """Generates appropriate integrator template."""
    # Kaula tides on the planet modifies the integrator defaults to avoid timestep issue.
    if effects["PLANET_TIDES_ENABLED"]:
        integrator = kaula_odex()
    else:
        integrator = base_odex()

    return integrator


def base_odex():
    """Generates basic Odex integrator template."""
    return {
        "Odex": {
            "step_size_reduction_factor": 0.66666666666666666666,
            "step_size_selection_b": 2.0,
            "max_step_size": 15778800000000.0,
            "max_integration_steps": 100000000,
        }
    }


def kaula_odex():
    """Generates Odex integrator template with values tuned for kaula simulations."""
    integrator = base_odex()
    integrator["Odex"]["max_step_size"] = 31557600000.0
    integrator["Odex"]["max_integration_steps"] = 500000000
    integrator["Odex"]["step_control_safety_a"] = 0.05
    integrator["Odex"]["step_control_safety_b"] = 0.2

    return integrator


def make_planets(planet_base, effects):
    """Generate all combinations of planets based on specified values of `planet_base` dictionary."""
    planets = []
    combis = [x for x in itertools.product(*planet_base.values())]

    for planet_vals in combis:
        (mass, radius, semi_major_axis, magnetic_field, is_destroyed) = planet_vals[:5]
        body = {}
        planet = {
            "mass": mass,
            "radius": radius,
            "semi_major_axis": semi_major_axis,
            "magnetic_field": magnetic_field,
            "is_destroyed": is_destroyed,
        }

        if effects["PLANET_TIDES_ENABLED"]:
            (
                inclination,
                eccentricity,
                spin,
                longitude_ascending_node,
                pericentre_omega,
                spin_inclination,
                radius_of_gyration,
                (particle_type, kaula_solid_file),
            ) = planet_vals[5:]
            planet.update(
                {
                    "inclination": inclination,
                    "eccentricity": eccentricity,
                    "spin": spin,
                    "longitude_ascending_node": longitude_ascending_node,
                    "pericentre_omega": pericentre_omega,
                    "spin_inclination": spin_inclination,
                    "radius_of_gyration_2": radius_of_gyration,
                }
            )

            body["tides"] = {
                "KaulaTides": {
                    "kaula": {
                        "particle_type": {
                            particle_type: {"solid_file": kaula_solid_file}
                        },
                    }
                }
            }

        body["kind"] = {"Planet": planet}
        planets.append(body)
    return planets


def make_stars(star_base, effects):
    """Generate all combinations of stars based on specified values of `star_base` dictionary."""
    stars = []

    combis = [x for x in itertools.product(*star_base.values())]
    for star_vals in combis:
        (
            mass,
            spin,
            core_envelope_coupling_constant,
            footpoint_conductance,
            star_file_path,
            sigma_bar,
        ) = star_vals[:6]

        body = {}
        star = {
            "mass": mass,
            "spin": spin,
            "core_envelope_coupling_constant": core_envelope_coupling_constant,
            "footpoint_conductance": footpoint_conductance,
            "evolution": "Disabled",
        }

        if effects["STAR_TIDES_ENABLED"]:
            body["tides"] = {"ConstantTimeLag": sigma_bar}

        if effects["MAGNETIC_EFFECT_ENABLED"]:
            body["magnetism"] = {"Wind": {}}

        if effects["STAR_EVOLUTION_ENABLED"]:
            star["evolution"] = {"Interpolated": {"star_file_path": star_file_path}}
        else:
            star["radiative_moment_of_inertia"] = star_vals[6]
            star["convective_moment_of_inertia"] = star_vals[7]
        body["kind"] = {"Star": star}
        stars.append(body)

    return stars


def generate_all_configs(
    start_time, final_time, disk_lifetime, planet_base, star_base, effects
):
    """Generates a simulation configuration file for each combination of planets and stars."""
    integrator = make_integrator(effects)
    planets = make_planets(planet_base, effects)
    stars = make_stars(star_base, effects)

    # Generate a simulation input config for all combinations
    # of the star and planet values.
    return (
        make_config(planet, star, disk_lifetime, integrator, start_time, final_time)
        for (planet, star) in itertools.product(planets, stars)
    )


def generate_all_effect_combinations(input_dict):
    """Generate all possible combinations of enabled effects for the simulations."""
    import itertools

    tags = {
        "MAGNETIC_EFFECT_ENABLED": "magnetism",
        "STAR_EVOLUTION_ENABLED": "star_evolution",
        "STAR_TIDES_ENABLED": "star_ctl_tides",
        "PLANET_TIDES_ENABLED": "planet_kaula_tides",
    }
    # Get keys and values from the input dictionary
    keys = input_dict.keys()
    values = input_dict.values()
    # Generate all combinations of enabled effects.
    combinations = [x for x in itertools.product(*values)]
    # Create a list of tuples (dictionary, label) where label is the enabled effects.
    result = []
    for combo in combinations:
        # Create the dictionary from the combination
        combo_dict = dict(zip(keys, combo))
        # Concatenate labels from keys with True values (enabled effects) for the filename.
        true_keys = "-".join(tags[key] for key, value in combo_dict.items() if value)
        if true_keys == "":
            true_keys = "no_effects"
        # Append the tuple (dictionary, concatenated string) to the result
        result.append((combo_dict, true_keys))

    return result
