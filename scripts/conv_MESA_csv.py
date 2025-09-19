import numpy as np
import mesa_reader as mr
import pandas as pd
from tqdm import tqdm
import os

from STEDA.MESA import loadMesa
from STEDA.constants import Rsun, Msun, Lsun

mass = 1.0  # in Msun
mesa_dir = "path_to_your_mesa_directory"  # Change this to your MESA directory
mesa_file = os.path.join(mesa_dir, "LOGS")

log = mr.MesaLogDir(mesa_file, memoize_profiles=False)

properties = {"age":                              np.zeros(len(log.profile_numbers)),
            "radius":                           np.zeros(len(log.profile_numbers)),
            "mass":                             np.zeros(len(log.profile_numbers)),
            "convective_radius":                np.zeros(len(log.profile_numbers)),
            "radiative_mass":                   np.zeros(len(log.profile_numbers)),
            "radiative_moment_of_inertia":      np.zeros(len(log.profile_numbers)),
            "convective_moment_of_inertia":     np.zeros(len(log.profile_numbers)),
            "luminosity":                       np.zeros(len(log.profile_numbers)),
            "convective_turnover_time":         np.zeros(len(log.profile_numbers)),
            "mass_loss_rate":                   np.zeros(len(log.profile_numbers)),
            }

for profile in tqdm(log.profile_numbers):
    Ms, Rs, Ls, age, fase, r, rho, drho, rho_mean, mass, g, dg, dN2, N, dN, Kt, i_int_enve, i_int_core, lc, vc = loadMesa(log, profile)
    
    h = log.history_data
    Mdot = -h.data_at_model_number("star_mdot",log.model_with_profile_number(profile))

    # Calculate the convective turnover time throughout the convective envelope
    if i_int_core < len(r) - 1 and i_int_core != 0:

        # # Calculate the convective turnover time at the interface between the radiative core and convective envelope
        # option 1
        # filter = [i >= i_int_core and i < i_int_core+30 for i in range(len(r))]
        # filter[-1] = False
        # tc_out = abs(np.percentile(lc[filter]/vc[filter],50))
        # option 2
        # filter = [i >= i_int_core and N[i] == 0. for i in range(len(r))]
        # filter[-1] = False
        # tc_out = abs(np.percentile(lc[filter]/vc[filter],50))

        # Calculate the convective turnover time at the center of the convective envelope
        i_cent = np.where(r>=0.5*(Rs+r[i_int_core]))[0][0]
        while (vc[i_cent] == 0.0 or N[i_cent] > 0.0) and i_cent > i_int_core: i_cent -= 1 # go down until we find a convective point
        tc_out = lc[i_cent]/vc[i_cent]
    else:
        tc_out = 1e99
    
    # adjusted_convective_mass = mass[i_int_core]/Ms
    # t_c_base_out = 10**(8.79 - 2. * abs(np.log10(adjusted_convective_mass))**(0.349) - 0.0194 * abs(np.log10(adjusted_convective_mass))**2 - 1.62 * min(np.log10(adjusted_convective_mass) + 8.55, 0.))

    properties["age"                              ][profile-1] = age # yr
    properties["radius"                           ][profile-1] = Rs / Rsun
    properties["mass"                             ][profile-1] = Ms / Msun
    properties["convective_radius"                ][profile-1] = r[i_int_core] / Rsun
    properties["radiative_mass"                   ][profile-1] = (mass[i_int_core]) / Msun # MESA radiative mass
    integrand = rho * r**4
    properties["radiative_moment_of_inertia"      ][profile-1] = max(8 * np.pi / 3 * np.trapz(integrand[:i_int_core], r[:i_int_core]), 1e44 * 1e7) / (Ms * Rs**2)
    properties["convective_moment_of_inertia"     ][profile-1] = max(8 * np.pi / 3 * np.trapz(integrand[i_int_core:], r[i_int_core:]), 1e44 * 1e7) / (Ms * Rs**2)
    properties["luminosity"                       ][profile-1] = Ls / Lsun
    properties["convective_turnover_time"         ][profile-1] = tc_out
    properties["mass_loss_rate"                   ][profile-1] = Mdot

output_csv = f"../examples/data/star/evolution/mesa_{10*mass:02.0f}.csv"
df = pd.DataFrame(properties)
df.to_csv(output_csv, index=False)
