import pandas as pd
import numpy as np

# This script processes stellar evolution data from MESA output CSV files.
# It reads the data, computes numerical derivatives for physical quantities,
# and filters out rows where changes are insignificant.
# The goal is to reduce the dataset size by keeping only rows with meaningful evolution,
# reducing computational cost for models that use this data.
# The filtered data is then saved back to the original CSV file location.

mass = 1.0  # in Msun
file_path = f"../examples/data/star/evolution/mesa_{10*mass:02.0f}.csv"

df = pd.read_csv(file_path)
periods = len(df) // 100  # Adjust the period based on the length of the DataFrame

# Compute numerical derivatives for all columns except the first, with respect to the index
derivatives = df.iloc[:, 1:].diff(periods=periods)
derivatives.columns = [f"d_{col}/d_index" for col in df.columns[1:]]

# Exclude convective_turnover_time from filtering criteria
columns = [col for col in derivatives.columns if col[2:-8] != "convective_turnover_time"]

# Calculate the maximum relative derivative for each row
maxderiv = np.zeros(len(derivatives))
for i in range(len(derivatives)):
    maxderiv[i] = max([abs(derivatives[col][i]/df[col[2:-8]][i]) for col in columns if df[col[2:-8]][i] != 0])

keep_indices = maxderiv > 0.001
keep_indices[:periods] = True  # Always keep the first 'periods' rows
for i in range(periods, len(derivatives)):
    if np.sum(keep_indices[i-periods:i]) == 0: # Ensure at least one row is kept in each segment
        keep_indices[i] = True

# Filter the DataFrame and recompute derivatives for the filtered data
filtered_df = df.iloc[keep_indices].reset_index(drop=True)
filtered_derivatives = filtered_df.iloc[:, 1:].diff(periods=periods)
filtered_derivatives.columns = [f"d_{col}/d_index" for col in filtered_df.columns[1:]]
filtered_maxderiv = np.zeros(len(filtered_derivatives))
for i in range(len(filtered_derivatives)):
    filtered_maxderiv[i] = max([abs(filtered_derivatives[col][i]/filtered_df[col[2:-8]][i]) for col in columns if filtered_df[col[2:-8]][i] != 0])

# Save the filtered DataFrame back to the original CSV file location
output_csv = f"../examples/data/star/evolution/mesa_{10*mass:02.0f}.csv"
filtered_df.to_csv(output_csv, index=False)
