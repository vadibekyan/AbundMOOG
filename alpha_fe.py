import pandas as pd
import numpy as np

# Load the RDB (tab-separated) file
file_path = "alpha_fe.rdb"  # update this path if needed
df = pd.read_csv(file_path, sep="\t")

# Define the elements and star columns
alpha_elements = ["MgI", "SiI", "TiI"]
stars = {
    "epsInd": ("epsInd_abund", "epsInd_abund_err"),
    "HD131977": ("HD131977_abund", "HD131977_abund_err"),
    "HD191408": ("HD191408_abund", "HD191408_abund_err")
}

# Filter to keep only the necessary elements
df_filtered = df[df['element'].isin(alpha_elements + ["FeI"])]
df_filtered.set_index("element", inplace=True)

# Weighted average function
def weighted_avg_and_uncertainty(values, errors):
    weights = 1 / np.square(errors)
    avg = np.sum(values * weights) / np.sum(weights)
    unc = np.sqrt(1 / np.sum(weights))
    return avg, unc

# Compute alpha - FeI for each star
results = {}
for star, (ab_col, err_col) in stars.items():
    alpha_vals = df_filtered.loc[alpha_elements, ab_col].values
    alpha_errs = df_filtered.loc[alpha_elements, err_col].values
    fei_val = df_filtered.loc["FeI", ab_col]
    fei_err = df_filtered.loc["FeI", err_col]

    alpha_avg, alpha_unc = weighted_avg_and_uncertainty(alpha_vals, alpha_errs)
    alpha_fei = alpha_avg - fei_val
    alpha_fei_unc = np.sqrt(alpha_unc**2 + fei_err**2)

    results[star] = {
        "alpha-FeI": alpha_fei,
        "uncertainty": alpha_fei_unc
    }

# Convert and print the results
results_df = pd.DataFrame(results).T
print(results_df)
