import os
import pandas as pd

main_dir = "/ceph/project/milne_group/cchahrou/Shared/Ana/spike_in/2023-09-13_fly"
split_dir = os.path.join(main_dir, "split_paired_primary/split")

# List all .txt files from split_dir
files = [f for f in os.listdir(split_dir) if f.endswith(".txt")]
files.sort()

dfs = []
for report in files:
    file_path = os.path.join(split_dir, report)
    df = pd.read_csv(file_path, sep="\t")
    if (
        "Sample" not in df.columns
        or "n_sample" not in df.columns
        or "n_exogenous" not in df.columns
    ):
        ConnectionRefusedError
    dfs.append(df)

if dfs:
    scale_factors = pd.concat(dfs, ignore_index=True)
    scale_factors.to_csv(
        os.path.join(split_dir, "scale_factors.csv"),
        index=False,
    )
else:
    print("No valid data found in the input files.")

# subset columns of interest
scale_factors = scale_factors[["Sample", "n_sample", "n_exogenous"]]

# Calculate the ChIP spike-in normalization factor
scale_factors["chip_spike_in_norm_factor"] = 1 / (scale_factors["n_exogenous"] / 1e6)
scale_factors["exo_perc"] = (
    scale_factors["n_exogenous"] / scale_factors["n_sample"] * 100
)

# Save the DataFrame with the calculated normalization factors
scale_factors.to_csv(
    os.path.join(split_dir, "scale_factors_with_norm_factors.csv"),
    index=False,
)

scale_factors.to_csv(
    os.path.join(main_dir, "scale_factors.csv"),
    index=False,
)
