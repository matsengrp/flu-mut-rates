import pandas as pd


file_paths = snakemake.input
outpath = snakemake.output.path

print(f"outpath is: {outpath}")

dfs = {}
roots = {}

for path in file_paths:
    segment = path.split("/")[-1].split(".")[0]
    the_df = pd.read_csv(path, sep="\t")
    the_df.query("~date.isna()", inplace=True)
    dfs[segment] = the_df

print("dataframe lengths:" + ", ".join(map(str, map(len, dfs.values()))))

strains = [set(df.strain) for df in dfs.values()]
common_strains = strains[0].intersection(*strains)
print("common strains:")
print(common_strains)
# No common strains

for segment, df in dfs.items():
    old_strain = df.query("date.str.startswith('1968')").strain.iat[0]
    roots[segment] = old_strain

with open(outpath, "w") as the_file:
    for segment, root in roots.items():
        if segment == "HA":
            # Hard code this root to check against Hugh's version.
            the_file.write(f"HA,A/Aichi/2/68_H3N2/1968|EF614251.1|1968\n")
        else:
            the_file.write(f"{segment},{root}\n")


# Based on the output, the roots are similar.
