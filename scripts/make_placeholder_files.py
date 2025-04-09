from itertools import permutations as perm


bad_sites_paths = snakemake.output.bad_sites_paths
secondary_structs = snakemake.output.secondary_structs
rates_tables = snakemake.output.rates_tables


def write_files(bad_sites_path, secondary_struct, rates_table):
    """
    Prepare temporary files used by ExpectedCounts.PossibleMutations. We will have a
    rates table after implementing context aware rates, we might have bad_sites, but
    we probably will not have secondary structure.
    """
    with open(bad_sites_path, "w") as the_file:
        the_file.write("-1\n")
    with open(secondary_struct, "w") as the_file:
        the_file.write("site\n")
    with open(rates_table, "w") as the_file:
        to_write = "mut_type,nt_site_boundary\n"
        to_write += "".join((f"{mt},0\n" for mt in map("".join, perm("ACGT", 2))))
        the_file.write(to_write)
    return None


any(write_files(*p) for p in zip(bad_sites_paths, secondary_structs, rates_tables))
