from os import makedirs
import numpy as np
import pandas as pd
from itertools import product, permutations
from sklearn.metrics import mean_squared_error

letters = ["A", "C", "G", "T"]


def load_synonymous_muts(subtype="H3N2"):
    # Read in curated dataset
    the_df = pd.read_csv("../results/rates.csv", keep_default_na=False)
    query = "mut_class=='synonymous' and subtype==@subtype and motif.str.len()==3"
    the_df.query(query, inplace=True)
    return the_df


class GeneralLinearModel:
    """
    ...

    Attributes:
        included_factors (???): ....local count, secondary structure, global positon, ...
        W (???): parameter values at training
        tau_squared (???): mse to log_counts at training

    """

    def __init__(self, included_factors):

        self.included_factors = included_factors
        self.W = {}
        self.tau_squared = {}

    def train(self, df_train):

        # Initialise dictionaries for model parameters and remaining variances
        self.W = {}
        self.tau_squared = {}

        # Check which mutation types are present in the training data
        present_mut_types = df_train.mut_type.unique()

        # Make a separate fit for every mutation type
        for mut_type in present_mut_types:
            # Select the current mutation type
            df_mut_type = df_train.query("mut_type==@mut_type")

            # Prepare log-counts and log-rates, dimensions: (# of sites, 1)
            log_counts = np.log(df_mut_type["actual_count"].values + 0.5).reshape(-1, 1)

            # Create data matrix X, dimensions: (# of sites, # of parameters in model)
            X = self.create_data_matrix(df_mut_type.copy(), mut_type)

            # Fit a linear model as log_counts = w @ X using a mean squared loss with a l2-regularization
            regularization_strength = 0.1
            regularization_matrix = regularization_strength * np.identity(X.shape[1])
            regularization_matrix[0, 0] = 0  # (don't regularize the offset term)
            w = np.linalg.inv(X.T @ X + regularization_matrix) @ X.T @ log_counts

            # Store model parameters
            self.W[mut_type] = w
            # Store remaining variance on log counts
            self.tau_squared[mut_type] = mean_squared_error(log_counts, X @ w)

        return None

    def create_data_matrix(self, mut_counts_df, mut_type=None):
        """
        ...mut_counts_df must have coluymns "motif"
        """

        factor_cols = {}

        # Get global context
        # if mut_type in ['AT', 'CG', 'GC']:
        #    factor_cols['global_context'] = [(mut_counts_df['nt_site'] > 21562).values.astype(int)]
        # elif mut_type in ['CT']:
        #    factor_cols['global_context'] = [(mut_counts_df['nt_site'] > 13467).values.astype(int)]
        # else:
        #    factor_cols['global_context'] = None

        # Get RNA structure
        # factor_cols['rna_structure'] = [mut_counts_df['unpaired'].values]

        # Get left and right context
        left_context = mut_counts_df["motif"].apply(lambda x: x[0]).values
        right_context = mut_counts_df["motif"].apply(lambda x: x[2]).values
        factor_cols["local_context"] = self.one_hot_l_r(left_context, right_context)

        # Offset term
        base = np.full(len(mut_counts_df), 1)

        # Add columns depending on which model is fitted and which factors are included
        columns = [base]
        for factor in self.included_factors:
            if factor_cols[factor] is not None:
                columns += factor_cols[factor]

        # Compile data matrix
        X = np.column_stack(columns)

        return X

    def test(self, df_test):

        mean_sq_errs = {}

        # Check which mutation types are present in the training data
        present_mut_types = df_test.mut_type.unique()

        # Compute mean squared error separately for every mutation type
        for mut_type in present_mut_types:
            # Select the current mutation type
            df_mut_type = df_test.query("mut_type==@mut_type")

            # Prepare log-counts, dimensions: (# of sites, 1)
            log_counts = np.log(df_mut_type["actual_count"].values + 0.5).reshape(-1, 1)

            # Create data matrix X, dimensions: (# of sites, # of parameters in model)
            X = self.create_data_matrix(df_mut_type.copy(), mut_type)

            # Compute the mean squared error of the fitted model on the training data
            mean_sq_errs[mut_type] = mean_squared_error(
                log_counts, (X @ self.W[mut_type]).flatten()
            )

        return mean_sq_errs

    def add_predicted_counts(self, df):

        # Check which mutation types are present in the data
        present_mut_types = df.mut_type.unique()

        # Prepare array for predicted counts and remaining variance
        predicted_counts = np.empty(len(df))
        remaining_variance = np.empty(len(df))

        # Compute mean squared error separately for every mutation type
        for mut_type in present_mut_types:
            # Select the current mutation type
            mask_mut_type = (df.mut_type == mut_type).values
            df_mut_type = df[mask_mut_type]

            # Create data matrix X, dimensions: (# of sites, # of parameters in model)
            X = self.create_data_matrix(df_mut_type.copy(), mut_type)

            # Compute the mean squared error of the fitted model on the training data
            predicted_counts[mask_mut_type] = (X @ self.W[mut_type]).flatten()
            remaining_variance[mask_mut_type] = self.tau_squared[mut_type]

        df["predicted_count"] = np.exp(predicted_counts) - 0.5
        df["tau_squared"] = remaining_variance

    @staticmethod
    def one_hot_l_r(context_l, context_r):
        sigma = {}
        for nt in ["C", "G", "T"]:
            sigma[nt + "_l"] = (context_l == nt).astype(int)
            sigma[nt + "_r"] = (context_r == nt).astype(int)

        return [
            sigma["C_l"],
            sigma["G_l"],
            sigma["T_l"],
            sigma["C_r"],
            sigma["G_r"],
            sigma["T_r"],
        ]

    def per_site_expected_counts(self):
        cols = ["mut_type", "motif"]
        rows = [
            (mut_type, f"{left}{mut_type[0]}{right}")
            for mut_type in map("".join, permutations(letters, 2))
            for left, right in product(letters, repeat=2)
        ]
        counts_df = pd.DataFrame(data=rows, columns=cols)
        self.add_predicted_counts(counts_df)
        return counts_df


def gather_and_write_counts(subtype, df_outpath, counts_outpath, summary_outpath):
    """..."""
    the_df = load_synonymous_muts(subtype=subtype)
    the_model = GeneralLinearModel(included_factors=["local_context"])
    the_model.train(df_train=the_df)

    the_model.add_predicted_counts(the_df)
    the_df.to_csv(df_outpath)
    counts_df = the_model.per_site_expected_counts()
    counts_df.to_csv(counts_outpath)

    header = "mut_type,mse,var,r2\n"
    lines_to_write = []
    mse_on_train = the_model.test(df_test=the_df)
    for mut_type, mse in mse_on_train.items():
        log_mut_count = np.log(the_df.query("mut_type==@mut_type").actual_count + 0.5)
        var = np.var(log_mut_count)
        r2 = 1 - mse / var
        lines_to_write.append(f"{mut_type},{mse},{var},{r2}\n")
    with open(summary_outpath, "w") as the_file:
        the_file.write(header)
        for line in lines_to_write:
            the_file.write(line)

    return None


if __name__ == "__main__":
    subtype = "H3N2"
    dir = f"../results/{subtype}"
    makedirs(dir, exist_ok=True)
    df_outpath = f"{dir}/full_df.csv"
    counts_outpath = f"{dir}/table.csv"
    summary_outpath = f"{dir}/summary.csv"
    gather_and_write_counts(subtype, df_outpath, counts_outpath, summary_outpath)
