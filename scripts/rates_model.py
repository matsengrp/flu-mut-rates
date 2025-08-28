from os import makedirs
import numpy as np
import pandas as pd
from itertools import product, permutations
from sklearn.metrics import mean_squared_error

letters = ["A", "C", "G", "T"]

# This script predicts rates without a pseudorate, so the rates at training and testing
# need to be positive. The model trains in log-space and the reported MSE and R^2 are
# for log rates.


def load_synonymous_muts(subtype=None):
    """
    Return a pandas dataframe for synonymous mutations for the given flu subtype at
    nucleotide sites with full length motifs (exclude the nucleotides sites at the start
    and end of a coding region). Check flu-syn-rates/results/syn_rates.csv for possible
    subtypes, use None for no restriction on subtype.
    """
    # Read in curated dataset. Unfortunately "NA" is a segment.
    the_df = pd.read_csv("../results/syn_rates.csv", keep_default_na=False)
    if subtype is not None:
        query = "mut_class=='synonymous' and subtype==@subtype and motif.str.len()==3"
    else:
        query = "mut_class=='synonymous' and motif.str.len()==3"
    the_df.query(query, inplace=True)
    return the_df


class GeneralLinearModel:
    """
    A log-linear model for predicting rates by least squares with ridge regression.

    Attributes:
        included_factors (list): A list of model factors, which should be some
            combination of "local_context", "global_context", and "rna_structure" (not
            yet implemented).
        W (dict): A dictionary mapping a mutation type to the model coefficients.
        tau_squared (dict): A dictionary mapping a mutation type to the MSE at training for log rates.
    """

    def __init__(self, included_factors):
        self.included_factors = included_factors
        self.W = {}
        self.tau_squared = {}

    def train(self, df_train):
        """
        Fit the model to the training data and fill the dictionaries self.W (model
        parameter values by mutation type) and self.tau_squareed (MSE by mutation type).

        Parameters:
            df_train  (pd.DatafFrame): A dataframe with columns "mut_type" and "rate".
                Additional columns are required based on self.included factors:
                "motif" for "local_context", "segment" for "global_context", and ???
                for "rna_structure".
        """
        # Initialise dictionaries for model parameters and remaining variances
        self.W = {}
        self.tau_squared = {}

        # Check which mutation types are present in the training data
        present_mut_types = df_train.mut_type.unique()

        # Make a separate fit for every mutation type
        for mut_type in present_mut_types:
            # Select the current mutation type
            df_mut_type = df_train.query("mut_type==@mut_type")

            # Prepare log-rates, dimensions: (# of sites, 1)
            log_rates = np.log(df_mut_type["rate"].values).reshape(-1, 1)

            # Create data matrix X, dimensions: (# of sites, # of parameters in model)
            X = self.create_data_matrix(df_mut_type.copy(), mut_type)

            # Fit a linear model as log_rates = w @ X using a mean squared loss with a
            # l2-regularization.
            regularization_strength = 0.1
            regularization_matrix = regularization_strength * np.identity(X.shape[1])
            regularization_matrix[0, 0] = 0  # (don't regularize the offset term)
            w = np.linalg.inv(X.T @ X + regularization_matrix) @ X.T @ log_rates

            # Store model parameters
            self.W[mut_type] = w
            # Store remaining variance on log counts
            self.tau_squared[mut_type] = mean_squared_error(log_rates, X @ w)

        return None

    def create_data_matrix(self, mut_counts_df):
        """
        Construct a numpy matrix, based on mut_counts_df, with the predictor variables
        of the model.

        Parameters:
            df_train  (pd.DatafFrame): Columns are required based on self.included factors:
                "motif" for "local_context", "segment" for "global_context", and ???
                for "rna_structure".
        """

        factor_cols = {
            "local_context": None,
            "global_context": None,
            "rna_structure": None,
        }

        if "local_context" in self.included_factors:
            left_context = mut_counts_df["motif"].apply(lambda x: x[0]).values
            right_context = mut_counts_df["motif"].apply(lambda x: x[2]).values
            factor_cols["local_context"] = self.one_hot_l_r(left_context, right_context)
        if "global_context" in self.included_factors:
            factor_cols["global_context"] = self.one_hot_segment(mut_counts_df.segment)
        if "rna_structure" in self.included_factors:
            raise NotImplementedError("We don't yet have rna secondary structure.")

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
        """
        Return a dictionary mapping mutation type to MSE for log rates.

        Parameters:
            df_test (pd.DatafFrame): A dataframe with columns "mut_type" and "rate".
                Additional columns are required based on self.included factors:
                "motif" for "local_context", "segment" for "global_context", and ???
                for "rna_structure".


        """
        mean_sq_errs = {}

        # Check which mutation types are present in the training data
        present_mut_types = df_test.mut_type.unique()

        # Compute mean squared error separately for every mutation type
        for mut_type in present_mut_types:
            # Select the current mutation type
            df_mut_type = df_test.query("mut_type==@mut_type")

            # Prepare log-counts, dimensions: (# of sites, 1)
            log_rates = np.log(df_mut_type["rate"].values).reshape(-1, 1)

            # Create data matrix X, dimensions: (# of sites, # of parameters in model)
            X = self.create_data_matrix(df_mut_type.copy(), mut_type)

            # Compute the mean squared error of the fitted model on the training data
            mean_sq_errs[mut_type] = mean_squared_error(
                log_rates, (X @ self.W[mut_type]).flatten()
            )

        return mean_sq_errs

    def add_predicted_rates(self, df):
        """
        Add columns "predicted_rate" (rate not log-rate) and "tau_squared" (MSE for
        log-rate, not rate) to the DataFrame.

        Parameters:
            df_train  (pd.DatafFrame): A dataframe with column "mut_type". Additional
                columns are required based on self.included factors: "motif" for
                "local_context", "segment" for "global_context", and ??? for
                "rna_structure".
        """
        # Check which mutation types are present in the data
        present_mut_types = df.mut_type.unique()

        # Prepare array for predicted counts and remaining variance
        predicted_rates = np.empty(len(df))
        remaining_variance = np.empty(len(df))

        # Compute mean squared error separately for every mutation type
        for mut_type in present_mut_types:
            # Select the current mutation type
            mask_mut_type = (df.mut_type == mut_type).values
            df_mut_type = df[mask_mut_type]
            df_mut_type = df.query("mut_type==@mut_type")

            # Create data matrix X, dimensions: (# of sites, # of parameters in model)
            X = self.create_data_matrix(df_mut_type.copy(), mut_type)

            # Compute the mean squared error of the fitted model on the training data
            predicted_rates[mask_mut_type] = (X @ self.W[mut_type]).flatten()
            remaining_variance[mask_mut_type] = self.tau_squared[mut_type]

        df["predicted_rate"] = np.exp(predicted_rates)
        df["tau_squared"] = remaining_variance

    def per_site_expected_rates(self):
        """
        Return a DataFrame listing the predicted rate for each combination of predictor
        variable values (depending on self.included_factors).
        """
        segments = ["PB1", "NS", "NA", "HA", "PB2", "MP", "PA", "NP"]
        mut_types = list(self.W.keys())

        cols = ["mut_type"]
        thing = [mut_types]
        if "local_context" in self.included_factors:
            cols.append("left")
            cols.append("right")
            thing.append(letters)
            thing.append(letters)
        if "global_context" in self.included_factors:
            cols.append("segment")
            thing.append(segments)

        rows = list(product(*thing))
        rates_df = pd.DataFrame(data=rows, columns=cols)
        if "local_context" in self.included_factors:
            rates_df["motif"] = (
                rates_df.left + rates_df.mut_type.str.get(0) + rates_df.right
            )
            rates_df.drop(columns=["left", "right"], inplace=True)

        self.add_predicted_rates(rates_df)
        return rates_df

    @staticmethod
    def one_hot_l_r(context_l, context_r):
        """
        Return a one hot encoding of the left neighboring nucleotide and the right
        neighboring nucleotide.
        """
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

    @staticmethod
    def one_hot_segment(segments):
        """
        Return a one hot encoding of the segements.
        """
        segment_order = ["PB1", "NS", "NA", "HA", "PB2", "MP", "PA", "NP"]
        return [(segments == seg).values.astype(int) for seg in segment_order[1:]]


def gather_and_write_rates(
    subtype, model_factors, df_outpath, rates_outpath, summary_outpath
):
    """
    Report rates and summary statistics for a GeneralizedLinearModel fit to data for the
    given subtype.

    Paraemeters:
        subtype (str): The flu subtype. Use None to include all data.
        model_factors (list): A list of model factors, which should be some combination
            of "local_context", "global_context", and "rna_structure" (not yet
            implemented).
        df_outpath (str): The path to write a csv of the training data and redicted rates.
        rates_outpath (str): The path to write a csv of the per site rate for each
            combination of predictor variable values.
        summary_outpath (str): The path to write a csv with MSE, variance, and R^2 for
            each mutation type.
    """
    the_df = load_synonymous_muts(subtype=subtype)
    the_model = GeneralLinearModel(included_factors=model_factors)
    the_model.train(df_train=the_df)

    the_model.add_predicted_rates(the_df)
    the_df.to_csv(df_outpath)
    rates_df = the_model.per_site_expected_rates()
    rates_df.to_csv(rates_outpath)

    header = "mut_type,mse,var,r2\n"
    lines_to_write = []
    mse_on_train = the_model.test(df_test=the_df)
    for mut_type, mse in mse_on_train.items():
        log_mut_rate = np.log(the_df.query("mut_type==@mut_type").rate)
        var = np.var(log_mut_rate)
        r2 = 1 - mse / var
        lines_to_write.append(f"{mut_type},{mse},{var},{r2}\n")
    with open(summary_outpath, "w") as the_file:
        the_file.write(header)
        for line in lines_to_write:
            the_file.write(line)

    return None


if __name__ == "__main__":
    model_hierarchy = [[], ["local_context"], ["local_context", "global_context"]]

    for model_factors in model_hierarchy:
        dir = f"../results/all_rows/"
        dir += "base" if len(model_factors) == 0 else "+".join(model_factors)
        makedirs(dir, exist_ok=True)
        df_outpath = f"{dir}/full_df.csv"
        counts_outpath = f"{dir}/table.csv"
        summary_outpath = f"{dir}/summary.csv"
        gather_and_write_rates(
            None, model_factors, df_outpath, counts_outpath, summary_outpath
        )
