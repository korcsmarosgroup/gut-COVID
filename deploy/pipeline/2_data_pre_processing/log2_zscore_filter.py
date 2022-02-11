# Gene expression filtration based on z-score
# Based on Hart et al, 2013 and a discussion here: https://www.biostars.org/p/94680/
# Aim: discarding gene expressions with low value
# Input: (Unfiltered) Gene expression table
# Output: Filtered gene expression table (expressed genes only)

import pandas as pd
from scipy.stats import gaussian_kde
import numpy as np
import argparse
import sys
import logging
import os

if os.path.isfile('log2_zscore_filter.log'):
    os.remove("log2_zscore_filter.log")

logging.basicConfig(filename = 'log2_zscore_filter.log', format = '%(asctime)s :: %(levelname)s :: %(funcName)s :: %(lineno)d :: %(message)s', level = logging.INFO)


def parse_args(args):
    help_text = \
        """
        === Log2 Zscore Fitlering ===
        
        Gene expression filtration based on z-score.
        """

    parser = argparse.ArgumentParser(description = help_text)

    parser.add_argument("-i", "--input-file",
                        help="<path to the input file> [mandatory]",
                        type=str,
                        dest="input_file",
                        action="store",
                        required=True)

    parser.add_argument("-o", "--output-file",
                        help="<path to the output file> [mandatory]",
                        type=str,
                        dest="output_file",
                        action="store",
                        required=True)

    results = parser.parse_args(args)
    return results.input_file, results.output_file


def tpm_log2_converter(num):
    try:
        num = float(num)
    except:
        return 'NaN'
    if num > 0:
        tpm = float(num) * 100
        normalised_value = np.log2(tpm)
        return normalised_value

logging.info("Log2 Zscore filter is starting")
input_file, output_file = parse_args(sys.argv[1:])

df_cells = pd.read_csv(input_file, sep = '\t')
logging.info("The input file is fine")

columns = df_cells.columns
rows = df_cells.index.values

df_cells_log2 = pd.DataFrame(columns = columns, index = rows)
logging.info("Log2 the data is done")
for column in columns:
    df_cells_log2[column] = df_cells[column].apply(tpm_log2_converter)

df_cells_log2_filtered = pd.DataFrame(columns = columns, index = rows)
logging.info("Filtered the log2 the data is done")

for i in columns:

    fpkm = df_cells_log2[i].tolist()

    fpkm = np.array(fpkm)
    fpkm_filtered = fpkm[np.logical_not(np.isnan(fpkm))]

    # creating the Gauss-curve
    kernel = gaussian_kde(fpkm_filtered)

    # creating the X axis -> divide the list for 100 units, xi numpy lists contains that 100 values
    # expected value = most oftest value, middle of Gaus-curve
    xi = np.linspace(fpkm_filtered.min(), fpkm_filtered.max(), 100)

    # calculate y for each x point
    yi = kernel.evaluate(xi)

    # expected value calculation, which x is by the max y value? (np.argmax(yi) = position)
    mu = xi[np.argmax(yi)]

    # fpkm > mu  = list of boolean values; mean of values right from the expected value
    U = fpkm_filtered[fpkm_filtered > mu].mean()

    # calculation of standard deviation
    sigma = (U - mu) * np.sqrt(np.pi / 2)

    # new score: deviancy from the mean divided by sigma (standard deviation)
    # z-value: relative value - deviation from the mean in the st.dev in the data -> Gaus-curve:  0.001 - 3 deviations
    zFPKM = (fpkm - mu) / sigma

    score_list = [fpkm[list(zFPKM).index(x)] if x > -3 else 'NaN' for x in zFPKM]

    s = pd.Series(score_list, index=rows)
    logging.info(f"{i} is done")
    df_cells_log2_filtered[i] = s

df_cells_log2_filtered.to_csv(output_file)
logging.info("Log2 Zscore filter is done")
