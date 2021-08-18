# Gene expression filtration based on z-score
# Based on Hart et al, 2013 and a discussion here: https://www.biostars.org/p/94680/
# Aim: discarding gene expressions with low value
# Input: (Unfiltered) Gene expression table
# Output: Filtered gene expression table (expressed genes only)

import pandas as pd
from scipy.stats import gaussian_kde
import numpy as np

df_cells = pd.read_csv('gene_expression_data')

columns = df_cells.columns
rows = df_cells.index.values


def tpm_log2_converter(num):
    try:
        num=float(num)
    except:
        return 'NaN'
    if num > 0:
        tpm = float(num) * 100
        normalised_value = np.log2(tpm)
        return normalised_value


df_cells_log2 = pd.DataFrame(columns=columns, index=rows)
for column in columns:
    df_cells_log2[column] = df_cells[column].apply(tpm_log2_converter)


df_cells_log2_filtered = pd.DataFrame(columns=columns, index=rows)


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
    print(s)
    df_cells_log2_filtered[i] = s

df_cells_log2_filtered.to_csv('output_file')
