import numpy as np
import pandas as pd
import time
import sys
import scipy.sparse as sp

np.random.seed(42)

file_path = sys.argv[1]

num_samples = 10000

tcgamutsfile = "../analyses/" + file_path + "/tcga-muts.csv"
quanumfile = "../analyses/" + file_path + "/quanum.csv"
outputfile = "../analyses/" + file_path + "/main-ci.csv"

tcga_muts = pd.read_csv(tcgamutsfile)
quanum = pd.read_csv(quanumfile)

A_excluded = tcga_muts.iloc[:, 1:].to_numpy()
A_sparse = sp.csr_matrix(A_excluded)

rows, cols = A_sparse.nonzero()  
nonzero_values = A_sparse.data

resampled_values = np.random.poisson(lam=nonzero_values[np.newaxis, :], size=(num_samples, len(nonzero_values)))

resampled_sparse_matrices = [
    sp.csr_matrix((resampled_values[i], (rows, cols)), shape=A_sparse.shape)
    for i in range(num_samples)
]

quanum_excluded = quanum.iloc[:, 1:].to_numpy()

results = [A.dot(quanum_excluded) for A in resampled_sparse_matrices]

results = np.hstack(results)

quantiles = np.quantile(results, [0.025, 0.975], axis=1)

quantiles = quantiles.T

first_column = tcga_muts["id"].to_numpy().reshape(-1,1)

result_matrix = np.hstack((first_column, quantiles))

np.savetxt(outputfile, result_matrix.astype(str), delimiter=",", fmt = "%s")


































