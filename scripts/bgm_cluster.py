import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn import mixture
import argparse

def read_gc_nmer_output(file_path):
    data = pd.read_csv(file_path, sep='\t', header=None)
    reads_names = data.iloc[:, 0].values
    return data, reads_names

def preprocess_data(data):
    kmer_data = data.iloc[:, 3:].values
    return kmer_data

def doPCA(X, n_components=16):
    pca = PCA(n_components=n_components, whiten=True)
    X2 = pca.fit_transform(X)
    return X2, pca.explained_variance_ratio_

def doBGM(X2, n_components=20):
    dpgmm = mixture.BayesianGaussianMixture(n_components=n_components, covariance_type='full').fit(X2)
    predict_Y = dpgmm.predict(X2)
    return predict_Y

def main(file_path, output_path, n_pca_components=16, n_bgm_components=20):
    data, reads_names = read_gc_nmer_output(file_path)
    kmer_data = preprocess_data(data)
    X2, explained_variance_ratio = doPCA(kmer_data, n_components=n_pca_components)
    predict_Y = doBGM(X2, n_components=n_bgm_components)

    output_df = pd.DataFrame({'ReadName': reads_names, 'Cluster': predict_Y})
    output_df.to_csv(output_path, sep='\t', index=False, header=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PCA and Bayesian Gaussian Mixture Model Clustering")
    parser.add_argument('-i', '--input', type=str, required=True, help="Path to the input file")
    parser.add_argument('-o', '--output', type=str, required=True, help="Path to the output file")
    parser.add_argument('-x', '--pca_n', type=int, default=16, help="Number of PCA components")
    parser.add_argument('-n', '--bgm_n', type=int, default=20, help="Number of BGM components")
    
    args = parser.parse_args()
    main(args.input, args.output, args.pca_n, args.bgm_n)
