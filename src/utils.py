import numpy as np
import pandas as pd
import networkx as nx
from statsmodels.stats.multitest import fdrcorrection

def seq_KNN(x, k=10):
    x = x.copy()

    is_miss = x.isna().sum(axis=1).sort_values() > 0
    x_bad = x.loc[is_miss, :]
    
    x_hat = x.loc[~is_miss, :].copy()
    for i in x_bad.index:
        miss_entry = x_bad.loc[i, :].copy()
        is_miss = x_bad.loc[i].isna()
        dist = ((x_hat.loc[:, ~is_miss] - miss_entry[~is_miss])**2).sum(axis=1).sort_values().iloc[:k]
        weight = 1/(dist+0.000000000000001)
        weight /= weight.sum()
        miss_entry[is_miss] = (x_hat.loc[weight.index, is_miss].T * weight).sum(axis=1)
        x_hat = pd.concat([x_hat, pd.DataFrame(miss_entry).T])
    
    return x_hat

def calc_strength_fraction(adj, region_group):
    W = 0
    for set_S in region_group:
        W += adj[set_S, :][:, set_S].sum()
    return W/(adj.sum()-W)

def calc_conductance(nx_adj, region_group):
    cond = 0
    for set_S in region_group:
        cond += nx.algorithms.conductance(nx_adj, set_S, weight="weight")
    return -cond

def calc_normalized_cut(nx_adj, region_group):
    ncut = 0
    for set_S in region_group:
        ncut += nx.algorithms.normalized_cut_size(nx_adj, set_S, weight="weight")
    return -ncut

def significant_test_for_index(corr_mtx, region_group, n_perm=1000):
    # 0. data preparation
    corr_array = np.array(corr_mtx)
    np.fill_diagonal(corr_array, 0)
    corr_array_nx = nx.from_numpy_array(corr_array).to_undirected()

    # 1. calc index
    seg_index = np.array([calc_strength_fraction(corr_array, region_group), 
                          calc_conductance(corr_array_nx, region_group), 
                          calc_normalized_cut(corr_array_nx, region_group)])

    # 2. permutation
    perm_index = np.empty(shape=(n_perm, 3))
    for k in range(n_perm):
        tmp_idx = np.random.permutation(corr_array.shape[0])
        perm_corr = corr_array[tmp_idx, :][:, tmp_idx]
        perm_corr_g = nx.from_numpy_array(perm_corr).to_undirected()

        perm_index[k, 0] = calc_strength_fraction(perm_corr, region_group)
        perm_index[k, 1] = calc_conductance(perm_corr_g, region_group)
        perm_index[k, 2] = calc_normalized_cut(perm_corr_g, region_group)

    pval = ((perm_index>seg_index).sum(axis=0)+1)/n_perm

    # 3. multiple test correction
    _, fdr_bh = fdrcorrection(pval, alpha=0.05, method="indep", is_sorted=False)

    # 4. results
    seg_index[1:] = -seg_index[1:]
    perm_index[:, 1:] = -perm_index[:, 1:]

    return seg_index, perm_index, pval, fdr_bh
