import pandas as pd

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




