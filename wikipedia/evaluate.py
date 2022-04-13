import pandas as pd
import numpy as np
from pyhts.reconciliation import wls
from pyhts.hierarchy import Hierarchy
import sys
import pathlib

current_path = pathlib.Path(__file__).parent
method = sys.argv[1]
clevel = int(sys.argv[2])
# method = 'ets'
# clevel = -1

output_file = current_path / "evaluate.csv"
input_ts = current_path / "../wikimedia2016.csv"
input_base = current_path / "base.csv"
measure = ['rmse']


def reconcile(method="ols", weighting=None, clevel=clevel):
    G = wls(ht, error.T, method, weighting, constraint_level=clevel)
    return G.dot(forecast.T)


def summary(acc):
    acc['level'] = ht.node_level
    summ = acc.groupby('level').apply(np.mean)
    summ['level'] = ht.level_name
    return summ


dataset = pd.read_csv(input_ts)
ts = dataset.iloc[:, 11:].values.T
ht = Hierarchy.from_long(dataset, list(dataset.columns[:11]), period=7)
f = pd.read_csv(input_base)
f['index'] = f['index'].astype('int')

if method in ['arima', 'ets']:
    fitted = f[f['method']==method].sort_values(by='index').iloc[:, :366].values.T
    forecast = f[f['method']==method].sort_values(by='index').iloc[:, 366:(366+28)].values.T
elif method == 'ets_arima':
    fitted = pd.concat([f[(f['method']=='ets') & (f['index']==0)], f[(f['method']=='arima') & (f['index'] > 0)]],axis=0).iloc[:, :366].values.T
    forecast = pd.concat([f[(f['method']=='ets') & (f['index']==0)], f[(f['method']=='arima') & (f['index'] > 0)]],axis=0).iloc[:, 366:(366+28)].values.T
elif method == 'arima_ets':
    fitted = pd.concat([f[(f['method']=='arima') & (f['index']==0)], f[(f['method']=='ets') & (f['index'] > 0)]],axis=0).iloc[:, :366].values.T
    forecast = pd.concat([f[(f['method']=='arima') & (f['index']==0)], f[(f['method']=='ets') & (f['index'] > 0)]],axis=0).iloc[:, 366:(366+28)].values.T
else:
    raise ValueError("method is not supported!")



train = ts[:(-28), :]
test = ts[(-28):, :]
train_all = ht.aggregate_ts(train)
test_all = ht.aggregate_ts(test)


error = fitted - train_all

ols = reconcile("ols").T # h * m
ols_all = ht.aggregate_ts(ols) # h * n

accs = summary(ht.accuracy(test, ols, hist=train, measure=measure))
accs['method'] = 'ols'
mint_shrink = reconcile("mint", "shrinkage").T
shrink_all = ht.aggregate_ts(mint_shrink)

acc = summary(ht.accuracy(test, mint_shrink, hist=train, measure=measure))
acc['method'] = 'shrinkage'
accs = pd.concat([accs, acc], axis=0)
wlss = reconcile("wls", "structural").T
wlss_all = ht.aggregate_ts(wlss)
acc = summary(ht.accuracy(test, wlss, train, measure=measure))
acc['method'] = 'structural'
accs = pd.concat([accs, acc], axis=0)
mintvar = reconcile("mint", "variance").T
var_all = ht.aggregate_ts(mintvar)


acc = summary(ht.accuracy(test, mintvar, train, measure=measure))
acc['method'] = 'variance'
accs = pd.concat([accs, acc], axis=0)

acc = summary(ht.accuracy_base(test, forecast, train, measure=measure))
acc['method'] = 'base'
accs = pd.concat([accs, acc], axis=0)
accs['clevel'] = clevel
accs['base'] = method

accs = accs[accs.level.isin(['Total', 'language', 'access', 'agent', 'purpose', 'network','bottom'])]

accs.to_csv(output_file, mode='a')






