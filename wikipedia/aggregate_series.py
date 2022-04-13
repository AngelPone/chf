import pandas as pd
from pyhts.hierarchy import Hierarchy
data = pd.read_csv('wikimedia2016.csv')

ht = Hierarchy.from_long(data, list(data.columns[:11]), period=7)

hist_bts = data.iloc[:, 11:(11+394)].values.T
hist_hierarchy = ht.aggregate_ts(hist_bts)
pd.DataFrame(hist_hierarchy).to_csv('wikimedia2016_hierarchy.csv', index=False)
