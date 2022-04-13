#!/usr/bin/env bash

# aggregate series
python aggregate_series.py
# produce base forecasts
Rscript base.R
# reconcile and evaluate forecasts
python evaluate.py ets -1
python evaluate.py ets 0
python evaluate.py ets_arima 0
python evaluate.py ets_arima -1
