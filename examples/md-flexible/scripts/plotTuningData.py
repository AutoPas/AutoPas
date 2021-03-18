#!/usr/bin/python3

import sys
import os.path
import pandas as pd
import plotly.express as px

if len(sys.argv) < 2:
    sys.exit('Usage: ' + os.path.basename(sys.argv[0]) + ' IterationPerformance.csv')
csvFile = sys.argv[1]
if not os.path.isfile(csvFile) :
    sys.exit('File not found: ' + csvFile)

data = pd.read_csv(csvFile)

fig = px.scatter(data,
        x = 'Iteration',
        y = 'Smoothed',
        labels = {'Smoothed' : 'Smoothed iteratePairwise() [ns]'},
        title = 'Smoothed Tuning Samples',
        hover_name = 'Traversal',
        hover_data = ['Data Layout', 'Newton 3', 'Load Estimator'],
        color = 'Container',
        )
fig.update_traces(marker=dict(size=6))
fig.show()
