#!/bin/python3

import re
import plotly.graph_objects as go

datafile='testTuning_2020-02-22_23-44-40/MDFlex_Case1.out'

regex='.* Collected times for +({.*}).*\[(.*)\].*: *([0-9]+)'

# lists to gather data series
configs = []
times = []
timesErrorPlus = []
timesErrorMinus = []
values = []

# parse file
with open(datafile) as f:
    for line in f.readlines():
        match = re.search(regex, line)
        if match:
            configs.append(match.group(1))
            thisTimes=[int(s) for s in match.group(2).split()]
            value=int(match.group(3))
            times.append(thisTimes)
            timesErrorPlus.append(max(thisTimes)-value)
            timesErrorMinus.append(value-min(thisTimes))
            values.append(value)

fig = go.Figure(data=go.Scatter(
        x=[n*3 for n in range(0, len(configs))],
        y=values,
        error_y=dict(
            type='data',
            symmetric=False,
            array=timesErrorPlus,
            arrayminus=timesErrorMinus,
            )
        ))
fig.show()
