#!/usr/bin/python3

import re
import plotly.graph_objects as go

datafile='testTuning_2020-02-24_14-11-02/MDFlex_Case1.out'

regex='.* Collected times for +({.*}).*\[(.*)\].*: *([0-9]+)'

# ------------------------------------------------- Functions -------------------------------------------------

def parseConfigToDict(confStr):
    return {key.strip():val.strip() for key,val in [pair.split(":") for pair in confStr[1:-1].split(',')]}

# --------------------------------------------------- Script --------------------------------------------------

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
            configs.append(parseConfigToDict(match.group(1)))
            thisTimes=[int(s) for s in match.group(2).split()]
            value=int(match.group(3))
            times.append(thisTimes)
            timesErrorPlus.append(max(thisTimes)-value)
            timesErrorMinus.append(value-min(thisTimes))
            values.append(value)

# create figure and define layout
fig = go.Figure(
        layout=dict(
                showlegend=True,
            ),
        )

# add data
fig.add_trace(go.Scatter(
            x=[n*3 for n in range(0, len(configs))],
            y=values,
            mode="lines+markers",
            error_y=dict(
                type='data',
                symmetric=False,
                array=timesErrorPlus,
                arrayminus=timesErrorMinus,
                ),
            hovertext=configs,
            marker=dict(
                color=["#%0.6X" % (hash(c["Container"]) % (256**3)) for c in configs],
                showscale=True,
                ),
            line=dict(
                color='#AAAAAA',
                dash='dash',
            ),
        ))
fig.show()
