#!/usr/bin/python3

import os
import sys
import csv
import re
from datetime import datetime

# -------------------------------------------- Functions --------------------------------------------

def printHelp():
    print("Usage: ./extractClusterGraph.py path/To/mdFlex/std.out [path/To/mdFelx/fullsearch/std.out]")
    print("If FullSearch is provided, measured runtimes will be added to the nodes")
    exit(0)
    
# extracts first capturing group of first line matching regex
def getStringFromLines(lines, regex):
    for line in lines:
        match=re.search(regex, line)
        if match:
            return match.group(1)

# ---------------------------------------------- Input ----------------------------------------------

# help message
for arg in sys.argv[1:]:
    if "--help" in arg:
        printHelp()

# first file contains graph
if len(sys.argv) > 1:
    graphFile = sys.argv[1]
else:
    printHelp()

# second file contains optional full search
fullSearchFile = None
if len(sys.argv) > 2:
    fullSearchFile = sys.argv[2]

# directory for graph output
outputDir="clusterGraph_"+datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


# ---------------------------------------------- Script ---------------------------------------------

graphNodeStartStr = 'GaussianCluster Graph: Nodes'
graphEdgeStartStr = 'GaussianCluster Graph: Edges'
graphEndStr = 'GaussianCluster Graph: End'

graphNodeStartPos = 0
graphEdgeStartPos = 0
graphEndPos = 0

# parse file
with open(graphFile) as f:
    lines = f.readlines()

# get lines containing markers
for i in range(len(lines)):
    if graphNodeStartStr in lines[i]:
        graphNodeStartPos = i
    elif graphEdgeStartStr in lines[i]:
        graphEdgeStartPos = i
    elif graphEndStr in lines[i]:
        graphEndPos = i

if (not graphNodeStartPos < graphEdgeStartPos < graphEndPos):
    print('File did not contain a valid graph!')
    exit(1)

print('Graph found from line %i to %i' % (graphNodeStartPos, graphEndPos))

# create output directory
try:
    os.mkdir(outputDir)
except OSError:
    print("Could not create the output directory: " + outputDir)
    exit(2)
    
if (fullSearchFile):
    # if fullSearchFile given, append time found in this file to each row
    
    with open(fullSearchFile) as f:
        fullSearchLines = f.readlines()
    
    nodesReader = csv.reader(lines[graphNodeStartPos+1:graphEdgeStartPos], delimiter=',', quotechar='"')
    
    # find index of labels and append runtime column
    header = next(nodesReader)
    labelIndex = header.index('Label')
    header.append('runtime')
    
    with open(os.path.join(outputDir, 'graphNodes.csv'), 'w', newline='') as f:
        nodesWriter = csv.writer(f, delimiter=',', quotechar='"')
        nodesWriter.writerow(header)
        for row in nodesReader:
            # extract label from row and find corresponding runtime in fullSearchFile
            runtime = getStringFromLines(fullSearchLines, row[labelIndex] + ".* Reduced value: ([0-9]+)")
            
            # append runtime to row if found
            if runtime:
                row.append(runtime)
            else:
                row.append('No time found')
            
            nodesWriter.writerow(row)
else:
    # without fullSearchFile just copy lines into file
    with open(os.path.join(outputDir, 'graphNodes.csv'), 'w') as f:
        for line in lines[graphNodeStartPos+1:graphEdgeStartPos]:
            f.write(line)

# copy edges into file
with open(os.path.join(outputDir, 'graphEdge.csv'), 'w') as f:
    for i in lines[graphEdgeStartPos+1:graphEndPos]:
        f.write(lines[i])
