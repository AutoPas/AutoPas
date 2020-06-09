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
    
def getStringFromLines(lines, regex):
    """extracts first capturing group of first line matching regex"""
    for line in lines:
        match = re.search(regex, line)
        if match:
            return match.group(1)
        
def getLabelDict(label):
    """
    Given label of form: { key1 : value1 , key2 : value 2 , ... }
    This function returns all key:value pairs as dict.
    """
    result = {}
    keyValues = re.search('{(.*)}', label)[1].split(',') # remove outer brackets and split by comma
    for keyValue in keyValues:
        key, value = keyValue.split(':')
        result[key.strip()] = value.strip()
        
    return result

def getID(label):
    """get a unique id for a label"""
    # if label not seen before, generate new id
    if not label in getID.labelDict:
        getID.labelDict[label] = len(getID.labelDict)
        
    return getID.labelDict[label]

getID.labelDict = {}
    

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

# ---------------------------------------------- Headers ---------------------------------------------

idStr = 'id'
timesetStr = 'timeset'
labelStr = 'Label'
sourceStr = 'Source'
targetStr = 'Target'

# ---------------------------------------------- Classes ---------------------------------------------

class GraphNode:
    """Class to store a node"""
    def __init__(self, label):
        self.id = str(getID(label))
        self.label = label
        self.timeValues = {}
        
    def getRow(self, nodeHeaders, labelHeaders):
        """Convert node to csv row"""
        # initialize empty list for each header
        columns = {header:[] for header in nodeHeaders}
            
        # collect time-value-pairs for each header
        for time,values in self.timeValues.items():
            for header,value in values.items():
                columns[header].append((time,time+1,value))
    
        # first rows: id, Label
        row = [self.id, self.label]
        
        # timeset format: <[time1,time1+1);[time2,time2+1);...;[timeN,timeN+1)>
        row.append('<{}>'.format(';'.join(['[{:d},{:d})'.format(time,time+1) for time in self.timeValues.keys()])))
        
        # extract label information
        configurations = getLabelDict(self.label)
        for header in labelHeaders:
            row.append(configurations[header])
        
        for header in nodeHeaders:
            # print each column in format: <[time1,time1+1,value1);[time2,time2+1);...;[timeN,timeN+1)>
            row.append('<{}>'.format(';'.join(['[{:d},{:d},{})'.format(*interval) for interval in columns[header]])))
            
        return row
        
class GraphEdge:
    """Class to store a edge"""
    def __init__(self, src, target):
        self.src = src
        self.target = target
        self.timeValues = {}
        
    def getRow(self, headers):
        """Convert edge to csv row"""
        # initialize empty list for each header
        columns = {header:[] for header in headers}
            
        # collect time-value-pairs for each header
        for time,values in self.timeValues.items():
            for header,value in values.items():
                columns[header].append((time,time+1,value))
    
        # first rows: Source, Target
        row = [str(self.src), str(self.target)]
        
        # timeset format: <[time1,time1+1);[time2,time2+1);...;[timeN,timeN+1)>
        row.append('<{}>'.format(';'.join(['[{:d},{:d})'.format(time,time+1) for time in self.timeValues.keys()])))

        for header in headers:
            # print each column in format: <[time1,time1+1,value1);[time2,time2+1);...;[timeN,timeN+1)>
            row.append('<{}>'.format(';'.join(['[{:d},{:d},{})'.format(*interval) for interval in columns[header]])))
            
        return row

# ---------------------------------------------- Script ---------------------------------------------

graphNodeStartStr = 'GaussianCluster Graph: Nodes'
graphEdgeStartStr = 'GaussianCluster Graph: Edges'
graphEndStr = 'GaussianCluster Graph: End'

graphNodeStartPos = []
graphEdgeStartPos = []
graphEndPos = []

# parse file
with open(graphFile) as f:
    lines = f.readlines()

# get lines containing markers
for i in range(len(lines)):
    if graphNodeStartStr in lines[i]:
        graphNodeStartPos.append(i)
    elif graphEdgeStartStr in lines[i]:
        graphEdgeStartPos.append(i)
    elif graphEndStr in lines[i]:
        graphEndPos.append(i)

# check that all markers appear equally often
if (not len(graphNodeStartPos) == len(graphEdgeStartPos) == len(graphEndPos)):
    print('Error: number of markers not equal! Node {:d}, Edge {:d}, End {:d}.'.format(len(graphNodeStartPos), len(graphEdgeStartPos), len(graphEndPos)))
    exit(1)
    
numGraphs = len(graphNodeStartPos)
    
# check that all markers repeatedly appear in order: node, edge, end
if (not all([(graphNodeStartPos[i] < graphEdgeStartPos[i] < graphEndPos[i]) for i in range(numGraphs)])) or (not all([(graphEndPos[i-1] < graphEdgeStartPos[i]) for i in range(1,numGraphs)])):
    print('Error: unexpected order of markers!')
    exit(1)

# ------ read all graphs ------
nodes = {}
nodeHeaders = set()
edges = {}
edgeHeaders = set()
for i in range(numGraphs):
    print('Graph found from line {:d} to {:d}'.format(graphNodeStartPos[i], graphEndPos[i]))
    
    # ----- read nodes -----
    nodesReader = csv.reader(lines[graphNodeStartPos[i]+1:graphEdgeStartPos[i]], delimiter=',', quotechar='"')
    
    # find index labels
    header = next(nodesReader)
    for c in range(len(header)):
        if header[c] == labelStr:
            labelIndex = c
        else:
            nodeHeaders.add(header[c])
    
    for row in nodesReader:
        if (row):
            # get values of current row
            values = { }
            for c in range(len(header)):
                value = row[c]
                
                # get Label seperately. Rest in dict values.
                if c == labelIndex:
                    label = value
                else:
                    values[header[c]] = value
            
            nodeId = getID(label)
            
            # create new node if not already otherwise check that label is the same
            if nodeId in nodes:
                if label != nodes[nodeId].label:
                    print('Error: label {} expected {}.'.format(label, nodes[nodeId].label))
                    exit(1)
            else:
                nodes[nodeId] = GraphNode(label)
                
            nodes[nodeId].timeValues[i] = values
    
    # ----- read edges -----
    edgeReader = csv.reader(lines[graphEdgeStartPos[i]+1:graphEndPos[i]], delimiter=',', quotechar='"')
    
    # find index of src and target
    header = next(edgeReader)
    for c in range(len(header)):
        if header[c] == sourceStr:
            srcIndex = c
        elif header[c] == targetStr:
            targetIndex = c
        else:
            edgeHeaders.add(header[c])
    
    for row in edgeReader:
        if (row):
            # get values of current row
            values = { }
            for c in range(len(header)):
                value = row[c]
                
                # get Source and Target seperately and convert them to unique ids. Rest in dict values.
                if c == srcIndex:
                    src = getID(value)
                elif c == targetIndex:
                    target = getID(value)
                else:
                    values[header[c]] = value
            
            # create new edge if not already exists
            if not (src,target) in edges:
                edges[src,target] = GraphEdge(src, target)
                
            edges[src,target].timeValues[i] = values
                
                
# convert header sets to lists
nodeHeaders = list(nodeHeaders)
edgeHeaders = list(edgeHeaders)

# Labels are used configuration: e.g. {Container: LinkedCells , CellSizeFactor: 2.000000 , Traversal: c08 , Data Layout: AoS , Newton 3: disabled}
# Extract header names
labelHeaders = list(getLabelDict(nodes[0].label).keys())
    
# ---- Write files ----

# create output directory
try:
    os.mkdir(outputDir)
except OSError:
    print('Could not create the output directory: ' + outputDir)
    exit(2)
    
if (fullSearchFile):
    # if fullSearchFile given, append time found in this file to each row
    
    with open(fullSearchFile) as f:
        fullSearchLines = f.readlines()
    
    with open(os.path.join(outputDir, 'graphNodes.csv'), 'w', newline='') as f:
        nodesWriter = csv.writer(f, delimiter=',', quotechar='"')
        nodesWriter.writerow([idStr, labelStr, timesetStr] + labelHeaders + nodeHeaders +['runtime'])
        for node in nodes.values():
            row = node.getRow(nodeHeaders, labelHeaders)
            
            # extract label from row and find corresponding runtime in fullSearchFile
            runtime = getStringFromLines(fullSearchLines, node.label + '.* Reduced value: ([0-9]+)')

            # append runtime to row if found otherwise discard row
            if runtime:
                row.append(runtime)
                nodesWriter.writerow(row)
            
else:
    # without fullSearchFile just write lines into file
    with open(os.path.join(outputDir, 'graphNodes.csv'), 'w', newline='') as f:
        nodesWriter = csv.writer(f, delimiter=',', quotechar='"')
        nodesWriter.writerow([idStr, labelStr, timesetStr] + labelHeaders + nodeHeaders)
        for node in nodes.values():
            nodesWriter.writerow(node.getRow(nodeHeaders, labelHeaders))

# write edges into file
with open(os.path.join(outputDir, 'graphEdges.csv'), 'w', newline='') as f:
    edgeWriter = csv.writer(f, delimiter=',', quotechar='"')
    edgeWriter.writerow([sourceStr, targetStr, timesetStr] + edgeHeaders)
    for edge in edges.values():
        edgeWriter.writerow(edge.getRow(edgeHeaders))