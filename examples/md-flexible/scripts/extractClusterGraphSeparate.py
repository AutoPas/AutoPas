#!/usr/bin/python3

import os
import sys
import csv
import re
import math
from datetime import datetime

# THIS SCRIPT NEEDS AT LEAST PYTHON 3.8.
assert sys.version_info >= (3, 8)

# Extracts graphs generated by GaussianCluster.
# The graphs is formatted into two files: nodesFile and edgesFile containing data for nodes and edges respectively.
#
# Output of GaussianCluster should contain multiple graphs marked by graphNodeStartStr, graphEdgeStartStr and graphEndStr.
# Between the lines containing graphNodeStartStr and graphEdgeStartStr is a valid csv-file containg all nodes.
# Between the lines containing graphEdgeStartStr and graphEndStr is a valid csv-file containg all edges.
#
# For each graph in the output a pair of nodesFile and edgesFile is created.

# ---------------------------------------------- Constants ---------------------------------------------

# marks the start of nodes data
graphNodeStartStr = 'GaussianCluster Graph: Nodes'
# marks the end of nodes and start of edge data
graphEdgeStartStr = 'GaussianCluster Graph: Edges'
# marks end of edge data
graphEndStr = 'GaussianCluster Graph: End'

# special headers
idStr = 'id'
labelStr = 'Label'
sourceStr = 'Source'
targetStr = 'Target'
fullSearchTimeStr = 'runtime'
fullSearchWeightStr = 'runtimeWeight'

# output file names
nodesFile = 'graph{}_nodes.csv'
edgesFile = 'graph{}_edges.csv'


# -------------------------------------------- Functions --------------------------------------------

def printHelp():
    print("Usage: ./extractClusterGraph.py path/To/mdFlex/std.out [path/To/mdFelx/fullsearch/std.out]")
    print("If FullSearch is provided, measured runtimes will be added to the nodes")
    sys.exit(0)


def getStringFromLines(lines, regex):
    """extracts first capturing group of first line matching regex"""
    for line in lines:
        match = re.search(regex, line)
        if match:
            return match.group(1)


def getLabelDict(label):
    """
    Given label of form: { key1 : value1 , key2 : value 2 , ... }
    e.g. {Container: LinkedCells , CellSizeFactor: 2.000000 , Traversal: c08 , Data Layout: AoS , Newton 3: disabled}
    This function returns all key:value pairs as dict.
    """
    result = {}
    keyValues = re.search('{(.*)}', label)[1].split(',')  # remove outer brackets and split by comma
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
outputDir = "clusterGraph_" + datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


# ---------------------------------------------- Classes ---------------------------------------------

class GraphNode:
    """Class to store a node"""

    def __init__(self, label, values):
        self.id = str(getID(label))
        self.label = label
        self.values = values

    def getRow(self, nodeHeaders, labelHeaders):
        """
        Convert node to csv row.
        Columns: id, Label, configurations... , remaining headers...
        """

        # first rows: id, Label
        row = [self.id, self.label]

        # extract label information
        configurations = getLabelDict(self.label)
        for header in labelHeaders:
            row.append(configurations[header])

        # get value for each header
        for header in nodeHeaders:
            row.append(self.values[header])

        return row


class GraphEdge:
    """Class to store a edge"""

    def __init__(self, src, target, values):
        self.src = src
        self.target = target
        self.values = values

    def getRow(self, headers):
        """
        Convert edge to csv row
        Columns: Source, Target, remaining headers...
        """

        # first rows: Source, Target
        row = [str(self.src), str(self.target)]

        for header in headers:
            row.append(self.values[header])

        return row


class Graph:
    """Graph consists of a list of nodes and edges"""

    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges


# ---------------------------------------------- Script ---------------------------------------------

graphNodeStartPos = []
graphEdgeStartPos = []
graphEndPos = []

# parse file
with open(graphFile) as f:
    lines = f.readlines()

# get lines containing markers
for i, line in enumerate(lines):
    if graphNodeStartStr in line:
        graphNodeStartPos.append(i)
    elif graphEdgeStartStr in line:
        graphEdgeStartPos.append(i)
    elif graphEndStr in line:
        graphEndPos.append(i)

# check that all markers appear equally often
if (not len(graphNodeStartPos) == len(graphEdgeStartPos) == len(graphEndPos)):
    print('Error: number of markers not equal! Node {:d}, Edge {:d}, End {:d}.'.format(len(graphNodeStartPos),
                                                                                       len(graphEdgeStartPos),
                                                                                       len(graphEndPos)))
    sys.exit(1)

numGraphs = len(graphNodeStartPos)
graphPos = [(graphNodeStartPos[i], graphEdgeStartPos[i], graphEndPos[i]) for i in range(numGraphs)]

# check that all markers repeatedly appear in order: node, edge, end
if (not all([node < edge < end for (node, edge, end) in graphPos])) or (
not all([(graphEndPos[i - 1] < graphEdgeStartPos[i]) for i in range(1, numGraphs)])):
    print('Error: unexpected order of markers!')
    sys.exit(1)

# ------ read all graphs ------
graphs = []
nodeHeaders = set()
edgeHeaders = set()
for nodeStart, edgeStart, graphEnd in graphPos:
    print('Graph found from line {:d} to {:d}'.format(nodeStart, graphEnd))

    # ----- read nodes -----
    nodesReader = csv.reader(lines[nodeStart + 1:edgeStart], delimiter=',', quotechar='"')

    # find index labels
    header = next(nodesReader)
    for c, column in enumerate(header):
        if column == labelStr:
            labelIndex = c
        else:
            nodeHeaders.add(column)

    nodes = {}
    for row in nodesReader:
        if (row):
            # get values of current row
            values = {}
            for c, column in enumerate(header):
                value = row[c]

                # get Label seperately. Rest in dict values.
                if c == labelIndex:
                    label = value
                else:
                    values[column] = value

            nodeId = getID(label)
            nodes[nodeId] = GraphNode(label, values)

    # ----- read edges -----
    edgeReader = csv.reader(lines[edgeStart + 1:graphEnd], delimiter=',', quotechar='"')

    # find index of src and target
    header = next(edgeReader)
    for c, column in enumerate(header):
        if column == sourceStr:
            srcIndex = c
        elif column == targetStr:
            targetIndex = c
        else:
            edgeHeaders.add(column)

    edges = {}
    for row in edgeReader:
        if (row):
            # get values of current row
            values = {}
            for c, column in enumerate(header):
                value = row[c]

                # get Source and Target seperately and convert them to unique ids. Rest in dict values.
                if c == srcIndex:
                    src = getID(value)
                elif c == targetIndex:
                    target = getID(value)
                else:
                    values[column] = value

            # let src the node with the lower index
            if (target < src):
                src, target = target, src

            edges[src, target] = GraphEdge(src, target, values)

    graphs.append(Graph(nodes, edges))

# convert header sets to lists
nodeHeaders = list(nodeHeaders)
edgeHeaders = list(edgeHeaders)

# Labels are used configuration: e.g. {Container: LinkedCells , CellSizeFactor: 2.000000 , Traversal: c08 , Data Layout: AoS , Newton 3: disabled}
# Extract header names
labelHeaders = list(getLabelDict(nodes[0].label).keys())

# ---- Write files ----

# create output directory
try:
    os.makedirs(outputDir)
except OSError:
    print('Could not create the output directory: ' + outputDir)
    sys.exit(2)

if (fullSearchFile):
    # if fullSearchFile given, append time found in this file to each row
    nodeHeaders.append(fullSearchTimeStr)
    edgeHeaders.append(fullSearchWeightStr)

    with open(fullSearchFile) as f:
        fullSearchLines = f.readlines()

    for graph in graphs:
        # get runtime in fullsearch for each node
        for nodeId, node in list(graph.nodes.items()):
            # extract label from row and find corresponding runtime in fullSearchFile
            runtime = getStringFromLines(fullSearchLines, node.label + '.* Reduced value: ([0-9]+)')

            # append runtime to row if found otherwise discard row
            if runtime:
                node.values[fullSearchTimeStr] = runtime
            else:
                del graph.nodes[nodeId]

        # generate a weight for each edge using the runtime of src and target
        for edgeId, edge in list(graph.edges.items()):
            # if src or target were removed delete edge instead
            if edge.src in graph.nodes and edge.target in graph.nodes:
                # get runtime of source and target
                sourceTime = float(graph.nodes[edge.src].values[fullSearchTimeStr])
                targetTime = float(graph.nodes[edge.target].values[fullSearchTimeStr])

                # weight = exp(-relativeDifference)
                if sourceTime == targetTime:
                    edge.values[fullSearchWeightStr] = 1
                else:
                    edge.values[fullSearchWeightStr] = math.exp(
                        -abs(sourceTime - targetTime) / max(sourceTime, targetTime))
            else:
                del graph.edges[edgeId]

allNodeHeaders = [idStr, labelStr] + labelHeaders + nodeHeaders
allEdgeHeaders = [sourceStr, targetStr] + edgeHeaders

# write nodes into file
for i, graph in enumerate(graphs):
    with open(os.path.join(outputDir, nodesFile.format(i)), 'w', newline='') as f:
        nodesWriter = csv.writer(f, delimiter=',', quotechar='"')
        nodesWriter.writerow(allNodeHeaders)
        for node in graph.nodes.values():
            nodesWriter.writerow(node.getRow(nodeHeaders, labelHeaders))

# write edges into file
for i, graph in enumerate(graphs):
    with open(os.path.join(outputDir, edgesFile.format(i)), 'w', newline='') as f:
        edgeWriter = csv.writer(f, delimiter=',', quotechar='"')
        edgeWriter.writerow(allEdgeHeaders)
        for edge in graph.edges.values():
            edgeWriter.writerow(edge.getRow(edgeHeaders))
