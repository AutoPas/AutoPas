#!/usr/bin/python3
import itertools

import yaml
from yaml.loader import SafeLoader

verletSkins = [0.1, 0.2, 0.4]
functors = ['Lennard-Jones (12-6)', 'Lennard-Jones (12-6) AVX']
cutoffs = [1, 1.5, 2]
boxMaxs = [[5, 5, 5], [80, 80, 80]]
cellSizes = [[1], [2]]
objects = []

particleGenShapes = [[5, 5, 5], [20, 20, 20], [50, 50, 50], [4, 50, 4], [7.2, 3.4, 43.1]]

cubeClosestPackedTemplate = {'CubeClosestPacked': {0: {'box-length': [40, 40, 40],
                                                       'bottomLeftCorner': [0, 0, 0],
                                                       'particle-spacing': 0.4,
                                                       'velocity': [0, 0, 0],
                                                       'particle-type': 0,
                                                       'particle-epsilon': 1,
                                                       'particle-sigma': 1,
                                                       'particle-mass': 1}}}
for boxLength, particleSpacing in itertools.product(particleGenShapes, [0.4, 0.6, 0.8]):
    cubeClosestPackedTemplate['CubeClosestPacked'][0]['box-length'] = boxLength
    cubeClosestPackedTemplate['CubeClosestPacked'][0]['particle-spacing'] = particleSpacing
    objects.append({'content': cubeClosestPackedTemplate,
                    'name': 'CubeClosestPacked' + str(boxLength) + str(particleSpacing)})

cubeUniformTemplate = {'CubeUniform': {0: {'numberOfParticles': 100,
                                           'box-length': [4, 4, 4],
                                           'bottomLeftCorner': [0, 0, 0],
                                           'velocity': [0, 0, 0],
                                           'particle-type': 0,
                                           'particle-epsilon': 1,
                                           'particle-sigma': 1,
                                           'particle-mass': 1}}}

for boxLength, numberOfParticles in itertools.product([[5, 5, 5], [50, 50, 50], [10, 60, 60]], [5000, 50000, 500000,
                                                                                                1000000]):
    if numberOfParticles > 1000000 and boxLength != [50, 50, 50]:
        continue
    cubeUniformTemplate['CubeUniform'][0]['numberOfParticles'] = numberOfParticles
    cubeUniformTemplate['CubeUniform'][0]['box-length'] = boxLength
    objects.append({'content': cubeUniformTemplate,
                    'name': 'CubeUniform' + str(boxLength) + str(numberOfParticles)})


cubeGaussTemplate = {'CubeGauss': {0: {'distribution-mean': [20, 20, 20],
                                       'distribution-stddeviation': [10, 10, 10],
                                       'numberOfParticles': 100,
                                       'box-length': [40, 40, 40],
                                       'bottomLeftCorner': [0, 0, 0],
                                       'velocity': [0, 0, 0],
                                       'particle-type': 0,
                                       'particle-epsilon': 1,
                                       'particle-sigma': 1,
                                       'particle-mass': 1}}}

for numberOfParticles, stddev in itertools.product([150000, 750000], [[10, 10, 10], [16, 8, 8]]):
    cubeGaussTemplate['CubeGauss'][0]['numberOfParticles'] = numberOfParticles
    cubeGaussTemplate['CubeGauss'][0]['distribution-stddeviation'] = stddev
    objects.append({'content': cubeGaussTemplate,
                    'name': 'CubeGauss' + str(stddev) + str(numberOfParticles)})


for verletSkin, functor, cutoff, boxMax, cellSize, object in \
        itertools.product(verletSkins, functors, cutoffs, boxMaxs, cellSizes, objects):
    with open('/home/tobias/tmp/tuning-data/template.yaml') as fileInput:
        data = yaml.load(fileInput, Loader=SafeLoader)
        data['verlet-skin-radius'] = verletSkin
        data['functor'] = functor
        data['cutoff'] = cutoff
        data['box-max'] = boxMax
        data['cell-size'] = cellSize
        data['Objects'] = object['content']

        outFilename = str(verletSkin) + '-' + str(functor) + '-' + str(cutoff) + str(boxMax) + '-' + str(cellSize) \
                      + '-' + object['name'] + '.yaml'
        with open(outFilename, 'w') as fileOutput:
            yaml.dump(data, fileOutput, sort_keys=False)
