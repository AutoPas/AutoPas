#!/usr/bin/python3
import itertools

import yaml
from yaml.loader import SafeLoader

template = """
container                        :  [LinkedCells, LinkedCellsReferences, VarVerletListsAsBuild, VerletClusterLists, VerletLists, VerletListsCells, PairwiseVerletLists]
verlet-rebuild-frequency         :  20
verlet-skin-radius               :  0.15
verlet-cluster-size              :  4
selector-strategy                :  Fastest-Absolute-Value
data-layout                      :  [AoS, SoA]
traversal                        :  [lc_sliced, lc_sliced_balanced, lc_sliced_c02, lc_c01, lc_c01_combined_SoA, lc_c01_cuda, lc_c04, lc_c04_HCP, lc_c04_combined_SoA, lc_c08, lc_c18, vcc_cluster_iteration_cuda, vcl_cluster_iteration, vcl_c06, vcl_c01_balanced, vcl_sliced, vcl_sliced_balanced, vcl_sliced_c02, vl_list_iteration, vlc_c01, vlc_c18, vlc_sliced, vlc_sliced_balanced, vlc_sliced_c02, vvl_as_built, vlp_c01, vlp_c18, vlp_sliced, vlp_sliced_balanced, vlp_sliced_c02]
tuning-strategy                  :  full-Search
mpi-strategy                     :  no-mpi
tuning-interval                  :  5000
tuning-samples                   :  3
tuning-max-evidence              :  10
functor                          :  Lennard-Jones (12-6) AVX
newton3                          :  [disabled, enabled]
cutoff                           :  1
box-min                          :  [-1.75, -1.75, -1.75]
box-max                          :  [7.25, 7.25, 7.25]
cell-size                        :  [1]
deltaT                           :  0.0
iterations                       :  10
periodic-boundaries              :  true
Objects:                         
  CubeClosestPacked:
    0:
      box-length                 :  [40, 40, 40]
      bottomLeftCorner           :  [0, 0, 0]
      particle-spacing           :  0.4
      velocity                   :  [0, 0, 0]
      particle-type              :  0
      particle-epsilon           :  1
      particle-sigma             :  1
      particle-mass              :  1
log-level                        :  info
no-flops                         :  true
no-end-config                    :  true
no-progress-bar                  :  true
"""

def make_object(domainSize, numParticles, distribution):
    if distribution[0] == 'uniform':
        objectSize = [dSize * distr for (dSize, distr) in zip(domainSize, distribution[1:4])]
        cubeUniform = {'CubeUniform': {0: {'numberOfParticles': numParticles,
                                                   'box-length': objectSize,
                                                   'bottomLeftCorner': [0, 0, 0],
                                                   'velocity': [0, 0, 0],
                                                   'particle-type': 0,
                                                   'particle-epsilon': 1,
                                                   'particle-sigma': 1,
                                                   'particle-mass': 1}}}
        return cubeUniform

    elif distribution[0] == 'gauss':
        distrMean = [size * distribution[1] for size in domainSize]
        distrStddev = [size * distribution[2] for size in domainSize]
        cubeGauss = {'CubeGauss': {0: {'distribution-mean': distrMean,
                                               'distribution-stddeviation': distrStddev,
                                               'numberOfParticles': numParticles,
                                               'box-length': domainSize.copy(),
                                               'bottomLeftCorner': [0, 0, 0],
                                               'velocity': [0, 0, 0],
                                               'particle-type': 0,
                                               'particle-epsilon': 1,
                                               'particle-sigma': 1,
                                               'particle-mass': 1}}}
        return cubeGauss
    elif distribution[0] == 'closest-packed':
        objectSize = [dSize * distribution[1] for dSize in domainSize]
        objectOffset = [dSize * distribution[2] for dSize in domainSize]
        particleSpacing = ((objectSize[0] * objectSize[1] * objectSize[2] * ((3/4)**(1/2)) * ((2/3)**(1/2))) /
                           (numParticles / 2))**(1/3)
        cubeClosestPacked = {'CubeClosestPacked': {0: {'box-length': objectSize,
                                                               'bottomLeftCorner': objectOffset,
                                                               'particle-spacing': particleSpacing,
                                                               'velocity': [0, 0, 0],
                                                               'particle-type': 0,
                                                               'particle-epsilon': 1,
                                                               'particle-sigma': 1,
                                                               'particle-mass': 1}}}
        return cubeClosestPacked
    else:
        print('Error making object')

def generate(domainSize, numParticles, distribution, cutoff, verletSkinToCutoffFactor, rebuildFrequencySkinFactorFactor,
             functor, cellSizeFactor):
    data = yaml.load(template, Loader=SafeLoader)
    data['box-max'] = domainSize
    data['cutoff'] = cutoff
    skin = cutoff * verletSkinToCutoffFactor
    data['verlet-skin-radius'] = skin
    data['verlet-rebuild-frequency'] = int(rebuildFrequencySkinFactorFactor * verletSkinToCutoffFactor)
    data['functor'] = functor
    data['cell-size'] = [cellSizeFactor]
    data['Objects'] = make_object(domainSize, numParticles, distribution)

    outFilename = str(domainSize) + '-' + str(numParticles) + '-' + str(distribution) + '-' + str(cutoff) + '-' \
                  + str(verletSkinToCutoffFactor) + '-' + str(rebuildFrequencySkinFactorFactor) + '-' \
                  + functor + '-' + str(cellSizeFactor) + '.yaml'
    with open(outFilename, 'w') as fileOutput:
        yaml.dump(data, fileOutput, sort_keys=False)

domainSizes = {'small': [5, 5, 5], 'middle': [25, 25, 25], 'big': [80, 80, 80], 'long': [10, 20, 60]}
particleCounts = {'very-few': 534, 'few': 12193, 'low-normal': 71034, 'high-normal': 230909, 'many': 1092804,
                  'huge': 8923403}
distributions = {'uniform-whole': ('uniform', 1, 1, 1),
                 'uniform-strip': ('uniform', 0.4, 0.6, 1),
                 'concentrated-gauss': ('gauss', 0.75, 0.25),
                 'concentrated-closest-packed': ('closest-packed', 0.3, 0.25)}
cutoffs = {'normal': 1, 'big': 2.5}
verletSkinToCutoffFactors = {'small': 0.1, 'normal': 0.2, 'big': 0.4}
rebuildFrequencySkinFactorFactor = 100
functors = {'lj-no-avx': 'Lennard-Jones (12-6)', 'lj-avx': 'Lennard-Jones (12-6) AVX'}
cellSizeFactors = {'normal': 1.0, 'big': 2.0}
# num threads could also be interesting

def isInteresting(domainSize, numParticles, distribution, cutoff, verletSkinToCutoffFactor,
                  rebuildFrequencySkinFactorFactor, functor, cellSizeFactor):
    tooDenseThreshold = 400 # particles per cell
    actualDomainSize = domainSizes[domainSize]
    numCells = actualDomainSize[0] * actualDomainSize[1] * actualDomainSize[2]
    uniformAvgParticlesPerCell = particleCounts[numParticles] / numCells

    return not (uniformAvgParticlesPerCell > tooDenseThreshold
            or (numParticles == 'very-few' and functor == 'lj-avx')
            #or (cellSizeFactor == 'big')
            #or (cutoff == 'big')
            #or (functor == 'lj-no-avx')
        )

numScenarios = 0
for domainSize, numParticles, distribution, cutoff, verletSkinToCutoffFactor, functor, cellSizeFactor in \
        itertools.product(domainSizes.items(), particleCounts.items(), distributions.items(), cutoffs.items(),
                          verletSkinToCutoffFactors.items(), functors.items(), cellSizeFactors.items()):
    if isInteresting(domainSize[0], numParticles[0], distribution[0], cutoff[0], verletSkinToCutoffFactor[0],
                     rebuildFrequencySkinFactorFactor, functor[0], cellSizeFactor[0]):
        numScenarios += 1
        generate(domainSize[1], numParticles[1], distribution[1], cutoff[1], verletSkinToCutoffFactor[1],
                rebuildFrequencySkinFactorFactor, functor[1], cellSizeFactor[1])

print(str(numScenarios) + ' scenarios generated')