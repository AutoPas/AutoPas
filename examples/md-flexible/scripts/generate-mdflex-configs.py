#!/usr/bin/python3

"""Generates a zoo of configurations for md-flexible

The script generates hundreds of yaml configuration files,
that aim to cover a wide range of possible simulation states.

It is not simply the cross product of the range of parameters below,
but also tries to exclude unfeasible combinations, e.g. huge particle
numbers in tiny domains.

These exclusion criteria may not be appropriate for your test cases.
E.g. they assume a sigma of no smaller than one and a cutoff of no 
larger than 3.

The files are generated in the folder where the script is executed.

"""

# Imports
import sys
# Python version check
if sys.version_info < (3, 8):
    sys.exit("This script requires Python 3.8 or higher for 'math.prod'!")
import math
import yaml

# Make sure the library version is sufficient
pyyaml_version = tuple(map(int, yaml.__version__.split('.')))
if pyyaml_version < (5, 1):
    raise ImportError('You need PyYAML version 5.1 or higher.')
from yaml.loader import SafeLoader

# if available use tqdm for a progress bar otherwise fall back to itertools
try:
    from tqdm.contrib import itertools
except ImportError:
    import itertools


# Script to generate a wide range of yaml input files intended to be used as benchmarks
# to judge the adaptability of tuning strategies.
template = """
functor                          :  Lennard-Jones (12-6) AVX
# container                        :  [all]    # Currently unfeasible. See below.
# Containers that are excluded:
#   - DirectSum             : Takes too long
#   - LinkedCellsReferences : Quasi equivalent to LinkedCells
#   - Octree                : Takes too long
#   - PairwiseVerletLists   : Uses too much memory resulting in chrashes
container                        :  [LinkedCells, VerletLists, VarVerletLists, VerletListsCells, VerletClusterLists]
verlet-rebuild-frequency         :  20
verlet-skin-radius               :  0.15
verlet-cluster-size              :  4
selector-strategy                :  Fastest-Absolute-Value
data-layout                      :  [all]
traversal                        :  [all]
tuning-strategies                :  []
tuning-interval                  :  5000
tuning-samples                   :  3
tuning-phases                    :  1
newton3                          :  [all]
cutoff                           :  1
cell-size                        :  [1]
deltaT                           :  0.0
boundary-type                    :  [periodic, periodic, periodic]
Sites:
  0:
    epsilon                      :  1.
    sigma                        :  1.
    mass                         :  1
Objects:                         
  CubeClosestPacked:
    0:
      box-length                 :  [40, 40, 40]
      bottomLeftCorner           :  [0, 0, 0]
      particle-spacing           :  0.4
      velocity                   :  [0, 0, 0]
      particle-type-id           :  0
log-level                        :  info
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
                                           'particle-type-id': 0,
                                          }}}
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
                                       'particle-type-id': 0,
                                      }}}
        return cubeGauss
    elif distribution[0] == 'closestPacked':
        objectSize = [dSize * distribution[1] for dSize in domainSize]
        objectOffset = [dSize * distribution[2] for dSize in domainSize]
        particleSpacing = ((objectSize[0] * objectSize[1] * objectSize[2] * ((3 / 4) ** (1 / 2)) * (
                (2 / 3) ** (1 / 2))) / (numParticles / 2)) ** (1 / 3)
        cubeClosestPacked = {'CubeClosestPacked': {0: {'box-length': objectSize,
                                                       'bottomLeftCorner': objectOffset,
                                                       'particle-spacing': particleSpacing,
                                                       'velocity': [0, 0, 0],
                                                       'particle-type-id': 0,
                                                      }}}
        return cubeClosestPacked
    else:
        print('Error making object')


def generate(domainSize,
             numParticles,
             distribution,
             cutoff,
             skinFactor,
             invMaxParticleSpeed,
             functor,
             cellSizeFactor):

    """
    Generator for config files.

    :param domainSize: Size of the domain .
    :param numParticles: Number of particles.
    :param distribution: Distribution of particles.
    :param cutoff: Cutoff distance .
    :param skinFactor: Total verlet skin as factor of the cutoff.
    :param invMaxParticleSpeed: Inverse of the (expected) maximum particle speed (for rebuild frequency calculation).
    :param functor: AutoPas Functor.
    :param cellSizeFactor: CSF.
    """

    data = yaml.load(template, Loader=SafeLoader)
    data['box-max'] = domainSize
    data['cutoff'] = cutoff
    data['verlet-skin-radius'] = cutoff * skinFactor
    data['verlet-rebuild-frequency'] = int(invMaxParticleSpeed * skinFactor)
    data['functor'] = functor
    data['cell-size'] = [cellSizeFactor]
    data['Objects'] = make_object(domainSize, numParticles, distribution)

    outFilename = ((str(domainSize)
                    + '-' + str(numParticles)
                    + '-' + str(distribution)
                    + '-' + str(cutoff)
                    + '-' + str(skinFactor)
                    + '-' + str(invMaxParticleSpeed)
                    + '-' + functor
                    + '-' + str(cellSizeFactor)
                    + '.yaml')
                   .replace('\'', '')
                   .replace(', ', '_')
                   .replace('(', '')
                   .replace(')', '')
                   .replace('Lennard-Jones 12-6', 'LJ')
                   .replace(' ', '_'))

    # This writes all flow style values (lists enclosed in []) to block style (list of -).
    # While they are not fully equivalent they are mostly interchangeable and not a problem in our context.
    with open(outFilename, 'w') as fileOutput:
        yaml.dump(data, fileOutput, sort_keys=False)


def isInteresting(domainSizeName,
                  numParticlesName,
                  distributionName,
                  cutoffName,
                  skinFactorName,
                  invMaxParticleSpeed,
                  functor,
                  cellSizeFactorName):

    """
    Filter to judge if a config is worth to generate.
    """

    cellSizeFactor = cellSizeFactors[cellSizeFactorName]
    domainSize = domainSizes[domainSizeName]
    cellSize = cutoffs[cutoffName] * (1 + skinFactors[skinFactorName]) * cellSizeFactor
    cellsPerDimension = [math.floor(d / cellSize) for d in domainSize]
    # If cells are larger than the domain in any dimension the scenario doesn't make sense
    if any(item == 0 for item in cellsPerDimension):
        return False
    numCells = math.prod(cellsPerDimension)
    uniformAvgParticlesPerCell = particleCounts[numParticlesName] / numCells
    # Currently unused:
    # volume=math.prod(domainSize)
    # density=particleCounts[numParticlesName]/volume

    # Exclusion rules. Mostly to avoid scenarios that take forever to test. Warning: These parameters may not be appropriate for your test cases.
    return not (False
                # avoid insanely dense scenarios
                or uniformAvgParticlesPerCell > 150
                # be even stricter for gauss because of their dense core
                or (distributionName == 'concentrated-gauss' and uniformAvgParticlesPerCell > 75)
                # bigger CSFs are for less dense scenarios
                or (uniformAvgParticlesPerCell > 50 and cellSizeFactor > 1)
                # smaller CSFs are for non-sparse scenarios
                or (uniformAvgParticlesPerCell < 10 and cellSizeFactor < 1)
                # for huge particle counts only maximal spread is interesting otherwise there are too many interactions
                or (numParticlesName == 'huge' and distributionName != 'uniform-whole')
                # only avx functor is relevant
                or (not functor == 'lj-avx')
                )


##################################### Start of Script #####################################

if __name__ == "__main__":
    domainSizes = {
        'small': [5, 5, 5],
        'middle': [25, 25, 25],
        'big': [80, 80, 80],
        'long': [10, 20, 60],
    }
    particleCounts = {
        'very-few': 534,
        'few': 12193,
        'low-normal': 71034,
        'high-normal': 230909,
        'many': 1092804,
        'huge': 4923403,
    }
    distributions = {
        # Parameters:
        #   Generator
        #   Object size X/Y/Z as fraction of the domain size in the respective dimension
        'uniform-whole': ('uniform', 1, 1, 1),
        'uniform-strip': ('uniform', 0.4, 0.6, 1),
        # Parameters:
        #   Generator
        #   Distribution mean as fraction of the domain
        #   Distribution Stddev as fraction of the domain
        'concentrated-gauss': ('gauss', 0.75, 0.25),
        # Parameters:
        #   Generator
        #   Object size as fraction of the domain
        #   Object offset as fraction of the domain
        'concentrated-closest-packed': ('closestPacked', 0.3, 0.25),
    }
    cutoffs = {
        'normal': 2.5,
        'big': 3.5,
    }
    skinFactors = {
        'small': 0.1,
        'normal': 0.2,
        'big': 0.4,
    }
    invMaxParticleSpeed = 100
    functors = {
        'lj-no-avx': 'Lennard-Jones (12-6)',
        'lj-avx': 'Lennard-Jones (12-6) AVX',
    }
    cellSizeFactors = {
        'small': 0.5,
        'normal': 1.0,
        'big': 2.0,
    }
    # TODO: num threads could also be interesting

    numScenarios = 0
    # loop over cartesian product of all generator options
    for domainSize, numParticles, distribution, cutoff, skinFactor, functor, cellSizeFactor in \
            itertools.product(domainSizes.items(), particleCounts.items(), distributions.items(), cutoffs.items(),
                              skinFactors.items(), functors.items(), cellSizeFactors.items()):
        # check if this combination is actually interesting to include in the benchmark
        if isInteresting(domainSize[0], numParticles[0], distribution[0], cutoff[0], skinFactor[0],
                         invMaxParticleSpeed, functor[0], cellSizeFactor[0]):
            numScenarios += 1
            generate(domainSize[1], numParticles[1], distribution[1], cutoff[1], skinFactor[1],
                     invMaxParticleSpeed, functor[1], cellSizeFactor[1])

    print(str(numScenarios) + ' scenarios generated')
