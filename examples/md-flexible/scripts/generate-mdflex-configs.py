#!/usr/bin/python3
import yaml
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
container                        :  [all]
verlet-rebuild-frequency         :  20
verlet-skin-radius-per-timestep  :  0.0075
verlet-cluster-size              :  4
selector-strategy                :  Fastest-Absolute-Value
data-layout                      :  [all]
traversal                        :  [all]
tuning-strategies                :  []
tuning-interval                  :  5000
tuning-samples                   :  3
tuning-max-evidence              :  10
newton3                          :  [all]
cutoff                           :  1
cell-size                        :  [1]
deltaT                           :  0.0
iterations                       :  10
boundary-type                    :  [periodic, periodic, periodic]
Sites:
  0:
    epsilon                          :  1.
    sigma                            :  1.
    mass                             :  1
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
    elif distribution[0] == 'closest-packed':
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
             verletSkinToCutoffFactor,
             rebuildFrequencySkinFactorFactor,
             functor,
             cellSizeFactor):
    data = yaml.load(template, Loader=SafeLoader)
    data['box-max'] = domainSize
    data['cutoff'] = cutoff
    skin = cutoff * verletSkinToCutoffFactor
    rebuildFrequency = int(rebuildFrequencySkinFactorFactor * verletSkinToCutoffFactor)
    data['verlet-skin-radius-per-timestep'] = skin / rebuildFrequency
    data['verlet-rebuild-frequency'] = rebuildFrequency
    data['functor'] = functor
    data['cell-size'] = [cellSizeFactor]
    data['Objects'] = make_object(domainSize, numParticles, distribution)

    outFilename = ((str(domainSize)
                    + '-' + str(numParticles)
                    + '-' + str(distribution)
                    + '-' + str(cutoff)
                    + '-' + str(verletSkinToCutoffFactor)
                    + '-' + str(rebuildFrequencySkinFactorFactor)
                    + '-' + functor
                    + '-' + str(cellSizeFactor)
                    + '.yaml')
                   .replace('\'', '')
                   .replace(', ', '_')
                   .replace('(', '')
                   .replace(')', '')
                   .replace('Lennard-Jones 12-6', 'LJ')
                   .replace(' ', '_'))

    with open(outFilename, 'w') as fileOutput:
        yaml.dump(data, fileOutput, sort_keys=False)


def isInteresting(domainSize,
                  numParticles,
                  distribution,
                  cutoff,
                  verletSkinToCutoffFactor,
                  rebuildFrequencySkinFactorFactor,
                  functor,
                  cellSizeFactor):
    tooDenseThreshold = 400  # particles per cell
    actualDomainSize = domainSizes[domainSize]
    numCells = actualDomainSize[0] * actualDomainSize[1] * actualDomainSize[2]
    uniformAvgParticlesPerCell = particleCounts[numParticles] / numCells

    return not (uniformAvgParticlesPerCell > tooDenseThreshold
                or (numParticles == 'very-few' and functor == 'lj-avx')
                or (numParticles == 'huge' and distribution != 'uniform-whole')
                # or (cellSizeFactor == 'big')
                # or (cutoff == 'big')
                # or (functor == 'lj-no-avx')
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
        'uniform-whole': ('uniform', 1, 1, 1),
        'uniform-strip': ('uniform', 0.4, 0.6, 1),
        'concentrated-gauss': ('gauss', 0.75, 0.25),
        'concentrated-closest-packed': ('closest-packed', 0.3, 0.25),
    }
    cutoffs = {
        'normal': 1,
        'big': 2.5,
    }
    verletSkinToCutoffFactors = {
        'small': 0.1,
        'normal': 0.2,
        'big': 0.4,
    }
    rebuildFrequencySkinFactorFactor = 100
    functors = {
        'lj-no-avx': 'Lennard-Jones (12-6)',
        'lj-avx': 'Lennard-Jones (12-6) AVX',
    }
    cellSizeFactors = {
        'normal': 1.0,
        'big': 2.0,
    }
    # TODO: num threads could also be interesting

    numScenarios = 0
    # loop over cartesian product of all generator options
    for domainSize, numParticles, distribution, cutoff, verletSkinToCutoffFactor, functor, cellSizeFactor in \
            itertools.product(domainSizes.items(), particleCounts.items(), distributions.items(), cutoffs.items(),
                              verletSkinToCutoffFactors.items(), functors.items(), cellSizeFactors.items()):
        # check if this combination is actually interesting to include in the benchmark
        if isInteresting(domainSize[0], numParticles[0], distribution[0], cutoff[0], verletSkinToCutoffFactor[0],
                         rebuildFrequencySkinFactorFactor, functor[0], cellSizeFactor[0]):
            numScenarios += 1
            generate(domainSize[1], numParticles[1], distribution[1], cutoff[1], verletSkinToCutoffFactor[1],
                     rebuildFrequencySkinFactorFactor, functor[1], cellSizeFactor[1])

    print(str(numScenarios) + ' scenarios generated')
