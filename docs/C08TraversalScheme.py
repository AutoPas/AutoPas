'''
 * @file C08TraversalScheme.cpp
 * @author C. Menges
 * @date 14.04.2019
'''

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
plots = [fig.add_subplot(1,2,1, projection='3d'),
         fig.add_subplot(1,2,2, projection='3d')]

####################### Old scheme #######################
# o = 0
# x = 1
# y = 2
# z = 3
# xy = 4
# yz = 5
# xz = 6
# xyz = 7
points1 = [[0,0,0], [1,0,0], [0,1,0], [0,0,1], [1,1,0],
         [0,1,1], [1,0,1], [1,1,1]]

dep1 = [[0,0], [0,2], [2,3], [0,3], [0,5], 
        [1,5], [1,2], [1,3], [0,1], [0,4], [4,3], [2,6], [0,6], [0,7]]

###################### New algorithm ######################
overlap = 4
points2 = []
dep2 = []

# generation of 3d coord representing single cells
for x in range(0,overlap+1):
    for y in range(0,overlap+1):
        for z in range(0,overlap+1):
            points2.append([x,y,z])

# generation of interactions between cells (3d coord)
for x in range(0,overlap+1):
    for y in range(0,overlap+1):
        for z in range(0,overlap+1):
            # index of current position
            index = ((overlap+1)**2)*x+y*(overlap+1)
            # origin (front left)
            dep2.append([z, index])
            # back left
            if y != overlap and z!= 0:
                dep2.append([(overlap+1)**2-overlap-1+z, index])
            # front right
            if x != overlap and (y != 0 or z!= 0):
                dep2.append([(overlap+1)**2*overlap+z, index])
            # back right
            if y != overlap and x != overlap and z!= 0:
                dep2.append([(overlap+1)**3-overlap-1+z, index])

###########################################################

for (p, points, dep) in zip(plots, [points1, points2], [dep1, dep2]):
    for (i, l) in enumerate(dep):
        line = [points[i] for i in l]
        [x,y,z] = zip(*line)
        p.plot(x, y, z, label=str(i)+": "+str(line))
    #p.legend()
    p.set_xlabel('x')
    p.set_ylabel('y')
    p.set_zlabel('z')

    locator = plticker.MultipleLocator(1.0)
    p.xaxis.set_major_locator(locator)
    p.yaxis.set_major_locator(locator)
    p.zaxis.set_major_locator(locator)

plt.show()