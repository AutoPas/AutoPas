#!/usr/bin/gnuplot -p

datafiles = "\
runtimes_DirectSum_AoS_N3on.csv                          \
runtimes_DirectSum_SoA_N3on.csv                          \
runtimes_LinkedCells_AoS_N3on.csv                        \
runtimes_LinkedCells_SoA_N3on.csv                        \
runtimes_VerletLists_AoS_20_0.3_N3on.csv                 \
runtimes_VerletLists_SoA_20_0.3_N3on.csv                 \
runtimes_VerletCells_AoS_20_0.3_N3on.csv                 \
runtimes_VerletCells_SoA_20_0.3_N3on.csv                 \
runtimes_VerletCluster_AoS_20_0.3_N3on.csv               \
runtimes_VerletCluster_SoA_20_0.3_N3on.csv               \
"

titles = "\
'Direct Sum AoS'                                    \
'Direct Sum SoA'                                    \
'Linked Cells AoS'                                  \
'Linked Cells SoA'                                  \
'Verlet Lists AoS'                                  \
'Verlet Lists SoA'                                  \
'Verlet Cells AoS'                                  \
'Verlet Cells SoA'                                  \
'Verlet Cluster AoS'                                \
'Verlet Cluster SoA'                                \
"

# set title "CooLMUC3 256 Threads"

# list of keywords used for coloring (same keyword in title = same color)
colorSelectors = "\
Linked \
Direct \
'Verlet Lists' \
'Verlet Cells' \
"

# list of keywords used for point type (same keyword in title = same point type)
pointSelectors = "\
SoA \
AoS \
"

# Names of the columns to be used for plotting and axis labels
xData = "NumParticles"
yData = "MFUPs/s"
#yData = "Time[micros]"
#yData = "SingleIteration[micros]"

# if your data file has no header (shame on you) this would be 0 else the number of header blocks
dataBlock = 1

# use this if you are not satisfied with the automatically chosen size
#set xrange [30:11000]
#set yrange [0.001:100]

set logscale x 2
set logscale y 2
set grid

set size square

set xlabel xData
set ylabel yData

set key autotitle columnheader
#set key bottom right
set key outside
set key center right
# set key bottom left
#set key top left
set key Left reverse

set xtics rotate by -45

fontsize="12"
set tics font ",". fontsize
set key font ",". fontsize
set xlabel font ",". fontsize
set ylabel font ",". fontsize
set title font "," . fontsize

set xrange [] writeback
set yrange [] writeback

set format x '2^{%L}'
set format y '2^{%L}'

plot for [i=1:words(datafiles)] \
    word(datafiles, i) \
    index dataBlock \
    using xData:yData \
    with linespoints \
    lc rgbcolor system("case '".word(titles, i)."' in \
                            *'".word(colorSelectors, 1)."'*) \
                                echo orange \
                            ;; \
                            *'".word(colorSelectors, 2)."'*) \
                                echo blue \
                            ;; \
                            *'".word(colorSelectors, 3)."'*) \
                                echo dark-green \
                            ;; \
                            *'".word(colorSelectors, 4)."'*) \
                                echo dark-violet \
                            ;; \
                            *) \
                                echo black \
                            ;; \
                         esac") \
    pt (system("case '".word(titles, i)."' in \
                            *'".word(pointSelectors, 1)."'*) \
                                echo 2 \
                            ;; \
                            *'".word(pointSelectors, 2)."'*) \
                                echo 4 \
                            ;; \
                            *'".word(pointSelectors, 3)."'*) \
                                echo 3 \
                            ;; \
                            *'".word(pointSelectors, 4)."'*) \
                                echo 8 \
                            ;; \
                            *) \
                                echo 11 \
                            ;; \
                         esac") + 0) \
    title word(titles, i)

# save margins
#set xrange restore
#set yrange restore
#replot 28000 / x linecolor 'gray'
