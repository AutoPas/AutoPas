#!/usr/bin/gnuplot -p

datafiles = "\
output-0-soa.txt \
output-0.txt \
output-1-soa.txt \
output-1.txt \
output-2-verlet-10-0.2-soa.txt \
output-2-verlet-10-0.2.txt \
"

titles = "\
'Linked Cells SoA' \
'Linked Cells AoS' \
'Direct Sum SoA' \
'Direct Sum AoS' \
'Verlet Lists SoA rebuild rate 10, skin 0.2*cutoff' \
'Verlet Lists AoS rebuild rate 10, skin 0.2*cutoff' \
"

# list of keywords used for coloring (same keyword in title = same color)
colorSelectors = "\
Linked \
Direct \
Verlet \
"

# list of keywords used for point type (same keyword in title = same point type)
pointSelectors = "\
SoA \
AoS \
"

# Names of the columns to be used for plotting and axis labels
xData = "Number of Molecules"
yData = "MFUPS"

# if your data file has no header (shame on you) this would be 0 else the number of header blocks
dataBlock = 0

# use this if you are not satisfied with the automatically chosen size
#set xrange [30:11000]
#set yrange [0.001:100]

set logscale x 10
set logscale y 10
set grid


set xlabel xData
set ylabel yData

set key autotitle columnheader
set key bottom right
set key bottom left
set key Left reverse

fontsize="12"
set tics font ",". fontsize
set key font ",". fontsize
set xlabel font ",". fontsize
set ylabel font ",". fontsize

set xrange [] writeback
set yrange [] writeback

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
set xrange restore
set yrange restore
replot 28000 / x linecolor 'gray'
