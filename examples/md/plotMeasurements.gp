#!/usr/bin/gnuplot -p

####################### All files: #######################

datafiles  = "\
output-0-soa.txt \
output-0.txt \
output-1-soa.txt \
output-1.txt \
output-2-verlet-10-0.2-soa.txt \
output-2-verlet-10-0.2.txt \
output-2-verlet-1-0.0-soa.txt \
output-2-verlet-1-0.0.txt \
output-2-verlet-20-0.3-soa.txt \
output-2-verlet-20-0.3.txt \
output-2-verlet-5-0.1-soa.txt \
output-2-verlet-5-0.1.txt \
"

titles = "\
'Linked Cells AoS' \
'Linked Cells SoA' \
'Direct Sum AoS' \
'Direct Sum SoA' \
'Verlet Lists SoA rebuild rate 10, skin 0.2*cutoff' \
'Verlet Lists AoS rebuild rate 10, skin 0.2*cutoff' \
'Verlet Lists SoA rebuild rate  1, skin 0.0*cutoff' \
'Verlet Lists AoS rebuild rate  1, skin 0.0*cutoff' \
'Verlet Lists SoA rebuild rate 20, skin 0.3*cutoff' \
'Verlet Lists AoS rebuild rate 20, skin 0.3*cutoff' \
'Verlet Lists SoA rebuild rate 5, skin 0.1*cutoff' \
'Verlet Lists AoS rebuild rate 5, skin 0.1*cutoff' \
"

##########################################################

# datafiles  = "\
# output-0-soa.txt \
# output-0.txt \
# output-1-soa.txt \
# output-1.txt \
# output-2-verlet-10-0.2-soa.txt \
# output-2-verlet-10-0.2.txt \
# "

# titles = "\
# 'Linked Cells SoA' \
# 'Linked Cells AoS' \
# 'Direct Sum SoA' \
# 'Direct Sum AoS' \
# 'Verlet Lists SoA rebuild rate 10, skin 0.2*cutoff' \
# 'Verlet Lists AoS rebuild rate 10, skin 0.2*cutoff' \
# "

# set key outside

set xlabel 'NumParticles'
set ylabel 'MFUPS'
set logscale x 10
set logscale y 10

set xrange [30:11000]
set yrange [0.001:100]

set key bottom left
set key Left reverse

fontsize="12"
set tics font ",". fontsize
set key font ",". fontsize
set xlabel font ",". fontsize
set ylabel font ",". fontsize

plot for [i=1:words(datafiles)] \
    word(datafiles, i) \
    using 1:4 \
    with linespoints \
    title word(titles, i)
