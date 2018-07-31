#!/usr/bin/gnuplot -p

datafiles = "\
linked_soa_c08 \
linked_soa_sliced \
verlet_soa \
"
titles = "\
'Linked Cells c08' \
'Linked Cells sliced' \
'Verlet Lists' \
"


set xlabel 'Threads'
set ylabel 'MFUPs/s'
set logscale x 2
set logscale y 2

set key autotitle columnheader
set key bottom right
# set key bottom left
# set key Left reverse

fontsize="12"
set tics font ",". fontsize
set key font ",". fontsize
set xlabel font ",". fontsize
set ylabel font ",". fontsize

plot for [i=1:words(datafiles)] \
    word(datafiles, i) \
    using 1:3 \
    with linespoints \
    title word(titles, i)
